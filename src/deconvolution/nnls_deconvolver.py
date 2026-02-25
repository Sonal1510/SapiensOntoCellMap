#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2026-02-22
Description :
NNLS-based cell type deconvolution for spatial (Visium, Xenium) and bulk RNA-seq.

Two complementary modes:

  MarkerDBDeconvolver  — Reference-free. Uses the SapiensOntoCellMap 14-database
                         marker knowledge base as the signature matrix. No scRNA-seq
                         reference required. Novel: no existing tool does this.

  ReferenceDeconvolver — Reference-based. Builds pseudobulk mean-expression signatures
                         from a provided scRNA-seq reference (AnnData / DataFrame).
                         Mirrors RCTD's core algorithm in pure Python/scipy.

Both classes share:
  - scipy.optimize.nnls solver per spot
  - log1p + L2 normalization of input expression (opt-in)
  - Shared gene alignment via SignatureMatrix.subset_genes()
  - Identical output format: pd.DataFrame (spots x cell_types), proportions sum to 1

Usage example:
    from src.deconvolution import MarkerDBDeconvolver, ReferenceDeconvolver

    # Reference-free
    m = MarkerDBDeconvolver(marker_db_path='data/processed_combined_db/master_cell_marker_db.csv',
                            tissue_filter='UBERON_0002097')
    m.build_signature()
    props = m.deconvolve(expr_df)   # expr_df: spots x genes DataFrame

    # Reference-based
    r = ReferenceDeconvolver(reference=ref_df, cell_type_col='cell_type')
    r.build_signature()
    props = r.deconvolve(expr_df)
"""

import logging

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from tqdm import tqdm

from .signature_builder import SignatureMatrix, build_signature_matrix

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _nnls_proportions(sig_matrix, expr_vector):
    """
    Solve NNLS for a single spot and return proportions normalized to sum to 1.

    Args:
        sig_matrix  : np.ndarray (n_genes, n_cell_types)
        expr_vector : np.ndarray (n_genes,)

    Returns:
        np.ndarray (n_cell_types,) — non-negative proportions summing to 1.
        All-zero vector if NNLS finds no solution (expression all zero).
    """
    props, _ = nnls(sig_matrix, expr_vector)
    total = props.sum()
    if total > 0:
        props = props / total
    return props.astype(np.float32)


def _to_dense(arr):
    """Convert sparse matrix to dense numpy array if needed."""
    try:
        from scipy.sparse import issparse
        if issparse(arr):
            return arr.toarray()
    except ImportError:
        pass
    return np.asarray(arr)


def _align_expression(expr, gene_names, sig_sub):
    """
    Extract and align expression columns to the genes in sig_sub.

    Returns:
        np.ndarray (n_spots, n_shared_genes) aligned to sig_sub.genes
    """
    gene_names_upper = [g.upper() for g in gene_names]
    shared_idx = [i for i, g in enumerate(gene_names_upper) if g in sig_sub.gene_index]
    if not shared_idx:
        raise ValueError(
            f"No shared genes between expression matrix and signature. "
            f"Check gene symbol format (should be HGNC approved, e.g. 'CD3D' not 'ENSG00000...')."
        )
    return expr[:, shared_idx]


def _normalize_expression(expr):
    """log1p + L2-normalize each row (spot) of the expression matrix."""
    expr = np.log1p(expr)
    norms = np.linalg.norm(expr, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return expr / norms


def _parse_expression_input(expression_matrix, spot_ids, gene_names):
    """
    Normalize expression input to (dense ndarray, spot_ids list, gene_names list).

    Accepts:
      - pd.DataFrame (spots x genes) — index = spot_ids, columns = gene_names
      - np.ndarray + explicit spot_ids + gene_names args
    """
    if isinstance(expression_matrix, pd.DataFrame):
        spot_ids = list(expression_matrix.index)
        gene_names = list(expression_matrix.columns)
        expr = _to_dense(expression_matrix.values).astype(np.float32)
    else:
        expr = _to_dense(expression_matrix).astype(np.float32)
        if gene_names is None:
            raise ValueError(
                "gene_names must be provided when expression_matrix is not a DataFrame."
            )
        gene_names = list(gene_names)
        spot_ids = list(spot_ids) if spot_ids is not None else list(range(expr.shape[0]))

    return expr, spot_ids, gene_names


# ---------------------------------------------------------------------------
# Public classes
# ---------------------------------------------------------------------------

class MarkerDBDeconvolver:
    """
    Reference-free cell type deconvolution using the SapiensOntoCellMap marker DB.

    Builds a gene x cell_type evidence-weighted signature matrix from
    master_cell_marker_db.csv, optionally filtered by UBERON tissue ID.
    Applies per-spot NNLS to estimate cell type proportions.

    No scRNA-seq reference data required — a key novel advantage over all
    existing deconvolution tools (RCTD, Cell2location, Spotlight, CARD).

    Args:
        marker_db_path    : str   — path to master_cell_marker_db.csv
        tissue_filter     : str or list[str] — UBERON ID(s) to restrict cell types.
                            E.g. 'UBERON_0002097' (skin of body). None = all tissues.
        min_source_weight : float — minimum evidence weight to include.
                            Use 2.0 to exclude Computational (0.5) and Company (1.0).
        normalize_cols    : bool  — L1-normalize signature columns (default True).
    """

    def __init__(
        self,
        marker_db_path,
        tissue_filter=None,
        min_source_weight=None,
        normalize_cols=True,
    ):
        self.marker_db_path = marker_db_path
        self.tissue_filter = tissue_filter
        self.min_source_weight = min_source_weight
        self.normalize_cols = normalize_cols
        self.signature_ = None
        self.proportions_ = None

    def build_signature(self):
        """
        Load master DB and build the gene x cell_type signature matrix.
        Must be called before deconvolve().

        Returns:
            self (for method chaining)
        """
        self.signature_ = build_signature_matrix(
            marker_db_path=self.marker_db_path,
            tissue_filter=self.tissue_filter,
            min_source_weight=self.min_source_weight,
            normalize_cols=self.normalize_cols,
        )
        return self

    def deconvolve(
        self,
        expression_matrix,
        spot_ids=None,
        gene_names=None,
        normalize_input=True,
    ):
        """
        Deconvolve spot expression into cell type proportions.

        Args:
            expression_matrix : pd.DataFrame (spots x genes) or np.ndarray.
                                Rows = spots/samples, columns = gene symbols.
            spot_ids          : list — row identifiers (auto-set from DataFrame index).
            gene_names        : list — gene symbols (auto-set from DataFrame columns).
            normalize_input   : bool — apply log1p + L2 normalization per spot (default True).

        Returns:
            pd.DataFrame (n_spots x n_cell_types) — proportions summing to 1 per row.
            Also stored in self.proportions_.
        """
        if self.signature_ is None:
            raise RuntimeError("Call build_signature() before deconvolve().")

        expr, spot_ids, gene_names = _parse_expression_input(
            expression_matrix, spot_ids, gene_names
        )
        n_spots = expr.shape[0]
        logger.info(f"Input: {n_spots} spots x {len(gene_names)} genes")

        sig_sub = self.signature_.subset_genes(gene_names)
        expr_aligned = _align_expression(expr, gene_names, sig_sub)

        if normalize_input:
            expr_aligned = _normalize_expression(expr_aligned)

        logger.info(
            f"Running NNLS deconvolution: "
            f"{n_spots} spots x {len(sig_sub.cell_types)} cell types "
            f"({len(sig_sub.genes)} shared genes)"
        )

        proportions = np.zeros((n_spots, len(sig_sub.cell_types)), dtype=np.float32)
        for i in tqdm(range(n_spots), desc="Deconvolving spots", unit="spot"):
            proportions[i] = _nnls_proportions(sig_sub.matrix, expr_aligned[i])

        self.proportions_ = pd.DataFrame(
            proportions,
            index=spot_ids,
            columns=sig_sub.cell_types,
        )
        logger.info(f"Done. Output: {self.proportions_.shape}")
        return self.proportions_

    def top_cell_types(self, n=5):
        """
        Return the top-n cell types by mean proportion across all spots.

        Returns:
            pd.Series — cell type label -> mean proportion, descending.
        """
        if self.proportions_ is None:
            raise RuntimeError("Call deconvolve() first.")
        return self.proportions_.mean(axis=0).sort_values(ascending=False).head(n)


class ReferenceDeconvolver:
    """
    Reference-based cell type deconvolution using pseudobulk signatures.

    Accepts a scRNA-seq reference as a pd.DataFrame (cells x genes) with a
    cell type label column, or as a pre-computed dict {cell_type: mean_vector}.
    Builds pseudobulk mean-expression signatures per cell type and applies NNLS.

    This implements the core logic of RCTD in pure Python/scipy — directly
    comparable to the reference-based gold standard without requiring R.

    Args:
        reference      : pd.DataFrame (cells x genes + label column) or
                         dict {cell_type_str: pd.Series/np.ndarray of mean expr}
        cell_type_col  : str  — column in DataFrame with cell type labels.
        normalize_cols : bool — L1-normalize pseudobulk columns (default True).
    """

    def __init__(self, reference, cell_type_col="cell_type", normalize_cols=True):
        self.reference = reference
        self.cell_type_col = cell_type_col
        self.normalize_cols = normalize_cols
        self.signature_ = None
        self.proportions_ = None

    def build_signature(self):
        """
        Build pseudobulk signature matrix from the reference.
        Must be called before deconvolve().

        Returns:
            self (for method chaining)
        """
        if isinstance(self.reference, dict):
            cell_types = list(self.reference.keys())
            first = next(iter(self.reference.values()))
            if isinstance(first, pd.Series):
                genes = [g.upper() for g in first.index]
                matrix = np.column_stack(
                    [self.reference[ct].values for ct in cell_types]
                ).astype(np.float32)
            else:
                genes = [f"gene_{i}" for i in range(len(first))]
                matrix = np.column_stack(
                    [np.asarray(self.reference[ct]) for ct in cell_types]
                ).astype(np.float32)

        elif isinstance(self.reference, pd.DataFrame):
            if self.cell_type_col not in self.reference.columns:
                raise ValueError(
                    f"cell_type_col '{self.cell_type_col}' not found in reference. "
                    f"Available columns: {list(self.reference.columns)}"
                )
            logger.info(
                f"Building pseudobulk signatures from {len(self.reference)} reference cells"
            )
            labels = self.reference[self.cell_type_col]
            expr_cols = [c for c in self.reference.columns if c != self.cell_type_col]
            genes = [g.upper() for g in expr_cols]

            ref_expr = _to_dense(self.reference[expr_cols].values).astype(np.float32)
            cell_types = sorted(labels.unique().tolist())

            matrix = np.zeros((len(genes), len(cell_types)), dtype=np.float32)
            for j, ct in enumerate(cell_types):
                mask = (labels == ct).values
                n_cells = mask.sum()
                if n_cells == 0:
                    logger.warning(f"Cell type '{ct}' has 0 cells in reference — skipping.")
                    continue
                matrix[:, j] = ref_expr[mask].mean(axis=0)
                logger.debug(f"  {ct}: {n_cells} cells -> pseudobulk computed")

            logger.info(
                f"Pseudobulk signature: {len(genes):,} genes x {len(cell_types)} cell types"
            )
        else:
            raise TypeError(
                "reference must be a pd.DataFrame (cells x genes + label col) "
                "or dict {cell_type: mean_expression_vector}."
            )

        if self.normalize_cols:
            col_sums = matrix.sum(axis=0, keepdims=True)
            col_sums[col_sums == 0] = 1.0
            matrix = matrix / col_sums

        self.signature_ = SignatureMatrix(
            matrix=matrix,
            cell_types=cell_types,
            genes=genes,
        )
        logger.info(f"Reference signature built: {self.signature_}")
        return self

    def deconvolve(
        self,
        expression_matrix,
        spot_ids=None,
        gene_names=None,
        normalize_input=True,
    ):
        """
        Deconvolve spot expression into cell type proportions.
        Identical interface to MarkerDBDeconvolver.deconvolve().

        Args:
            expression_matrix : pd.DataFrame (spots x genes) or np.ndarray.
            spot_ids          : list — row identifiers (auto-set from DataFrame index).
            gene_names        : list — gene symbols (auto-set from DataFrame columns).
            normalize_input   : bool — apply log1p + L2 normalization per spot.

        Returns:
            pd.DataFrame (n_spots x n_cell_types) — proportions summing to 1 per row.
        """
        if self.signature_ is None:
            raise RuntimeError("Call build_signature() before deconvolve().")

        expr, spot_ids, gene_names = _parse_expression_input(
            expression_matrix, spot_ids, gene_names
        )
        n_spots = expr.shape[0]
        logger.info(f"Input: {n_spots} spots x {len(gene_names)} genes")

        sig_sub = self.signature_.subset_genes(gene_names)
        expr_aligned = _align_expression(expr, gene_names, sig_sub)

        if normalize_input:
            expr_aligned = _normalize_expression(expr_aligned)

        logger.info(
            f"Running NNLS deconvolution: "
            f"{n_spots} spots x {len(sig_sub.cell_types)} cell types "
            f"({len(sig_sub.genes)} shared genes)"
        )

        proportions = np.zeros((n_spots, len(sig_sub.cell_types)), dtype=np.float32)
        for i in tqdm(range(n_spots), desc="Deconvolving spots", unit="spot"):
            proportions[i] = _nnls_proportions(sig_sub.matrix, expr_aligned[i])

        self.proportions_ = pd.DataFrame(
            proportions,
            index=spot_ids,
            columns=sig_sub.cell_types,
        )
        logger.info(f"Done. Output: {self.proportions_.shape}")
        return self.proportions_

    def top_cell_types(self, n=5):
        """Return top-n cell types by mean proportion across all spots."""
        if self.proportions_ is None:
            raise RuntimeError("Call deconvolve() first.")
        return self.proportions_.mean(axis=0).sort_values(ascending=False).head(n)
