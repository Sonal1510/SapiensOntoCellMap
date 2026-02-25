#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2026-02-22
Description :
Build a gene x cell_type evidence-weighted signature matrix from the
SapiensOntoCellMap master marker database, optionally filtered by tissue (UBERON).

Evidence weighting per (gene, cell_type) row:
  source_type_weight : Experiment=4, Single-Cell Sequencing=3,
                       Literature=2, Review=2, Company=1, Computational=0.5
  cross_db_multiplier: 1 + log2(n_unique_databases reporting this gene-cell pair)
  final_weight       : source_type_weight x cross_db_multiplier

Aggregation: max evidence weight across all rows for each (gene, cell_type) pair.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Evidence weights by source_type (controlled vocabulary from master DB)
SOURCE_TYPE_WEIGHTS = {
    "Experiment": 4.0,
    "Single-Cell Sequencing": 3.0,
    "Literature": 2.0,
    "Review": 2.0,
    "Company": 1.0,
    "Computational": 0.5,
}
DEFAULT_SOURCE_WEIGHT = 1.0


class SignatureMatrix:
    """
    Container for a gene x cell_type signature matrix with associated metadata.

    Attributes:
        matrix      : np.ndarray, shape (n_genes, n_cell_types), float32
        cell_types  : list[str]  — cell type labels (e.g. "macrophage (CL_0000235)")
        genes       : list[str]  — gene symbols, uppercase
        tissue_filter: UBERON filter used when building this matrix (or None)
        n_raw_rows  : int        — number of rows in master DB before filtering
        gene_index  : dict       — {gene_upper: row_index}
        cell_index  : dict       — {cell_label: col_index}
    """

    def __init__(self, matrix, cell_types, genes, tissue_filter=None, n_raw_rows=None):
        self.matrix = matrix
        self.cell_types = list(cell_types)
        self.genes = list(genes)
        self.tissue_filter = tissue_filter
        self.n_raw_rows = n_raw_rows
        self.gene_index = {g: i for i, g in enumerate(genes)}
        self.cell_index = {c: i for i, c in enumerate(cell_types)}

    def __repr__(self):
        return (
            f"SignatureMatrix({len(self.genes)} genes x {len(self.cell_types)} cell types, "
            f"tissue_filter={self.tissue_filter!r})"
        )

    def subset_to_cell_type_labels(self, labels):
        """
        Return a new SignatureMatrix restricted to the specified cell type column labels.

        Args:
            labels: list[str] — cell type labels to keep (matched against self.cell_types).

        Returns:
            SignatureMatrix with only the matching columns (genes unchanged).

        Raises:
            ValueError if no label matches.
        """
        label_set = set(labels)
        keep_idx = [i for i, ct in enumerate(self.cell_types) if ct in label_set]
        if not keep_idx:
            raise ValueError(
                f"None of the {len(labels)} requested labels match the signature's "
                f"{len(self.cell_types)} cell types."
            )
        logger.info(
            f"Cell type subset: {len(keep_idx)}/{len(labels)} requested labels matched "
            f"(of {len(self.cell_types)} total)."
        )
        new_cell_types = [self.cell_types[i] for i in keep_idx]
        new_matrix = self.matrix[:, keep_idx]
        return SignatureMatrix(
            matrix=new_matrix,
            cell_types=new_cell_types,
            genes=self.genes,
            tissue_filter=self.tissue_filter,
        )

    def subset_genes(self, target_genes):
        """
        Return a new SignatureMatrix restricted to genes present in target_genes.

        Args:
            target_genes: list[str] — gene names (case-insensitive)

        Returns:
            SignatureMatrix with only the overlapping genes (rows).

        Raises:
            ValueError if there is no gene overlap.
        """
        target_upper = [g.upper() for g in target_genes]
        keep_pairs = [(i, g) for i, g in enumerate(target_upper) if g in self.gene_index]

        if not keep_pairs:
            raise ValueError(
                f"No overlap between the {len(target_upper)} input genes and the "
                f"{len(self.genes)} genes in the signature matrix. "
                f"Check that gene symbols are HGNC-approved and not Ensembl IDs."
            )

        logger.info(
            f"Gene intersection: {len(keep_pairs)}/{len(target_upper)} input genes "
            f"found in signature ({len(self.genes)} total)."
        )

        sig_row_idx = [self.gene_index[g] for _, g in keep_pairs]
        new_genes = [g for _, g in keep_pairs]
        new_matrix = self.matrix[sig_row_idx, :]

        return SignatureMatrix(
            matrix=new_matrix,
            cell_types=self.cell_types,
            genes=new_genes,
            tissue_filter=self.tissue_filter,
        )


def build_signature_matrix(
    marker_db_path,
    tissue_filter=None,
    min_source_weight=None,
    normalize_cols=True,
    use_cell_names=True,
):
    """
    Build a gene x cell_type evidence-weighted signature matrix from the master DB.

    Args:
        marker_db_path   : str  — path to master_cell_marker_db.csv
        tissue_filter    : str or list[str] — UBERON ID(s) to restrict cell types
                           (e.g. 'UBERON_0002097'). None = all tissues.
        min_source_weight: float or None — drop entries below this threshold.
                           Use 2.0 to exclude Computational (0.5) and Company (1.0).
        normalize_cols   : bool — L1-normalize each cell type column to sum to 1.
        use_cell_names   : bool — label columns as "cell_name (cell_id)" if True,
                           else use cell_id only.

    Returns:
        SignatureMatrix
    """
    logger.info(f"Loading marker database: {marker_db_path}")
    usecols = ["gene", "cell_id", "cell_name", "tissue_id", "source_type", "database"]
    df = pd.read_csv(marker_db_path, usecols=usecols, dtype=str)
    n_raw = len(df)
    logger.info(f"Loaded {n_raw:,} rows.")

    # Normalize gene symbols
    df["gene"] = df["gene"].str.strip().str.upper()
    df = df.dropna(subset=["gene", "cell_id"])

    # Optional tissue filter
    if tissue_filter is not None:
        if isinstance(tissue_filter, str):
            tissue_filter = [tissue_filter]
        # Normalize separators (UBERON_xxx or UBERON:xxx both accepted)
        norm_filter = {t.replace(":", "_") for t in tissue_filter}
        df["_tissue_norm"] = df["tissue_id"].str.replace(":", "_")
        df = df[df["_tissue_norm"].isin(norm_filter)].drop(columns=["_tissue_norm"])
        logger.info(f"After tissue filter {tissue_filter}: {len(df):,} rows")
        if df.empty:
            raise ValueError(
                f"No rows remain after tissue filter {tissue_filter}. "
                f"Verify UBERON IDs against the database tissue_id column."
            )

    # Map source_type to evidence weight
    df["source_weight"] = (
        df["source_type"].map(SOURCE_TYPE_WEIGHTS).fillna(DEFAULT_SOURCE_WEIGHT)
    )

    if min_source_weight is not None:
        before = len(df)
        df = df[df["source_weight"] >= min_source_weight]
        logger.info(
            f"min_source_weight={min_source_weight}: kept {len(df):,}/{before:,} rows"
        )

    # Cross-DB agreement multiplier: 1 + log2(n_unique_databases per gene-cell pair)
    db_counts = (
        df.groupby(["gene", "cell_id"])["database"]
        .nunique()
        .rename("n_db")
        .reset_index()
    )
    df = df.merge(db_counts, on=["gene", "cell_id"], how="left")
    df["cross_db_mult"] = 1.0 + np.log2(df["n_db"].clip(lower=1))
    df["evidence_weight"] = df["source_weight"] * df["cross_db_mult"]

    # Aggregate: max evidence weight per (gene, cell_id)
    agg = (
        df.groupby(["gene", "cell_id", "cell_name"])["evidence_weight"]
        .max()
        .reset_index()
    )

    # Build cell type labels
    if use_cell_names:
        agg["cell_label"] = (
            agg["cell_name"].str.strip() + " (" + agg["cell_id"].str.strip() + ")"
        )
    else:
        agg["cell_label"] = agg["cell_id"].str.strip()

    # Pivot to gene x cell_type matrix
    pivot = agg.pivot_table(
        index="gene", columns="cell_label", values="evidence_weight", aggfunc="max"
    )
    pivot = pivot.fillna(0.0)
    pivot.columns.name = None

    genes = list(pivot.index)
    cell_types = list(pivot.columns)
    matrix = pivot.values.astype(np.float32)

    # L1-normalize columns (each cell type signature sums to 1)
    if normalize_cols:
        col_sums = matrix.sum(axis=0, keepdims=True)
        col_sums[col_sums == 0] = 1.0
        matrix = matrix / col_sums

    logger.info(f"Signature matrix built: {len(genes):,} genes x {len(cell_types):,} cell types")
    return SignatureMatrix(
        matrix=matrix,
        cell_types=cell_types,
        genes=genes,
        tissue_filter=tissue_filter,
        n_raw_rows=n_raw,
    )
