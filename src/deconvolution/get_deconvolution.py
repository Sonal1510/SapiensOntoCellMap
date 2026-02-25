#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2026-02-22
Description :
CLI entry point for SapiensOntoCellMap cell type deconvolution.

Supports two modes (selectable via --mode):
  marker_db  — Reference-free NNLS using the SapiensOntoCellMap marker DB (novel).
  reference  — Reference-based pseudobulk NNLS from a scRNA-seq reference.
  both       — Run both and save both outputs for direct comparison.

Usage:
  python -m src.deconvolution.get_deconvolution \\
      --expression_matrix spots_x_genes.csv \\
      --output_dir ./deconv_out/ \\
      --tissue UBERON_0002097 \\
      --mode marker_db

  python -m src.deconvolution.get_deconvolution \\
      --expression_matrix spots_x_genes.csv \\
      --output_dir ./deconv_out/ \\
      --reference scrnaseq_ref.csv \\
      --cell_type_col cell_type \\
      --mode both

Input CSV formats:
  expression_matrix : rows = spot/sample IDs (index col 0), columns = gene symbols
  reference         : rows = cell IDs (index col 0), columns = gene symbols + one
                      column matching --cell_type_col

Output files:
  proportions_marker_db.csv  — spots x cell_types, values = proportions (sum to 1)
  proportions_reference.csv  — same format for reference-based mode
  summary_marker_db.csv      — mean proportion per cell type, sorted descending
  summary_reference.csv      — same for reference-based mode
"""

import argparse
import logging
import os
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _find_default_db():
    """Auto-locate master_cell_marker_db.csv relative to this file."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.normpath(
        os.path.join(here, "..", "..", "data", "processed_combined_db", "master_cell_marker_db.csv")
    )
    return candidate if os.path.exists(candidate) else None


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="SapiensOntoCellMap: Cell type deconvolution",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--expression_matrix", required=True,
        help="CSV: spots x genes. Row index = spot IDs, columns = gene symbols.",
    )
    p.add_argument(
        "--output_dir", required=True,
        help="Directory to write output CSV files.",
    )
    p.add_argument(
        "--marker_db", default=None,
        help="Path to master_cell_marker_db.csv. Auto-detected if omitted.",
    )
    p.add_argument(
        "--tissue", default=None,
        help=(
            "UBERON tissue ID(s) to restrict cell types. "
            "E.g. 'UBERON_0002097' or comma-separated list. "
            "Omit to use all tissues (slower, more cell types)."
        ),
    )
    p.add_argument(
        "--min_source_weight", type=float, default=None,
        help=(
            "Minimum evidence weight. "
            "Experiment=4, Single-Cell Sequencing=3, Literature=2, "
            "Company=1, Computational=0.5. "
            "Use 2.0 to exclude Computational and Company sources."
        ),
    )
    p.add_argument(
        "--mode", choices=["marker_db", "reference", "both"], default="marker_db",
        help="Deconvolution mode.",
    )
    p.add_argument(
        "--reference", default=None,
        help=(
            "scRNA-seq reference CSV: cells x genes + cell type label column. "
            "Required for --mode reference or both."
        ),
    )
    p.add_argument(
        "--cell_type_col", default="cell_type",
        help="Column name in reference CSV containing cell type labels.",
    )
    p.add_argument(
        "--no_normalize_input", action="store_true",
        help="Skip log1p + L2 normalization of input expression (not recommended).",
    )
    return p.parse_args(argv)


def _save_outputs(props, prefix, output_dir):
    """Save full proportion matrix and per-cell-type summary."""
    props_path = os.path.join(output_dir, f"proportions_{prefix}.csv")
    props.to_csv(props_path)
    logger.info(f"  Proportions saved -> {props_path}")

    summary = props.mean(axis=0).sort_values(ascending=False).rename("mean_proportion")
    summary_path = os.path.join(output_dir, f"summary_{prefix}.csv")
    summary.to_csv(summary_path, header=True)
    logger.info(f"  Summary saved     -> {summary_path}")

    # Print top-10 to console
    top10 = summary.head(10)
    logger.info(f"  Top 10 cell types ({prefix}):")
    for ct, prop in top10.items():
        logger.info(f"    {prop:.4f}  {ct}")


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)

    # Load expression matrix
    logger.info(f"Loading expression matrix: {args.expression_matrix}")
    expr_df = pd.read_csv(args.expression_matrix, index_col=0)
    logger.info(f"Expression matrix: {expr_df.shape[0]} spots x {expr_df.shape[1]} genes")

    normalize_input = not args.no_normalize_input
    tissue_filter = (
        [t.strip() for t in args.tissue.split(",")]
        if args.tissue else None
    )

    # --- Mode: marker_db (reference-free, novel) ---
    if args.mode in ("marker_db", "both"):
        from .nnls_deconvolver import MarkerDBDeconvolver

        marker_db = args.marker_db or _find_default_db()
        if not marker_db or not os.path.exists(marker_db):
            logger.error(
                "Cannot find master_cell_marker_db.csv. "
                "Specify --marker_db path."
            )
            sys.exit(1)

        logger.info("=== Reference-free marker-DB deconvolution ===")
        deconv = MarkerDBDeconvolver(
            marker_db_path=marker_db,
            tissue_filter=tissue_filter,
            min_source_weight=args.min_source_weight,
        )
        deconv.build_signature()
        props = deconv.deconvolve(expr_df, normalize_input=normalize_input)
        _save_outputs(props, "marker_db", args.output_dir)

    # --- Mode: reference (pseudobulk NNLS) ---
    if args.mode in ("reference", "both"):
        from .nnls_deconvolver import ReferenceDeconvolver

        if not args.reference:
            logger.error("--reference is required for --mode reference or both.")
            sys.exit(1)
        if not os.path.exists(args.reference):
            logger.error(f"Reference file not found: {args.reference}")
            sys.exit(1)

        logger.info(f"=== Reference-based deconvolution (ref: {args.reference}) ===")
        ref_df = pd.read_csv(args.reference, index_col=0)
        deconv_ref = ReferenceDeconvolver(
            reference=ref_df,
            cell_type_col=args.cell_type_col,
        )
        deconv_ref.build_signature()
        props_ref = deconv_ref.deconvolve(expr_df, normalize_input=normalize_input)
        _save_outputs(props_ref, "reference", args.output_dir)

    logger.info(f"Deconvolution complete. Results in: {args.output_dir}")


if __name__ == "__main__":
    main()
