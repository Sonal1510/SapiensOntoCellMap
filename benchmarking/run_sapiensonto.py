#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Annotation Runner
=====================================================
Runs SapiensOntoCellMap on the 3 benchmark datasets downloaded by
download_datasets.py. Handles DEG computation for datasets that need it
(Xenium) and passes the correct input format for each platform.

Datasets
--------
1. PBMC3k          — scRNA-seq  | DEG: SpaceRanger graphclust CSV (in analysis.tar)
2. Xenium Skin     — Xenium     | DEG: computed via scanpy Wilcoxon from cell_feature_matrix.h5
3. Visium Melanoma — Visium     | DEG: SpaceRanger graphclust CSV (in analysis.tar)

Prerequisites
-------------
    python benchmarking/download_datasets.py   # download raw data first

Usage
-----
    python benchmarking/run_sapiensonto.py
    python benchmarking/run_sapiensonto.py --datasets pbmc3k xenium_skin
    python benchmarking/run_sapiensonto.py --skip_existing
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

_BENCH_DIR   = Path(__file__).parent.resolve()
_PROJECT_DIR = _BENCH_DIR.parent
_DATA_DIR    = _BENCH_DIR / "data"
_RESULTS_DIR = _BENCH_DIR / "results"
_ANNOTATE    = _PROJECT_DIR / "src" / "cluster_annotation" / "get_cluster_annotation.py"
_MARKER_DB   = _PROJECT_DIR / "data" / "processed_combined_db" / "master_cell_marker_db.csv"


# ---------------------------------------------------------------------------
# DEG computation for Xenium (scanpy Wilcoxon)
# ---------------------------------------------------------------------------

def compute_xenium_degs(data_dir: Path, output_dir: Path, skip_existing: bool) -> Path:
    """
    Compute Wilcoxon rank-sum DEGs per cluster from Xenium cell_feature_matrix.h5.
    Clusters are read from analysis.zarr.zip (Leiden graphclust).
    Returns path to the saved DEG CSV (scanpy format).
    """
    deg_path = output_dir / "xenium_skin_degs.csv"
    if deg_path.exists() and skip_existing:
        logger.info(f"Using cached DEG file: {deg_path}")
        return deg_path

    try:
        import scanpy as sc
        import pandas as pd
        import numpy as np
        import zarr
    except ImportError as e:
        logger.error(f"Missing dependency: {e}. Install: pip install scanpy zarr")
        sys.exit(1)

    h5_path   = data_dir / "cell_feature_matrix.h5"
    zarr_path = data_dir / "analysis.zarr"   # extracted from analysis.zarr.zip

    if not h5_path.exists():
        logger.error(f"cell_feature_matrix.h5 not found at {h5_path}. Run download_datasets.py first.")
        sys.exit(1)

    logger.info("Loading Xenium cell_feature_matrix.h5 ...")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    # Load cluster assignments from zarr
    cluster_labels = None
    if zarr_path.exists():
        logger.info("Loading cluster assignments from analysis.zarr ...")
        try:
            z = zarr.open(str(zarr_path), mode="r")
            # Xenium zarr: clustering/gene_expression_graphclust/clusters
            clusters_arr = z["clustering"]["gene_expression_graphclust"]["clusters"][:]
            barcodes_arr = z["clustering"]["gene_expression_graphclust"]["barcodes"][:]
            barcodes = [b.decode() if isinstance(b, bytes) else b for b in barcodes_arr]
            cluster_labels = dict(zip(barcodes, clusters_arr.astype(str)))
        except Exception as e:
            logger.warning(f"Could not read zarr clusters: {e} — falling back to leiden clustering")

    if cluster_labels:
        # Align cluster labels to adata order
        adata.obs["cluster"] = [cluster_labels.get(bc, "NA") for bc in adata.obs_names]
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")
        adata = adata[adata.obs["cluster"] != "NA"].copy()
        groupby_col = "cluster"
    else:
        logger.info("Running leiden clustering (no pre-computed clusters found) ...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.pp.pca(adata, n_comps=30)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=0.5)
        groupby_col = "leiden"

    logger.info("Normalising counts ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    logger.info(f"Computing Wilcoxon DEGs per cluster ({adata.obs[groupby_col].nunique()} clusters) ...")
    sc.tl.rank_genes_groups(adata, groupby=groupby_col, method="wilcoxon", use_raw=False)

    # Export in scanpy format (names, scores, pvals_adj, logfoldchanges per cluster)
    logger.info("Exporting DEG CSV ...")
    results = adata.uns["rank_genes_groups"]
    groups  = results["names"].dtype.names
    rows = []
    for grp in groups:
        for gene, score, pval, lfc in zip(
            results["names"][grp],
            results["scores"][grp],
            results["pvals_adj"][grp],
            results["logfoldchanges"][grp],
        ):
            rows.append({
                "names":          gene,
                "scores":         score,
                "pvals_adj":      pval,
                "logfoldchanges": lfc,
                "group":          grp,
            })

    deg_df = pd.DataFrame(rows)
    deg_df.to_csv(deg_path, index=False)
    logger.info(f"DEGs saved: {deg_path}  ({len(deg_df):,} rows)")
    return deg_path


# ---------------------------------------------------------------------------
# SpaceRanger / CellRanger graphclust DEG path
# ---------------------------------------------------------------------------

def find_graphclust_deg(data_dir: Path, dataset_key: str) -> Path:
    """
    Locate the graphclust differential_expression.csv inside an extracted
    analysis directory. Works for both CellRanger (PBMC3k) and SpaceRanger (Visium).
    """
    for root, dirs, files in os.walk(data_dir):
        if "differential_expression.csv" in files and "_graphclust" in root:
            return Path(root) / "differential_expression.csv"
    logger.error(
        f"Could not find graphclust/differential_expression.csv under {data_dir}.\n"
        f"Make sure you ran download_datasets.py (which extracts the analysis tar)."
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Annotation runner
# ---------------------------------------------------------------------------

def run_annotation(
    sample_name: str,
    input_path: Path,
    deg_type: str,
    output_dir: Path,
    tissue: str = None,
    deg_format: str = None,
    background_n: int = None,
    skip_existing: bool = False,
) -> Path:
    """
    Call get_cluster_annotation.py CLI and return the output directory.
    """
    sample_out = output_dir / sample_name
    summary_csv = sample_out / f"{sample_name}_top_annotation_summary.csv"

    if summary_csv.exists() and skip_existing:
        logger.info(f"Skipping {sample_name} (results exist)")
        return sample_out

    sample_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(_ANNOTATE),
        str(input_path),
        sample_name,
        str(sample_out),
        f"--deg_type={deg_type}",
        f"--marker_db={_MARKER_DB}",
        "--tissue_priority_ratio=0.3",
    ]

    if tissue:
        cmd.append(f"--tissue={tissue}")
    if deg_format:
        cmd.append(f"--deg_format={deg_format}")
    if background_n:
        cmd.append(f"--background_gene_count={background_n}")

    logger.info(f"\n{'='*60}")
    logger.info(f"Running SapiensOntoCellMap: {sample_name}")
    logger.info(f"Input : {input_path}")
    logger.info(f"Output: {sample_out}")
    logger.info(f"{'='*60}")

    result = subprocess.run(cmd, cwd=str(_PROJECT_DIR))
    if result.returncode != 0:
        logger.error(f"Annotation failed for {sample_name} (exit code {result.returncode})")
        sys.exit(result.returncode)

    logger.info(f"Annotation complete: {sample_name}")
    return sample_out


# ---------------------------------------------------------------------------
# Per-dataset logic
# ---------------------------------------------------------------------------

def run_pbmc3k(skip_existing: bool) -> None:
    data_dir = _DATA_DIR / "pbmc3k"
    if not data_dir.exists():
        logger.error("PBMC3k data not found. Run: python benchmarking/download_datasets.py --datasets pbmc3k")
        sys.exit(1)

    deg_csv = find_graphclust_deg(data_dir, "pbmc3k")
    logger.info(f"PBMC3k DEG file: {deg_csv}")

    run_annotation(
        sample_name="pbmc3k",
        input_path=deg_csv,
        deg_type="scrna",
        output_dir=_RESULTS_DIR,
        tissue=None,                  # pan-tissue — PBMCs are blood, use all-tissue
        deg_format="generic",         # CellRanger 1.x graphclust format
        background_n=20000,           # full transcriptome background
        skip_existing=skip_existing,
    )


def run_xenium_skin(skip_existing: bool) -> None:
    data_dir = _DATA_DIR / "xenium_skin"
    if not data_dir.exists():
        logger.error("Xenium skin data not found. Run: python benchmarking/download_datasets.py --datasets xenium_skin")
        sys.exit(1)

    deg_csv = compute_xenium_degs(data_dir, _RESULTS_DIR / "xenium_skin", skip_existing)

    run_annotation(
        sample_name="xenium_skin",
        input_path=deg_csv,
        deg_type="spatial",
        output_dir=_RESULTS_DIR,
        tissue="skin",
        deg_format="scanpy",
        skip_existing=skip_existing,
    )


def run_visium_melanoma(skip_existing: bool) -> None:
    data_dir = _DATA_DIR / "visium_melanoma"
    if not data_dir.exists():
        logger.error("Visium melanoma data not found. Run: python benchmarking/download_datasets.py --datasets visium_melanoma")
        sys.exit(1)

    deg_csv = find_graphclust_deg(data_dir, "visium_melanoma")
    logger.info(f"Visium Melanoma DEG file: {deg_csv}")

    run_annotation(
        sample_name="visium_melanoma",
        input_path=deg_csv,
        deg_type="spatial",
        output_dir=_RESULTS_DIR,
        tissue="skin",
        deg_format="generic",
        skip_existing=skip_existing,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

RUNNERS = {
    "pbmc3k":          run_pbmc3k,
    "xenium_skin":     run_xenium_skin,
    "visium_melanoma": run_visium_melanoma,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run SapiensOntoCellMap annotation on all benchmark datasets."
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=list(RUNNERS.keys()),
        default=list(RUNNERS.keys()),
        help="Which datasets to annotate (default: all 3)",
    )
    parser.add_argument(
        "--skip_existing",
        action="store_true",
        help="Skip datasets where results already exist",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info(f"Datasets : {', '.join(args.datasets)}")
    logger.info(f"Results  : {_RESULTS_DIR}")

    for key in args.datasets:
        RUNNERS[key](args.skip_existing)

    logger.info("\nAll annotations complete.")
    logger.info(f"Results: {_RESULTS_DIR}")


if __name__ == "__main__":
    main()
