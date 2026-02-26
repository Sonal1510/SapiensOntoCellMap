#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Tabula Sapiens
==================================================
Author : Sonal Rashmi
Date   : 2026-02-25

Benchmark on skin and blood subsets of the Tabula Sapiens atlas.
Reference: The Tabula Sapiens Consortium, Science 2022 (doi:10.1126/science.abl4896)

Requires:
  pip install cellxgene-census scanpy celltypist

Usage
-----
  python benchmarking/benchmark_tabula_sapiens.py

  # Skip download (use cached h5ad)
  python benchmarking/benchmark_tabula_sapiens.py --skip_download

  # Skip R tools
  python benchmarking/benchmark_tabula_sapiens.py --no_r_tools
"""

import argparse
import logging
import os
import sys

import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

_BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(_BENCH_DIR)
sys.path.insert(0, _PROJECT_DIR)

from benchmarking.annotation_tool_exec.celltypist_runner import CellTypistRunner
from benchmarking.annotation_tool_exec.sapiensonto_runner import SapiensOntoRunner
from benchmarking.annotation_tool_exec.sctype_runner import ScTypeRunner
from benchmarking.annotation_tool_exec.singler_runner import SingleRRunner
from benchmarking.download.datasets.tabula_sapiens_downloader import TabulaSapiensDownloader
from benchmarking.figures.benchmark_figures import BenchmarkFigures
from benchmarking.ground_truth.tabula_sapiens_gt import TABULA_SAPIENS_GROUND_TRUTH
from benchmarking.metrics.benchmark_metrics import BenchmarkMetrics

try:
    from config.config import PROCESSED_COMBINED_DATABASE_FILE
    MARKER_DB = PROCESSED_COMBINED_DATABASE_FILE
except ImportError:
    MARKER_DB = os.path.join(_PROJECT_DIR, "data", "processed_combined_db", "master_cell_marker_db.csv")


# ---------------------------------------------------------------------------
# DEG computation (Tabula Sapiens: Wilcoxon on cell_type obs column)
# ---------------------------------------------------------------------------

def compute_degs(adata, output_dir: str, groupby: str = "cell_type") -> str:
    """
    Run Wilcoxon rank-sum test per cell type and export in Scanpy CSV format.
    Returns path to the saved DEG CSV. Idempotent — skips if already present.
    """
    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy required. pip install scanpy")
        sys.exit(1)

    deg_path = os.path.join(output_dir, "degs.csv")
    if os.path.exists(deg_path):
        logger.info(f"Using cached DEG file: {deg_path}")
        return deg_path

    logger.info(f"Computing Wilcoxon DEGs grouped by '{groupby}'...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method="wilcoxon",
        n_genes=min(200, adata.n_vars),
        use_raw=False,
    )
    records = []
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    for grp in groups:
        names     = adata.uns["rank_genes_groups"]["names"][grp]
        pvals_adj = adata.uns["rank_genes_groups"]["pvals_adj"][grp]
        logfc     = adata.uns["rank_genes_groups"]["logfoldchanges"][grp]
        for g, p, lfc in zip(names, pvals_adj, logfc):
            records.append({
                "names":          g,
                "group":          grp,
                "pvals_adj":      float(p),
                "logfoldchanges": float(lfc),
            })
    pd.DataFrame(records).to_csv(deg_path, index=False)
    logger.info(f"DEGs: {len(records):,} rows → {deg_path}")
    return deg_path


# ---------------------------------------------------------------------------
# Main benchmark function
# ---------------------------------------------------------------------------

def run_tabula_sapiens_benchmark(
    output_dir: str,
    skip_download: bool = False,
    background_n: int = 20000,
    no_celltypist: bool = False,
    no_r_tools: bool = False,
) -> pd.DataFrame:
    """
    Run the full Tabula Sapiens (skin + blood) benchmark.

    Returns
    -------
    pd.DataFrame
        Per-tool accuracy summary.
    """
    if not os.path.exists(MARKER_DB):
        logger.error(f"Marker DB not found: {MARKER_DB}")
        logger.error("Run scripts/build_marker_db.py to build the database first.")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # 1. Download / load Tabula Sapiens subset
    downloader = TabulaSapiensDownloader(output_dir=output_dir)
    h5ad_path = str(downloader.download(force=False))

    try:
        import scanpy as sc
        adata = sc.read_h5ad(h5ad_path)
    except ImportError:
        logger.error("scanpy required. pip install scanpy")
        sys.exit(1)

    deg_csv = compute_degs(adata, output_dir, groupby="cell_type")

    # 2. Run annotation tools
    metrics = BenchmarkMetrics(ground_truth=TABULA_SAPIENS_GROUND_TRUTH)

    # SapiensOntoCellMap
    sapiensonto = SapiensOntoRunner()
    predictions = sapiensonto.run(
        deg_csv, output_dir,
        sample_name="tabula_sapiens_benchmark",
        background_n=background_n,
    )
    metrics.add_predictions("SapiensOntoCellMap", predictions)

    # CellTypist
    if not no_celltypist:
        celltypist = CellTypistRunner()
        ct_preds = celltypist.run(
            deg_csv, output_dir,
            adata_path=h5ad_path,
            model_name="Healthy_COVID19_PBMC",  # broad immune model
            cluster_key="cell_type",
        )
        metrics.add_predictions("CellTypist", ct_preds)

    # scType (R)
    if not no_r_tools:
        sctype = ScTypeRunner()
        st_preds = sctype.run(deg_csv, output_dir, tissue_type="Skin")
        if st_preds:
            metrics.add_predictions("scType", st_preds)

    # SingleR (R)
    if not no_r_tools:
        singler = SingleRRunner()
        sr_preds = singler.run(
            deg_csv, output_dir,
            adata_path=h5ad_path,
            reference="HumanPrimaryCellAtlasData",
            cluster_key="cell_type",
        )
        if sr_preds:
            metrics.add_predictions("SingleR", sr_preds)

    # 3. Save metrics
    summary = metrics.summary()
    detail_df = metrics.to_dataframe()

    metrics_dir = os.path.join(output_dir, "metrics")
    os.makedirs(metrics_dir, exist_ok=True)
    summary.to_csv(os.path.join(metrics_dir, "accuracy_summary.csv"), index=False)
    detail_df.to_csv(os.path.join(metrics_dir, "per_cluster_results.csv"), index=False)
    logger.info(f"Metrics saved to {metrics_dir}")

    # 4. Generate figures
    figs = BenchmarkFigures(output_dir=os.path.join(output_dir, "figures"))
    figs.plot_accuracy_bar(summary, dataset_name="Tabula Sapiens", metric="top1_accuracy")

    print("\n=== Tabula Sapiens Benchmark Summary ===")
    print(summary.to_string(index=False))
    return summary


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="SapiensOntoCellMap — Tabula Sapiens benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--skip_download", action="store_true",
                   help="Use cached h5ad if present")
    p.add_argument("--background_n", type=int, default=20000)
    p.add_argument("--no_celltypist", action="store_true")
    p.add_argument("--no_r_tools", action="store_true",
                   help="Skip scType and SingleR")
    p.add_argument("--output_dir",
                   default=os.path.join(_BENCH_DIR, "results", "tabula_sapiens"),
                   help="Output directory")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_tabula_sapiens_benchmark(
        output_dir=args.output_dir,
        skip_download=args.skip_download,
        background_n=args.background_n,
        no_celltypist=args.no_celltypist,
        no_r_tools=args.no_r_tools,
    )
