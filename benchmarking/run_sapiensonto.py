#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Annotation Runner
=====================================================
Runs SapiensOntoCellMap on the benchmark datasets downloaded by
download_datasets.py.

Datasets
--------
1. PBMC3k              — scRNA-seq  | DEG: CellRanger kmeans/8_clusters CSV (no graphclust in v1.1.0)
2. Atera Breast Cancer — Atera WTA  | DEG: computed via scanpy Wilcoxon from cell_feature_matrix.h5

Prerequisites
-------------
    python benchmarking/download_datasets.py   # download raw data first

Usage
-----
    python benchmarking/run_sapiensonto.py
    python benchmarking/run_sapiensonto.py --datasets pbmc3k
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
# CellRanger kmeans DEG path (PBMC3k)
# ---------------------------------------------------------------------------

def find_kmeans_deg(data_dir: Path, n_clusters: int = 8) -> Path:
    """
    Locate the kmeans differential_expression.csv for a specific cluster count.
    Used for PBMC3k (CellRanger 1.1.0), which has no graphclust output.
    Prefers n_clusters; falls back to the closest available count.
    """
    target = f"{n_clusters}_clusters"
    for root, dirs, files in os.walk(data_dir):
        if "differential_expression.csv" in files and target in root:
            return Path(root) / "differential_expression.csv"

    # Fallback: any kmeans differential_expression.csv
    candidates = []
    for root, dirs, files in os.walk(data_dir):
        if "differential_expression.csv" in files and "kmeans" in root:
            candidates.append(Path(root) / "differential_expression.csv")
    if candidates:
        candidates.sort()
        chosen = candidates[len(candidates) // 2]
        logger.warning(f"kmeans/{n_clusters}_clusters not found; using {chosen}")
        return chosen

    logger.error(
        f"Could not find kmeans/differential_expression.csv under {data_dir}.\n"
        f"Make sure you ran download_datasets.py."
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# DEG computation for Atera (scanpy Wilcoxon from cell_feature_matrix.h5)
# ---------------------------------------------------------------------------

def compute_atera_degs(data_dir: Path, output_dir: Path, skip_existing: bool,
                       dataset_name: str = "atera_breast_cancer") -> Path:
    """
    Compute Wilcoxon rank-sum DEGs per cluster from Atera cell_feature_matrix.h5.

    Clusters are read from cell_groups.csv (per-barcode labels from 10x CF).
    Saves DEG CSV (scanpy format) and barcode→cluster CSV for GT comparison.

    Returns path to the saved DEG CSV.
    """
    deg_path = output_dir / f"{dataset_name}_degs.csv"
    if deg_path.exists() and skip_existing:
        logger.info(f"Using cached DEG file: {deg_path}")
        return deg_path

    try:
        import scanpy as sc
        import pandas as pd
    except ImportError as e:
        logger.error(f"Missing dependency: {e}. Install: pip install scanpy")
        sys.exit(1)

    h5_path        = data_dir / "cell_feature_matrix.h5"
    cell_groups_path = data_dir / "cell_groups.csv"

    for p in [h5_path, cell_groups_path]:
        if not p.exists():
            logger.error(f"{p.name} not found at {p}.")
            sys.exit(1)

    logger.info("Loading Atera cell_feature_matrix.h5 ...")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    logger.info(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    logger.info("Loading cluster assignments from cell_groups.csv ...")
    cg = pd.read_csv(cell_groups_path)
    # Expected columns: Barcode, Group (cell type label)
    barcode_col = next((c for c in cg.columns if "barcode" in c.lower()), cg.columns[0])
    group_col   = next((c for c in cg.columns if "group" in c.lower() or "cell" in c.lower()
                        and c != barcode_col), cg.columns[1])
    cg = cg.rename(columns={barcode_col: "barcode", group_col: "cluster"})
    cg = cg.set_index("barcode")["cluster"]

    adata.obs["cluster"] = cg.reindex(adata.obs_names).values
    n_missing = adata.obs["cluster"].isna().sum()
    if n_missing > 0:
        logger.warning(f"  {n_missing:,} cells not in cell_groups.csv — dropping")
        adata = adata[adata.obs["cluster"].notna()].copy()
    adata.obs["cluster"] = adata.obs["cluster"].astype("category")
    logger.info(f"  {adata.n_obs:,} cells across {adata.obs['cluster'].nunique()} clusters")

    logger.info("Normalising counts ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    logger.info(f"Computing Wilcoxon DEGs ({adata.obs['cluster'].nunique()} clusters) ...")
    sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", use_raw=False)

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

    output_dir.mkdir(parents=True, exist_ok=True)
    deg_df = pd.DataFrame(rows)
    deg_df.to_csv(deg_path, index=False)
    logger.info(f"DEGs saved: {deg_path}  ({len(deg_df):,} rows)")

    cell_assignments_path = output_dir / f"{dataset_name}_cell_clusters.csv"
    pd.DataFrame({
        "barcode": adata.obs_names,
        "cluster": adata.obs["cluster"].astype(str),
    }).to_csv(cell_assignments_path, index=False)
    logger.info(f"Cell cluster assignments saved: {cell_assignments_path}")

    return deg_path


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
    """Call get_cluster_annotation.py CLI and return the output directory."""
    sample_out  = output_dir / sample_name
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

def convert_cellranger_v1_to_standard(src_csv: Path, out_csv: Path) -> Path:
    """
    Convert CellRanger 1.x wide-format DEG CSV to standard Visium-style format.

    CellRanger 1.x columns: Gene ID, Gene Name, Cluster N Weight, Cluster N UMI counts/cell
    Output columns:         Feature ID, Feature Name, Cluster N Log2 fold change,
                            Cluster N Adjusted p value, Cluster N Mean Counts

    No adjusted p-values exist in CellRanger 1.x — p=0.0001 for positive Weight, p=1.0 otherwise.
    """
    if out_csv.exists():
        logger.info(f"Using cached converted DEG file: {out_csv}")
        return out_csv

    import pandas as pd
    import re

    df = pd.read_csv(src_csv)
    df.rename(columns={"Gene ID": "Feature ID", "Gene Name": "Feature Name"}, inplace=True)

    weight_cols  = [c for c in df.columns if re.match(r"Cluster \d+ Weight", c)]
    cluster_nums = [re.search(r"Cluster (\d+) Weight", c).group(1) for c in weight_cols]

    result = df[["Feature ID", "Feature Name"]].copy()
    for n in cluster_nums:
        w = df[f"Cluster {n} Weight"]
        umi_col = f"Cluster {n} UMI counts/cell"
        result[f"Cluster {n} Log2 fold change"] = w
        result[f"Cluster {n} Adjusted p value"] = w.apply(lambda x: 0.0001 if x > 0 else 1.0)
        result[f"Cluster {n} Mean Counts"]      = df[umi_col] if umi_col in df.columns else 0.0

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_csv, index=False)
    logger.info(f"Converted CellRanger v1 DEG to standard format: {out_csv}  ({len(result):,} genes)")
    return out_csv


def run_pbmc3k(skip_existing: bool) -> None:
    data_dir = _DATA_DIR / "pbmc3k"
    if not data_dir.exists():
        logger.error("PBMC3k data not found. Run: python benchmarking/download_datasets.py --datasets pbmc3k")
        sys.exit(1)

    raw_deg_csv   = find_kmeans_deg(data_dir, n_clusters=8)
    logger.info(f"PBMC3k DEG file: {raw_deg_csv}")

    converted_csv = _RESULTS_DIR / "pbmc3k" / "pbmc3k_deg_converted.csv"
    deg_csv       = convert_cellranger_v1_to_standard(raw_deg_csv, converted_csv)

    run_annotation(
        sample_name="pbmc3k",
        input_path=deg_csv,
        deg_type="spatial",
        output_dir=_RESULTS_DIR,
        tissue=None,
        deg_format=None,
        background_n=20000,
        skip_existing=skip_existing,
    )


def run_atera_breast_cancer(skip_existing: bool) -> None:
    data_dir = _DATA_DIR / "atera_breast_cancer"
    if not data_dir.exists():
        logger.error(
            "Atera breast cancer data not found. "
            "Extract WTA_Preview_FFPE_Breast_Cancer_outs.zip to benchmarking/data/atera_breast_cancer/"
        )
        sys.exit(1)

    out_dir = _RESULTS_DIR / "atera_breast_cancer"
    deg_csv = compute_atera_degs(data_dir, out_dir, skip_existing, dataset_name="atera_breast_cancer")

    run_annotation(
        sample_name="atera_breast_cancer",
        input_path=deg_csv,
        deg_type="scrna",
        output_dir=_RESULTS_DIR,
        tissue="breast",
        deg_format="scanpy",
        skip_existing=skip_existing,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

RUNNERS = {
    "pbmc3k":              run_pbmc3k,
    "atera_breast_cancer": run_atera_breast_cancer,
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
        help="Which datasets to annotate (default: all)",
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
