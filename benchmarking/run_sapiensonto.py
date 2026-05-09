#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Annotation Runner
=====================================================
Runs SapiensOntoCellMap on the benchmark datasets.

Datasets
--------
1. PBMC3k              — scRNA-seq  | DEG: CellRanger kmeans/8_clusters CSV (no graphclust in v1.1.0)
2. Atera Breast Cancer — Atera WTA  | DEG: pre-computed graphclust CSV from NAS analysis.tar.gz
                                       GT:  cell_groups.csv from NAS (per-barcode, 20 cell types)

Prerequisites
-------------
    python benchmarking/download_datasets.py --datasets pbmc3k   # PBMC3k only
    # Atera reads directly from NAS — /Volumes/shainlab/Sonal/sapiensontocellmap_atera/

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
# Atera NAS helpers — extract pre-computed files from analysis.tar.gz
# ---------------------------------------------------------------------------

_ATERA_NAS = Path("/Volumes/shainlab/Sonal/sapiensontocellmap_atera")

_ATERA_DEG_TAR_PATH      = "analysis/diffexp/gene_expression_graphclust/differential_expression.csv"
_ATERA_CLUSTERS_TAR_PATH = "analysis/clustering/gene_expression_graphclust/clusters.csv"
_ATERA_CELL_GROUPS_URL   = (
    "https://cf.10xgenomics.com/samples/atera/dev/"
    "WTA_Preview_FFPE_Breast_Cancer/WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv"
)
_HEADERS = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"}


def extract_atera_nas_files(out_dir: Path, skip_existing: bool) -> tuple[Path, Path, Path]:
    """
    Prepare the three files needed for Atera benchmarking:
      1. graphclust differential_expression.csv — extracted from NAS analysis.tar.gz
      2. graphclust clusters.csv               — extracted from NAS analysis.tar.gz
      3. cell_groups.csv (GT labels)           — downloaded from 10x CF if not cached

    All outputs are written to out_dir (local results directory).
    Returns (deg_csv, clusters_csv, cell_groups_csv) as local Paths.
    """
    import tarfile
    import urllib.request

    nas_tar       = _ATERA_NAS / "analysis.tar.gz"
    deg_dest      = out_dir / "atera_differential_expression.csv"
    clusters_dest = out_dir / "atera_clusters.csv"
    cg_dest       = out_dir / "atera_cell_groups.csv"

    if not nas_tar.exists():
        logger.error(f"NAS file not found: {nas_tar}. Mount /Volumes/shainlab first.")
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    if deg_dest.exists() and skip_existing:
        logger.info(f"Using cached DEG CSV: {deg_dest}")
    else:
        logger.info(f"Extracting DEG CSV from {nas_tar.name} ...")
        with tarfile.open(nas_tar, "r:gz") as tf:
            member = tf.getmember(_ATERA_DEG_TAR_PATH)
            with tf.extractfile(member) as src, open(deg_dest, "wb") as dst:
                dst.write(src.read())
        logger.info(f"  Saved: {deg_dest}")

    if clusters_dest.exists() and skip_existing:
        logger.info(f"Using cached clusters CSV: {clusters_dest}")
    else:
        logger.info(f"Extracting clusters CSV from {nas_tar.name} ...")
        with tarfile.open(nas_tar, "r:gz") as tf:
            member = tf.getmember(_ATERA_CLUSTERS_TAR_PATH)
            with tf.extractfile(member) as src, open(clusters_dest, "wb") as dst:
                dst.write(src.read())
        logger.info(f"  Saved: {clusters_dest}")

    if cg_dest.exists() and skip_existing:
        logger.info(f"Using cached cell_groups CSV: {cg_dest}")
    else:
        logger.info(f"Downloading cell_groups.csv from 10x CF ...")
        req = urllib.request.Request(_ATERA_CELL_GROUPS_URL, headers=_HEADERS)
        with urllib.request.urlopen(req) as resp, open(cg_dest, "wb") as dst:
            dst.write(resp.read())
        logger.info(f"  Saved: {cg_dest}")

    return deg_dest, clusters_dest, cg_dest


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
    out_dir = _RESULTS_DIR / "atera_breast_cancer"
    deg_csv, clusters_csv, cell_groups_csv = extract_atera_nas_files(out_dir, skip_existing)

    run_annotation(
        sample_name="atera_breast_cancer",
        input_path=deg_csv,
        deg_type="spatial",
        output_dir=_RESULTS_DIR,
        tissue="breast",
        deg_format="generic",
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
