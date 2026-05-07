#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Dataset Downloader
=====================================================
Downloads processed (CellRanger/SpaceRanger/Xenium) output files for the
3 public benchmark datasets. All files are freely accessible without login.

Datasets
--------
1. PBMC3k          — 3k PBMCs from a Healthy Donor (10x Chromium, CellRanger 1.1.0)
2. Xenium Skin     — Human Skin, non-diseased, section 1 (Xenium 1.9.0)
3. Visium Melanoma — Human Melanoma IF Stained FFPE (CytAssist SpaceRanger 2.0.0)

Usage
-----
    python benchmarking/download_datasets.py
    python benchmarking/download_datasets.py --datasets pbmc3k xenium_skin
    python benchmarking/download_datasets.py --output_dir /path/to/data
    python benchmarking/download_datasets.py --skip_existing      # skip already-downloaded files
"""

import argparse
import logging
import os
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

_DEFAULT_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# ---------------------------------------------------------------------------
# Dataset definitions — all CF URLs verified 200 on 2026-05-07
# ---------------------------------------------------------------------------

DATASETS = {
    "pbmc3k": {
        "name": "3k PBMCs from a Healthy Donor",
        "platform": "scRNA-seq (CellRanger 1.1.0)",
        "tissue": "Blood (normal)",
        "source": "https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0",
        "files": [
            {
                "url": "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                "filename": "pbmc3k_filtered_gene_bc_matrices.tar.gz",
                "extract": True,
                "description": "Filtered gene-barcode matrix (MTX format)",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_analysis.tar.gz",
                "filename": "pbmc3k_analysis.tar.gz",
                "extract": True,
                "description": "Cluster assignments + differential expression",
            },
        ],
    },

    "xenium_skin": {
        "name": "Human Skin, Non-Diseased, Section 1 (Xenium Multi-Tissue Panel)",
        "platform": "Xenium (XOA 1.9.0)",
        "tissue": "Skin (normal FFPE)",
        "source": "https://www.10xgenomics.com/datasets/human-skin-data-xenium-human-multi-tissue-and-cancer-panel-1-standard",
        "files": [
            {
                "url": "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hSkin_nondiseased_section_1_FFPE/Xenium_V1_hSkin_nondiseased_section_1_FFPE_cell_feature_matrix.h5",
                "filename": "cell_feature_matrix.h5",
                "extract": False,
                "description": "Cell x gene expression matrix (HDF5)",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hSkin_nondiseased_section_1_FFPE/Xenium_V1_hSkin_nondiseased_section_1_FFPE_cells.csv.gz",
                "filename": "cells.csv.gz",
                "extract": False,
                "description": "Cell coordinates and metadata",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hSkin_nondiseased_section_1_FFPE/Xenium_V1_hSkin_nondiseased_section_1_FFPE_analysis.zarr.zip",
                "filename": "analysis.zarr.zip",
                "extract": True,
                "description": "Cluster assignments (Leiden/graph-based)",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/xenium/1.9.0/Xenium_V1_hSkin_nondiseased_section_1_FFPE/Xenium_V1_hSkin_nondiseased_section_1_FFPE_experiment.xenium",
                "filename": "experiment.xenium",
                "extract": False,
                "description": "Experiment metadata (JSON)",
            },
        ],
    },

    "visium_melanoma": {
        "name": "Human Melanoma, IF Stained (FFPE)",
        "platform": "Visium CytAssist (SpaceRanger 2.0.0)",
        "tissue": "Skin melanoma (FFPE)",
        "source": "https://www.10xgenomics.com/datasets/human-melanoma-if-stained-ffpe-2-standard",
        "files": [
            {
                "url": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5",
                "filename": "filtered_feature_bc_matrix.h5",
                "extract": False,
                "description": "Filtered spot x gene matrix (HDF5)",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz",
                "filename": "spatial.tar.gz",
                "extract": True,
                "description": "Spatial coordinates + tissue image",
            },
            {
                "url": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_analysis.tar.gz",
                "filename": "analysis.tar.gz",
                "extract": True,
                "description": "Cluster assignments + differential expression",
            },
        ],
    },
}


# ---------------------------------------------------------------------------
# Downloader
# ---------------------------------------------------------------------------

_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"
}


def _progress_hook(block_num: int, block_size: int, total_size: int) -> None:
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(100, downloaded * 100 // total_size)
        mb = downloaded / 1_048_576
        total_mb = total_size / 1_048_576
        print(f"\r    {mb:.1f} / {total_mb:.1f} MB ({pct}%)", end="", flush=True)


def _download_file(url: str, dest: Path, skip_existing: bool) -> bool:
    if dest.exists() and skip_existing:
        logger.info(f"  Skipping (exists): {dest.name}")
        return False
    logger.info(f"  Downloading: {dest.name}")
    try:
        req = urllib.request.Request(url, headers=_HEADERS)
        with urllib.request.urlopen(req) as response, open(dest, "wb") as out:
            total = int(response.headers.get("Content-Length", 0))
            downloaded = 0
            block = 1 << 16  # 64 KB
            while True:
                chunk = response.read(block)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = min(100, downloaded * 100 // total)
                    print(f"\r    {downloaded/1e6:.1f} / {total/1e6:.1f} MB ({pct}%)", end="", flush=True)
        print()
        return True
    except Exception as e:
        print()
        logger.error(f"  Failed: {e}")
        if dest.exists():
            dest.unlink()
        return False


def _extract(archive: Path, dest_dir: Path) -> None:
    logger.info(f"  Extracting: {archive.name}")
    if archive.suffix == ".gz" and archive.stem.endswith(".tar"):
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(dest_dir)
    elif archive.suffix == ".zip":
        with zipfile.ZipFile(archive, "r") as zf:
            zf.extractall(dest_dir)
    else:
        logger.warning(f"  Unknown archive format: {archive.name} — skipping extract")


def download_dataset(key: str, dataset: dict, output_dir: Path, skip_existing: bool) -> None:
    dataset_dir = output_dir / key
    dataset_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"\n{'='*60}")
    logger.info(f"Dataset : {dataset['name']}")
    logger.info(f"Platform: {dataset['platform']}")
    logger.info(f"Tissue  : {dataset['tissue']}")
    logger.info(f"Output  : {dataset_dir}")
    logger.info(f"{'='*60}")

    for file_spec in dataset["files"]:
        dest = dataset_dir / file_spec["filename"]
        logger.info(f"\n  [{file_spec['description']}]")
        downloaded = _download_file(file_spec["url"], dest, skip_existing)
        if (downloaded or dest.exists()) and file_spec.get("extract"):
            _extract(dest, dataset_dir)

    logger.info(f"\n  Done: {dataset['name']}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download processed benchmark datasets for SapiensOntoCellMap."
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=list(DATASETS.keys()),
        default=list(DATASETS.keys()),
        help="Which datasets to download (default: all 3)",
    )
    parser.add_argument(
        "--output_dir",
        default=_DEFAULT_OUTPUT_DIR,
        help=f"Root output directory (default: {_DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--skip_existing",
        action="store_true",
        help="Skip files that already exist in the output directory",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory : {output_dir.resolve()}")
    logger.info(f"Datasets         : {', '.join(args.datasets)}")
    logger.info(f"Skip existing    : {args.skip_existing}")

    for key in args.datasets:
        download_dataset(key, DATASETS[key], output_dir, args.skip_existing)

    logger.info("\nAll downloads complete.")
    logger.info(f"Data directory: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
