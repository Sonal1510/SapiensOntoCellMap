"""
SapiensOntoCellMap — Integration Test Suite
============================================
Runs the full database build pipeline (download → parse → combine → validate)
and prints a validation summary with row counts and quality metrics.

Usage
-----
  # Full pipeline (download + build + validate)
  python3 test/test_classes.py

  # Skip download (re-use cached raw files)
  python3 test/test_classes.py --skip_download

  # Build only (no download)
  python3 test/test_classes.py --build_only
"""
import argparse
import logging
import os
import sys

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    from src.download.bio_database_downloader import BioDataDownloader
    from src.db_manager.database_creator import DatabaseCreate
    from config.config import (
        PROCESSED_DATA_DIR,
        RECOVER_ID_DATA_DIR,
        PROCESSED_COMBINED_DATA_DIR,
        DATABASE_CONFIG,
    )
except ImportError as e:
    logger.error(f"Import error: {e}")
    logger.error("Run from the project root: python3 test/test_classes.py")
    sys.exit(1)


def run_downloader():
    """Download all 14 source databases + HGNC reference."""
    logger.info("=== Step 1: Downloading databases ===")
    downloader = BioDataDownloader()
    downloader.download_all_databases()
    downloader.download_reference_data()
    logger.info("Download complete.")


def run_build_and_validate():
    """Parse, combine, and validate all databases. Print quality summary."""
    import pandas as pd

    logger.info("=== Step 2: Building combined database ===")
    db_creator = DatabaseCreate()
    db_creator.run()

    logger.info("=== Step 3: Validating output files ===")
    missing = []
    for db_name in DATABASE_CONFIG:
        for label, path in [
            (f"processed/{db_name}", os.path.join(PROCESSED_DATA_DIR, f"{db_name}_processed.csv")),
            (f"recovery/{db_name}", os.path.join(RECOVER_ID_DATA_DIR, f"{db_name}_recovery_log.csv")),
        ]:
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                missing.append(label)

    combined_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
    if not os.path.exists(combined_path) or os.path.getsize(combined_path) == 0:
        missing.append("combined_db/master_cell_marker_db.csv")

    if missing:
        logger.error(f"FAILED — {len(missing)} expected files are missing or empty:")
        for f in missing:
            logger.error(f"  - {f}")
        sys.exit(1)

    # Summary stats
    df = pd.read_csv(combined_path, low_memory=False)
    quarantine_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "quarantine_log.csv")
    n_quarantine = len(pd.read_csv(quarantine_path)) if os.path.exists(quarantine_path) else "n/a"

    print()
    print("=" * 60)
    print("  SapiensOntoCellMap — Database Build Summary")
    print("=" * 60)
    print(f"  Total rows         : {len(df):,}")
    print(f"  Unique genes       : {df['gene'].nunique():,}")
    print(f"  Unique cell types  : {df['cell_name'].nunique():,}")
    print(f"  Unique tissues     : {df['tissue_name'].nunique():,}")
    print(f"  Source databases   : {df['database_name'].nunique()}")
    print(f"  CL IDs assigned    : {df['cell_id'].str.startswith('CL:').sum():,}")
    print(f"  UBERON IDs assigned: {df['tissue_id'].str.startswith('UBERON:').sum():,}")
    print(f"  Quarantined rows   : {n_quarantine}")
    print("=" * 60)
    print()
    print("  Per-database row counts:")
    for db, n in df["database_name"].value_counts().items():
        print(f"    {db:<45s} {n:>7,}")
    print()
    print("  RESULT: SUCCESS — all expected files present.")
    print("=" * 60)


def parse_args():
    p = argparse.ArgumentParser(description="SapiensOntoCellMap integration tests")
    p.add_argument("--skip_download", action="store_true",
                   help="Skip download step (use cached raw files)")
    p.add_argument("--build_only", action="store_true",
                   help="Run build + validate only (same as --skip_download)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not args.skip_download and not args.build_only:
        run_downloader()
    run_build_and_validate()
