#!/usr/bin/env python3
"""
SapiensOntoCellMap — Combined Benchmark Runner
===============================================
Runs all configured datasets and comparator tools, then generates
the combined publication figure.

Usage
-----
  # Run all datasets
  python benchmarking/run_all_benchmarks.py

  # Run specific datasets
  python benchmarking/run_all_benchmarks.py --datasets pbmc3k tabula_sapiens

  # Check tool availability without running benchmarks
  python benchmarking/run_all_benchmarks.py --check_tools

  # Skip download (use cached data)
  python benchmarking/run_all_benchmarks.py --skip_download
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

from benchmarking.download.tools.tool_installer import ToolInstaller
from benchmarking.figures.benchmark_figures import BenchmarkFigures

AVAILABLE_DATASETS = ["pbmc3k", "tabula_sapiens"]
RESULTS_BASE = os.path.join(_BENCH_DIR, "results")


def check_tools() -> None:
    installer = ToolInstaller(["celltypist", "scsa", "sctype", "singler"])
    availability = installer.check_all(raise_on_missing=False)
    print("\n=== Tool Availability ===")
    for tool, available in availability.items():
        status = "✅  available" if available else "❌  NOT found"
        print(f"  {tool:<20s} {status}")
    print()


def run_pbmc3k(skip_download: bool = False) -> pd.DataFrame:
    """Run PBMC3k benchmark and return summary DataFrame."""
    from benchmarking.benchmark_pbmc3k import run_pbmc3k_benchmark
    return run_pbmc3k_benchmark(
        output_dir=os.path.join(RESULTS_BASE, "pbmc3k"),
        skip_download=skip_download,
    )


def run_tabula_sapiens(skip_download: bool = False) -> pd.DataFrame:
    """Run Tabula Sapiens benchmark and return summary DataFrame."""
    from benchmarking.benchmark_tabula_sapiens import run_tabula_sapiens_benchmark
    return run_tabula_sapiens_benchmark(
        output_dir=os.path.join(RESULTS_BASE, "tabula_sapiens"),
        skip_download=skip_download,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="SapiensOntoCellMap — combined benchmark runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--datasets", nargs="+", choices=AVAILABLE_DATASETS,
        default=AVAILABLE_DATASETS,
        help="Datasets to benchmark (default: all)",
    )
    parser.add_argument(
        "--skip_download", action="store_true",
        help="Use cached data files (skip download step)",
    )
    parser.add_argument(
        "--check_tools", action="store_true",
        help="Check tool availability and exit",
    )
    args = parser.parse_args()

    if args.check_tools:
        check_tools()
        return

    dataset_summaries: dict[str, pd.DataFrame] = {}

    runners = {
        "pbmc3k":          run_pbmc3k,
        "tabula_sapiens":  run_tabula_sapiens,
    }

    for dataset in args.datasets:
        logger.info(f"=== Running benchmark: {dataset} ===")
        try:
            summary = runners[dataset](skip_download=args.skip_download)
            dataset_summaries[dataset] = summary
            print(f"\n--- {dataset} summary ---")
            print(summary.to_string(index=False))
        except Exception as exc:
            logger.error(f"Benchmark {dataset} failed: {exc}")

    if len(dataset_summaries) >= 2:
        logger.info("Generating combined publication figure...")
        fig_dir = os.path.join(RESULTS_BASE, "figures")
        figs = BenchmarkFigures(output_dir=fig_dir)
        path = figs.plot_combined_figure(
            dataset_summaries={
                k.replace("_", " ").title(): v
                for k, v in dataset_summaries.items()
            }
        )
        logger.info(f"Combined figure saved: {path}")


if __name__ == "__main__":
    main()
