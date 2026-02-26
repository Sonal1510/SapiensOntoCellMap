"""
SapiensOntoCellMap Runner
==========================
Wraps the get_cluster_annotation.py CLI as an annotation tool runner,
producing predictions in the standard {cluster: cell_type} dict format.
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys

import pandas as pd

from benchmarking.annotation_tool_exec.base_runner import BaseAnnotationRunner

logger = logging.getLogger(__name__)

_PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_ANNOTATION_SCRIPT = os.path.join(
    _PROJECT_DIR, "src", "cluster_annotation", "get_cluster_annotation.py"
)


class SapiensOntoRunner(BaseAnnotationRunner):
    """
    Runs SapiensOntoCellMap annotation via the CLI entry point.

    Reads the top annotation summary CSV to produce cluster → cell type predictions.

    Additional kwargs accepted by run():
      tissue : str       — tissue filter (e.g. 'skin', 'blood')
      background_n : int — background gene count for hypergeometric test
      no_hierarchy : bool — skip CL ontology hierarchical annotation
    """

    @property
    def tool_name(self) -> str:
        return "sapiensonto"

    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict[str, str]:
        tool_dir = self.tool_output_dir(output_dir)
        sample_name = kwargs.get("sample_name", "benchmark")
        tissue = kwargs.get("tissue", None)
        background_n = kwargs.get("background_n", 20000)
        no_hierarchy = kwargs.get("no_hierarchy", False)

        cmd = [
            sys.executable, _ANNOTATION_SCRIPT,
            deg_csv, sample_name, tool_dir,
            "--deg_type", "scrna",
            "--background_gene_count", str(background_n),
        ]
        if tissue:
            cmd += ["--tissue", tissue]
        if no_hierarchy:
            cmd.append("--no_hierarchy")

        logger.info(f"[{self.tool_name}] Running: {' '.join(cmd[-6:])}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"[{self.tool_name}] FAILED:\n{result.stderr[-1000:]}")
            return {}

        # Parse top annotation summary CSV
        summary_pattern = os.path.join(tool_dir, f"{sample_name}_top_annotation_summary.csv")
        if not os.path.exists(summary_pattern):
            logger.error(f"[{self.tool_name}] Output not found: {summary_pattern}")
            return {}

        df = pd.read_csv(summary_pattern)
        predictions: dict[str, str] = {}
        for _, row in df.iterrows():
            cluster = str(row.get("Cluster", ""))
            cell_type = str(row.get("Cell_type", ""))
            if cluster and cell_type:
                predictions[cluster] = cell_type

        logger.info(f"[{self.tool_name}] Annotated {len(predictions)} clusters.")
        return predictions
