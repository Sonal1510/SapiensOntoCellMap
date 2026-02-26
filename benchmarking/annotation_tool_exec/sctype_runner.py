"""
scType Runner
=============
Wraps scType (Ianevski et al. Nat Commun 2022) via an Rscript subprocess.

scType has no CRAN or PyPI package. It is loaded directly from GitHub
in the R wrapper script (benchmarking/r_scripts/run_sctype.R).

The Python–R boundary is a CSV file:
  Input:  DEG CSV (Seurat FindAllMarkers format)
  Output: <tool_dir>/sctype_predictions.csv  (columns: cluster, predicted_cell_type)

Requires: Rscript in PATH, R packages: HGNChelper, openxlsx
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess

import pandas as pd

from benchmarking.annotation_tool_exec.base_runner import BaseAnnotationRunner

logger = logging.getLogger(__name__)

_RSCRIPT = os.path.join(
    os.path.dirname(__file__), "..", "r_scripts", "run_sctype.R"
)


class ScTypeRunner(BaseAnnotationRunner):
    """
    Runs scType via Rscript and parses the output CSV.

    Additional kwargs accepted by run():
      tissue_type : str  — scType tissue type string, e.g. 'Immune system', 'Skin'
    """

    @property
    def tool_name(self) -> str:
        return "sctype"

    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict[str, str]:
        if shutil.which("Rscript") is None:
            logger.error(f"[{self.tool_name}] Rscript not found in PATH.")
            return {}

        tissue_type = kwargs.get("tissue_type", "Immune system")
        tool_dir = self.tool_output_dir(output_dir)
        out_csv = os.path.join(tool_dir, "sctype_predictions.csv")
        r_script = os.path.abspath(_RSCRIPT)

        if not os.path.exists(r_script):
            logger.error(f"[{self.tool_name}] R script not found: {r_script}")
            return {}

        cmd = ["Rscript", r_script, deg_csv, out_csv, tissue_type]
        logger.info(f"[{self.tool_name}] Running Rscript run_sctype.R...")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"[{self.tool_name}] FAILED:\n{result.stderr[-1000:]}")
            return {}

        if not os.path.exists(out_csv):
            logger.error(f"[{self.tool_name}] Output CSV not produced: {out_csv}")
            return {}

        df = pd.read_csv(out_csv)
        predictions = dict(zip(df["cluster"].astype(str), df["predicted_cell_type"]))
        logger.info(f"[{self.tool_name}] Annotated {len(predictions)} clusters.")
        return predictions
