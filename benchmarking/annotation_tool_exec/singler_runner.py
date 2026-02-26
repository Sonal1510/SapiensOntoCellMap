"""
SingleR Runner
==============
Wraps SingleR (Aran et al. Nat Immunol 2019) via an Rscript subprocess.

The Python–R boundary is CSV-based:
  Input:  .h5ad path (AnnData exported to disk)
  Output: <tool_dir>/singler_predictions.csv  (columns: cluster, predicted_cell_type)

Requires: Rscript in PATH, R packages: SingleR, celldex (BioConductor)
  Install: BiocManager::install(c("SingleR", "celldex"))
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
    os.path.dirname(__file__), "..", "r_scripts", "run_singler.R"
)


class SingleRRunner(BaseAnnotationRunner):
    """
    Runs SingleR via Rscript and parses the output CSV.

    Additional kwargs accepted by run():
      adata_path  : str — path to .h5ad file (required; deg_csv is ignored)
      reference   : str — celldex reference name
                          e.g. 'HumanPrimaryCellAtlasData', 'MonacoImmuneData'
      cluster_key : str — adata.obs column for cluster (default: 'louvain')
    """

    @property
    def tool_name(self) -> str:
        return "singler"

    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict[str, str]:
        if shutil.which("Rscript") is None:
            logger.error(f"[{self.tool_name}] Rscript not found in PATH.")
            return {}

        adata_path = kwargs.get("adata_path")
        if not adata_path or not os.path.exists(adata_path):
            logger.error(f"[{self.tool_name}] adata_path is required and must exist.")
            return {}

        reference = kwargs.get("reference", "HumanPrimaryCellAtlasData")
        cluster_key = kwargs.get("cluster_key", "louvain")
        tool_dir = self.tool_output_dir(output_dir)
        out_csv = os.path.join(tool_dir, "singler_predictions.csv")
        r_script = os.path.abspath(_RSCRIPT)

        if not os.path.exists(r_script):
            logger.error(f"[{self.tool_name}] R script not found: {r_script}")
            return {}

        cmd = ["Rscript", r_script, adata_path, out_csv, reference, cluster_key]
        logger.info(f"[{self.tool_name}] Running Rscript run_singler.R...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

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
