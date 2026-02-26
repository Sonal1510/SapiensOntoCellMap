"""
CellTypist Runner
==================
Wraps CellTypist (Dominguez Conde et al. Science 2022) for benchmark comparison.

CellTypist uses logistic regression trained on the Human Cell Atlas reference.
It operates on AnnData objects (not DEG CSVs), so this runner accepts an
optional adata_path parameter alongside deg_csv.

Requires: pip install celltypist
"""
from __future__ import annotations

import logging
import os

import pandas as pd

from benchmarking.annotation_tool_exec.base_runner import BaseAnnotationRunner

logger = logging.getLogger(__name__)

# Default models: Immune_All_High for blood/PBMC; Human_Lung_Atlas for lung
DEFAULT_MODEL = "Immune_All_High"


class CellTypistRunner(BaseAnnotationRunner):
    """
    Runs CellTypist annotation and converts per-cell predictions to
    per-cluster majority-vote predictions.

    Additional kwargs accepted by run():
      adata_path : str     — path to .h5ad file (required; deg_csv is ignored)
      model_name : str     — CellTypist model name (default: Immune_All_High)
      cluster_key : str    — adata.obs column for cluster assignment
                             (default: 'louvain')
    """

    @property
    def tool_name(self) -> str:
        return "celltypist"

    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict[str, str]:
        try:
            import celltypist
            import scanpy as sc
        except ImportError as exc:
            logger.error(f"[{self.tool_name}] Not installed: {exc}. pip install celltypist")
            return {}

        adata_path = kwargs.get("adata_path")
        if not adata_path or not os.path.exists(adata_path):
            logger.error(f"[{self.tool_name}] adata_path is required and must exist.")
            return {}

        model_name = kwargs.get("model_name", DEFAULT_MODEL)
        cluster_key = kwargs.get("cluster_key", "louvain")
        tool_dir = self.tool_output_dir(output_dir)

        logger.info(f"[{self.tool_name}] Loading {adata_path}...")
        adata = sc.read_h5ad(adata_path)

        logger.info(f"[{self.tool_name}] Downloading model '{model_name}' if needed...")
        celltypist.models.download_models(model=model_name, force_update=False)

        logger.info(f"[{self.tool_name}] Running annotation...")
        predictions = celltypist.annotate(
            adata, model=model_name, majority_voting=True
        )
        result_adata = predictions.to_adata()

        # Per-cell predictions → per-cluster majority vote
        if cluster_key not in result_adata.obs.columns:
            logger.error(
                f"[{self.tool_name}] Cluster key '{cluster_key}' not in adata.obs. "
                f"Available: {list(result_adata.obs.columns)}"
            )
            return {}

        per_cell = result_adata.obs[[cluster_key, "majority_voting"]].copy()
        per_cell.columns = ["cluster", "predicted"]
        majority = (
            per_cell.groupby("cluster")["predicted"]
            .agg(lambda x: x.value_counts().idxmax())
            .to_dict()
        )

        out_csv = os.path.join(tool_dir, "celltypist_predictions.csv")
        pd.DataFrame(
            [{"cluster": k, "predicted_cell_type": v} for k, v in majority.items()]
        ).to_csv(out_csv, index=False)

        logger.info(f"[{self.tool_name}] Annotated {len(majority)} clusters → {out_csv}")
        return majority
