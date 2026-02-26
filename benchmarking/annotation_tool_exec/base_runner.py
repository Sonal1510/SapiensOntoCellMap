"""
Base Annotation Tool Runner
============================
Abstract base class that all annotation tool wrappers must implement.

Interface contract:
  run(deg_csv, output_dir, **kwargs) → dict[cluster_id, predicted_cell_type]

Every runner writes its tool-specific output to:
  <output_dir>/<tool_name>/

No runner may write outside its designated subdirectory.
"""
from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)


class BaseAnnotationRunner(ABC):
    """
    Abstract base for annotation tool runners.

    Subclasses wrap one external tool (Python package or R script) and
    expose a uniform interface. The caller never needs to know the
    tool's internal API — only the DEG CSV path and output directory.
    """

    @property
    @abstractmethod
    def tool_name(self) -> str:
        """
        Short identifier used for output subdirectory and logging.
        Examples: 'sapiensonto', 'celltypist', 'sctype', 'singler'.
        """

    @abstractmethod
    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict[str, str]:
        """
        Run cell type annotation on the provided DEG CSV.

        Parameters
        ----------
        deg_csv : str
            Path to Seurat FindAllMarkers or Scanpy-format DEG CSV.
        output_dir : str
            Base output directory. All output is written to
            <output_dir>/<tool_name>/.
        **kwargs
            Tool-specific parameters (tissue, model_name, background_n, etc.)

        Returns
        -------
        dict[str, str]
            Mapping of {cluster_id: predicted_cell_type_name}.
            Returns an empty dict if the tool fails or produces no predictions.
        """

    def tool_output_dir(self, base_output_dir: str) -> str:
        """Create and return the tool-specific subdirectory."""
        d = os.path.join(base_output_dir, self.tool_name)
        os.makedirs(d, exist_ok=True)
        return d
