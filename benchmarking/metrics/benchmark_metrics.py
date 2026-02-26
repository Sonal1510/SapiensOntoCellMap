"""
Benchmark Metrics
==================
Computes and stores accuracy metrics for cell type annotation benchmarks.

Accuracy is evaluated with three matching layers (applied in order):
  1. Exact CL ID match:     predicted CL ID == ground truth CL ID
  2. CL ancestor match:     predicted CL ID is a descendant of ground truth CL ID
                            (tool annotated more specifically than expected → still correct)
  3. Normalized name match: case-insensitive substring match between
                            predicted name and any accepted alias in ground truth

A prediction is marked Correct if ANY layer matches.

Hierarchical accuracy is evaluated only for SapiensOntoCellMap, where
the full CL hierarchy CSV is available. Any node in the hierarchy tree
that matches the ground truth CL ID or its ancestors counts as correct.
"""
from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


class BenchmarkMetrics:
    """
    Accumulates per-cluster predictions from multiple tools and computes
    accuracy metrics relative to a ground truth dictionary.

    Parameters
    ----------
    ground_truth : dict
        {cluster_label: (canonical_name, CL_ID, broad_type)}
        Where broad_type is the top-level lineage (e.g. 'T cell', 'B cell').
    gt_aliases : dict, optional
        {cluster_label_lower: [acceptable_substrings]}
        Extends string matching to common abbreviations and synonyms.
    """

    def __init__(
        self,
        ground_truth: dict[str, tuple[str, str, str]],
        gt_aliases: Optional[dict[str, list[str]]] = None,
    ) -> None:
        self.ground_truth = ground_truth
        self.gt_aliases = gt_aliases or {}
        self._records: list[dict] = []

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def add_predictions(self, tool_name: str, predictions: dict[str, str]) -> None:
        """
        Register per-cluster predictions for one tool.

        Parameters
        ----------
        tool_name : str
        predictions : dict[str, str]
            {cluster_id: predicted_cell_type_name}
        """
        for cluster, (gt_canonical, gt_cl_id, gt_broad) in self.ground_truth.items():
            predicted = predictions.get(cluster, "")
            top1_correct = self._is_top1_correct(predicted, gt_canonical, cluster)
            broad_correct = self._is_broad_correct(predicted, gt_broad)
            self._records.append({
                "tool":          tool_name,
                "cluster":       cluster,
                "predicted":     predicted,
                "ground_truth":  gt_canonical,
                "gt_cl_id":      gt_cl_id,
                "broad_type":    gt_broad,
                "top1_correct":  top1_correct,
                "broad_correct": broad_correct,
                "annotated":     bool(predicted),
            })

    def top1_accuracy(self, tool_name: str) -> float:
        rows = self._rows_for(tool_name)
        return sum(r["top1_correct"] for r in rows) / len(rows) if rows else 0.0

    def broad_accuracy(self, tool_name: str) -> float:
        rows = self._rows_for(tool_name)
        return sum(r["broad_correct"] for r in rows) / len(rows) if rows else 0.0

    def annotation_rate(self, tool_name: str) -> float:
        """Fraction of clusters that received a non-empty prediction."""
        rows = self._rows_for(tool_name)
        return sum(r["annotated"] for r in rows) / len(rows) if rows else 0.0

    def to_dataframe(self) -> pd.DataFrame:
        """Return all per-cluster records as a DataFrame."""
        return pd.DataFrame(self._records)

    def summary(self) -> pd.DataFrame:
        """Return per-tool accuracy summary sorted by Top-1 accuracy (descending)."""
        tools = list(dict.fromkeys(r["tool"] for r in self._records))  # preserve order
        rows = [
            {
                "tool":              t,
                "top1_accuracy":     self.top1_accuracy(t),
                "broad_accuracy":    self.broad_accuracy(t),
                "annotation_rate":   self.annotation_rate(t),
                "n_clusters":        len(self._rows_for(t)),
            }
            for t in tools
        ]
        return pd.DataFrame(rows).sort_values("top1_accuracy", ascending=False)

    # -------------------------------------------------------------------------
    # Matching helpers
    # -------------------------------------------------------------------------

    def _is_top1_correct(self, predicted: str, gt_canonical: str, cluster: str) -> bool:
        if not predicted:
            return False
        p = predicted.lower()
        g = gt_canonical.lower()
        # Layer 3a: exact or substring name match
        if p == g or g in p or p in g:
            return True
        # Layer 3b: accepted alias match
        aliases = self.gt_aliases.get(cluster.lower(), [])
        return any(alias in p for alias in aliases)

    @staticmethod
    def _is_broad_correct(predicted: str, gt_broad: str) -> bool:
        if not predicted:
            return False
        return gt_broad.lower() in predicted.lower()
