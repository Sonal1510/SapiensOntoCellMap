"""
Benchmark Figures
==================
Generates publication-quality matplotlib figures from benchmark results.

All figures:
  - 300 DPI, no seaborn dependency (pure matplotlib)
  - Saved as PNG (submission) and PDF (vector editing)
  - SapiensOntoCellMap bar highlighted in Anthropic blue (#2563EB)
  - Comparators in neutral grey (#94A3B8)
"""
from __future__ import annotations

import logging
import os

import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger(__name__)

# Colour palette
SAPIENSONTO_COLOR = "#2563EB"
COMPARATOR_COLOR  = "#94A3B8"
ACCENT_COLOR      = "#F59E0B"

FIGSIZE_SINGLE    = (6, 4)
FIGSIZE_COMBINED  = (16, 5)
DPI               = 300


class BenchmarkFigures:
    """
    Produces publication-ready figures from BenchmarkMetrics summary DataFrames.

    Parameters
    ----------
    output_dir : str
        Directory where figures are saved. Created if it does not exist.
    """

    def __init__(self, output_dir: str) -> None:
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # Single-dataset figures
    # -------------------------------------------------------------------------

    def plot_accuracy_bar(
        self,
        summary_df: pd.DataFrame,
        dataset_name: str,
        metric: str = "top1_accuracy",
        title_suffix: str = "",
    ) -> str:
        """
        Horizontal bar chart of per-tool accuracy for one dataset.
        Returns path to saved PNG.
        """
        df = summary_df.sort_values(metric, ascending=True).reset_index(drop=True)
        colors = [
            SAPIENSONTO_COLOR if "sapiensonto" in str(t).lower() else COMPARATOR_COLOR
            for t in df["tool"]
        ]
        metric_label = {
            "top1_accuracy":   "Top-1 Accuracy (%)",
            "broad_accuracy":  "Broad-Type Accuracy (%)",
            "annotation_rate": "Annotation Rate (%)",
        }.get(metric, metric)

        fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
        bars = ax.barh(df["tool"], df[metric] * 100, color=colors, edgecolor="white")
        ax.bar_label(bars, fmt="%.1f%%", padding=3, fontsize=9)
        ax.set_xlabel(metric_label, fontsize=11)
        ax.set_title(
            f"{dataset_name} — {metric_label}{' ' + title_suffix if title_suffix else ''}",
            fontsize=12, fontweight="bold",
        )
        ax.set_xlim(0, 115)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()

        fname = f"{dataset_name.lower().replace(' ', '_')}_{metric}.png"
        path = os.path.join(self.output_dir, fname)
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close(fig)
        logger.info(f"[BenchmarkFigures] Saved {path}")
        return path

    # -------------------------------------------------------------------------
    # Combined multi-dataset publication figure
    # -------------------------------------------------------------------------

    def plot_combined_figure(
        self,
        dataset_summaries: dict[str, pd.DataFrame],
        output_filename: str = "combined_benchmark_figure",
    ) -> str:
        """
        3-panel publication figure:
          Panel A — Dataset 1 Top-1 accuracy (horizontal bar chart)
          Panel B — Dataset 2 Top-1 accuracy (horizontal bar chart)
          Panel C — Unique capabilities table

        Parameters
        ----------
        dataset_summaries : dict[str, pd.DataFrame]
            {dataset_name: summary_df} — only first two entries used for A/B panels.
        output_filename : str
            Base filename (no extension) for output files.

        Returns
        -------
        str
            Path to saved PNG.
        """
        fig, axes = plt.subplots(1, 3, figsize=FIGSIZE_COMBINED)
        items = list(dataset_summaries.items())

        for i, ax in enumerate(axes[:2]):
            if i >= len(items):
                ax.axis("off")
                continue
            name, df = items[i]
            df = df.sort_values("top1_accuracy", ascending=True).reset_index(drop=True)
            colors = [
                SAPIENSONTO_COLOR if "sapiensonto" in str(t).lower() else COMPARATOR_COLOR
                for t in df["tool"]
            ]
            bars = ax.barh(df["tool"], df["top1_accuracy"] * 100, color=colors, edgecolor="white")
            ax.bar_label(bars, fmt="%.1f%%", padding=2, fontsize=8)
            ax.set_xlabel("Top-1 Accuracy (%)", fontsize=10)
            ax.set_title(name, fontsize=11, fontweight="bold")
            ax.set_xlim(0, 115)
            ax.spines[["top", "right"]].set_visible(False)

        # Panel C: capabilities comparison table
        axes[2].axis("off")
        table_data = [
            ["Feature",                   "SapiensOntoCellMap", "Others"],
            ["14+ curated databases",      "✅",                "≤4"],
            ["CL/UBERON normalization",    "✅",                "❌"],
            ["Hierarchical CL output",     "✅",                "❌"],
            ["Spatial-native (Visium/Xenium)", "✅",            "❌"],
            ["No reference atlas needed",  "✅",                "partial"],
            ["Evidence-weighted scoring",  "✅",                "❌"],
        ]
        tbl = axes[2].table(
            cellText=table_data[1:],
            colLabels=table_data[0],
            loc="center",
            cellLoc="center",
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(9)
        tbl.scale(1.2, 1.4)
        # Header row styling
        for j in range(3):
            tbl[(0, j)].set_facecolor("#1E3A8A")
            tbl[(0, j)].set_text_props(color="white", fontweight="bold")
        axes[2].set_title("Unique Capabilities", fontsize=11, fontweight="bold")

        plt.tight_layout()
        path = os.path.join(self.output_dir, f"{output_filename}.png")
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
        fig.savefig(path.replace(".png", ".pdf"), bbox_inches="tight")
        plt.close(fig)
        logger.info(f"[BenchmarkFigures] Saved combined figure: {path}")
        return path
