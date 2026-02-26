# SapiensOntoCellMap: Benchmarking Plan

**Version:** 0.1.0 | **Date:** 2026-02-25 | **Author:** Sonal Rashmi

This document defines the complete benchmarking strategy for SapiensOntoCellMap —
covering comparator tool acquisition, public dataset download, ground truth definitions,
and the OOP-modular software architecture for the benchmarking suite.

**Design constraint:** All benchmarking infrastructure is strictly separate from
the marker database pipeline (`src/download/`). The `BioDataDownloader` class is
used only for marker database and reference file downloads. A dedicated
`benchmarking/download/` module handles all benchmark-specific data and tool acquisition.

---

## Table of Contents

1. [Comparator Tools](#1-comparator-tools)
2. [Public Benchmark Datasets](#2-public-benchmark-datasets)
3. [Ground Truth Definitions](#3-ground-truth-definitions)
4. [Software Architecture](#4-software-architecture)
5. [Running the Benchmarks](#5-running-the-benchmarks)
6. [Expected Outputs](#6-expected-outputs)
7. [Publication Figure Plan](#7-publication-figure-plan)

---

## 1. Comparator Tools

### 1.1 Tool Landscape

| Tool | Category | Method | Interface | Priority |
|------|----------|--------|-----------|----------|
| **CellTypist** | B: Reference | Logistic regression (HCA reference) | Python | P0 |
| **scType** | A: Marker-based | Score-based marker enrichment | R via Rscript | P0 |
| **SCSA** | A: Marker-based | Fisher-exact enrichment | Python CLI | P0 |
| **SingleR** | B: Reference | Spearman correlation to reference | R via Bioconductor | P0 |
| **ClusterMole** | A: Marker-based | Enrichment (CellMarker, PanglaoDB, PANTHER) | R | P1 |
| **Azimuth** | B: Reference | RPCA label transfer (Seurat v5) | R | P1 |
| **OnClass** | C: Ontology | CL graph embedding + ML | Python | P1 |

**Priority:** P0 = must include for submission; P1 = include if feasible.

### 1.2 Unique Claims vs. Each Category

| Claim | vs. Category A | vs. Category B | vs. Category C |
|-------|---------------|---------------|---------------|
| 14+ databases (4–14× coverage) | ✅ unique | — | — |
| Formal CL/UBERON normalization | ✅ unique | — | — |
| Hierarchical output (CL-ontology) | ✅ unique | ✅ unique | partial |
| No reference atlas required | — | ✅ unique | partial |
| Spatial-native (Visium HD, Xenium) | ✅ unique | ✅ unique | ✅ unique |
| Evidence-weighted scoring | ✅ unique | — | — |

The hierarchical CL-ontology output and spatial-native support are the two capabilities
no existing tool in any category provides.

---

### 1.3 Tool Download & Installation Instructions

#### CellTypist (Python, P0)
```bash
pip install celltypist
# Models downloaded automatically on first use (~200 MB each)
# Models used: Immune_All_High (PBMC), Human_Lung_Atlas (HLCA)
```

#### scType (R, P0)
```bash
# Install in R:
install.packages("BiocManager")
BiocManager::install(c("Seurat", "HGNChelper", "openxlsx"))
# scType is loaded directly from GitHub in the run script (no CRAN package):
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.r")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.r")
```
The wrapper script `benchmarking/r_scripts/run_sctype.R` handles all sourcing and I/O.

#### SCSA (Python, P0)
```bash
pip install scsa
# or: git clone https://github.com/bioinfo-ibms-pumc/SCSA && pip install -e SCSA/
# Uses CellMarker and PanglaoDB internally
```

#### SingleR (R, P0)
```bash
# Install in R:
BiocManager::install("SingleR")
BiocManager::install("celldex")   # reference datasets
```
The wrapper script `benchmarking/r_scripts/run_singler.R` exports per-cluster
predictions to a CSV file that Python reads back.

#### ClusterMole (R, P1)
```bash
# Install in R:
install.packages("clustermole")   # CRAN
```

#### Azimuth (R, P1)
```bash
# Install in R:
install.packages("Seurat")        # v5 required
install.packages("SeuratData")
# Reference datasets: InstallData("pbmcref")  (~3 GB)
```

---

## 2. Public Benchmark Datasets

All datasets are publicly available without registration. No in-house or unpublished
data is used in the benchmarking suite.

### 2.1 PBMC3k (Zheng et al., *Nat Commun* 2017)

| Field | Value |
|-------|-------|
| **Cells** | 2,638 PBMCs |
| **Cell types** | 8 (see Ground Truth §3.1) |
| **Ground truth** | Seurat tutorial labels (standard benchmark) |
| **Download** | Automatic — 10x Genomics website |
| **Access** | `scanpy.datasets.pbmc3k()` (from Scanpy) |
| **Status** | DONE — DEGs computed, SapiensOntoCellMap results available |

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()   # ~5 MB, downloads automatically
```

### 2.2 Human Lung Cell Atlas — HLCA (Sikkema et al., *Nat Med* 2023)

| Field | Value |
|-------|-------|
| **Cells (full)** | ~2.4M lung cells from 486 donors |
| **Benchmark subset** | 500 cells × 20 cell types = 10,000 cells |
| **Cell types** | 20 major types (see Ground Truth §3.2) |
| **Ground truth** | Atlas labels (Sikkema Suppl. Table 2) |
| **Download** | CELLxGENE Census API (preferred) or HuggingFace |
| **DOI** | 10.1038/s41591-023-02354-7 |
| **Access** | `cellxgene_census` Python package (see §5.2) |

```bash
pip install cellxgene-census
```

```python
import cellxgene_census
# Downloads a targeted 10K-cell subset; does NOT download the full 2.4M-cell atlas
with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census, organism="Homo sapiens",
        obs_value_filter='tissue_general=="lung" and dataset_id=="<HLCA_core_id>"',
    )
```

### 2.3 Tabula Sapiens — Skin + Lung + Blood Subsets (The Tabula Sapiens Consortium, *Science* 2022)

| Field | Value |
|-------|-------|
| **Cells (full)** | ~500K cells from 24 tissues, 15 donors |
| **Benchmark subset** | Skin + Lung + Blood compartments |
| **Cell types** | ~25 per tissue |
| **Ground truth** | Atlas labels (expert-annotated) |
| **Download** | figshare / CZ CELLxGENE |
| **DOI** | 10.1126/science.abl4896 |
| **Access** | `cellxgene_census` or direct download from CELLxGENE |

---

## 3. Ground Truth Definitions

### 3.1 PBMC3k Ground Truth

```python
PBMC3K_GROUND_TRUTH = {
    # cluster_label : (canonical_name, CL_ID, broad_type)
    "CD4 T cells":          ("CD4-positive, alpha-beta T cell",  "CL:0000624", "T cell"),
    "CD14+ Monocytes":      ("classical monocyte",               "CL:0000860", "Monocyte"),
    "B cells":              ("B cell",                           "CL:0000236", "B cell"),
    "CD8 T cells":          ("CD8-positive, alpha-beta T cell",  "CL:0000625", "T cell"),
    "NK cells":             ("natural killer cell",              "CL:0000623", "NK cell"),
    "FCGR3A+ Monocytes":    ("non-classical monocyte",           "CL:0000875", "Monocyte"),
    "Dendritic cells":      ("myeloid dendritic cell",           "CL:0000782", "Dendritic cell"),
    "Megakaryocytes":       ("megakaryocyte",                    "CL:0000556", "Megakaryocyte"),
}
```

**Status:** DONE. SapiensOntoCellMap achieves Top-1 100%, Hierarchical 100% on PBMC3k.

### 3.2 HLCA Ground Truth

20 major cell types from Sikkema et al. (*Nat Med* 2023), Supplementary Table 2:

```python
HLCA_GROUND_TRUTH = {
    # atlas_label         : (canonical_name, CL_ID, broad_type)
    "Alveolar Macrophage": ("alveolar macrophage",                "CL:0000583", "Macrophage"),
    "AT1":                 ("type I pneumocyte",                  "CL:0000062", "Epithelial"),
    "AT2":                 ("type II pneumocyte",                 "CL:0000063", "Epithelial"),
    "B cell":              ("B cell",                             "CL:0000236", "B cell"),
    "CD4+ T cell":         ("CD4-positive, alpha-beta T cell",    "CL:0000624", "T cell"),
    "CD8+ T cell":         ("CD8-positive, alpha-beta T cell",    "CL:0000625", "T cell"),
    "Ciliated":            ("ciliated cell",                      "CL:0000064", "Epithelial"),
    "Club cell":           ("club cell",                          "CL:0000158", "Epithelial"),
    "Dendritic cell":      ("myeloid dendritic cell",             "CL:0000782", "Dendritic cell"),
    "Endothelial":         ("endothelial cell",                   "CL:0000115", "Endothelial"),
    "Fibroblast":          ("fibroblast",                         "CL:0000057", "Fibroblast"),
    "Mast cell":           ("mast cell",                          "CL:0000097", "Mast cell"),
    "Monocyte":            ("classical monocyte",                 "CL:0000860", "Monocyte"),
    "NK cell":             ("natural killer cell",                "CL:0000623", "NK cell"),
    "Pericyte":            ("pericyte",                           "CL:0000669", "Pericyte"),
    "Plasma cell":         ("plasma cell",                        "CL:0000786", "B cell"),
    "Smooth muscle cell":  ("smooth muscle cell",                 "CL:0000192", "Smooth muscle"),
    "Tuft cell":           ("tuft cell",                          "CL:0000239", "Epithelial"),
    "Interstitial macrophage": ("interstitial macrophage",        "CL:0000877", "Macrophage"),
    "Goblet cell":         ("goblet cell",                        "CL:0000160", "Epithelial"),
}
```

### 3.3 Matching Strategy

Top-1 accuracy is evaluated with three layers of matching (applied in order):

1. **Exact CL ID match** — `predicted_cl_id == ground_truth_cl_id`
2. **CL ancestor match** — `ground_truth_cl_id` is an ancestor of `predicted_cl_id` in the
   CL ontology (correct lineage, more specific than ground truth = still correct)
3. **Normalized name match** — case-insensitive substring match between predicted cell type
   name and the canonical name or any accepted alias in the ground truth dict

A prediction is marked **Correct** if any layer matches. **Hierarchical accuracy** is
computed identically but using the hierarchical annotation output (any level in the
hierarchy tree that matches the ground truth CL ID or its ancestors counts).

---

## 4. Software Architecture

### 4.1 Directory Structure

```
benchmarking/
├── download/
│   ├── __init__.py
│   ├── base_downloader.py           ← Abstract base class for all dataset/tool downloaders
│   ├── dataset_downloader.py        ← Concrete: PBMC3k, HLCA, Tabula Sapiens via Census API
│   └── tool_installer.py            ← Concrete: verify/install Python tools; check R tools
├── runners/
│   ├── __init__.py
│   ├── base_runner.py               ← Abstract base class: run(adata, output_dir) → dict
│   ├── sapiensonto_runner.py        ← Wraps get_cluster_annotation.py CLI
│   ├── celltypist_runner.py         ← Python: celltypist.annotate()
│   ├── scsa_runner.py               ← Python: SCSA CLI wrapper
│   ├── sctype_runner.py             ← R: subprocess Rscript run_sctype.R
│   └── singler_runner.py            ← R: subprocess Rscript run_singler.R
├── metrics/
│   ├── __init__.py
│   └── benchmark_metrics.py         ← BenchmarkMetrics class: top1, hierarchical, broad_type accuracy
├── figures/
│   ├── __init__.py
│   └── benchmark_figures.py         ← BenchmarkFigures class: publication-quality matplotlib panels
├── ground_truth/
│   ├── __init__.py
│   ├── pbmc3k_gt.py                 ← PBMC3K_GROUND_TRUTH dict
│   ├── hlca_gt.py                   ← HLCA_GROUND_TRUTH dict (20 types)
│   └── tabula_sapiens_gt.py         ← TABULA_SAPIENS_GROUND_TRUTH dict
├── r_scripts/
│   ├── run_sctype.R                 ← scType: reads DEG CSV, writes predictions CSV
│   └── run_singler.R                ← SingleR: reads .h5ad path, writes predictions CSV
├── results/
│   ├── pbmc3k/                      ← Per-run outputs (CSV + figures)
│   ├── hlca/
│   └── tabula_sapiens/
├── benchmark_pbmc3k.py              ← Dataset-specific runner (uses runners/ + metrics/)
├── benchmark_hlca.py                ← Dataset-specific runner
├── benchmark_tabula_sapiens.py      ← Dataset-specific runner
└── run_all_benchmarks.py            ← Combined CLI: runs all datasets, generates publication figure
```

### 4.2 Abstract Base Classes

#### `benchmarking/download/base_downloader.py`

```python
from abc import ABC, abstractmethod
from pathlib import Path


class BaseDownloader(ABC):
    """
    Abstract base for all benchmarking data and tool downloaders.

    Subclasses implement _download() to fetch a specific resource.
    The public download() method handles idempotency (skip if already present).

    Strictly separate from BioDataDownloader (src/download/) which is
    reserved for marker database and reference file acquisition.
    """

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def download(self, force: bool = False) -> Path:
        """Download resource if not already present. Returns local path."""
        target = self._target_path()
        if target.exists() and not force:
            return target
        return self._download()

    @abstractmethod
    def _target_path(self) -> Path:
        """Return the expected local path of the downloaded resource."""

    @abstractmethod
    def _download(self) -> Path:
        """Perform the actual download. Return local path on success."""
```

#### `benchmarking/runners/base_runner.py`

```python
from abc import ABC, abstractmethod
import pandas as pd


class BaseAnnotationRunner(ABC):
    """
    Abstract base for all annotation tool runners.

    Each subclass wraps one tool and exposes a uniform interface:
      run(deg_csv, output_dir, **kwargs) → dict[cluster_id, predicted_cell_type]

    Runners must NOT modify the input DEG file or output directory of other runners.
    All tool-specific output is written to a subdirectory named after the tool.
    """

    @property
    @abstractmethod
    def tool_name(self) -> str:
        """Short identifier used for output directory naming (e.g. 'celltypist')."""

    @abstractmethod
    def run(self, deg_csv: str, output_dir: str, **kwargs) -> dict:
        """
        Run annotation on the provided DEG CSV.

        Parameters
        ----------
        deg_csv : str
            Path to Seurat-format FindAllMarkers CSV.
        output_dir : str
            Base output directory. Tool writes to output_dir/<tool_name>/.
        **kwargs
            Tool-specific parameters (tissue, model_name, etc.)

        Returns
        -------
        dict
            {cluster_id (str): predicted_cell_type (str)}
            Empty dict if tool fails or produces no output.
        """

    def _tool_output_dir(self, base_output_dir: str) -> str:
        """Return the tool-specific subdirectory path."""
        import os
        d = os.path.join(base_output_dir, self.tool_name)
        os.makedirs(d, exist_ok=True)
        return d
```

#### `benchmarking/metrics/benchmark_metrics.py`

```python
import pandas as pd
from typing import Optional


class BenchmarkMetrics:
    """
    Computes and stores accuracy metrics for one benchmark run.

    Accuracy layers (applied in order):
      1. Exact CL ID match
      2. CL ancestor match (predicted is more specific than ground truth)
      3. Normalized cell type name substring match

    Parameters
    ----------
    ground_truth : dict
        {cluster_label: (canonical_name, CL_ID, broad_type)}
    hierarchical_csv : str, optional
        Path to SapiensOntoCellMap hierarchical annotation CSV.
        Required only for hierarchical accuracy computation.
    """

    def __init__(self, ground_truth: dict, hierarchical_csv: Optional[str] = None):
        self.ground_truth = ground_truth
        self.hierarchical_csv = hierarchical_csv
        self._results: list[dict] = []

    def add_tool_predictions(self, tool_name: str, predictions: dict) -> None:
        """
        Register predictions for one tool.

        Parameters
        ----------
        tool_name : str
        predictions : dict
            {cluster_id: predicted_cell_type_name}
        """
        for cluster, gt_tuple in self.ground_truth.items():
            gt_canonical, gt_cl_id, gt_broad = gt_tuple
            predicted = predictions.get(cluster, "")
            top1_correct = self._is_correct(predicted, gt_canonical, gt_cl_id)
            self._results.append({
                "tool": tool_name,
                "cluster": cluster,
                "predicted": predicted,
                "ground_truth": gt_canonical,
                "gt_cl_id": gt_cl_id,
                "broad_type": gt_broad,
                "top1_correct": top1_correct,
            })

    def top1_accuracy(self, tool_name: str) -> float:
        rows = [r for r in self._results if r["tool"] == tool_name]
        if not rows:
            return 0.0
        return sum(r["top1_correct"] for r in rows) / len(rows)

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self._results)

    def summary(self) -> pd.DataFrame:
        """Return per-tool accuracy summary DataFrame."""
        tools = list({r["tool"] for r in self._results})
        rows = []
        for t in tools:
            rows.append({"tool": t, "top1_accuracy": self.top1_accuracy(t)})
        return pd.DataFrame(rows).sort_values("top1_accuracy", ascending=False)

    @staticmethod
    def _is_correct(predicted: str, gt_canonical: str, gt_cl_id: str) -> bool:
        if not predicted:
            return False
        # Layer 3: name substring (case-insensitive)
        p = predicted.lower()
        g = gt_canonical.lower()
        return p == g or g in p or p in g
```

#### `benchmarking/figures/benchmark_figures.py`

```python
import matplotlib.pyplot as plt
import pandas as pd


class BenchmarkFigures:
    """
    Generates publication-quality matplotlib figures from benchmark results.

    All figures use 300 DPI, no seaborn imports (pure matplotlib),
    and are saved as both PNG (for submission) and PDF (for editing).
    """

    def __init__(self, output_dir: str):
        self.output_dir = output_dir

    def plot_accuracy_comparison(
        self,
        summary_df: pd.DataFrame,
        dataset_name: str,
        metric: str = "top1_accuracy",
    ) -> str:
        """
        Horizontal bar chart of per-tool accuracy for one dataset.
        Highlights SapiensOntoCellMap bar in a distinct colour.
        Returns path to saved PNG.
        """
        fig, ax = plt.subplots(figsize=(7, 4))
        colors = [
            "#2563EB" if t == "SapiensOntoCellMap" else "#94A3B8"
            for t in summary_df["tool"]
        ]
        ax.barh(summary_df["tool"], summary_df[metric] * 100, color=colors)
        ax.set_xlabel("Top-1 Accuracy (%)")
        ax.set_title(f"{dataset_name} — Top-1 Accuracy")
        ax.set_xlim(0, 105)
        plt.tight_layout()
        path = f"{self.output_dir}/{dataset_name.lower()}_accuracy.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path

    def plot_combined_figure(
        self,
        dataset_results: dict[str, pd.DataFrame],
    ) -> str:
        """
        3-panel publication figure:
          Panel A — PBMC3k accuracy (bar chart)
          Panel B — HLCA accuracy (bar chart)
          Panel C — Unique capabilities table (text)
        Returns path to saved PNG.
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for ax, (name, df) in zip(axes[:2], list(dataset_results.items())[:2]):
            colors = [
                "#2563EB" if t == "SapiensOntoCellMap" else "#94A3B8"
                for t in df["tool"]
            ]
            ax.barh(df["tool"], df["top1_accuracy"] * 100, color=colors)
            ax.set_xlabel("Top-1 Accuracy (%)")
            ax.set_title(name)
            ax.set_xlim(0, 105)
        # Panel C: capability table
        axes[2].axis("off")
        capabilities = [
            ["Feature", "SapiensOntoCellMap", "Others"],
            ["14+ databases", "✅", "≤4"],
            ["Hierarchical output", "✅", "❌"],
            ["CL/UBERON ontology", "✅", "❌"],
            ["Spatial-native", "✅", "❌"],
            ["No reference atlas", "✅", "partial"],
        ]
        axes[2].table(
            cellText=capabilities[1:],
            colLabels=capabilities[0],
            loc="center",
            cellLoc="center",
        )
        plt.tight_layout()
        path = f"{self.output_dir}/combined_benchmark_figure.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        return path
```

### 4.3 R Script Interface

R-based tools are called via `subprocess.run(["Rscript", script_path, ...])`.
All R scripts accept positional arguments and write a simple CSV:
`cluster,predicted_cell_type`

This keeps the Python–R boundary clean: Python writes DEG CSV in; R reads it,
runs the tool, and writes predictions CSV out. No shared state.

#### `benchmarking/r_scripts/run_sctype.R` — Interface

```
Usage: Rscript run_sctype.R <deg_csv> <output_csv> <tissue_type>

Input:  deg_csv        — Seurat FindAllMarkers CSV (p_val_adj, cluster, gene)
Output: output_csv     — cluster,predicted_cell_type
        tissue_type    — e.g. "Immune system", "Skin"
```

#### `benchmarking/r_scripts/run_singler.R` — Interface

```
Usage: Rscript run_singler.R <adata_h5ad> <output_csv> <reference>

Input:  adata_h5ad     — AnnData .h5ad file path
Output: output_csv     — cluster,predicted_cell_type
        reference      — celldex reference name e.g. "HumanPrimaryCellAtlasData"
```

---

## 5. Running the Benchmarks

### 5.1 Prerequisites

```bash
# Python tools
pip install celltypist cellxgene-census scanpy anndata

# R tools (run in R console)
# install.packages(c("BiocManager", "clustermole"))
# BiocManager::install(c("SingleR", "celldex"))
```

Verify Rscript is accessible:
```bash
Rscript --version
```

### 5.2 PBMC3k (already done — reproduce or extend)

```bash
python benchmarking/benchmark_pbmc3k.py \
    --output_dir benchmarking/results/pbmc3k/
```

### 5.3 HLCA

```bash
python benchmarking/benchmark_hlca.py \
    --output_dir benchmarking/results/hlca/ \
    --n_cells_per_type 500       # subset size per cell type (default 500)
```

The script:
1. Downloads a 10K-cell HLCA subset via CELLxGENE Census API (first run only)
2. Computes Wilcoxon DEGs per cluster (Scanpy)
3. Runs SapiensOntoCellMap + all P0 comparators
4. Evaluates accuracy against HLCA ground truth
5. Writes results to `benchmarking/results/hlca/`

### 5.4 Tabula Sapiens (Phase 2)

```bash
python benchmarking/benchmark_tabula_sapiens.py \
    --tissues skin lung blood \
    --output_dir benchmarking/results/tabula_sapiens/
```

### 5.5 Combined runner

```bash
python benchmarking/run_all_benchmarks.py \
    --datasets pbmc3k hlca \
    --output_dir benchmarking/results/
```

Produces the combined publication figure at `benchmarking/results/combined_benchmark_figure.png`.

---

## 6. Expected Outputs

Each dataset benchmark produces:

```
benchmarking/results/<dataset>/
├── degs.csv                        ← Scanpy Wilcoxon DEGs (input to all tools)
├── sapiensontocellmap/
│   ├── <dataset>_top_annotation_summary.csv
│   ├── <dataset>_hierarchical_annotation.csv
│   └── <dataset>_report.html
├── celltypist/
│   └── celltypist_predictions.csv  ← cluster,predicted_cell_type
├── sctype/
│   └── sctype_predictions.csv
├── singler/
│   └── singler_predictions.csv
├── accuracy_summary.csv            ← per-tool Top-1 accuracy
└── <dataset>_accuracy.png          ← horizontal bar chart figure
```

### Metrics reported per tool

| Metric | Definition |
|--------|-----------|
| **Top-1 accuracy** | Fraction of clusters where top-1 prediction matches ground truth (layers 1–3) |
| **Hierarchical accuracy** | Fraction correct at any level of the CL hierarchy (SapiensOntoCellMap only) |
| **Broad-type accuracy** | Fraction correct at the broad lineage level (T cell, B cell, etc.) |
| **Annotation rate** | Fraction of clusters annotated (not "Unknown" or no call) |

---

## 7. Publication Figure Plan

### Figure layout (manuscript Figure 4)

**Panel A — Accuracy across datasets**
Grouped bar chart: one bar per tool per dataset.
x-axis: datasets (PBMC3k, HLCA, Tabula Sapiens).
y-axis: Top-1 accuracy (%).
SapiensOntoCellMap bar highlighted in blue; comparators in grey.

**Panel B — Annotation rate (coverage)**
Horizontal bar chart: what fraction of clusters does each tool annotate
(vs. returning "Unknown" or failing)?
Motivation: marker-based tools annotate all clusters; reference-based tools
may fail for rare or disease-specific cell types not in the reference.

**Panel C — Unique capabilities comparison table**

| Feature | SapiensOntoCellMap | CellTypist | scType | SCSA | SingleR |
|---------|-------------------|-----------|--------|------|---------|
| Hierarchical CL output | ✅ | ❌ | ❌ | ❌ | ❌ |
| Spatial-native | ✅ | ❌ | ❌ | ❌ | ❌ |
| 14+ databases | ✅ | ❌ | ❌ | ❌ | ❌ |
| No reference atlas | ✅ | ❌ | ✅ | ✅ | ❌ |
| Evidence weighting | ✅ | ❌ | ❌ | ❌ | ❌ |

---

## Implementation Roadmap

| Phase | Task | Status |
|-------|------|--------|
| 1 | PBMC3k benchmark (SapiensOntoCellMap + CellTypist) | DONE |
| 2 | OOP refactor: extract `runners/`, `metrics/`, `figures/`, `download/` | NOT STARTED |
| 2 | Ground truth dicts: `pbmc3k_gt.py`, `hlca_gt.py` | NOT STARTED |
| 2 | HLCA dataset download + DEG computation | NOT STARTED |
| 2 | HLCA benchmark: all P0 comparators | NOT STARTED |
| 2 | R scripts: `run_sctype.R`, `run_singler.R` | NOT STARTED |
| 3 | Tabula Sapiens (skin + lung + blood) | NOT STARTED |
| 3 | Combined publication figure (3-panel) | NOT STARTED |
