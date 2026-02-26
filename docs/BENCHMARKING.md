# SapiensOntoCellMap: Benchmarking Reference

**Version:** 0.2.0 | **Date:** 2026-02-25 | **Author:** Sonal Rashmi

This document describes the benchmarking strategy, datasets, comparator tools,
and evaluation metrics used to validate SapiensOntoCellMap for publication.

**Separation of concerns:** All benchmark data and tool acquisition is handled by
`benchmarking/download/` — completely separate from `src/download/`, which is
reserved exclusively for marker database and reference file management.

---

## Table of Contents

1. [Comparator Tools](#1-comparator-tools)
2. [Benchmark Datasets](#2-benchmark-datasets)
3. [Ground Truth Definitions](#3-ground-truth-definitions)
4. [Evaluation Metrics](#4-evaluation-metrics)
5. [Module Architecture](#5-module-architecture)
6. [Running the Benchmarks](#6-running-the-benchmarks)
7. [Expected Outputs](#7-expected-outputs)
8. [Publication Figure Plan](#8-publication-figure-plan)

---

## 1. Comparator Tools

### 1.1 Tool Landscape

| Tool | Category | Method | Reference | Priority |
|------|----------|--------|-----------|----------|
| **CellTypist** | B: Reference | Logistic regression (HCA reference) | Dominguez Conde et al. *Science* 2022 | P0 |
| **scType** | A: Marker-based | Score-based marker enrichment | Ianevski et al. *Nat Commun* 2022 | P0 |
| **SCSA** | A: Marker-based | Fisher-exact enrichment | Cao et al. *Front Genet* 2020 | P0 |
| **SingleR** | B: Reference | Spearman correlation to reference | Aran et al. *Nat Immunol* 2019 | P0 |
| **ClusterMole** | A: Marker-based | Enrichment (CellMarker, PanglaoDB) | Brügger et al. *Bioinformatics* 2024 | P1 |
| **Azimuth** | B: Reference | RPCA label transfer (Seurat v5) | Hao et al. *Nat Biotech* 2024 | P1 |

**Priority:** P0 = included in submission; P1 = included if feasible.

### 1.2 Competitive Positioning

| Claim | vs. Category A | vs. Category B | vs. Category C |
|-------|---------------|---------------|---------------|
| 14+ curated databases | ✅ unique | — | — |
| Formal CL/UBERON normalization | ✅ unique | — | — |
| Hierarchical CL-ontology output | ✅ unique | ✅ unique | partial |
| No reference atlas required | — | ✅ unique | partial |
| Spatial-native (Visium HD, Xenium) | ✅ unique | ✅ unique | ✅ unique |
| Evidence-weighted scoring (6 tiers) | ✅ unique | — | — |

The **hierarchical CL-ontology output** and **spatial-native support** are the two
capabilities no existing published tool in any category provides.

### 1.3 Tool Installation

Installation instructions for all comparator tools are documented in
`benchmarking/download/tools/tool_installer.py` and reproduced here for reference.

**Python tools** (install into the project virtualenv):

| Tool | Command |
|------|---------|
| CellTypist | `pip install celltypist` |
| SCSA | `pip install scsa` |
| cellxgene-census | `pip install cellxgene-census` (dataset access) |
| scanpy | `pip install scanpy` |

**R tools** (run in an R console):

| Tool | Command |
|------|---------|
| scType | `install.packages(c("HGNChelper", "openxlsx"))` — sourced from GitHub at runtime |
| SingleR | `BiocManager::install(c("SingleR", "celldex", "zellkonverter"))` |

Verify Rscript is accessible before running R-based tools:
```bash
Rscript --version
python benchmarking/run_all_benchmarks.py --check_tools
```

---

## 2. Benchmark Datasets

All datasets are publicly available without registration. No in-house or
unpublished data is used in the benchmarking suite.

### 2.1 PBMC3k (Zheng et al., *Nat Commun* 2017)

| Field | Value |
|-------|-------|
| **Cells** | 2,638 PBMCs (peripheral blood mononuclear cells) |
| **Cell types** | 8 (see §3.1) |
| **Ground truth** | Seurat tutorial labels — the de-facto standard for scRNA-seq annotation benchmarks |
| **Access** | `scanpy.datasets.pbmc3k_processed()` — downloads automatically (~5 MB) |
| **DOI** | 10.1038/ncomms14049 |
| **Status** | DONE — SapiensOntoCellMap achieves Top-1 100%, Hierarchical 100% |

### 2.2 Tabula Sapiens — Skin + Blood Subsets (The Tabula Sapiens Consortium, *Science* 2022)

| Field | Value |
|-------|-------|
| **Full atlas** | ~500K cells from 24 tissues, 15 donors |
| **Benchmark subset** | Skin + Blood compartments, ≤300 cells per cell type |
| **Cell types** | ~20 per tissue (see §3.2) |
| **Ground truth** | Expert-curated atlas labels |
| **Access** | CELLxGENE Census API (`cellxgene-census` package) — downloads targeted subset only; full atlas is NOT downloaded |
| **DOI** | 10.1126/science.abl4896 |
| **Status** | NOT STARTED |

---

## 3. Ground Truth Definitions

Ground truth dictionaries are defined in `benchmarking/ground_truth/`.

### 3.1 PBMC3k Ground Truth (8 cell types)

Cluster labels from the Seurat PBMC3k tutorial (Hao et al. 2021), with their
canonical Cell Ontology identifiers:

| Cluster Label | Canonical Name | CL ID | Broad Lineage |
|---------------|----------------|--------|---------------|
| CD4 T cells | CD4-positive, alpha-beta T cell | CL:0000624 | T cell |
| CD14+ Monocytes | classical monocyte | CL:0000860 | Monocyte |
| B cells | B cell | CL:0000236 | B cell |
| CD8 T cells | CD8-positive, alpha-beta T cell | CL:0000625 | T cell |
| NK cells | natural killer cell | CL:0000623 | NK cell |
| FCGR3A+ Monocytes | non-classical monocyte | CL:0000875 | Monocyte |
| Dendritic cells | myeloid dendritic cell | CL:0000451 | Dendritic cell |
| Megakaryocytes | megakaryocyte | CL:0000556 | Megakaryocyte |

### 3.2 Tabula Sapiens Ground Truth

Selected cell types from the Tabula Sapiens atlas (blood + skin):

**Blood compartment:** CD4+ T cell, CD8+ T cell, B cell, classical monocyte,
non-classical monocyte, NK cell, plasmablast, platelet, erythrocyte, neutrophil.

**Skin compartment:** keratinocyte, dermal fibroblast, endothelial cell,
CD4+ T cell, CD8+ T cell, macrophage, mast cell, melanocyte, pericyte, Schwann cell.

Full dictionaries with CL IDs are in `benchmarking/ground_truth/tabula_sapiens_gt.py`.

### 3.3 Matching Strategy

Top-1 accuracy is evaluated with three layers of matching (applied in order):

1. **Exact CL ID match** — predicted CL ID equals ground truth CL ID
2. **CL ancestor match** — predicted CL ID is a descendant of the ground truth CL ID
   (more specific annotation → still counted as correct)
3. **Normalized name match** — case-insensitive substring match between the predicted
   cell type name and the canonical name or any accepted alias

A prediction is marked **Correct** if any layer matches.

---

## 4. Evaluation Metrics

| Metric | Definition |
|--------|-----------|
| **Top-1 accuracy** | Fraction of clusters where the top-1 prediction matches ground truth via 3-layer matching |
| **Broad-type accuracy** | Fraction where the predicted cell type belongs to the correct broad lineage (T cell, B cell, etc.) |
| **Annotation rate** | Fraction of clusters that received a non-empty prediction (vs. "Unknown" / no call) |
| **Hierarchical accuracy** | SapiensOntoCellMap only: fraction correct at any level of the CL hierarchy tree |

Metrics are computed by `benchmarking/metrics/benchmark_metrics.py` and saved to
`benchmarking/results/<dataset>/metrics/accuracy_summary.csv`.

---

## 5. Module Architecture

```
benchmarking/
├── download/
│   ├── datasets/
│   │   ├── base_dataset_downloader.py   — abstract base (idempotent download())
│   │   ├── pbmc3k_downloader.py         — PBMC3k via scanpy.datasets
│   │   └── tabula_sapiens_downloader.py — targeted subset via CELLxGENE Census API
│   └── tools/
│       └── tool_installer.py            — availability checker for all comparator tools
├── annotation_tool_exec/
│   ├── base_runner.py                   — abstract base: run(deg_csv, output_dir) → dict
│   ├── sapiensonto_runner.py            — wraps get_cluster_annotation.py CLI
│   ├── celltypist_runner.py             — Python: celltypist.annotate()
│   ├── sctype_runner.py                 — R: Rscript run_sctype.R
│   └── singler_runner.py               — R: Rscript run_singler.R
├── metrics/
│   └── benchmark_metrics.py            — BenchmarkMetrics: 3-layer accuracy evaluation
├── figures/
│   └── benchmark_figures.py            — BenchmarkFigures: publication-quality matplotlib
├── ground_truth/
│   ├── pbmc3k_gt.py                     — PBMC3K_GROUND_TRUTH + GT_ALIASES + CL_LINEAGE_MAP
│   └── tabula_sapiens_gt.py             — TABULA_SAPIENS_GROUND_TRUTH (blood + skin)
├── r_scripts/
│   ├── run_sctype.R                     — scType wrapper (reads DEG CSV, writes predictions CSV)
│   └── run_singler.R                    — SingleR wrapper (reads .h5ad, writes predictions CSV)
├── results/
│   ├── pbmc3k/                          — per-tool outputs, metrics/, figures/
│   └── tabula_sapiens/                  — per-tool outputs, metrics/, figures/
├── benchmark_pbmc3k.py                  — PBMC3k benchmark runner (uses modules above)
├── benchmark_tabula_sapiens.py          — Tabula Sapiens benchmark runner
└── run_all_benchmarks.py               — combined CLI + publication figure
```

**R script interface** — all R tools use a CSV-based boundary:
- Input: DEG CSV or .h5ad path (passed as CLI arguments to Rscript)
- Output: `cluster,predicted_cell_type` CSV written to the tool's output subdirectory
- No shared state between Python and R processes

---

## 6. Running the Benchmarks

```bash
# Verify all tools are installed
python benchmarking/run_all_benchmarks.py --check_tools

# PBMC3k only (all P0 tools)
python benchmarking/benchmark_pbmc3k.py

# PBMC3k, Python tools only (no R required)
python benchmarking/benchmark_pbmc3k.py --no_r_tools

# Tabula Sapiens (requires cellxgene-census)
python benchmarking/benchmark_tabula_sapiens.py

# All datasets + combined publication figure
python benchmarking/run_all_benchmarks.py

# Skip download if data already cached
python benchmarking/run_all_benchmarks.py --skip_download
```

---

## 7. Expected Outputs

Each dataset benchmark produces:

```
benchmarking/results/<dataset>/
├── degs.csv                      — Wilcoxon DEG CSV (Scanpy format)
├── <dataset>_processed.h5ad      — preprocessed AnnData (gitignored — large file)
├── sapiensonto/                  — SapiensOntoCellMap outputs (HTML report, CSVs)
├── celltypist/
│   └── celltypist_predictions.csv
├── sctype/
│   └── sctype_predictions.csv
├── singler/
│   └── singler_predictions.csv
├── metrics/
│   ├── accuracy_summary.csv      — per-tool Top-1, broad-type, annotation rate
│   └── per_cluster_results.csv   — per-cluster detail (all tools)
└── figures/
    ├── <dataset>_top1_accuracy.png
    ├── <dataset>_top1_accuracy.pdf
    └── combined_benchmark_figure.png  (run_all_benchmarks.py only)
```

---

## 8. Publication Figure Plan

### Figure layout (manuscript Figure 4)

**Panel A — PBMC3k accuracy**
Horizontal bar chart: per-tool Top-1 accuracy. SapiensOntoCellMap bar in blue (#2563EB),
comparators in grey (#94A3B8). Bar labels show exact percentages.

**Panel B — Tabula Sapiens accuracy**
Same layout as Panel A for the skin + blood multi-tissue benchmark.

**Panel C — Unique capabilities table**

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
| 1 | PBMC3k: SapiensOntoCellMap + CellTypist | DONE |
| 2 | OOP module refactor: `annotation_tool_exec/`, `metrics/`, `figures/`, `download/` | DONE |
| 2 | Ground truth dicts: PBMC3k, Tabula Sapiens | DONE |
| 2 | R scripts: `run_sctype.R`, `run_singler.R` | DONE |
| 3 | Tabula Sapiens: download + DEG + all P0 tools | NOT STARTED |
| 3 | Combined publication figure | NOT STARTED |
