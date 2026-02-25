# SapiensOntoCellMap

**Ontology-Aware Cell Type Annotation and Deconvolution for Single-Cell and Spatial Transcriptomics**

---

## Overview

SapiensOntoCellMap is a Python toolkit for cell type annotation of scRNA-seq and spatial transcriptomics data. It integrates **14 curated marker databases** into a single ontology-normalized resource and uses it for (1) statistically rigorous cell type enrichment, (2) hierarchical annotation through the Cell Ontology graph, and (3) reference-free pseudobulk deconvolution of spatial data.

**Key differentiator:** Annotation is performed at two resolutions — fine-grained per-database (Level 1) and CL-ontology-normalized evidence-weighted (Level 2) — with hierarchical traversal of the Cell Ontology to report confidence at every level of resolution from broad cell lineage down to specific cell subtypes.

---

## Features

### Database
- **14 integrated sources:** CellMarker 2.0, PanglaoDB, CellxGene, HuBMAP HRA, WIMMS, Human Skin SCC (Ji et al. Cell 2020), Skin Fibroblast Atlas, Epidermal Cluster Atlas (GSE147482), and 4 curated skin scRNA-seq databases (GSE130973/163973/138669/156326) + full-thickness skin atlas (MERFISH + scRNAseq)
- **442,000 marker–cell type associations**, 22,469 unique genes, 1,208 cell types, 266 tissues
- UBERON tissue IDs and Cell Ontology (CL) cell type IDs throughout
- 6-tier evidence weighting: Experiment (1.0) > Single-Cell Sequencing (0.9) > Company (0.8) > Literature (0.7) > Review (0.6) > Computational (0.5)
- Interactive self-contained database visualizer (`sapiens_visualizer.html`, ~20 MB)

### Annotation Engine
- Hypergeometric enrichment test with global Benjamini–Hochberg FDR correction across all (cluster × cell type) pairs
- Evidence-weighted scoring (Weighted Enrichment, Weighted Recall, Combined Score)
- Hierarchical annotation via CL ontology graph traversal — reports confidence at every ancestor node from root to the significant hit
- HGNC gene alias resolution (maps historical and synonym gene symbols to current approved names)
- Auto-detection of DEG format: Seurat, Scanpy, or generic CSV (by column name)
- Two contexts: tissue-specific and all-tissue, with priority-ranked top annotation summary

### Cell Type Composition
- Annotation-derived composition scores: per-cluster `Combined_Score` (= Weighted_Enrichment × −log₁₀(adj_p_value)) normalised to sum to 1.0 across all significantly enriched Level 2 cell types
- Works for both scRNA-seq and spatial data — no expression matrix required
- Uses the marker database as intended (via enrichment testing), avoiding L1-norm bias inherent to NNLS on sparse marker lists

### Outputs
- Interactive HTML report (6 tabs): Cell Type Summary, DEG Browser, Enrichment Visuals, Hypergeometric Results, Hierarchy (icicle chart), Composition
- Top annotation summary CSV with cell type, confidence, broad type, N databases, overlapping genes
- Per-context/level all-results and significant-results CSVs
- Hierarchical annotation CSV

---

## Installation

### Requirements

- Python ≥ 3.10
- Recommended: [uv](https://docs.astral.sh/uv/) for fast environment management

### Quick install

```bash
git clone https://github.com/Sonal1510/SapiensOntoCellMap.git
cd SapiensOntoCellMap

# Create environment and install (uv recommended)
uv venv
uv pip install -e .

# Or with standard pip
python3 -m venv SapiensOntoCellMap_env
source SapiensOntoCellMap_env/bin/activate
pip install -e .
```

For HDF5 support (required for Xenium expression matrices):
```bash
pip install h5py
```

### Build the marker database

Download and parse all 14 source databases (~10–20 minutes, internet required):

```bash
python3 test/test_classes.py
```

This creates `data/processed_combined_db/master_cell_marker_db.csv` and the HGNC alias map at `data/reference/hgnc_complete_set.txt`.

---

## Usage

### scRNA-seq annotation

```bash
python3 src/cluster_annotation/get_cluster_annotation.py \
    /path/to/cluster_markers.csv  SAMPLE_NAME  /path/to/output_dir \
    --deg_type scrna \
    --deg_format seurat \
    --marker_db data/processed_combined_db/master_cell_marker_db.csv \
    --tissue skin \
    --hgnc_map data/reference/hgnc_complete_set.txt
```

**Supported DEG formats** (auto-detected or forced with `--deg_format`):
- `seurat` — `FindAllMarkers` output (`p_val_adj`, `avg_log2FC`, `cluster`, `gene`)
- `scanpy` — `rank_genes_groups` export (`pvals_adj`, `logfoldchanges`, `group`, `names`)
- `generic` — any CSV with gene symbols; p-values/FC columns auto-matched by name

### Spatial annotation

```bash
python3 src/cluster_annotation/get_cluster_annotation.py \
    /path/to/spaceranger_outs/  SAMPLE_NAME  /path/to/output_dir \
    --deg_type spatial \
    --marker_db data/processed_combined_db/master_cell_marker_db.csv \
    --tissue skin \
    --hgnc_map data/reference/hgnc_complete_set.txt
```

The UMAP and cluster assignments are auto-detected from the Space Ranger output directory.

### CLI reference

| Argument | Default | Description |
|----------|---------|-------------|
| `--deg_type` | required | `scrna` or `spatial` |
| `--deg_format` | auto | `seurat`, `scanpy`, or `generic` |
| `--marker_db` | required | Path to `master_cell_marker_db.csv` |
| `--tissue` | None | Tissue name substring for marker DB filter (e.g. `skin`) |
| `--pval` | 0.05 | Adjusted p-value threshold |
| `--log2fc` | 1.0 | log2 fold-change threshold |
| `--min_overlap` | 2 | Minimum gene overlap to report |
| `--background_gene_count` | auto | Override hypergeometric N (see §N for scRNA-seq) |
| `--hgnc_map` | None | HGNC complete set file for gene alias resolution |
| `--no_hierarchy` | — | Skip CL ontology hierarchical annotation |
| `--umap_csv` | None | scRNA-seq UMAP coordinates CSV (Barcode, UMAP-1, UMAP-2) |
| `--cell_cluster_csv` | None | scRNA-seq cluster assignments CSV (Barcode, Cluster) |
| `--no_deconvolution` | — | Disable the Cell Type Composition tab in the HTML report |
| `--no_auto_spatial_filter` | — | Disable auto-calibration of mean_counts threshold |

---

## Output Report Tabs

| Tab | Contents |
|-----|----------|
| **Cell Type Summary** | Top annotation per cluster table, UMAP (spatial auto-detected), marker gene heatmap |
| **DEG Browser** | Interactive per-cluster DEG table with p-value/FC/mean-counts filters |
| **Enrichment Visuals** | Clustered heatmap (adj p < 0.05), p-value violin, log2FC violin, mean-counts box |
| **Hypergeometric Result** | Full and significant results tables (Level 1 and Level 2) |
| **Hierarchy** | Icicle chart of CL ontology traversal with confidence scores |
| **Composition** | Annotation-derived cell type composition scores per cluster (stacked bar + table) |

---

## Output Columns (Level 2 Results)

| Column | Definition |
|--------|-----------|
| `adj_p_value` | Benjamini–Hochberg corrected hypergeometric p-value (global FDR) |
| `Enrichment_ratio` | (k/n) / (K/N) — observed vs. expected overlap fraction |
| `Weighted_Recall` | W_overlap / W_ref — fraction of marker evidence weight captured |
| `Weighted_Enrichment` | (W_overlap/n) / (W_ref/N) — weighted enrichment ratio |
| `Combined_Score` | Weighted_Enrichment × −log10(adj_p_value) — ranking metric |
| `N_Databases` | Number of independent databases contributing at least one overlapping marker gene |

See [`docs/STATISTICAL_METHODS.md`](docs/STATISTICAL_METHODS.md) for complete statistical documentation.

---

## Project Structure

```
SapiensOntoCellMap/
├── config/
│   └── config.py                    # DATABASE_CONFIG, paths, constants
├── data/
│   ├── raw/                         # Auto-downloaded source files
│   ├── reference/                   # HGNC alias map
│   ├── recovered_ids_dfs/           # Ontology normalization QC logs
│   └── processed_combined_db/       # master_cell_marker_db.csv + visualizer HTML
├── docs/
│   └── STATISTICAL_METHODS.md       # Complete statistical reference
├── src/
│   ├── parser/                      # Per-database parsers + ontology utilities
│   ├── db_manager/                  # DatabaseCreate orchestrator, DatabaseValidator
│   ├── download/                    # BioDataDownloader
│   ├── cluster_annotation/          # Annotation engine + HTML report
│   │   ├── get_cluster_annotation.py    # CLI entry point
│   │   ├── get_marker_enrichment_test.py # Hypergeometric enrichment engine
│   │   ├── get_html_report.py           # Jinja2 HTML report (6 tabs)
│   │   └── hierarchical_annotation.py   # CL ontology traversal
│   ├── deconvolution/               # (Legacy: NNLS solver, unused)
│   │   ├── signature_builder.py         # Evidence-weighted signature matrix builder
│   │   └── nnls_deconvolver.py          # NNLS solver (not used in current pipeline)
│   └── visualization/               # Interactive database explorer
│       └── db_process.py                # CSV → JSON → sapiens_visualizer.html
├── test/
│   └── test_classes.py              # Full pipeline integration test (downloads + builds DB)
├── pyproject.toml                   # Package metadata and dependencies
└── requirements.txt                 # Pinned environment snapshot
```

---

## Adding a New Database

1. Add an entry to `DATABASE_CONFIG` in `config/config.py`
2. Use `GenericFileParser` (`parser_key: "generic"`) when possible
3. Run `python3 test/test_classes.py` to rebuild
4. Check `data/processed_combined_db/quarantine_log.csv` for any rejected rows

---

## Tested Platforms

| Data type | Format | Expression matrix |
|-----------|--------|-------------------|
| scRNA-seq (10x Flex) | Seurat `FindAllMarkers` CSV | — |
| Visium HD | Space Ranger `outs/` | MEX (`filtered_feature_bc_matrix/`) |
| Xenium | Xenium Ranger `outs/` | HDF5 (`cell_feature_matrix.h5`) |

---

## Benchmarking

### PBMC3k (Zheng et al. 2017, n=2,638 cells, 8 cell types)

Standard scRNA-seq annotation benchmark comparing SapiensOntoCellMap vs. CellTypist (`Immune_All_High` model):

| Tool | Top-1 Accuracy | Hierarchical Accuracy |
|------|---------------|----------------------|
| **SapiensOntoCellMap** | **100%** | **100%** |
| CellTypist (Immune_All_High) | 62.5% | — |

Per-cluster annotations (SapiensOntoCellMap):

| Ground Truth | Predicted Cell Type | Match |
|-------------|--------------------|----|
| CD4 T cells | naive thymus-derived CD4-positive, alpha-beta T cell | ✅ |
| CD14+ Monocytes | CD14-positive monocyte | ✅ |
| B cells | B-2 B cell | ✅ |
| CD8 T cells | effector CD8-positive, alpha-beta T cell | ✅ |
| NK cells | CD16-positive, CD56-dim natural killer cell, human | ✅ |
| FCGR3A+ Monocytes | non-classical monocyte | ✅ |
| Dendritic cells | myeloid dendritic cell, human | ✅ |
| Megakaryocytes | platelet | ✅ |

Reproduce:
```bash
pip install scanpy celltypist
python3 benchmarking/benchmark_pbmc3k.py
```

---

## Citation

Manuscript in preparation. Please contact the authors before citing.

---

## License

MIT License. See `LICENSE` for details.
