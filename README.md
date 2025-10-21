# 🧬 SapiensOntoCellMap
### *An Ontology-Integrated Framework for Curated Human Cell Marker Databases and Cell Type Enrichment Analysis*

---

## Overview

**SapiensOntoCellMap** is a Python-based toolkit designed to address two major challenges in single-cell biology:

1. **Database Chaos** – It automatically integrates over 8 public and manually curated cell marker databases into a single, standardized, and ontology-aware resource.  
2. **Cell Annotation** – It provides a powerful enrichment tool to accurately annotate cell clusters from scRNA-seq and spatial transcriptomics data using this integrated database.

The entire pipeline is **configurable and extensible**, allowing new databases to be added with minimal effort.

---

## 🚀 Key Features

- **Automated Database Integration**  
  Downloads and parses 8+ databases (including *CellMarkerDB*, *PanglaoDB*, *CellxGene*, *HuBMAP*) into a unified SQLite/CSV database.

- **Ontology-Driven Normalization**  
  Leverages the *Cell Ontology (CL)* and *Uberon* to standardize cell and tissue names, using fuzzy-matching to resolve ambiguous terms.

- **Powerful Enrichment Analysis**  
  Implements a **hypergeometric test** to statistically determine the most likely cell types for a given list of marker genes.

- **Publication-Quality Outputs**  
  Generates a comprehensive, interactive **HTML report** with sortable tables, heatmaps, and DEG (Differentially Expressed Gene) browser plots.

- **Extensible Framework**  
  New databases can be added easily by updating a single configuration file.

---

## ⚙️ Quick Start: Installation & Usage

Get up and running in **3 simple steps**.

### 1. Installation

Clone the repository and install dependencies:

```bash
# Clone the repository
git clone https://github.com/Sonal1510/SapiensOntoCellMap.git
cd SapiensOntoCellMap

# Install dependencies
pip install -r requirements.txt
```

## 2. Build the Database

Run the test script to execute the full data integration pipeline. This will download all source files and build the final master cell marker database.

```bash
# This script downloads data and runs the parsing pipeline
python3 test/test_classes.py
```

This command will populate the `data/` directory and create the final database at:

```text
data/processed_combined_db/master_cell_marker_db.csv
```

## 3. Annotate Your Data (Example)

Once the database is built, use the annotation tool to identify cell types for your clusters.

```bash
python3 src/cluster_annotation/get_cluster_annotation.py \
    --deg_file_path /path/to/your/deg_file.csv \
    --marker_db_path data/processed_combined_db/master_cell_marker_db.csv \
    --output_dir /path/to/your/results
```

**Notes:**

- **Input DEG file formats supported:** CSVs exported from Seurat, Scanpy, or similar tools (ensure a column with gene symbols).  
- **Output:** An interactive HTML report (sortable tables, heatmaps, DEG browser plots) saved in `--output_dir`.

---

## 🔄 How It Works: The Data Pipeline

SapiensOntoCellMap operates in two main stages:

### 1. Database Integration Pipeline

The pipeline processes multiple public and curated data sources into a final standardized marker database.

**Workflow:**

- **Data Acquisition** — Downloads source files as defined in `config/config.py`.  
- **Parsing & Standardization** — Each source is read by a dedicated parser in `src/parser/`. The central `BaseParser` normalizes cell/tissue names to official ontology IDs.  
- **Logging & QC** — All normalization attempts (exact, fuzzy, or failed) are logged in `data/recovered_ids_dfs/` for transparency and manual curation.  
- **Assembly** — Standardized data are combined into a single master CSV file.

**Data Directory Structure:**

```text
data/
├── raw/                            # Original source data (auto-downloaded)
├── pre_manually_downloaded_files/  # Optional manually added sources
├── recovered_ids_dfs/              # Normalization QC logs
└── processed_combined_db/          # Final master database output
```

## 2. Cell Type Annotation

The annotation tool identifies likely cell types based on DEG lists for each cluster.

**Input:**

- DEG file (gene lists per cluster)

**Enrichment Method:**

- For each cluster's gene list, the tool performs a hypergeometric test against marker sets of every cell type in the master database.

**Output:**

- Ranked list of enriched cell types with raw p-values, FDR-corrected p-values, and enrichment scores.  
- Interactive HTML report containing sortable tables, summary heatmaps, and DEG browser visualizations.

---

## 📊 Project Status (as of Oct 2025)

| Phase   | Description                                                      | Status         |
|---------|------------------------------------------------------------------|:--------------:|
| Phase 1 | Data Integration — downloads, normalizes, and combines sources   | ✅ Complete     |
| Phase 2 | Annotation Tool — enrichment & HTML report generation            | ✅ Complete     |

---

## 📁 Project Structure

```text
SapiensOntoCellMap/
├── README.md
├── config/
│   └── config.py                 # Central configuration for all databases
├── data/
│   ├── raw/                      # Auto-downloaded source files
│   ├── pre_manually_downloaded_files/
│   ├── recovered_ids_dfs/        # Normalization QC logs
│   └── processed_combined_db/    # Final master database
├── requirements.txt
├── src/
│   ├── download/                 # Scripts for data download
│   ├── parser/                   # Scripts for parsing and normalization
│   └── cluster_annotation/       # Enrichment analysis and reporting
└── test/
    └── test_classes.py           # Main script to run the full pipeline
```

## 🛠️ Configuration

- `config/config.py` contains source URLs, local filenames, parser mappings, and any database-specific options.  
- To add a new database: add an entry to `config/config.py` and implement a parser in `src/parser/` (or reuse `BaseParser` where possible).

---

## ✅ Quality Control & Logging

- All normalization attempts are saved under `data/recovered_ids_dfs/` with details on whether the mapping to Cell Ontology / Uberon was **exact**, **fuzzy**, or **failed**.  
- These logs are intended for manual review to improve mappings and extend the fuzzy-matching rules.

---

## 🧩 Extensibility

- Parsers follow a common interface (see `src/parser/base_parser.py`) so new sources can be integrated with minimal code.  
- The final master table (`master_cell_marker_db.csv`) is a flattened CSV that includes standardized ontology IDs and source provenance for every marker.

---

## 🔬 Usage Tips

- Use gene symbols consistent with HGNC-approved symbols for best matching.  
- If your DEG lists contain multiple identifier types (ENSEMBL, RefSeq), convert to gene symbols first.  
- Review `data/recovered_ids_dfs/` after building the DB to inspect any failed mappings that may affect annotation.
