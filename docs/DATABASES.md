# SapiensOntoCellMap: Database & Reference File Registry

**Version:** 0.1.0 | **Date:** 2026-02-25 | **Author:** Sonal Rashmi

This document is the authoritative reference for every external data source
ingested by SapiensOntoCellMap — marker databases, ontology files, and
reference gene sets. It also provides a step-by-step guide for adding new
databases, which is the primary extensibility mechanism of the tool.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Marker Databases (14 sources)](#2-marker-databases-14-sources)
3. [Reference Files](#3-reference-files)
4. [How to Add a New Database](#4-how-to-add-a-new-database)
5. [Source Type Vocabulary](#5-source-type-vocabulary)
6. [Keeping the Database Current](#6-keeping-the-database-current)

---

## 1. Architecture Overview

```
config/config.py
  └── DATABASE_CONFIG          ← single source of truth for all 14 databases
  └── MSIGDB_*, HGNC_*         ← reference file constants

src/download/bio_database_downloader.py
  └── BioDataDownloader
        ├── download_all_databases()      ← fetches marker DB raw files
        └── download_reference_data()     ← fetches HGNC + MSigDB gene sets

src/parser/
  ├── generic_parser.py        ← handles most new databases (recommended)
  ├── base_parser.py           ← ontology normalisation (UBERON + CL IDs)
  └── <db>_parser.py           ← dedicated parsers for complex formats

src/db_manager/
  ├── database_creator.py      ← orchestrates all parsers → master_cell_marker_db.csv
  └── database_validator.py    ← enforces 14-column schema; quarantines bad rows
```

The master output is `data/processed_combined_db/master_cell_marker_db.csv`
(14 columns; see schema in `src/db_manager/database_validator.py`).

---

## 2. Marker Databases (14 sources)

### 2.1 CellMarkerDB
| Field | Value |
|-------|-------|
| **Key** | `cellmarkerdb` |
| **Parser** | `CellMarkerDBParser` (`parser_key: "cellmarker"`) |
| **Source type** | Literature / Experiment (mixed; per-row in source) |
| **Tissue scope** | Pan-human |
| **Download** | Automatic — `http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx` |
| **File** | `data/raw/Cell_marker_Human.xlsx` |
| **Citation** | Zhang et al., *Nucleic Acids Res* 2023 (doi:10.1093/nar/gkac900) |
| **Notes** | Curated from >100,000 publications. Broadest tissue coverage of any single source. Updated annually. |

### 2.2 HuBMAP (Human BioMolecular Atlas Program)
| Field | Value |
|-------|-------|
| **Key** | `hubmap` |
| **Parser** | `HuBMapDBParser` (`parser_key: "hubmap"`) |
| **Source type** | Experiment |
| **Tissue scope** | Pan-human (27 organs) |
| **Download** | Manual — see note below |
| **File** | `data/pre_manually_downloaded_files/hubmap_summary.csv` |
| **Citation** | HuBMAP Consortium, *Nature* 2019 (doi:10.1038/s41586-019-1629-x) |
| **Notes** | Download the latest ASCT+B table from https://humanatlas.io/asctb-tables, select each organ, and use the provided parsing script to produce `hubmap_summary.csv`. The URL includes a session ID that changes on each visit, making fully automated download impractical. |

### 2.3 CellxGene (CZ CELLxGENE Census)
| Field | Value |
|-------|-------|
| **Key** | `cellxgene` |
| **Parser** | `CellxGeneDBParser` (`parser_key: "cellxgene"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Pan-human |
| **Download** | Automatic — `https://cellguide.cellxgene.cziscience.com/.../marker_gene_data.json.gz` |
| **File** | `data/raw/marker_gene_data.json.gz` |
| **Citation** | CZI Single-Cell Biology Program (doi:10.1101/2023.10.30.563174) |
| **Notes** | Computational marker genes derived from the CELLxGENE Census (~60M cells). Updated quarterly. Re-run the downloader to get the latest snapshot. |

### 2.4 PanglaoDB
| Field | Value |
|-------|-------|
| **Key** | `panglaodb` |
| **Parser** | `PanglaoParser` (`parser_key: "panglao"`) |
| **Source type** | Computational |
| **Tissue scope** | Pan-human (mouse genes excluded at parse time) |
| **Download** | Automatic — `https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz` |
| **File** | `data/raw/PanglaoDB_markers_27_Mar_2020.tsv.gz` |
| **Citation** | Franzen et al., *Database* 2019 (doi:10.1093/database/baz046) |
| **Notes** | **Frozen at March 2020** — no updates since publication. Gene content reflects the 2020 HGNC namespace. Alias resolution via HGNC partially compensates for gene symbol drift. Consider supplementing with newer sources (see Section 6). |

### 2.5 WIMMS (Wound & Inflammatory Melanocyte Markers)
| Field | Value |
|-------|-------|
| **Key** | `wimms` |
| **Parser** | `WimmsMelanocyteParser` (`parser_key: "wimms"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (melanocytes, wound healing) |
| **Download** | Manual — https://wimms.tanlab.org/ → Signature and Datasets → Download All |
| **File** | `data/pre_manually_downloaded_files/wimms_signature.csv` |
| **Notes** | The download URL contains a session ID that changes each visit. Save the CSV manually and place it in `pre_manually_downloaded_files/`. |

### 2.6 Human SCC 2020 (Ji et al., *Cell* 2020)
| Field | Value |
|-------|-------|
| **Key** | `human_scc_cell_2020` |
| **Parser** | `HumanSccCell2020Parser` (`parser_key: "human_scc_cell_2020"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (squamous cell carcinoma, tumour microenvironment) |
| **Download** | Manual — Supplementary Table S2 from https://www.cell.com/cell/fulltext/S0092-8674(20)30672-3 |
| **File** | `data/pre_manually_downloaded_files/mmc2.xlsx` |
| **Citation** | Ji et al., *Cell* 2020 (doi:10.1016/j.cell.2020.06.060) |
| **Notes** | scRNA-seq of 10 SCC tumours. Resolves tumour-specific keratinocyte, immune, and stromal states at single-cell resolution. |

### 2.7 Epithelial Clusters GSE147482 (Vu et al., *Nat Commun* 2020)
| Field | Value |
|-------|-------|
| **Key** | `epi_cluster_gse147482` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Computational |
| **Tissue scope** | Skin (epithelial compartment) |
| **Download** | Automatic — Springer supplementary table |
| **File** | `data/raw/Epi_Cluster_gse147482.xlsx` |
| **Citation** | Vu et al., *Nat Commun* 2020 (doi:10.1038/s41467-020-18075-7) |
| **parser_config** | `gene_col: "Row"`, `base_cell_name: "Epithelial"`, `cell_subtype_col: "Cluster"` |

### 2.8 Skin Fibroblast Atlas
| Field | Value |
|-------|-------|
| **Key** | `skin_fibroblast_atlas` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Computational |
| **Tissue scope** | Skin (fibroblasts) |
| **Download** | Automatic — GitHub raw content |
| **File** | `data/raw/Skin_fibroblast_atlas.csv` |
| **Citation** | Solé-Boldo et al., *Nat Commun* 2020 (doi:10.1038/s41467-020-15900-x) |
| **parser_config** | `gene_col: "Gene"`, `base_cell_name: "Fibroblast"`, `cell_subtype_col: "Group"` |

### 2.9 ScarCellMarker GSE130973 (Ji et al., *Cell* 2020)
| Field | Value |
|-------|-------|
| **Key** | `scar_cell_marker_gse130973` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (primary SCC, 11 patients) |
| **Download** | Manual — http://124.220.48.30:3838/ScarCellMarker/#tab-8377-7 |
| **File** | `data/pre_manually_downloaded_files/GSE130973_2020_sc_skin_Cell_anno_V2_marker_genes.csv` |
| **Citation** | Ji et al., *Cell* 2020 (doi:10.1016/j.cell.2020.06.060) |
| **parser_config** | `gene_col: "gene"`, `cell_name_col: "cluster"`, Seurat FindAllMarkers format |

### 2.10 ScarCellMarker GSE163973 (Chen et al., *Nat Commun* 2021)
| Field | Value |
|-------|-------|
| **Key** | `scar_cell_marker_gse163973` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (keloid, 16 patients) |
| **Download** | Manual — ScarCellMarker portal above |
| **File** | `data/pre_manually_downloaded_files/GSE163973_2021_NC_sc_keloid_Cell_anno_V2_marker_genes.csv` |
| **Citation** | Chen et al., *Nat Commun* 2021 (doi:10.1038/s41467-021-26849-2) |

### 2.11 ScarCellMarker GSE138669 (Theocharidis et al., *JID* 2022)
| Field | Value |
|-------|-------|
| **Key** | `scar_cell_marker_gse138669` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (wound healing, fibroblasts + keratinocytes) |
| **Download** | Manual — ScarCellMarker portal above |
| **File** | `data/pre_manually_downloaded_files/GSE138669_2018_JID_sc_skin_Cell_anno_V2_marker_genes.csv` |
| **Citation** | Theocharidis et al., *JID* 2022 (doi:10.1016/j.jid.2021.07.178) |

### 2.12 ScarCellMarker GSE156326 (Zhu et al., *Nat Commun* 2021)
| Field | Value |
|-------|-------|
| **Key** | `scar_cell_marker_gse156326` |
| **Parser** | `GenericFileParser` (`parser_key: "generic"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (hypertrophic scar, myofibroblasts + macrophages, 6 patients) |
| **Download** | Manual — ScarCellMarker portal above |
| **File** | `data/pre_manually_downloaded_files/GSE156326_2021_NC_sc_hyper_Cell_anno_V2_marker_genes.csv` |
| **Citation** | Zhu et al., *Nat Commun* 2021 (doi:10.1038/s41467-021-26386-y) |

### 2.13 Skin Atlas MERFISH (Joost et al., *Cell Systems* 2025?)
| Field | Value |
|-------|-------|
| **Key** | `skin_atlas_MERFISH` |
| **Parser** | `SkinAtlasParser` (`parser_key: "skin_atlas"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (spatial MERFISH panel, full-thickness) |
| **Download** | Manual — Supplementary Table S3 from the associated publication |
| **File** | `data/pre_manually_downloaded_files/Table_S3_MERFISHClusterMarkers.xlsx` |
| **Notes** | Multi-sheet Excel; each sheet is one cell-type cluster. Parsed by `SkinAtlasParser`. |

### 2.14 Skin Atlas scRNA-seq (Joost et al., same study)
| Field | Value |
|-------|-------|
| **Key** | `skin_atlas_scRNAseq` |
| **Parser** | `SkinAtlasParser` (`parser_key: "skin_atlas"`) |
| **Source type** | Single-Cell Sequencing |
| **Tissue scope** | Skin (scRNA-seq companion to MERFISH atlas) |
| **Download** | Manual — Supplementary Table S4 |
| **File** | `data/pre_manually_downloaded_files/Table_S4_scRNAseqClusterMarkers.xlsx` |

---

## 3. Reference Files

These files are downloaded by `BioDataDownloader.download_reference_data()` and
stored in `data/reference/`. They are not marker databases; they support the
enrichment pipeline.

### 3.1 HGNC Complete Set (gene alias resolution)
| Field | Value |
|-------|-------|
| **Constant** | `HGNC_COMPLETE_SET_FILE` |
| **File** | `data/reference/hgnc_complete_set.txt` |
| **Download** | Automatic — Google Cloud Storage (HGNC public bucket) |
| **Format** | Tab-separated; columns include `symbol`, `alias_symbol`, `prev_symbol` |
| **Purpose** | Resolves gene name aliases and deprecated symbols during enrichment. ~42,000 approved HGNC entries, ~58,000 alias mappings. |
| **Citation** | HGNC, European Bioinformatics Institute (https://www.genenames.org) |
| **Update cadence** | Monthly. Re-run `download_reference_data()` to refresh. |

### 3.2 MSigDB HALLMARK_G2M_CHECKPOINT
| Field | Value |
|-------|-------|
| **Constant** | `MSIGDB_G2M_FILE` |
| **File** | `data/reference/HALLMARK_G2M_CHECKPOINT.grp` |
| **Download** | Automatic — MSigDB public endpoint |
| **Format** | Plain text `.grp`; one HGNC gene symbol per line; comment lines start with `#` |
| **Purpose** | G2/M cell cycle genes used for `Proliferative_Flag` detection (200 genes). |
| **Citation** | Liberzon et al., *Cell Systems* 2015 (doi:10.1016/j.cels.2015.12.004); Subramanian et al., *PNAS* 2005 (doi:10.1073/pnas.0506580102) |

### 3.3 MSigDB HALLMARK_E2F_TARGETS
| Field | Value |
|-------|-------|
| **Constant** | `MSIGDB_E2F_FILE` |
| **File** | `data/reference/HALLMARK_E2F_TARGETS.grp` |
| **Download** | Automatic — MSigDB public endpoint |
| **Format** | Plain text `.grp`; same as G2M |
| **Purpose** | E2F transcription factor target genes — the canonical S-phase / cell-cycle entry gene set in the MSigDB Hallmark collection (200 genes). Used together with G2M for `Proliferative_Flag` detection. |
| **Citation** | Same as G2M above. |
| **Notes** | The MSigDB Hallmark collection does not contain a separate `HALLMARK_G1S_CHECKPOINT` gene set. `HALLMARK_E2F_TARGETS` is the established Hallmark equivalent covering G1/S cell-cycle genes. If either MSigDB file is absent, `Proliferative_Flag` is set to `None` (not computed). No hardcoded fallback gene list is used. Lines not matching HGNC gene symbol format are rejected with a warning (guards against HTML redirect responses from the MSigDB endpoint). |

---

## 4. How to Add a New Database

This is the primary extensibility mechanism of SapiensOntoCellMap. The tool was
designed from the outset to be database-agnostic — no code changes are needed
for databases that fit the generic format.

### Decision tree

```
New database
    │
    ├── Is it a simple TSV/CSV/Excel with one gene per row,
    │   a tissue column (or fixed tissue), and a cell type column?
    │   └── YES → use GenericFileParser (Case A — no new code)
    │
    └── NO → does it have a complex format (JSON, multi-sheet, nested)?
            └── Write a dedicated parser class (Case B — ~50 lines of code)
```

---

### Case A: Add via `GenericFileParser` (most common, no new code)

**Step 1 — Obtain the file**

Place it in `data/pre_manually_downloaded_files/` (manual) or use a direct URL
(automatic download).

**Step 2 — Add an entry to `DATABASE_CONFIG` in `config/config.py`**

```python
"my_new_database": {
    "source": [
        "https://example.com/path/to/markers.csv",  # URL (or "" for manual)
        "my_markers.csv",                            # filename in data/raw/
        "csv"                                        # "csv", "tsv", or "xlsx"
    ],
    "parser_key": "generic",
    "parser_config": {
        # Required:
        "database_name": "My New Database",          # label in master DB
        "gene_col": "Gene",                          # column with gene symbols

        # Tissue — provide exactly ONE of these two:
        "tissue_name": "Skin",                       # fixed tissue for all rows
        # "tissue_name_col": "Tissue",               # OR column that contains it

        # Cell type — provide exactly ONE of these (or combine):
        "cell_name_col": "CellType",                 # full cell name column
        # "base_cell_name": "Fibroblast",            # fixed prefix
        # "cell_subtype_col": "Subtype",             # appended to base_cell_name

        # Evidence tier (see Section 5):
        "source_type": "Single-Cell Sequencing",

        # Optional — columns to include in source_info:
        "source_info_cols": ["p_val_adj", "avg_log2FC"]
    }
}
```

**Step 3 — Rebuild the master database**

```bash
python3 scripts/build_marker_db.py
```

This downloads the file (if URL provided), parses it, runs
`DatabaseValidator`, and writes the updated `master_cell_marker_db.csv`.
Check `data/processed_combined_db/quarantine_log.csv` for rejected rows.

**Step 4 — Verify**

```bash
python3 -c "
import pandas as pd
df = pd.read_csv('data/processed_combined_db/master_cell_marker_db.csv')
print(df[df['database_name'] == 'My New Database'].head())
print('Rows:', len(df[df['database_name'] == 'My New Database']))
"
```

---

### Case B: Dedicated parser (complex format)

Required when the source has a non-tabular format (JSON, multi-sheet Excel,
nested structure) or needs custom cell-type name construction logic.

**Step 1 — Create `src/parser/mydb_parser.py`**

```python
from src.parser.base_parser import BaseParser

class MyDBParser(BaseParser):
    """
    Parser for My New Database.

    Input: DataFrame loaded by DatabaseFileParser from <file>.
    Output: Standardised DataFrame with columns expected by DatabaseValidator.
    """
    def __init__(self, df: pd.DataFrame, database_name: str):
        super().__init__()          # inherits ontology normalisation
        self.df = df
        self.database_name = database_name

    def parse(self) -> pd.DataFrame:
        # 1. Rename/extract relevant columns
        # 2. Set database_name, source_type
        # 3. Call self.normalise_tissue(name) → UBERON ID  (from BaseParser)
        #    Call self.normalise_cell(name)   → CL ID      (from BaseParser)
        # 4. Return DataFrame with all 14 schema columns
        ...
```

**Step 2 — Register the parser key**

In `src/db_manager/database_creator.py`, add to `PARSER_MAPPING`:

```python
from src.parser.mydb_parser import MyDBParser

PARSER_MAPPING = {
    ...
    "mydb": MyDBParser,   # your new key
}
```

**Step 3 — Add to `DATABASE_CONFIG`**

```python
"my_new_database": {
    "source": ["https://...", "myfile.xlsx", "xlsx"],
    "parser_key": "mydb",
    "parser_config": {
        "database_name": "My New Database",
    }
}
```

**Step 4 — Rebuild and verify** (same as Case A, Step 3–4)

---

### Schema reference (14 columns)

All parsers must produce exactly these columns. `DatabaseValidator` rejects
rows that violate any constraint.

| Column | Constraint |
|--------|-----------|
| `database_name` | Non-null string |
| `tissue_name` | Non-null string |
| `tissue_id` | Must start with `UBERON:` |
| `cell_name` | Non-null string |
| `cell_id` | Must start with `CL:` |
| `db_cell_name` | Non-null (Level 1 lookup key) |
| `gene` | Non-null string |
| `source_type` | Must be in controlled vocabulary (Section 5) |
| `source_info` | String (p-values, logFC, etc.) |
| `species` | Default `"Human"` |
| `organ` | Derived from tissue |
| `organ_id` | UBERON ID of organ |
| `year` | Integer or null |
| `other_info` | Free text |

---

### Common mistakes

| Mistake | Fix |
|---------|-----|
| `tissue_id` fails validator | Check that `BaseParser.normalise_tissue()` found a match; inspect `data/recovered_ids_dfs/` |
| `cell_id` is `CL:0000000` (root) | Cell type name too generic; add more specific aliases to ontology_utils |
| Many rows in quarantine | Set `source_info_cols` only to numeric quality columns; avoid free-text columns that contain commas |
| New DB overrepresents one tissue | Expected — document bias in the publication methods section |

---

## 5. Source Type Vocabulary

The `source_type` field controls the evidence weight applied during enrichment
scoring. Use exactly these strings (case-sensitive):

| source_type | Weight | When to use |
|------------|--------|-------------|
| `Experiment` | 1.0 | Flow cytometry, IHC, FACS-sorted bulk RNA-seq with validated cell populations |
| `Single-Cell Sequencing` | 0.9 | Seurat/Scanpy FindMarkers on annotated scRNA-seq or spatial clusters |
| `Company` | 0.8 | Commercial antibody datasheets (e.g. BioLegend, R&D Systems) |
| `Literature` | 0.7 | Review articles or textbook statements about canonical markers |
| `Review` | 0.6 | Systematic reviews or meta-analyses without primary experimental data |
| `Computational` | 0.5 | Computational predictions, NLP extraction, or co-expression modules without experimental validation |

**Key distinction:** The `source_type` reflects the *biological evidence tier*,
not the computational method used to identify markers. Seurat FindAllMarkers
applied to a peer-reviewed scRNA-seq study is `Single-Cell Sequencing`, not
`Computational`.

---

## 6. Keeping the Database Current

The most important limitation of any static marker database is staleness.
SapiensOntoCellMap mitigates this through:

### Short-term (any time)
- **Re-download automatic sources:** `python3 src/download/bio_database_downloader.py`
  refreshes CellMarkerDB, CellxGene, PanglaoDB, and reference files.
- **Add a new GenericFileParser entry** for any recently published scRNA-seq atlas
  using the Case A workflow above (15–30 minutes per database).

### Medium-term (per publication cycle)
- **HuBMAP ASCT+B tables** are updated annually. Download the latest version and
  replace `hubmap_summary.csv`.
- **HGNC gene aliases** drift continuously. Re-run `download_reference_data()` to
  refresh `hgnc_complete_set.txt` before each major analysis run.

### Long-term (architectural)
- **PanglaoDB is frozen at 2020.** CellxGene Census (Section 2.3) already supplements
  it with a continuously updated pan-human source (~60M cells, updated quarterly).
  Additional pan-human atlases can be added via the Case A workflow (Section 4)
  as new resources are published.
