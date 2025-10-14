# config/config.py

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Centralized file for all URLs, paths, and parser configurations.
"""

import os
from src.parser.cellmarkerdb_parser import CellMarkerDBParser
from src.parser.hubmap_parser import HuBMapDBParser
from src.parser.cellxgene_parser import CellxGeneDBParser
from src.parser.panglaodb_parser import PanglaoParser
from src.parser.wimms_parser import WimmsMelanocyteParser
from src.parser.generic_parser import GenericFileParser

# --- Directory Paths --- (No changes here)
_CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DATA_DIR = os.path.join(_CURRENT_DIR, '..', 'data')
RAW_DATA_DIR = os.path.join(BASE_DATA_DIR, 'raw')
# ... (rest of the paths are the same)
PROCESSED_DATA_DIR = os.path.join(BASE_DATA_DIR, 'processed_db_dfs')
RECOVER_ID_DATA_DIR = os.path.join(BASE_DATA_DIR, 'recovered_ids_dfs')
PROCESSED_COMBINED_DATA_DIR = os.path.join(BASE_DATA_DIR, 'processed_combined_db')
SQLITE_DB_PATH = os.path.join(BASE_DATA_DIR, 'sapiens_ontocellmap.sqlite')

os.makedirs(RAW_DATA_DIR, exist_ok=True)
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
os.makedirs(RECOVER_ID_DATA_DIR, exist_ok=True)
os.makedirs(PROCESSED_COMBINED_DATA_DIR, exist_ok=True)


# --- NEW: Dynamic Database and Parser Configuration ---
DATABASE_CONFIG = {
    "cellmarkerdb": {
        "source": ["http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx", "Cell_marker_Human.xlsx", "xlsx"],
        "parser_class": CellMarkerDBParser,
        "parser_config": {} # This parser takes the df directly, no extra config needed
    },
    "hubmap": {
        "source": ["", "pre_manually_downloaded_files/hubmap_summary.csv", "csv"],
        "parser_class": HuBMapDBParser,
        "parser_config": {}
    },
    "cellxgene": {
        "source": ["https://cellguide.cellxgene.cziscience.com/1716401368/computational_marker_genes/marker_gene_data.json.gz", "marker_gene_data.json.gz", "json"],
        "parser_class": CellxGeneDBParser,
        "parser_config": {}
    },
    "panglaodb": {
        "source": ["https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", "PanglaoDB_markers_27_Mar_2020.tsv.gz", "tsv"],
        "parser_class": PanglaoParser,
        "parser_config": {}
    },
    "wimms": {
        "source": ["https://wimms.tanlab.org/session/9b2f3eed2686b8a8c7da5250a564a808/download/download_signatures?w=", "wimms_signature.csv", "csv"],
        "parser_class": WimmsMelanocyteParser,
        "parser_config": {}
    },
    "epi_cluster_gse147482": {
        "source": ["https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-18075-7/MediaObjects/41467_2020_18075_MOESM4_ESM.xlsx", "Epi_Cluster_gse147482.xlsx", "xlsx"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "Epi_clusters_GSE147482",
            "gene_col": "Row",
            "tissue_name": "Skin",
            "base_cell_name": "Epithelial",
            "cell_subtype_col": "Cluster",
            "source_type": "Computational",
            "source_info_cols": ["DE_Score"]
        }
    },
    "skin_fibroblast_atlas": {
        "source": ["https://raw.githubusercontent.com/haniffalab/skin_fibroblast_atlas/main/figures/SupplementaryTable3.csv", "Skin_fibroblast_atlas.csv", "csv"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "Fibroblast Skin Atlas",
            "gene_col": "Gene",
            "tissue_name": "Skin",
            "base_cell_name": "Fibroblast",
            "cell_subtype_col": "Group",
            "source_type": "Computational",
            "source_info_cols": ["p-value", "logFC"]
        }
    },
    "scar_cell_marker_gse130973": {
        "source": ["", "pre_manually_downloaded_files/GSE130973_2020_sc_skin_Cell_anno_V2_marker_genes.csv", "csv"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "ScarCellMarker_GSE130973",
            "gene_col": "gene",
            "tissue_name": "Skin",
            "cell_name_col": "cluster",
            "source_type": "Computational",
            "source_info_cols": ["p_val_adj", "avg_logFC"]
        }
    },
    "scar_cell_marker_gse163973": {
        "source": ["", "pre_manually_downloaded_files/GSE163973_2021_NC_sc_keloid_Cell_anno_V2_marker_genes.csv", "csv"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "ScarCellMarker_GSE163973",
            "gene_col": "gene",
            "tissue_name": "Skin",
            "cell_name_col": "cluster",
            "source_type": "Computational",
            "source_info_cols": ["p_val_adj", "avg_logFC"]
        }
    },
    "scar_cell_marker_gse138669": {
        "source": ["", "pre_manually_downloaded_files/GSE138669_2018_JID_sc_skin_Cell_anno_V2_marker_genes.csv", "csv"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "ScarCellMarker_GSE138669",
            "gene_col": "gene",
            "tissue_name": "Skin",
            "cell_name_col": "cluster",
            "source_type": "Computational",
            "source_info_cols": ["p_val_adj", "avg_logFC"]
        }
    },
    "scar_cell_marker_gse156326": {
        "source": ["", "pre_manually_downloaded_files/GSE156326_2021_NC_sc_hyper_Cell_anno_V2_marker_genes.csv", "csv"],
        "parser_class": GenericFileParser,
        "parser_config": {
            "database_name": "ScarCellMarker_GSE156326",
            "gene_col": "gene",
            "tissue_name": "Skin",
            "cell_name_col": "cluster",
            "source_type": "Computational",
            "source_info_cols": ["p_val_adj", "avg_logFC"]
        }
    }
}

# HuBMAP: https://humanatlas.io/asctb-tables manually download the latest version and further once parsed the summary file to get the organ csv url and download 
# ScarCellMarker: http://124.220.48.30:3838/ScarCellMarker/#tab-8377-7 manually downloaded by clicking on the ratio buttons of each study