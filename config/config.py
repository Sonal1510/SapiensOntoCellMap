# config/config.py

"""
Author        : Sonal Rashmi (expert review by Gemini)
Date          : 15/08/2025
Description   : Centralized file for all URLs and file paths used for databases
"""

import os

# --- Directory Paths ---
# Get the absolute path to the directory this script is in
_CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# Navigate up one level to the root and then into 'data'
BASE_DATA_DIR = os.path.join(_CURRENT_DIR, '..', 'data')
RAW_DATA_DIR = os.path.join(BASE_DATA_DIR, 'raw')
PROCESSED_DATA_DIR = os.path.join(BASE_DATA_DIR, 'processed_db_dfs')
RECOVER_ID_DATA_DIR = os.path.join(BASE_DATA_DIR, 'recovered_ids_dfs')
PROCESSED_COMBINED_DATA_DIR = os.path.join(BASE_DATA_DIR, 'processed_combined_db')
SQLITE_DB_PATH = os.path.join(BASE_DATA_DIR, 'sapiens_ontocellmap.sqlite')

# Ensure directories exist
os.makedirs(RAW_DATA_DIR, exist_ok=True)
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
os.makedirs(RECOVER_ID_DATA_DIR, exist_ok=True)
os.makedirs(PROCESSED_COMBINED_DATA_DIR, exist_ok=True)
# --- Databases sources ---

DATABASE_SOURCE_DICTIONARY = {
	"CELLMARKER_DB" : ["http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx", "Cell_marker_Human.xlsx", "xlsx"],
	"HUBMAP_DB" : ["", "hubmap_summary.csv", "csv"],
	"CELLXGENE_DB" : ["https://cellguide.cellxgene.cziscience.com/1716401368/computational_marker_genes/marker_gene_data.json.gz", "marker_gene_data.json.gz", "json"],
	"PANGLAO_DB" : ["https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", "PanglaoDB_markers_27_Mar_2020.tsv.gz", "tsv"]
}


# HuBMAP: https://humanatlas.io/asctb-tables manually download the latest version and further once parsed the summary file to get the organ csv url and download 