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
PROCESSED_DATA_DIR = os.path.join(BASE_DATA_DIR, 'processed')
SQLITE_DB_PATH = os.path.join(BASE_DATA_DIR, 'sapiens_ontocellmap.sqlite')

# Ensure directories exist
os.makedirs(RAW_DATA_DIR, exist_ok=True)
os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)

# --- URLs for data sources ---
UBERON_HUMAN_VIEW_URL = "http://purl.obolibrary.org/obo/uberon/human-view.json"
CELL_ONTOLOGY_HUMAN_VIEW_URL = "https://github.com/obophenotype/cell-ontology/releases/download/v2025-07-30/human-view.tsv"
CELLMARKER_HUMAN_EXCEL_URL = "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx"

# --- Local file names after download ---
UBERON_OBO_FILENAME = "human-view.obo"
CELL_ONTOLOGY_TSV_FILENAME = "human-view.tsv"
CELLMARKER_EXCEL_FILENAME = "Cell_marker_Human.xlsx"

# --- Processed file names ---
PROCESSED_UBERON_PKL = "uberon_parsed.pkl"
PROCESSED_CELL_ONTOLOGY_PKL = "cell_ontology_parsed.pkl"
PROCESSED_CELLMARKER_PKL = "cellmarker_parsed.pkl"
