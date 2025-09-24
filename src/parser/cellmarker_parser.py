#!/usr/bin/python3

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Parses cellmarker database and combine information from cell ontology and uberon 
"""

import os
import sys
import pandas as pd
from typing import Optional, List, Union
from cell_ontology_parser import CellOntologyParser

# Add the project root to sys.path to resolve absolute imports from the config directory.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from config.config import (
    CELLMARKER_EXCEL_FILENAME,
    RAW_DATA_DIR
)

class CellMarkerDBParser:
    """
    A class to parse and provide lookup functionalities for Cell Marker database and combine with the cell ontology using cellontology_id and uberon using uberonogology_id
    """
    def __init__(self):
        """
        Initializes the CellMarkerDBParser by loading and processing the file.
        It constructs file paths and prepares lookup structures for efficient querying.
        """
        self.cell_ontology_obj = CellOntologyParser()
        self.cell_marker_db_file = os.path.join(RAW_DATA_DIR, CELLMARKER_EXCEL_FILENAME)
        self.cell_marker_db = pd.read_excel(self.cell_marker_db_file)
        self.cell_marker_db = self.cell_marker_db[self.cell_marker_db['cancer_type'] == "Normal"]

        # get the cell ontology based cell names 
        self.cell_marker_db['cell_ontology_based_names'] =  [self.cell_ontology_obj.get_cell_name_given_cell_id(cell_id) for cell_id in self.cell_marker_db['cellontology_id']]

        # get uberon based dev layer, organs and tissues
        


