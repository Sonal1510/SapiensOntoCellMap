#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses Panglao database and combines information from Cell Ontology and Uberon.
"""
import os
import sys
import pandas as pd
from typing import Optional, List, Union

# Add the project root to sys.path to resolve absolute imports from the config directory.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from src.parser.ontology_utils import *
from src.download.bio_database_downloader import BioDataDownloader
from config.config import RAW_DATA_DIR, DATABASE_SOURCE_DICTIONARY

class PanglaoParser:
    """
    A class to parse and provide lookup functionalities for the Panglao database 
    and combine it with the Cell Ontology (using cell_id) and Uberon (using tissue_id).
    """
    def __init__(self, panglaodb_df: pd.DataFrame):
        """
        Initializes the PanglaoParser by processing the pre-loaded DataFrame.
        It maps names to ontology IDs and logs all mapping attempts.

        Args:
            panglaodb_df (pd.DataFrame): The input DataFrame from the Panglao database.
        """
        if not isinstance(panglaodb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Initialization and Setup ---
        self.p = CellxGeneOntologyParser()
        self.recovery_logs = []
        self.tissue_lookup = {k.lower(): v for k, v in self.p.uberon_name_to_id.items()}
        self.cell_lookup = {k.lower(): v for k, v in self.p.cl_name_to_id.items()}
        
        # --- Step 2: Filter Data and Create a Safe Copy ---
        # Filter for human and mixed-species markers and work on a copy.
        df = panglaodb_df[panglaodb_df['species'].isin(['Hs', 'Mm Hs'])].copy()

        # --- Step 3: Map to Standard Schema ---
        df['tissue_name'] = df['organ']
        df['cell_name'] = df['cell type']
        df['gene'] = df['official gene symbol']
        
        # --- Step 4: Map Names to Ontology IDs with Logging ---
        df['tissue_id'] = df['tissue_name'].apply(self._map_and_log_tissue)
        df['cell_id'] = df['cell_name'].apply(self._map_and_log_cell_by_name)

        # --- Step 5: Format Source Metadata ---
        df['source_type'] = "Literature"
        df['source_info'] = df.apply(
            lambda row: (
                f"canonical_marker_flag={row.get('canonical marker', '')}"
                f"|human_sensitivity={row.get('sensitivity_human', '')}"
                f"|human_specificity={row.get('specificity_human', '')}"
                f"|mouse_sensitivity={row.get('sensitivity_mouse', '')}"
                f"|mouse_specificity={row.get('specificity_mouse', '')}"
            ),
            axis=1
        )
        df['database'] = "PanglaoDB"

        # --- Step 6: Finalize Schema and Store Results ---
        final_cols = ["tissue_name", "tissue_id", "cell_name", "cell_id", "gene", "source_type", "source_info", "database"]
        self.processed_df = df[final_cols]
        self.recovery_df = pd.DataFrame(self.recovery_logs)

    def _map_and_log_tissue(self, tissue_name: str) -> Optional[str]:
        """Maps tissue name to ID, logs the process, and attempts fuzzy match as fallback."""
        if pd.isna(tissue_name):
            return None
        
        # Attempt 1: Direct dictionary lookup (case-insensitive)
        match_id = self.tissue_lookup.get(tissue_name.lower())
        if match_id:
            self.recovery_logs.append({
                "type": "tissue", "query": tissue_name, "match": tissue_name.lower(),
                "match_id": match_id, "match_score": 1.0, "match_type": "direct_match"
            })
            return match_id
        
        # Attempt 2: Fuzzy matching
        match = self.p.find_best_tissue_name_match(tissue_name)
        if match:
            self.recovery_logs.append({
                "type": "tissue", "query": match["query"], "match": match["match_label"],
                "match_id": match["match_id"], "match_score": match["score"], "match_type": "fuzzy_match"
            })
            return match["match_id"]

        # Log failure
        self.recovery_logs.append({
            "type": "tissue", "query": tissue_name, "match": None, "match_id": None, 
            "match_score": 0.0, "match_type": "no_match"
        })
        return None

    def _map_and_log_cell_by_name(self, cell_name: str) -> Optional[str]:
        """Maps cell name to ID, logs the process, and attempts fuzzy match as fallback."""
        if pd.isna(cell_name):
            return None

        # Attempt 1: Direct dictionary lookup (case-insensitive)
        match_id = self.cell_lookup.get(cell_name.lower())
        if match_id:
            self.recovery_logs.append({
                "type": "cell", "query": cell_name, "match": cell_name.lower(),
                "match_id": match_id, "match_score": 1.0, "match_type": "direct_match"
            })
            return match_id

        # Attempt 2: Fuzzy matching
        match = self.p.find_best_cell_name_match(cell_name)
        if match:
            self.recovery_logs.append({
                "type": "cell", "query": match["query"], "match": match["match_label"],
                "match_id": match["match_id"], "match_score": match["score"], "match_type": "fuzzy_match"
            })
            return match["match_id"]

        # Log failure
        self.recovery_logs.append({
            "type": "cell", "query": cell_name, "match": None, "match_id": None,
            "match_score": 0.0, "match_type": "no_match"
        })
        return None