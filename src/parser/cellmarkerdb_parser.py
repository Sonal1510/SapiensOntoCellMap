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

# Add the project root to sys.path to resolve absolute imports from the config directory.
# This assumes 'ontology_utils' can be resolved from this path.
# You might need to adjust the path based on your project structure.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Assuming 'ontology_utils' contains a class (e.g., 'OntologyParser') that provides the required methods.
from src.parser.ontology_utils import *

class CellMarkerDBParser:
    """
    A class to parse and provide lookup functionalities for the Cell Marker database 
    and combine it with the Cell Ontology (using cellontology_id) and Uberon (using uberonongology_id).
    """
    def __init__(self, cellmarkerdb_df: pd.DataFrame):
        """
        Initializes the CellMarkerDBParser by loading and processing the DataFrame.
        It constructs file paths and prepares lookup structures for efficient querying.

        Args:
            cellmarkerdb_df (pd.DataFrame): The input DataFrame from the CellMarker database.
        """
        if not isinstance(cellmarkerdb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Initialization and Setup ---
        self.p = CellxGeneOntologyParser() 
        self.recovery_logs = []
        df = cellmarkerdb_df.copy() # Work on a copy to avoid modifying the original DataFrame

        # --- Step 2: Schema Mapping ---
        df['tissue_name'] = df['tissue_class']
        df['tissue_id']   = df['uberonongology_id']
        df['cell_id']     = df['cellontology_id']
        df['gene']        = df['marker']
        df['source_type'] = df['marker_source']
        df['source_info'] = (
            df['Title'].astype(str) + ";" +
            df['journal'].astype(str) + ";" +
            df['year'].astype(str)
        )

        # --- Step 3: Identify Missing IDs ---
        missing_tissue_mask = df['tissue_id'].isna()
        missing_cell_mask   = df['cell_id'].isna()

        # --- Step 4: Recover Missing IDs using Fuzzy Matching ---
        # Apply recovery functions only on rows with missing IDs
        if missing_tissue_mask.any():
            df.loc[missing_tissue_mask, 'tissue_id'] = (
                df.loc[missing_tissue_mask, 'tissue_name'].apply(self.safe_get_tissue)
            )
        if missing_cell_mask.any():
            df.loc[missing_cell_mask, 'cell_id'] = (
                df.loc[missing_cell_mask, 'cell_name'].apply(self.safe_get_cell)
            )
        
        # --- Step 5: Update Names Based on Mapped/Recovered IDs ---
        df['cell_name'] = df.apply(self.get_cell_name, axis=1)
        df['tissue_name'] = df.apply(self.get_tissue_name, axis=1)

        # --- Step 6: Finalize Schema ---
        df["database"] = "cellmarkerdb"
        self.processed_df = df[
            ['tissue_name', 'tissue_id', 'cell_name', 'cell_id', 'gene', 'source_type', 'source_info', 'database']
        ]
        
        # --- Step 7: Store Recovery Logs ---
        self.recovery_df = pd.DataFrame(self.recovery_logs)

    def safe_get_tissue(self, x: str) -> Optional[str]:
        """Safely find the best tissue name match and log the result."""
        if pd.isna(x): 
            return None
        match = self.p.find_best_tissue_name_match(x)
        if match:
            self.recovery_logs.append({
                "type": "tissue",
                "query": match["query"],
                "match": match["match_label"],
                "match_id": match["match_id"],
                "match_score": match["score"],
                "match_type": match["type"],
            })
            return match["match_id"]
        return None

    def safe_get_cell(self, x: str) -> Optional[str]:
        """Safely find the best cell name match and log the result."""
        if pd.isna(x): 
            return None
        match = self.p.find_best_cell_name_match(x)
        if match:
            self.recovery_logs.append({
                "type": "cell",
                "query": match["query"],
                "match": match["match_label"],
                "match_id": match["match_id"],
                "match_score": match["score"],
                "match_type": match["type"],
            })
            return match["match_id"]
        return None

    def get_cell_name(self, row: pd.Series) -> str:
        """Get the official cell name from its ID; otherwise, format the existing name."""
        cancer_type = row.get("cancer_type", "")
        original_cell_name = row.get('cell_name', "")
        
        if pd.notna(row['cell_id']):
            # Look up the official name from the ontology ID
            mapped_name = self.p.cl_id_to_name.get(row['cell_id'].replace('_', ':'))
            # Use the official name if found, otherwise fall back to the original name
            base_name = mapped_name if mapped_name else original_cell_name
        else:
            base_name = original_cell_name
            
        # Append cancer type if it exists
        if cancer_type:
            return f"{base_name}_{cancer_type}"
        return base_name

    def get_tissue_name(self, row: pd.Series) -> Optional[str]:
        """Get the official tissue name from its ID; otherwise, return the existing name."""
        if pd.isna(row['tissue_id']):
            return row.get('tissue_name')
            
        # Look up the official name from the ontology ID
        mapped_name = self.p.uberon_id_to_name.get(row['tissue_id'].replace('_', ':'))
        # Return the official name if found, otherwise fall back to the original name
        return mapped_name if mapped_name else row.get('tissue_name')