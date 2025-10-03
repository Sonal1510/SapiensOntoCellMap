#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses CellxGene database and combines information from Cell Ontology and Uberon.
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

class CellxGeneDBParser:
    """
    A class to parse and provide lookup functionalities for the CellxGene database 
    and combine it with the Cell Ontology (using cell_id) and Uberon (using tissue_id).
    """
    def __init__(self, cellxgenedb_df: pd.DataFrame):
        """
        Initializes the CellxGeneDBParser by loading and processing the DataFrame.
        It unnests the complex dictionary structure, maps names to ontology IDs, and logs all mapping attempts.

        Args:
            cellxgenedb_df (pd.DataFrame): The input DataFrame from the CellxGene database.
        """
        if not isinstance(cellxgenedb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Initialization and Setup ---
        self.p = CellxGeneOntologyParser()  
        self.recovery_logs = []
        df_to_process = cellxgenedb_df.copy()

        # --- Step 2: Unnest the DataFrame ---
        # The source data has tissues in the index, organisms in columns, and dicts in cells.
        rows = []
        for tissue, row in df_to_process.iterrows():
            for organism in df_to_process.columns:
                cell_dict = row[organism]
                if isinstance(cell_dict, dict):  # Skip NaNs or other non-dict entries
                    for cell_id, markers in cell_dict.items():
                        for marker in markers:
                            rows.append({
                                "organism": organism,
                                "tissue_name": tissue,
                                "cell_id": cell_id,
                                "gene": marker.get("gene"),
                                "marker_score": marker.get("marker_score"),
                                "pc": marker.get("pc"),
                                "me": marker.get("me"),
                            })
        
        if not rows:
            # Handle case where no data was parsed
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()
            return

        cellxgene_df = pd.DataFrame(rows)
        # Filter for human data and create a safe copy to avoid SettingWithCopyWarning
        human_cellxgene_df = cellxgene_df[cellxgene_df["organism"] == "Homo sapiens"].copy()

        # --- Step 3: Map Tissue Name to Uberon ID with Logging ---
        self.tissue_lookup = {k.lower(): v for k, v in self.p.uberon_name_to_id.items()}
        human_cellxgene_df["tissue_id"] = human_cellxgene_df["tissue_name"].apply(self._map_and_log_tissue)
        
        # --- Step 4: Map Cell ID to Cell Name with Logging ---
        human_cellxgene_df["cell_name"] = human_cellxgene_df["cell_id"].apply(self._map_and_log_cell)
        
        # --- Step 5: Format Source Metadata ---
        human_cellxgene_df["source_type"] = "computational"
        human_cellxgene_df["source_info"] = (
            "marker_score: " + human_cellxgene_df["marker_score"].astype(str) +
            ", pc: " + human_cellxgene_df["pc"].astype(str) +
            ", me: " + human_cellxgene_df["me"].astype(str)
        )
        human_cellxgene_df["database"] = "cellxgene"

        # --- Step 6: Finalize Schema ---
        self.processed_df = human_cellxgene_df[
            ["tissue_name", "tissue_id", "cell_name", "cell_id", "gene", "source_type", "source_info", "database"]
        ]
        
        # --- Step 7: Store Recovery Logs ---
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
        
        # Attempt 2: Fuzzy matching for names that failed direct lookup
        match = self.p.find_best_tissue_name_match(tissue_name)
        if match:
            self.recovery_logs.append({
                "type": "tissue", "query": match["query"], "match": match["match_label"],
                "match_id": match["match_id"], "match_score": match["score"], "match_type": "fuzzy_match"
            })
            return match["match_id"]

        # Log failure if no match was found
        self.recovery_logs.append({
            "type": "tissue", "query": tissue_name, "match": None, "match_id": None, 
            "match_score": 0.0, "match_type": "no_match"
        })
        return None

    def _map_and_log_cell(self, cell_id: str) -> str:
        """Maps cell ID to a cell name and logs the process."""
        if pd.isna(cell_id):
            return "" # Return empty string for consistency
        
        # Attempt direct lookup
        match_name = self.p.cl_id_to_name.get(cell_id)
        if match_name:
            self.recovery_logs.append({
                "type": "cell", "query": cell_id, "match": match_name, "match_id": cell_id,
                "match_score": 1.0, "match_type": "direct_match"
            })
            return match_name
        
        # Log failure and fallback to the original ID as the name
        self.recovery_logs.append({
            "type": "cell", "query": cell_id, "match": None, "match_id": cell_id,
            "match_score": 0.0, "match_type": "no_match"
        })
        return cell_id