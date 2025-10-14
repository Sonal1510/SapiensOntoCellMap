#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses CellxGene database and uses BaseParser for normalization.
"""
import pandas as pd
from src.parser.base_parser import BaseParser # Import the new base class

class CellxGeneDBParser:
    """
    A class to parse the CellxGene database. It unnests the data structure
    and then uses the BaseParser for ontology mapping and normalization.
    """
    def __init__(self, cellxgenedb_df: pd.DataFrame):
        if not isinstance(cellxgenedb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Unnest the DataFrame (Parser-specific logic) ---
        rows = []
        for tissue, row in cellxgenedb_df.iterrows():
            for organism in cellxgenedb_df.columns:
                cell_dict = row[organism]
                if isinstance(cell_dict, dict):
                    for cell_id, markers in cell_dict.items():
                        for marker in markers:
                            rows.append({
                                "organism": organism,
                                "db_tissue_name": tissue, # Use db_ prefix
                                "db_cell_id": cell_id,     # Use db_ prefix
                                "gene": marker.get("gene"),
                                "marker_score": marker.get("marker_score"),
                                "pc": marker.get("pc"),
                                "me": marker.get("me"),
                            })
        
        if not rows:
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()
            return

        df = pd.DataFrame(rows)
        human_cellxgene_df = df[df["organism"] == "Homo sapiens"].copy()

        # --- Step 2: Format Source Metadata ---
        human_cellxgene_df["source_type"] = "computational"
        human_cellxgene_df["source_info"] = (
            "marker_score: " + human_cellxgene_df["marker_score"].astype(str) +
            ", pc: " + human_cellxgene_df["pc"].astype(str) +
            ", me: " + human_cellxgene_df["me"].astype(str)
        )
        human_cellxgene_df["database"] = "cellxgene"
        
        # CellxGene provides no cell names, only IDs
        human_cellxgene_df["db_cell_name"] = None 
        human_cellxgene_df["db_tissue_id"] = None

        # --- Step 3: Use the BaseParser for normalization ---
        normalizer = BaseParser()
        self.processed_df, self.recovery_df = normalizer.normalize_dataframe(human_cellxgene_df)