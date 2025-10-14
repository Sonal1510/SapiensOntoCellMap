#!/usr/bin/python3
"""
Author         : Gemini
Date           : 14/10/2025
Description    : A base parser to handle common ontology normalization and schema finalization.
                 Optimized to use vectorized mapping for performance and memory efficiency.
"""
import pandas as pd
from typing import Tuple, Dict
import sys
import os
try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.parser.ontology_utils import CellxGeneOntologyParser
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

class BaseParser:
    """
    Handles the common logic for normalizing tissue/cell names and IDs against ontologies.
    """
    def __init__(self):
        """Initializes the base parser with an ontology helper."""
        self.p = CellxGeneOntologyParser()
        self.recovery_logs = []

    def normalize_dataframe(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Takes a pre-formatted DataFrame and performs all normalization steps using
        an optimized, vectorized approach to handle large datasets efficiently.
        """
        self.recovery_logs = []
        if df.empty:
            return pd.DataFrame(), pd.DataFrame()

        df_proc = df.copy()

        # Steps 1-4 (No changes here)
        tissue_id_map, tissue_flag_map = self._create_lookup_maps(
            unique_names=df_proc['db_tissue_name'].dropna().unique(),
            finder_func=self.p.find_best_tissue_name_match,
            log_type='tissue'
        )
        mapped_tissue_ids = df_proc['db_tissue_name'].map(tissue_id_map)
        df_proc['tissue_id'] = df_proc['db_tissue_id'].fillna(mapped_tissue_ids)
        df_proc['nlp_tissue_flag'] = df_proc['db_tissue_name'].map(tissue_flag_map)

        cell_id_map, cell_flag_map = self._create_lookup_maps(
            unique_names=df_proc['db_cell_name'].dropna().unique(),
            finder_func=self.p.find_best_cell_name_match,
            log_type='cell'
        )
        mapped_cell_ids = df_proc['db_cell_name'].map(cell_id_map)
        df_proc['cell_id'] = df_proc['db_cell_id'].fillna(mapped_cell_ids)
        df_proc['nlp_cell_flag'] = df_proc['db_cell_name'].map(cell_flag_map)
        
        self._log_presupplied_ids(df_proc)

        # 5. Get official names from final IDs using ontology's internal map (very fast)
        # --- FIX: Standardize IDs to use colons (':') before mapping ---
        # Create temporary Series with corrected IDs to ensure a successful key lookup.
        # This is a vectorized operation and very fast.
        corrected_tissue_ids = df_proc['tissue_id'].str.replace('_', ':', n=1, regex=False)
        corrected_cell_ids = df_proc['cell_id'].str.replace('_', ':', n=1, regex=False)

        # Map using the corrected Series
        df_proc['tissue_name'] = corrected_tissue_ids.map(self.p.uberon_id_to_name)
        df_proc['cell_name'] = corrected_cell_ids.map(self.p.cl_id_to_name)

        # For cellmarkerdb specifically, handle cancer type suffix
        if 'cancer_type' in df_proc.columns:
             df_proc['cell_name'] = df_proc.apply(
                lambda row: f"{row['cell_name']}_{row['cancer_type']}" if pd.notna(row['cell_name']) and pd.notna(row['cancer_type']) and row['cancer_type'] else row['cell_name'],
                axis=1
            )
            
        # 6. Finalize Schema (No changes here)
        final_cols = [
            'db_tissue_name', 'db_tissue_id', 'db_cell_name', 'db_cell_id',
            'tissue_name', 'tissue_id', 'cell_name', 'cell_id',
            'nlp_tissue_flag', 'nlp_cell_flag',
            'gene', 'source_type', 'source_info', 'database'
        ]
        
        for col in final_cols:
            if col not in df_proc.columns:
                df_proc[col] = None
        
        recovery_df = pd.DataFrame(self.recovery_logs)
        return df_proc[final_cols], recovery_df

    def _create_lookup_maps(self, unique_names: list, finder_func, log_type: str) -> Tuple[Dict, Dict]:
        # ... (no changes to this method) ...
        id_map = {}
        flag_map = {}
        for name in unique_names:
            match = finder_func(name)
            if match:
                id_map[name] = match["match_id"]
                flag_map[name] = match.get("type", "fuzzy") # Corrected from match.get('type') to use "type" as key
                self.recovery_logs.append({**match, "type": log_type})
            else:
                self.recovery_logs.append({
                    "type": log_type, "query": name, "match": None, "match_id": None, 
                    "match_score": 0.0, "match_type": "no_match"
                })
        return id_map, flag_map


    def _log_presupplied_ids(self, df: pd.DataFrame):
        # ... (no changes to this method) ...
        tissue_presupplied = df[df['db_tissue_id'].notna() & (df['tissue_id'] == df['db_tissue_id'])]
        for _, row in tissue_presupplied.iterrows():
            self.recovery_logs.append({
                "type": "tissue", "query": row['db_tissue_name'] or row['db_tissue_id'],
                "match": self.p.get_tissue_name_given_id(row['tissue_id'].replace('_',':')),
                "match_id": row['tissue_id'], "match_score": 100, "match_type": "pre_supplied_id"
            })
            
        cell_presupplied = df[df['db_cell_id'].notna() & (df['cell_id'] == df['db_cell_id'])]
        for _, row in cell_presupplied.iterrows():
            self.recovery_logs.append({
                "type": "cell", "query": row['db_cell_name'] or row['db_cell_id'],
                "match": self.p.get_cell_name_given_id(row['cell_id'].replace('_',':')),
                "match_id": row['cell_id'], "match_score": 100, "match_type": "pre_supplied_id"
            })