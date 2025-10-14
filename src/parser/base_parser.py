#!/usr/bin/python3
"""
Author         : Gemini
Date           : 14/10/2025
Description    : A base parser to handle common ontology normalization and schema finalization.
"""
import pandas as pd
from typing import Tuple, List

from src.parser.ontology_utils import CellxGeneOntologyParser

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
        Takes a pre-formatted DataFrame and performs all normalization steps.

        The input DataFrame must contain a subset of the following columns:
        'db_tissue_name', 'db_tissue_id', 'db_cell_name', 'db_cell_id',
        and other required fields like 'gene', 'source_type', 'source_info', 'database'.

        Returns:
            A tuple containing the final processed DataFrame and the recovery log DataFrame.
        """
        self.recovery_logs = [] # Reset logs for each run

        # --- Step 1 & 2: Retrieve missing IDs using db_names ---
        df['tissue_id'] = df.apply(
            lambda row: self._find_tissue_id(row.get('db_tissue_id'), row.get('db_tissue_name')),
            axis=1
        )
        df['cell_id'] = df.apply(
            lambda row: self._find_cell_id(row.get('db_cell_id'), row.get('db_cell_name')),
            axis=1
        )

        # --- Step 3: Create NLP flags ---
        # The recovery log contains the match type, which we can use to create flags.
        # This step is implicitly handled when we build the recovery DataFrame.
        # We will merge the flags back later.

        # --- Step 5: Get official names from final IDs ---
        df['tissue_name'] = df['tissue_id'].apply(lambda x: self.p.get_tissue_name_given_id(x.replace('_', ':')) if pd.notna(x) else None)
        df['cell_name'] = df['cell_id'].apply(lambda x: self.p.get_cell_name_given_id(x.replace('_', ':')) if pd.notna(x) else None)

        # --- Step 4: Create and format the Recovery DataFrame ---
        recovery_df = pd.DataFrame(self.recovery_logs)
        
        # Add flags to the main dataframe
        if not recovery_df.empty:
            tissue_flags = recovery_df[recovery_df['type'] == 'tissue'][['query', 'match_type']].rename(columns={'query': 'db_tissue_name', 'match_type': 'nlp_tissue_flag'})
            cell_flags = recovery_df[recovery_df['type'] == 'cell'][['query', 'match_type']].rename(columns={'query': 'db_cell_name', 'match_type': 'nlp_cell_flag'})
            
            # Use pd.merge to add the flags
            if not tissue_flags.empty:
                df = pd.merge(df, tissue_flags, on='db_tissue_name', how='left')
            else:
                 df['nlp_tissue_flag'] = None

            if not cell_flags.empty:
                df = pd.merge(df, cell_flags, on='db_cell_name', how='left')
            else:
                df['nlp_cell_flag'] = None
        else:
            df['nlp_tissue_flag'] = None
            df['nlp_cell_flag'] = None
            
        # For cellmarkerdb specifically, handle cancer type suffix
        if 'cancer_type' in df.columns:
             df['cell_name'] = df.apply(
                lambda row: f"{row['cell_name']}_{row['cancer_type']}" if pd.notna(row['cancer_type']) and row['cancer_type'] else row['cell_name'],
                axis=1
            )


        # --- Step 7: Finalize Schema ---
        final_cols = [
            'db_tissue_name', 'db_tissue_id', 'db_cell_name', 'db_cell_id',
            'tissue_name', 'tissue_id', 'cell_name', 'cell_id',
            'nlp_tissue_flag', 'nlp_cell_flag',
            'gene', 'source_type', 'source_info', 'database'
        ]
        
        # Add columns that might be missing in some parsers
        for col in final_cols:
            if col not in df.columns:
                df[col] = None

        return df[final_cols], recovery_df

    def _find_tissue_id(self, tissue_id, tissue_name):
        """Helper to find tissue ID, preferring existing ID, then fuzzy matching name."""
        if pd.notna(tissue_id):
            self.recovery_logs.append({
                "type": "tissue", "query": tissue_name or tissue_id, "match": self.p.get_tissue_name_given_id(tissue_id),
                "match_id": tissue_id, "match_score": 100, "match_type": "pre_supplied_id"
            })
            return tissue_id

        if pd.notna(tissue_name):
            match = self.p.find_best_tissue_name_match(tissue_name)
            if match:
                self.recovery_logs.append({**match, "type": "tissue"}) # Add type key
                return match["match_id"]
        
        self.recovery_logs.append({
            "type": "tissue", "query": tissue_name, "match": None, "match_id": None, 
            "match_score": 0.0, "match_type": "no_match"
        })
        return None

    def _find_cell_id(self, cell_id, cell_name):
        """Helper to find cell ID, preferring existing ID, then fuzzy matching name."""
        if pd.notna(cell_id):
            self.recovery_logs.append({
                "type": "cell", "query": cell_name or cell_id, "match": self.p.get_cell_name_given_id(cell_id),
                "match_id": cell_id, "match_score": 100, "match_type": "pre_supplied_id"
            })
            return cell_id

        if pd.notna(cell_name):
            match = self.p.find_best_cell_name_match(cell_name)
            if match:
                self.recovery_logs.append({**match, "type": "cell"}) # Add type key
                return match["match_id"]
        
        self.recovery_logs.append({
            "type": "cell", "query": cell_name, "match": None, "match_id": None, 
            "match_score": 0.0, "match_type": "no_match"
        })
        return None