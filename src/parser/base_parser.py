#!/usr/bin/python3
"""
Author          : Gemini
Date            : 14/10/2025
Description     : A base parser to handle common ontology normalization and schema finalization.
                  Optimized to use vectorized mapping for performance and memory efficiency.
                  Includes logic for pre-search remapping and post-search keyword imputation.
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
    Implements a multi-stage imputation process:
    1. Use pre-supplied db_id if present.
    2. If not, apply pre-search remapping to db_name (e.g., "Undifferentiated X" -> "X").
    3. Search with the (remapped) name using the 3-stage hybrid search.
    4. If search fails, apply keyword-based imputation (e.g., "cancer" -> neoplastic cell ID).
    """
    def __init__(self):
        """Initializes the base parser with an ontology helper and imputation/remapping rules."""
        self.p = CellxGeneOntologyParser()
        self.recovery_logs = []
        
        # --- Rule: Handle non-specific tissues (from previous request) ---
        self.NON_SPECIFIC_TISSUES = ['all tissue', 'undetermined', 'undefined', 'All Tissues', 'Undefined']
        
        # --- NEW: Rule 4 - Pre-search remapping ---
        # Maps a lowercase query string to a new query string *before* searching
        self.PRE_SEARCH_REMAP = {
            'undifferentiated keratinocytes': 'keratinocyte'
            # Add more pre-search rules here, e.g.:
            # 'tumor-infiltrating t-cell': 't-cell'
        }
        
        # --- NEW: Rule 3 - Post-search keyword imputation ---
        # Maps keywords to a specific ID if the main search fails (result is NaN)
        
        # For Cells:
        cancer_cell_id = 'CL:0001063'
        cancer_cell_name = self.p.get_cell_name_given_id(cancer_cell_id) or 'neoplastic cell'
        
        self.CELL_KEYWORD_IMPUTATIONS = {
            # Key: The ID to impute
            # Value: A dict with the official name and a list of keywords
            cancer_cell_id: {
                'name': cancer_cell_name,
                'keywords': ['cancer', 'neoplastic', 'tumor', 'malignant', 'carcinoma']
            }
            # Add more cell keyword rules here, e.g.:
            # 'CL:0000037': {
            #     'name': self.p.get_cell_name_given_id('CL:0000037') or 'fibroblast',
            #     'keywords': ['fibroblast-like', 'fibroblastic']
            # }
        }

        # For Tissues:
        self.TISSUE_KEYWORD_IMPUTATIONS = {
            # Add tissue keyword rules here, e.g.:
            # 'UBERON:0002107': {
            #     'name': self.p.get_tissue_name_given_id('UBERON:0002107') or 'liver',
            #     'keywords': ['hepatic']
            # }
        }

    def normalize_dataframe(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Takes a pre-formatted DataFrame and performs all normalization steps using
        an optimized, vectorized approach to handle large datasets efficiently.
        """
        self.recovery_logs = []
        if df.empty:
            return pd.DataFrame(), pd.DataFrame()

        df_proc = df.copy()

        # --- Step 1: Create lookup maps (Implements Rule 2 & 4) ---
        # This step runs the 3-stage search on all unique names,
        # *after* applying the pre-search remapping rules.
        tissue_id_map, tissue_flag_map = self._create_lookup_maps(
            unique_names=df_proc['db_tissue_name'].dropna().unique(),
            finder_func=self.p.find_best_tissue_name_match,
            log_type='tissue'
        )
        
        cell_id_map, cell_flag_map = self._create_lookup_maps(
            unique_names=df_proc['db_cell_name'].dropna().unique(),
            finder_func=self.p.find_best_cell_name_match,
            log_type='cell'
        )

        # --- Step 2: Map Tissue IDs (Implements Rule 1) ---
        # Get IDs from the name-based search
        mapped_tissue_ids = df_proc['db_tissue_name'].map(tissue_id_map)
        # Use db_tissue_id if it exists, otherwise use the mapped ID
        df_proc['tissue_id'] = df_proc['db_tissue_id'].fillna(mapped_tissue_ids)
        df_proc['nlp_tissue_flag'] = df_proc['db_tissue_name'].map(tissue_flag_map)

        # --- Step 3: Map Cell IDs (Implements Rule 1) ---
        mapped_cell_ids = df_proc['db_cell_name'].map(cell_id_map)
        df_proc['cell_id'] = df_proc['db_cell_id'].fillna(mapped_cell_ids)
        df_proc['nlp_cell_flag'] = df_proc['db_cell_name'].map(cell_flag_map)
        
        # --- Step 4: Post-Search Keyword Imputation (Implements Rule 3) ---
        # This block runs *only* if cell_id is still NaN after Steps 1-3.
        
        # --- Step 4.5: Impute TISSUE IDs by keyword ---
        missing_tissue_mask = df_proc['tissue_id'].isna()
        if missing_tissue_mask.any():
            # Get just the rows where tissue_id is still NaN
            db_tissue_names_lower = df_proc.loc[missing_tissue_mask, 'db_tissue_name'].str.lower()
            
            for tissue_id, config in self.TISSUE_KEYWORD_IMPUTATIONS.items():
                keyword_regex = '|'.join(config['keywords'])
                # Find matches within the 'missing' subset
                impute_mask = db_tissue_names_lower.str.contains(keyword_regex, na=False, regex=True)
                
                # Get the original DataFrame indices of these rows
                rows_to_impute_idx = db_tissue_names_lower[impute_mask].index
                
                if not rows_to_impute_idx.empty:
                    # Apply imputation to the main DataFrame
                    df_proc.loc[rows_to_impute_idx, 'tissue_id'] = tissue_id
                    df_proc.loc[rows_to_impute_idx, 'nlp_tissue_flag'] = 'imputed_keyword'
                    
                    # Log this
                    for _, row in df_proc.loc[rows_to_impute_idx].iterrows():
                        self.recovery_logs.append({
                            "log_type": "tissue", "query": row['db_tissue_name'],
                            "match": config['name'], "match_id": tissue_id,
                            "score": 100, "match_type": "imputed_keyword"
                        })
                    
                    # Update the missing_tissue_mask to exclude rows we just imputed
                    missing_tissue_mask = df_proc['tissue_id'].isna()
                    if not missing_tissue_mask.any():
                        break # All tissues imputed, stop looping

        # --- Step 4.6: Impute CELL IDs by keyword ---
        missing_cell_mask = df_proc['cell_id'].isna()
        if missing_cell_mask.any():
            db_cell_names_lower = df_proc.loc[missing_cell_mask, 'db_cell_name'].str.lower()

            for cell_id, config in self.CELL_KEYWORD_IMPUTATIONS.items():
                keyword_regex = '|'.join(config['keywords'])
                impute_mask = db_cell_names_lower.str.contains(keyword_regex, na=False, regex=True)
                
                rows_to_impute_idx = db_cell_names_lower[impute_mask].index
                
                if not rows_to_impute_idx.empty:
                    df_proc.loc[rows_to_impute_idx, 'cell_id'] = cell_id
                    df_proc.loc[rows_to_impute_idx, 'nlp_cell_flag'] = 'imputed_keyword'
                    
                    # Log this
                    for _, row in df_proc.loc[rows_to_impute_idx].iterrows():
                        self.recovery_logs.append({
                            "log_type": "cell", "query": row['db_cell_name'],
                            "match": config['name'], "match_id": cell_id,
                            "score": 100, "match_type": "imputed_keyword"
                        })
                    
                    # Update mask
                    missing_cell_mask = df_proc['cell_id'].isna()
                    if not missing_cell_mask.any():
                        break # All cells imputed, stop looping

        # --- Step 5: Log any IDs that were provided in the original data ---
        self._log_presupplied_ids(df_proc)

        # --- Step 6: Get official names from final IDs ---
        # (This is where NaNs will remain if all imputation steps failed)
        corrected_tissue_ids = df_proc['tissue_id'].str.replace('_', ':', n=1, regex=False)
        corrected_cell_ids = df_proc['cell_id'].str.replace('_', ':', n=1, regex=False)

        df_proc['tissue_name'] = corrected_tissue_ids.map(self.p.uberon_id_to_name)
        df_proc['cell_name'] = corrected_cell_ids.map(self.p.cl_id_to_name)

        # --- Step 6.5 - Handle non-specific tissue names ---
        # This ensures 'tissue_name' (final) is the 'db_tissue_name' if it's non-specific
        nonspecific_mask = df_proc['nlp_tissue_flag'] == 'non_specific'
        
        if nonspecific_mask.any():
            df_proc.loc[nonspecific_mask, 'tissue_name'] = df_proc.loc[nonspecific_mask, 'db_tissue_name']
            
        # --- Step 7: Finalize Schema ---
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
        """
        (MODIFIED)
        Creates lookup maps for IDs and NLP flags based on a list of unique names.
        This now implements:
        1. Non-specific tissue check
        2. Pre-search remapping (Rule 4)
        3. Standard ontology search (Rule 2)
        """
        id_map = {}
        flag_map = {}
        for name in unique_names:
            if not name or pd.isna(name):
                continue
                
            # --- Step 1: Handle non-specific tissues ---
            if log_type == 'tissue' and name.lower() in self.NON_SPECIFIC_TISSUES:
                id_map[name] = None
                flag_map[name] = 'non_specific'
                self.recovery_logs.append({
                    "log_type": log_type,
                    "query": name,
                    "match": name, # Match is the name itself
                    "match_id": None,
                    "score": 100,
                    "match_type": "non_specific"
                })
                continue # Skip the normal search

            # --- NEW: Step 2: Apply Pre-Search Remapping (Rule 4) ---
            query_name = self.PRE_SEARCH_REMAP.get(name.lower(), name)
            was_remapped = (query_name != name)

            # --- Step 3: Standard ontology search (Rule 2) ---
            match = finder_func(query_name) # Use the (potentially) remapped name
            
            if match:
                # The ID is mapped from the *original* name for the pandas .map() function
                id_map[name] = match["match_id"] 
                match_type = match.get("type", "fuzzy")
                
                # Prepare log entry
                log_entry = {**match, "log_type": log_type}
                
                if was_remapped:
                    log_entry["original_query"] = name # Log the original query
                    log_entry["query"] = query_name # Log what was *actually* searched
                    log_entry["match_type"] = "remapped_search"
                else:
                    log_entry["match_type"] = match_type
                
                # Use the new, more descriptive type for the flag
                flag_map[name] = log_entry["match_type"] 

                # Rename 'type' key from ontology_utils to avoid confusion
                if "type" in log_entry:
                    del log_entry["type"]
                    
                self.recovery_logs.append(log_entry)
            else:
                # Log a "no_match" against the *original* name
                self.recovery_logs.append({
                    "log_type": log_type, "query": name, "match": None, "match_id": None, 
                    "score": 0, "match_type": "no_match"
                })
        return id_map, flag_map


    def _log_presupplied_ids(self, df: pd.DataFrame):
        """
        (MODIFIED)
        Logs items where the ID was pre-supplied.
        Adds the correct 'match_type' and 'log_type' flags.
        """
        # Find rows where db_tissue_id was present AND was used as the final tissue_id
        tissue_presupplied = df[df['db_tissue_id'].notna() & (df['tissue_id'] == df['db_tissue_id'])]
        for _, row in tissue_presupplied.iterrows():
            self.recovery_logs.append({
                "log_type": "tissue", "query": row['db_tissue_name'] or row['db_tissue_id'],
                "match": self.p.get_tissue_name_given_id(row['tissue_id'].replace('_',':')),
                "match_id": row['tissue_id'], "score": 100, "match_type": "pre_supplied_id"
            })
            
        cell_presupplied = df[df['db_cell_id'].notna() & (df['cell_id'] == df['db_cell_id'])]
        for _, row in cell_presupplied.iterrows():
            self.recovery_logs.append({
                "log_type": "cell", "query": row['db_cell_name'] or row['db_cell_id'],
                "match": self.p.get_cell_name_given_id(row['cell_id'].replace('_',':')),
                "match_id": row['cell_id'], "score": 100, "match_type": "pre_supplied_id"
            })