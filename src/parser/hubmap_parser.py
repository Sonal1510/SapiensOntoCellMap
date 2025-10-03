#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses HuBMAP database and combines information from Cell Ontology and Uberon.
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

class HuBMapDBParser:
    """
    A class to parse and provide lookup functionalities for the HuBMAP database 
    and combine it with the Cell Ontology (using cell_id) and Uberon (using tissue_id).
    """
    def __init__(self, hubmapdb_df: pd.DataFrame):
        """
        Initializes the HuBMapDBParser by downloading and processing data files listed in the input DataFrame.
        It unnests the complex dictionary structure, maps names to ontology IDs, and logs all mapping attempts.

        Args:
            hubmapdb_df (pd.DataFrame): The input DataFrame from the HuBMAP database containing URLs to data files.
        """
        if not isinstance(hubmapdb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Initialization and Setup ---
        self.p = CellxGeneOntologyParser() 
        self.recovery_logs = []
        self.downloader = BioDataDownloader()
        self.tissue_lookup = {k.lower(): v for k, v in self.p.uberon_name_to_id.items()}
        
        hubmap_organ_df_list = {}

        # --- Step 2: Download and Parse Data for Each Organ ---
        for organ in set(hubmapdb_df['Organ'].tolist()):
            if organ != "anatomical systems":
                url_file = hubmapdb_df.loc[hubmapdb_df['Organ'] == organ, 'csv'].iloc[0]
                file_name = f"hubmap_{organ.replace(' ', '_')}.csv"
                file_path = os.path.join(RAW_DATA_DIR, file_name)

                # Download the specific organ's data file
                # Note: Calling a "protected" method; replace with a public one if available.
                self.downloader._download_file(url_file, file_path)

                print(f"--- Parsing {file_name} ---")
                try:
                    temp_df = pd.read_csv(
                        file_path, header=10, on_bad_lines='skip', engine='python'
                    )
                    
                    if not temp_df.empty:
                        # --- Step 3: Process and Standardize the Organ Data ---
                        parsed_df = self._process_organ_dataframe(temp_df, organ)
                        hubmap_organ_df_list[organ] = parsed_df
                        print(f"✅ Successfully processed {file_name}")
                    else:
                        print(f"⚠️ Skipping {file_name} due to empty data.")

                except Exception as e:
                    print(f"❌ Error while parsing {file_name}: {e}")

        # --- Step 4: Combine and Finalize ---
        if hubmap_organ_df_list:
            hubmap_all_organs_df = pd.concat(hubmap_organ_df_list.values(), ignore_index=True)
            # Ensure final schema is consistent with other parsers
            final_cols = ["tissue_name", "tissue_id", "cell_name", "cell_id", "gene", "source_type", "source_info", "database"]
            self.processed_df = hubmap_all_organs_df[final_cols]
        else:
            self.processed_df = pd.DataFrame()
        
        self.recovery_df = pd.DataFrame(self.recovery_logs)

    def _process_organ_dataframe(self, current_df: pd.DataFrame, organ: str) -> pd.DataFrame:
        """Helper function to process the dataframe for a single organ."""
        parsed_df = pd.DataFrame()

        # Determine cell/anatomical structure columns
        if any(col.startswith("CT/") for col in current_df.columns):
            parsed_df['cell_id'] = current_df.get('CT/1/ID')
            parsed_df['cell_name'] = current_df.get('CT/1/LABEL')
        else:
            parsed_df['cell_id'] = current_df.get('AS/1/ID')
            parsed_df['cell_name'] = current_df.get('AS/1/LABEL')

        parsed_df['source_type'] = "Literature"

        # Gather all reference (REF/*) columns into source_info
        ref_cols = [col for col in current_df.columns if col.startswith("REF/")]
        if ref_cols:
            ref_info = current_df[ref_cols].fillna('').astype(str)
            parsed_df['source_info'] = ref_info.apply(lambda row: "; ".join(val for val in row if val), axis=1)
        else:
            parsed_df['source_info'] = ""

        # Gather and explode gene markers (one row per gene)
        gene_marker_cols = [c for c in current_df.columns if ('BGene' in c or 'BProtein' in c) and 'LABEL' in c]
        if gene_marker_cols:
            # Join all gene columns for a row, then split into a list
            parsed_df['gene'] = current_df[gene_marker_cols].apply(
                lambda row: ', '.join(row.dropna().astype(str)), axis=1
            ).str.split(', ')
            parsed_df = parsed_df.explode('gene')
            parsed_df['gene'] = parsed_df['gene'].str.strip()
        else:
            parsed_df['gene'] = ''

        # Add metadata
        parsed_df['tissue_name'] = organ
        parsed_df['tissue_id'] = self._map_and_log_tissue(organ)
        parsed_df['database'] = "HuBMAP"
        
        # Clean up and return
        parsed_df.dropna(subset=['cell_id', 'gene'], inplace=True)
        parsed_df = parsed_df[parsed_df['gene'] != '']
        return parsed_df

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