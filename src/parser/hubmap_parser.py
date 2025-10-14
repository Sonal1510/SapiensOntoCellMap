#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses HuBMAP database and uses BaseParser for normalization.
"""
import os
import sys
import pandas as pd
from typing import Optional

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.download.bio_database_downloader import BioDataDownloader
    from src.parser.base_parser import BaseParser
    from config.config import RAW_DATA_DIR

except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

class HuBMapDBParser:
    """
    A class to parse the HuBMAP database, preparing data for the BaseParser.
    """
    def __init__(self, df: pd.DataFrame):
        if not isinstance(df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        self.downloader = BioDataDownloader()
        hubmap_organ_df_list = {}

        for organ in set(df['Organ'].tolist()):
            if organ != "anatomical systems":
                # ... (file downloading logic remains the same) ...
                url_file = df.loc[df['Organ'] == organ, 'csv'].iloc[0]
                file_name = f"hubmap_{organ.replace(' ', '_')}.csv"
                file_path = os.path.join(RAW_DATA_DIR, file_name)
                self.downloader._download_file(url_file, file_path)
                
                try:
                    temp_df = pd.read_csv(file_path, header=10, on_bad_lines='skip', engine='python')
                    if not temp_df.empty:
                        parsed_df = self._process_organ_dataframe(temp_df, organ)
                        hubmap_organ_df_list[organ] = parsed_df
                except Exception as e:
                    print(f"❌ Error while parsing {file_name}: {e}")

        if hubmap_organ_df_list:
            hubmap_all_organs_df = pd.concat(hubmap_organ_df_list.values(), ignore_index=True)
            
            # Use the BaseParser for normalization
            normalizer = BaseParser()
            self.processed_df, self.recovery_df = normalizer.normalize_dataframe(hubmap_all_organs_df)
        else:
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()

    def _process_organ_dataframe(self, current_df: pd.DataFrame, organ: str) -> pd.DataFrame:
        """Helper function to extract data into the intermediate 'db_' schema."""
        parsed_df = pd.DataFrame()

        # Determine cell/anatomical structure columns and map to 'db_' schema
        if any(col.startswith("CT/") for col in current_df.columns):
            parsed_df['db_cell_id'] = current_df.get('CT/1/ID')
            parsed_df['db_cell_name'] = current_df.get('CT/1/LABEL')
        else:
            parsed_df['db_cell_id'] = current_df.get('AS/1/ID')
            parsed_df['db_cell_name'] = current_df.get('AS/1/LABEL')

        # Map tissue info
        parsed_df['db_tissue_name'] = organ
        parsed_df['db_tissue_id'] = None # HuBMAP does not provide tissue IDs directly

        # Gather all reference (REF/*) columns into source_info
        ref_cols = [col for col in current_df.columns if col.startswith("REF/")]
        ref_info = current_df[ref_cols].fillna('').astype(str)
        parsed_df['source_info'] = ref_info.apply(lambda row: "; ".join(val for val in row if val), axis=1)
        parsed_df['source_type'] = "Literature"
        parsed_df['database'] = "HuBMAP"

        # Gather and explode gene markers
        gene_marker_cols = [c for c in current_df.columns if ('BGene' in c or 'BProtein' in c) and 'LABEL' in c]
        if gene_marker_cols:
            parsed_df['gene'] = current_df[gene_marker_cols].apply(
                lambda row: ', '.join(row.dropna().astype(str)), axis=1
            ).str.split(', ')
            parsed_df = parsed_df.explode('gene')
            parsed_df['gene'] = parsed_df['gene'].str.strip()
        else:
            parsed_df['gene'] = ''

        # Clean up and return
        parsed_df.dropna(subset=['db_cell_id', 'gene'], inplace=True)
        parsed_df = parsed_df[parsed_df['gene'] != '']
        return parsed_df