#!/usr/bin/python3

"""
Author         : Sonal Rashmi (expert review by Gemini)
Date           : 15/08/2025
Description    : Manages the parsing of biological datasets.
"""

import sys
import os
import gzip
import json
import csv
import zipfile
import pandas as pd

# This try/except block is good for robust importing.
try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from config.config import RAW_DATA_DIR, DATABASE_SOURCE_DICTIONARY
except ImportError:
    print("Error: Could not import from 'config.config'.")
    print("Please ensure that this script is run from within the project structure,")
    print("and that 'config/config.py' exists at the project root.")
    sys.exit(1)


class BioCellularDatabaseParser:
    """
    A class to parse downloaded biological datasets from the raw data folder.
    """
    def __init__(self):
        """
        Initializes the parser. The __init__ is now lightweight.
        """
        self.raw_dir = RAW_DATA_DIR
        self.database_source_dict = DATABASE_SOURCE_DICTIONARY
        self.database_df_dict = {}

    def parse_all_sources(self):
        """
        Main execution method to iterate through, check, and parse all database files.
        """
        print("Starting database parsing process...")
        for db_name, (_, _filename, _file_type) in self.database_source_dict.items():
            db_file_path = os.path.join(self.raw_dir, _filename)

            # --- Centralized check for file existence ---
            if not os.path.exists(db_file_path):
                print(f"Warning: File not found for '{db_name}' at '{db_file_path}'. Skipping.")
                continue

            print(f"\n--- Processing: {db_name} ---")
            
            if file_type:
                db_df = self._parse_file(db_file_path, _file_type)
                if db_df is not None:
                    self.database_df_dict[db_name] = db_df
            else:
                print(f"Could not determine file type for {db_file_path}. Skipping.")
        
        print("\n--- Parsing complete ---")
        return self.database_df_dict

    def _parse_file(self, db_file, content_type):
        """
        Parses a single file into a Pandas DataFrame based on its determined file type.
        (Renamed, corrected, and now an internal method)
        """
        df = None
        try:
            if content_type == 'xlsx':
                df = pd.read_excel(db_file)
            elif content_type == 'csv':
                df = pd.read_csv(db_file)
            elif content_type == 'tsv':
                df = pd.read_csv(db_file, sep='\t')
            elif content_type == 'json':
                try:
                    df = pd.read_json(db_file)
                except ValueError as e:
                    if "Trailing data" in str(e):
                        print("Standard JSON parsing failed, trying line-delimited JSON (jsonl)...")
                        df = pd.read_json(db_file, lines=True)
            
            if df is not None:
                print(f"Successfully parsed into a DataFrame with shape: {df.shape}")
                print("First 5 rows:")
                print(df.head())
            else:
                print(f"Parser not implemented for content type '{content_type}' or file is empty.")

        except Exception as e:
            # Fixed variable name from 'filename' to 'db_file'
            print(f"An error occurred while parsing {os.path.basename(db_file)}: {e}")
            return None # Return None on failure
            
        return df

