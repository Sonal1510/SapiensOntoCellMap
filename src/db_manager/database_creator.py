#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Orchestrates the parsing of all biological databases, combines them, and saves the results.
"""
import os
import sys
import pandas as pd

# Add the project root to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Import all the corrected parsers and config
from src.parser.cellmarkerdb_parser import CellMarkerDBParser
from src.parser.cellxgene_parser import CellxGeneDBParser
from src.parser.hubmap_parser import HuBMapDBParser
from src.parser.panglaodb_parser import PanglaoParser
from src.parser.database_file_parser import DatabaseFileParser
from config.config import (
    RAW_DATA_DIR, 
    PROCESSED_DATA_DIR,
    RECOVER_ID_DATA_DIR,
    PROCESSED_COMBINED_DATA_DIR,
    DATABASE_SOURCE_DICTIONARY
)

class DatabaseCreate:
    """
    A class to orchestrate the processing of multiple biological databases.
    It loads raw data, runs each specific parser, combines the standardized outputs,
    and saves all resulting DataFrames to designated directories.
    """
    def __init__(self):
        """
        Initializes the DatabaseCreate orchestrator.
        """
        print("Initializing database creation process...")
        self.db_df_dict = {}
        self.processed_dfs = {}
        self.recovery_dfs = {}
        self.combined_df = pd.DataFrame()

        # --- 1. Load all raw database files into memory ---
        print("\n--- Step 1: Loading Raw Database Files ---")
        for db_name, (_, file_name, file_type) in DATABASE_SOURCE_DICTIONARY.items():
            file_path = os.path.join(RAW_DATA_DIR, file_name)
            print(f"Loading {db_name} from {file_path}...")
            self.db_df_dict[db_name] = DatabaseFileParser(file_path, file_type).dataframe
        
        # --- 2. Run each specific parser ---
        print("\n--- Step 2: Parsing and Standardizing Databases ---")
        
        # CellMarkerDB
        print("Parsing CellMarkerDB...")
        cellmarker_obj = CellMarkerDBParser(self.db_df_dict['CELLMARKER_DB'])
        self.processed_dfs['cellmarkerdb'] = cellmarker_obj.processed_df
        self.recovery_dfs['cellmarkerdb'] = cellmarker_obj.recovery_df

        # CellxGene
        print("Parsing CellxGene...")
        cellxgene_obj = CellxGeneDBParser(self.db_df_dict['CELLXGENE_DB'])
        self.processed_dfs['cellxgene'] = cellxgene_obj.processed_df
        self.recovery_dfs['cellxgene'] = cellxgene_obj.recovery_df
        
        # HuBMAP
        print("Parsing HuBMAP...")
        hubmap_obj = HuBMapDBParser(self.db_df_dict['HUBMAP_DB'])
        self.processed_dfs['hubmap'] = hubmap_obj.processed_df
        self.recovery_dfs['hubmap'] = hubmap_obj.recovery_df

        # PanglaoDB
        print("Parsing PanglaoDB...")
        panglao_obj = PanglaoParser(self.db_df_dict['PANGLAO_DB'])
        self.processed_dfs['panglaodb'] = panglao_obj.processed_df
        self.recovery_dfs['panglaodb'] = panglao_obj.recovery_df
        
        # --- 3. Combine all processed DataFrames ---
        print("\n--- Step 3: Combining All Parsed DataFrames ---")
        self.combined_df = pd.concat(self.processed_dfs.values(), ignore_index=True)
        print(f"Successfully combined data. Total records: {len(self.combined_df)}")

        # --- 4. Save all DataFrames to their respective folders ---
        print("\n--- Step 4: Saving All Output Files ---")
        
        # Save individual processed files
        for name, df in self.processed_dfs.items():
            output_path = os.path.join(PROCESSED_DATA_DIR, f"{name}_processed.csv")
            df.to_csv(output_path, index=False)
            print(f"Saved processed {name} data to {output_path}")

        # Save individual recovery logs
        for name, df in self.recovery_dfs.items():
            output_path = os.path.join(RECOVER_ID_DATA_DIR, f"{name}_recovery_log.csv")
            df.to_csv(output_path, index=False)
            print(f"Saved {name} recovery log to {output_path}")

        # Save the final combined database
        combined_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
        self.combined_df.to_csv(combined_path, index=False)
        print(f"Saved combined master database to {combined_path}")
        
        print("\nDatabase creation process finished successfully!")


