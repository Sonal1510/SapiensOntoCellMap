#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Orchestrates the parsing of all biological databases, combines them, and saves the results.
"""
import os
import sys
import pandas as pd

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    
    # --- IMPORTS MOVED HERE ---
    from src.parser.database_file_parser import DatabaseFileParser
    from config.config import (
        RAW_DATA_DIR, 
        PROCESSED_DATA_DIR,
        RECOVER_ID_DATA_DIR,
        PROCESSED_COMBINED_DATA_DIR,
        DATABASE_CONFIG
    )
    # Import all parser classes that will be used
    from src.parser.cellmarkerdb_parser import CellMarkerDBParser
    from src.parser.hubmap_parser import HuBMapDBParser
    from src.parser.cellxgene_parser import CellxGeneDBParser
    from src.parser.panglaodb_parser import PanglaoParser
    from src.parser.wimms_parser import WimmsMelanocyteParser
    from src.parser.generic_parser import GenericFileParser
    from src.parser.human_scc_cell_2020_parser import HumanSccCell2020Parser
    
except ImportError as e:
    print(f"❌ A critical import error occurred in database_creator.py: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

# --- NEW: Map string keys from config to actual Parser Classes ---
PARSER_MAPPING = {
    "cellmarker": CellMarkerDBParser,
    "hubmap": HuBMapDBParser,
    "cellxgene": CellxGeneDBParser,
    "panglao": PanglaoParser,
    "wimms": WimmsMelanocyteParser,
    "generic": GenericFileParser,
    "human_scc_cell_2020": HumanSccCell2020Parser
}

class DatabaseCreate:
    """
    A class to dynamically orchestrate the processing of multiple biological databases
    based on the central configuration file.
    """
    def __init__(self):
        self.db_df_dict = {}
        self.processed_dfs = {}
        self.recovery_dfs = {}
        self.combined_df = pd.DataFrame()

    def run(self):
        """Executes the entire database creation pipeline."""
        print("Initializing database creation process...")
        
        # --- Step 1: Loading Raw Database Files (No changes needed) ---
        print("\n--- Step 1: Loading Raw Database Files ---")
        for db_name, config in DATABASE_CONFIG.items():
            _, file_name, file_type = config["source"]
            file_path = os.path.join(RAW_DATA_DIR, file_name)
            if os.path.exists(file_path):
                print(f"Loading {db_name} from {file_path}...")
                self.db_df_dict[db_name] = DatabaseFileParser(file_path, file_type).dataframe
            else:
                print(f"⚠️ Warning: File not found for {db_name} at {file_path}. Skipping.")

        # --- Step 2: Dynamically run each parser based on config ---
        print("\n--- Step 2: Parsing and Standardizing Databases ---")
        for db_name, config in DATABASE_CONFIG.items():
            if db_name not in self.db_df_dict or self.db_df_dict[db_name] is None:
                continue

            print(f"Parsing {db_name}...")
            parser_key = config["parser_key"]
            parser_config = config["parser_config"]
            
            if parser_key not in PARSER_MAPPING:
                print(f"❌ Error: Parser key '{parser_key}' for database '{db_name}' not found in PARSER_MAPPING.")
                continue
            
            parser_class = PARSER_MAPPING[parser_key]
            
            # Instantiate the parser class with its specific configuration
            parser_obj = parser_class(df=self.db_df_dict[db_name], **parser_config)

            self.processed_dfs[db_name] = parser_obj.processed_df
            self.recovery_dfs[db_name] = parser_obj.recovery_df

        # --- Steps 3 & 4 (No changes needed) ---
        # ... (Combine and Save logic remains the same) ...
        print("\n--- Step 3: Combining All Parsed DataFrames ---")
        if self.processed_dfs:
            self.combined_df = pd.concat(self.processed_dfs.values(), ignore_index=True)
            print(f"Successfully combined data. Total records: {len(self.combined_df)}")
        else:
            print("No data was processed. Combined DataFrame is empty.")


        print("\n--- Step 4: Saving All Output Files ---")
        for name, df in self.processed_dfs.items():
             if df is not None and not df.empty:
                output_path = os.path.join(PROCESSED_DATA_DIR, f"{name}_processed.csv")
                df.to_csv(output_path, index=False)
                print(f"Saved processed {name} data to {output_path}")

        for name, df in self.recovery_dfs.items():
            if df is not None and not df.empty:
                output_path = os.path.join(RECOVER_ID_DATA_DIR, f"{name}_recovery_log.csv")
                df.to_csv(output_path, index=False)
                print(f"Saved {name} recovery log to {output_path}")

        if not self.combined_df.empty:
            combined_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
            self.combined_df.to_csv(combined_path, index=False)
            print(f"Saved combined master database to {combined_path}")
        
        print("\nDatabase creation process finished successfully!")