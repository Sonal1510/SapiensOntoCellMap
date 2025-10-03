import sys
import os
import pandas as pd

# Add the project root to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# --- Corrected Imports ---
# Import the downloader, the main creator/orchestrator, and the config paths
from src.download.bio_database_downloader import BioDataDownloader
from src.db_manager.database_creator import DatabaseCreate
from config.config import (
    PROCESSED_DATA_DIR,
    RECOVER_ID_DATA_DIR,
    PROCESSED_COMBINED_DATA_DIR,
)

def run_downloader_tests():
    """
    Tests the BioDataDownloader class by attempting to download all data.
    """
    print("--- Running BioDataDownloader Tests ---")
    downloader = BioDataDownloader()

    print("\nAttempting to download all databases...")
    # Assuming download_all_databases() returns True on success and handles errors internally
    success = downloader.download_all_databases()
    if success:
        print("🎉 Downloader test: SUCCESS!")
    else:
        print("😞 Downloader test: FAILED!")
    print("--- BioDataDownloader Tests Completed ---\n")

def run_creation_and_parser_tests():
    """
    Tests the entire database creation pipeline managed by the DatabaseCreate class.
    This serves as an integration test for all individual parsers.
    """
    print("--- Running Database Creation and Parser Tests ---")
    
    try:
        # Instantiate and run the main orchestrator
        db_creator = DatabaseCreate()
        db_creator.run()

        # --- Verification Step ---
        print("\n--- Verifying Output Files ---")
        
        # List of key output files we expect to be created
        expected_files = {
            "Combined DB": os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv"),
            "Processed CellMarker": os.path.join(PROCESSED_DATA_DIR, "cellmarkerdb_processed.csv"),
            "Recovery CellMarker": os.path.join(RECOVER_ID_DATA_DIR, "cellmarkerdb_recovery_log.csv"),
            "Processed CellxGene": os.path.join(PROCESSED_DATA_DIR, "cellxgene_processed.csv"),
            "Recovery CellxGene": os.path.join(RECOVER_ID_DATA_DIR, "cellxgene_recovery_log.csv"),
            "Processed HuBMAP": os.path.join(PROCESSED_DATA_DIR, "hubmap_processed.csv"),
            "Recovery HuBMAP": os.path.join(RECOVER_ID_DATA_DIR, "hubmap_recovery_log.csv"),
            "Processed PanglaoDB": os.path.join(PROCESSED_DATA_DIR, "panglaodb_processed.csv"),
            "Recovery PanglaoDB": os.path.join(RECOVER_ID_DATA_DIR, "panglaodb_recovery_log.csv"),
        }
        
        missing_files = []
        for name, path in expected_files.items():
            # Check if file exists and is not empty
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                missing_files.append(f"({name}: {path})")

        if not missing_files:
            print("🎉 Creation and Parser test: SUCCESS! All expected files were created.")
        else:
            print("😞 Creation and Parser test: FAILED! The following files are missing or empty:")
            for mf in missing_files:
                print(f"  - {mf}")

    except Exception as e:
        print(f"😞 Creation and Parser test: FAILED with an unexpected error: {e}")
    
    print("--- Database Creation and Parser Tests Completed ---")


if __name__ == "__main__":
    # First, test the downloader to ensure data is present
    #run_downloader_tests()
    
    # Second, test the entire creation and parsing pipeline
    run_creation_and_parser_tests()