import sys
import os
import pandas as pd

try:
    # This structure is robust for running from different locations.
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.download.bio_database_downloader import BioDataDownloader
    from src.db_manager.database_creator import DatabaseCreate
    from config.config import (
        PROCESSED_DATA_DIR,
        RECOVER_ID_DATA_DIR,
        PROCESSED_COMBINED_DATA_DIR,
        DATABASE_CONFIG # Import the central config
    )

except ImportError:
    print("Error: Could not import from 'config.config'.")
    print("Please ensure that this script is run from within the project structure,")
    print("and that 'config/config.py' exists at the project root.")
    sys.exit(1)


def run_downloader_tests():
    """Tests the BioDataDownloader class."""
    print("--- Running BioDataDownloader Tests ---")
    downloader = BioDataDownloader()
    downloader.download_all_databases() # Run the download process
    print("🎉 Downloader test completed. Check logs for successes or failures.")
    print("--- BioDataDownloader Tests Completed ---\n")

def run_creation_and_parser_tests():
    """Tests the entire database creation pipeline managed by DatabaseCreate."""
    print("--- Running Database Creation and Parser Tests ---")
    
    try:
        db_creator = DatabaseCreate()
        db_creator.run() # Correctly call the run() method

        # --- DYNAMIC Verification Step ---
        print("\n--- Verifying Output Files ---")
        
        expected_files = {}
        # Dynamically generate the list of expected files from the config
        for db_name in DATABASE_CONFIG.keys():
            expected_files[f"Processed {db_name}"] = os.path.join(PROCESSED_DATA_DIR, f"{db_name}_processed.csv")
            expected_files[f"Recovery {db_name}"] = os.path.join(RECOVER_ID_DATA_DIR, f"{db_name}_recovery_log.csv")
        
        # Add the final combined file
        expected_files["Combined DB"] = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
        
        missing_files = []
        for name, path in expected_files.items():
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
    #run_downloader_tests()
    run_creation_and_parser_tests()