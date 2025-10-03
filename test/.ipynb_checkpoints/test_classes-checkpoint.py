import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# --- FIX 1: Import both the downloader and the parser classes ---
from src.download.bio_database_downloader import BioDataDownloader
from src.parser.bio_database_parser import BioCellularDatabaseParser # <-- Added this import

def run_downloader_tests():
    """
    Tests the BioDataDownloader class by attempting to download all data.
    """
    print("--- Running BioDataDownloader Tests ---")
    downloader = BioDataDownloader()

    print("\nAttempting to download all databases...")
    # Assuming download_all_databases() returns True on success
    success = downloader.download_all_databases() 
    if success:
        print("🎉 Downloader test: SUCCESS!")
    else:
        print("😞 Downloader test: FAILED!")
    print("--- BioDataDownloader Tests Completed ---")

def run_parser_tests():
    """
    Tests the BioCellularDatabaseParser class by attempting to parse all data.
    """
    print("\n--- Running BioCellularDatabaseParser Tests ---")
    db_parser = BioCellularDatabaseParser()

    print("Attempting to parse all data sources...")
    db_df_dict = db_parser.parse_all_sources()
    
    # This is a great way to check for success!
    if len(db_parser.database_source_dict) == len(db_df_dict) and len(db_df_dict) > 0:
        print("🎉 Parser test: SUCCESS! All sources parsed correctly.")
        for name, df in db_df_dict.items():
            print(f"  - '{name}': DataFrame with shape {df.shape}")
    else:
        print("😞 Parser test: FAILED! Not all sources were parsed.")
    print("--- BioCellularDatabaseParser Tests Completed ---")


if __name__ == "__main__":
    run_downloader_tests()
    run_parser_tests()