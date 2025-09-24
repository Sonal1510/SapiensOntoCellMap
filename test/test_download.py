import sys
import os

# Add the project root to sys.path to allow absolute imports
# This assumes test_downloader.py is in the 'tests' directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now import the class from your src package
from src.download.bio_database_downloader import BioDataDownloader

def run_downloader_tests():
    """
    A simple function to test the BioDataDownloader class.
    This will attempt to download data and print success/failure messages.
    """
    print("--- Running BioDataDownloader Tests ---")
    downloader = BioDataDownloader()

    print("\nAttempting to download Uberon data...")
    uberon_path = downloader.download_uberon_data()
    if uberon_path:
        print(f"🎉 Uberon data test: SUCCESS! Saved to: {uberon_path}")
    else:
        print("😞 Uberon data test: FAILED!")

    print("\nAttempting to download CellMarker data...")
    cellmarker_path = downloader.download_cellmarker_data()
    if cellmarker_path:
        print(f"🎉 CellMarker data test: SUCCESS! Saved to: {cellmarker_path}")
    else:
        print("😞 CellMarker data test: FAILED!")

    print("\nAttempting to download Cell Ontology data...")
    cell_ontology_path = downloader.download_cell_ontology_data()
    if cell_ontology_path:
        print(f"🎉 Cell Ontology data test: SUCCESS! Saved to: {cell_ontology_path}")
    else:
        print("😞 Cell Ontology data test: FAILED!")

    print("\n--- BioDataDownloader Tests Completed ---")

if __name__ == "__main__":
    run_downloader_tests()