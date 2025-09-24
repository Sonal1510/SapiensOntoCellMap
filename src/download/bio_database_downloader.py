#!/usr/bin/python3

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Manages the download of biological datasets.
"""

import requests
import os
import gzip # Not currently used, but imported. If needed in future, it's here.
import shutil # Not currently used, but imported. If needed in future, it's here.
import sys

# Add the project root to the sys.path
# This assumes the script is run from the project root or the level above src.
# os.path.dirname(__file__) is 'src/download'
# os.path.join(..., '..') moves to 'src'
# os.path.join(..., '..') moves to the project root '.'
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Import specific constants from the config file
# This import now works because the project root has been added to sys.path
from config.config import (
    UBERON_HUMAN_VIEW_URL, UBERON_OBO_FILENAME,
    CELL_ONTOLOGY_HUMAN_VIEW_URL, CELL_ONTOLOGY_TSV_FILENAME,
    CELLMARKER_HUMAN_EXCEL_URL, CELLMARKER_EXCEL_FILENAME,
    RAW_DATA_DIR
)

class BioDataDownloader:
    """
    A class to manage the download and organization of biological datasets.
    """
    def __init__(self):
        """
        Initializes the BioDataDownloader with the specified output directory
        and ensures the directory exists.
        """
        self.output_dir = RAW_DATA_DIR
        # Ensure the raw data directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory set to: {self.output_dir}")

    def _download_file(self, url, file_name):
        """
        Private helper method to download a file from a URL to the output directory.
        Includes error handling and checks for HTTP status codes.

        Args:
            url (str): The URL of the file to download.
            file_name (str): The name to save the downloaded file as.

        Returns:
            str: The full path to the saved file if successful, None otherwise.
        """
        save_path = os.path.join(self.output_dir, file_name)
        print(f"Attempting to download from: {url}\nSaving to: {save_path}")
        try:
            # Define headers to mimic a web browser, often helps with preventing blocks
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            # Use stream=True to handle large files efficiently
            with requests.get(url, stream=True, headers=headers, timeout=10) as response:
                # Raise an HTTPError for bad responses (4xx or 5xx), effectively checking the URL.
                response.raise_for_status()

                # Write content in chunks to the file
                with open(save_path, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk: # Filter out keep-alive new chunks
                            file.write(chunk)
            print(f"✅ File downloaded successfully to: {save_path}")
            return save_path
        except requests.exceptions.RequestException as e:
            # Catch all request-related errors (connection, timeout, HTTP errors)
            print(f"❌ Error downloading file from {url}: {e}")
            return None
        except Exception as e:
            # Catch any other unexpected errors during the download process
            print(f"❌ An unexpected error occurred during download of {url}: {e}")
            return None

    def download_uberon_data(self):
        """
        Downloads the Uberon human-view OBO file.
        """
        print("\n--- Downloading Uberon Data ---")
        return self._download_file(UBERON_HUMAN_VIEW_URL, UBERON_OBO_FILENAME)

    def download_cell_ontology_data(self):
        """
        Downloads the Cell Ontology human-view TSV file.
        """
        print("\n--- Downloading Cell Ontology Data ---")
        return self._download_file(CELL_ONTOLOGY_HUMAN_VIEW_URL, CELL_ONTOLOGY_TSV_FILENAME)

    def download_cellmarker_data(self):
        """
        Downloads the CellMarker human Excel file.
        """
        print("\n--- Downloading CellMarker Data ---")
        return self._download_file(CELLMARKER_HUMAN_EXCEL_URL, CELLMARKER_EXCEL_FILENAME)
        
