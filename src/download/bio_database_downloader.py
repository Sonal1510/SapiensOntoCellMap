#!/usr/bin/python3

"""
Author         : Sonal Rashmi (expert review by Gemini)
Date           : 15/08/2025
Description    : Manages the download of biological datasets, skipping existing files.
"""

import requests
import os
import sys
try:
    # This structure is robust for running from different locations.
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from config.config import (
        RAW_DATA_DIR, DATABASE_CONFIG, REFERENCE_DATA_DIR,
        HGNC_COMPLETE_SET_URL, HGNC_COMPLETE_SET_FILE,
        MSIGDB_G2M_URL, MSIGDB_G2M_FILE,
        MSIGDB_E2F_URL, MSIGDB_E2F_FILE,
    )
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

class BioDataDownloader:
    """
    A class to manage the download and organization of biological datasets.
    Skips downloading files that already exist.
    """
    def __init__(self):
        """
        Initializes the BioDataDownloader with the specified output directory
        and ensures the directory exists.
        """
        self.output_dir = RAW_DATA_DIR
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory set to: {self.output_dir}")
        self.skipped_files = [] # Keep track of skipped files

    def _download_file(self, url, file_path): # Renamed file_name to file_path for clarity
        """
        Private helper method to download a file from a URL to the output directory.
        Checks if the file already exists before downloading.
        Includes error handling, streaming for large files, and status code checks.

        Args:
            url (str): The URL of the file to download.
            file_path (str): The full path where the file should be saved.

        Returns:
            str: The full path to the file (either existing or newly downloaded) if available, None otherwise.
        """
        # --- Check if file already exists ---
        if os.path.exists(file_path):
            print(f"⚠️ File already exists at: {file_path}. Skipping download.")
            self.skipped_files.append(file_path) # Add to skipped list
            return file_path # Return the existing path
        # --- End Check ---

        print(f"Attempting to download from: {url}\nSaving to: {file_path}")
        try:
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            with requests.get(url, stream=True, headers=headers, timeout=60) as response: # Increased timeout slightly
                response.raise_for_status()

                # Write content to the file in chunks.
                with open(file_path, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)

                print(f"✅ File downloaded successfully to: {file_path}")
                return file_path
        except requests.exceptions.Timeout:
             print(f"❌ Timeout error downloading file from {url}")
             return None
        except requests.exceptions.RequestException as e:
            print(f"❌ Error downloading file from {url}: {e}")
            # Optionally remove partially downloaded file
            if os.path.exists(file_path):
                 try:
                     os.remove(file_path)
                     print(f"🧹 Removed partially downloaded file: {file_path}")
                 except OSError as oe:
                      print(f"❌ Error removing partial file {file_path}: {oe}")
            return None
        except Exception as e:
            print(f"❌ An unexpected error occurred during download of {url}: {e}")
            return None

    def download_all_databases(self):
        """
        Iterates through the DATABASE_CONFIG and downloads each file if it doesn't exist.

        Returns:
            list: A list of file paths for all available files (existing or newly downloaded).
        """
        available_files = []
        self.skipped_files = [] # Reset skipped list for this run

        for db_name, config in DATABASE_CONFIG.items():
            source_info = config.get('source')
            if not source_info or len(source_info) != 3:
                print(f"⚠️ Warning: Skipping '{db_name}' due to malformed 'source' configuration.")
                continue

            url, file_name, file_type = source_info
            _file_path = os.path.join(self.output_dir, file_name)

            if url: # If URL is provided, attempt download (which includes the existence check)
                print(f"\n--- Processing {db_name} Data ---")
                file_path_result = self._download_file(url, _file_path)
                if file_path_result:
                    available_files.append(file_path_result)
            else: # If no URL, just check if the file exists manually
                print(f"\n--- Checking for manually placed {db_name} Data ---")
                if os.path.exists(_file_path):
                    print(f"✅ File exists at: {_file_path}")
                    available_files.append(_file_path)
                    # No need to add to skipped_files here, as it wasn't an attempted download
                else:
                    print(f"❌ File not found at {_file_path}. Please place it manually.")

        print("\n--- Processing Summary ---")
        newly_downloaded_count = len(available_files) - len(self.skipped_files)

        if newly_downloaded_count > 0:
             print(f"Successfully downloaded {newly_downloaded_count} new file(s).")
        if self.skipped_files:
             print(f"Skipped download for {len(self.skipped_files)} existing file(s).")
        if not available_files:
             print("No database files are available in the target directory.")
        elif newly_downloaded_count == 0 and not self.skipped_files:
             # This case covers when only manually placed files were found
             print(f"Found {len(available_files)} manually placed file(s). No downloads attempted.")


        print(f"\nTotal available files: {len(available_files)}")
        # Optionally list all available files
        # for path in available_files:
        #      print(f" - {os.path.basename(path)}") # Just show filename

        return available_files

    def download_reference_data(self):
        """
        Downloads reference data files (e.g., HGNC gene alias map) needed by
        the enrichment pipeline. These are separate from the marker databases
        and stored in data/reference/.

        Returns:
            dict: {name: file_path} for all available reference files.
        """
        os.makedirs(REFERENCE_DATA_DIR, exist_ok=True)
        available = {}

        # --- HGNC Complete Set (gene alias resolution) ---
        print("\n--- Processing HGNC Complete Set (gene alias map) ---")
        result = self._download_file(HGNC_COMPLETE_SET_URL, HGNC_COMPLETE_SET_FILE)
        if result:
            available['hgnc_complete_set'] = result
        else:
            print("❌ HGNC download failed. Gene alias resolution will be unavailable.")

        # --- MSigDB Hallmark cell-cycle gene sets (Proliferative_Flag) ---
        # Liberzon et al., Cell Systems 2015; Subramanian et al., PNAS 2005
        # G2M: HALLMARK_G2M_CHECKPOINT (200 genes, G2/M transition)
        # E2F: HALLMARK_E2F_TARGETS    (200 genes, S-phase / cell-cycle entry)
        print("\n--- Processing MSigDB HALLMARK_G2M_CHECKPOINT ---")
        result = self._download_file(MSIGDB_G2M_URL, MSIGDB_G2M_FILE)
        if result:
            available['msigdb_g2m'] = result
        else:
            print("❌ MSigDB G2M download failed. Proliferative_Flag will not be computed.")

        print("\n--- Processing MSigDB HALLMARK_E2F_TARGETS ---")
        result = self._download_file(MSIGDB_E2F_URL, MSIGDB_E2F_FILE)
        if result:
            available['msigdb_e2f'] = result
        else:
            print("❌ MSigDB E2F download failed. Proliferative_Flag will not be computed.")

        print(f"\nTotal reference files available: {len(available)}")
        return available