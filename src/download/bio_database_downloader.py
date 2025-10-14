#!/usr/bin/python3

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Manages the download of biological datasets.
"""

import requests
import os
import sys
try:
    # This structure is robust for running from different locations.
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from config.config import RAW_DATA_DIR, DATABASE_CONFIG
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

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
        # The config file already creates this, but it's safe to have it here too.
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory set to: {self.output_dir}")

    def _download_file(self, url, file_name):
        """
        Private helper method to download a file from a URL to the output directory.
        Includes error handling, streaming for large files, and status code checks.

        Args:
            url (str): The URL of the file to download.
            file_name (str): The name to save the downloaded file as.

        Returns:
            str: The full path to the saved file if successful, None otherwise.
        """
        print(f"Attempting to download from: {url}\nSaving to: {file_name}")
        try:
            # Define headers to mimic a web browser, which can prevent being blocked.
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
            }
            # Use stream=True to handle large files efficiently without loading them into memory.
            with requests.get(url, stream=True, headers=headers, timeout=30) as response:
                # Raise an HTTPError for bad responses (4xx client or 5xx server errors).
                response.raise_for_status()

                # Write content to the file in chunks.
                with open(file_name, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)
                        
            print(f"✅ File downloaded successfully to: {file_name}")
            return file_name
        except requests.exceptions.RequestException as e:
            # Catch all request-related errors (connection, timeout, HTTP errors).
            print(f"❌ Error downloading file from {url}: {e}")
            return None
        except Exception as e:
            # Catch any other unexpected errors during the download process.
            print(f"❌ An unexpected error occurred during download of {url}: {e}")
            return None

    def download_all_databases(self):
        """
        Iterates through the DATABASE_CONFIG and downloads each file.

        Returns:
            list: A list of file paths for all successfully downloaded files.
        """
        downloaded_files = []
        # CORRECTED: Iterate over key-value pairs directly from .items()
        for db_name, config in DATABASE_CONFIG.items():
            # Use .get() for safer access in case 'source' key is ever missing
            source_info = config.get('source')
            if not source_info or len(source_info) != 3:
                print(f"⚠️ Warning: Skipping '{db_name}' due to malformed 'source' configuration.")
                continue

            url, file_name, file_type = source_info
            _file_path = os.path.join(self.output_dir, file_name)
            
            if url:
                print(f"\n--- Downloading {db_name} Data ---")
                file_path = self._download_file(url, _file_path)
                if file_path:
                    downloaded_files.append(file_path)
            else:
                print(f"\n--- Checking for manually downloaded {db_name} Data ---")
                if os.path.exists(_file_path):
                    print(f"--- ✅ {db_name} exists at {_file_path} ---")
                else:
                    print(f"--- ❌ {db_name} not found at {_file_path}. Please download it manually. ---")
        
        print("\n--- Download Summary ---")
        if downloaded_files:
            print(f"Successfully downloaded {len(downloaded_files)} file(s):")
            for path in downloaded_files:
                print(f" - {path}")
        else:
            print("No new files were downloaded.")
            
        return downloaded_files