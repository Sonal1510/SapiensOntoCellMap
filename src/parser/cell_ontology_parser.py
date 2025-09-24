#!/usr/bin/python3

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Parses Cell Ontology (CL) IDs to cell names and vice versa.
                  It reads the CL tsv file, processes it, and provides lookup methods.
"""

import os
import sys
import pandas as pd
from typing import Optional, List, Union

# Add the project root to sys.path to resolve absolute imports from the config directory.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from config.config import (
    CELL_ONTOLOGY_TSV_FILENAME,
    RAW_DATA_DIR
)

class CellOntologyParser:
    """
    A class to parse and provide lookup functionalities for Cell Ontology data.
    It loads the Cell Ontology TSV file, extracts cell IDs and labels,
    and provides methods to query cell names by ID and IDs by cell name.
    """
    def __init__(self):
        """
        Initializes the CellOntologyParser by loading and processing the TSV file.
        It constructs file paths and prepares lookup structures for efficient querying.
        """
        self.cell_ontology_filepath = os.path.join(RAW_DATA_DIR, CELL_ONTOLOGY_TSV_FILENAME)
        self.cell_ontology_df: Optional[pd.DataFrame] = None
        self._id_to_label_map: Optional[dict] = None
        self._label_to_id_map: Optional[dict] = None

        self._load_and_process_data()

    def _load_and_process_data(self):
        """
        Private helper method to load the TSV file and prepare the DataFrame
        and lookup dictionaries. Includes error handling for file operations.
        """
        if not os.path.exists(self.cell_ontology_filepath):
            print(f"Error: Cell Ontology TSV file not found at {self.cell_ontology_filepath}")
            print("Please ensure the file is downloaded and present in the RAW_DATA_DIR.")
            return

        try:
            # Read the TSV file, using the tab delimiter
            self.cell_ontology_df = pd.read_csv(self.cell_ontology_filepath, sep="\t")

            # Strip any leading/trailing whitespace from column names for robust access
            self.cell_ontology_df.columns = self.cell_ontology_df.columns.str.strip()

            # Ensure the required columns exist
            required_columns = ['ID [IRI]', 'LABEL']
            if not all(col in self.cell_ontology_df.columns for col in required_columns):
                print(f"Error: Missing required columns in {self.cell_ontology_filepath}. "
                      f"Expected: {required_columns}")
                self.cell_ontology_df = None # Invalidate the DataFrame
                return

            # Split the 'ID [IRI]' column to get the short ID and create a new 'ID' column
            self.cell_ontology_df['ID'] = self.cell_ontology_df['ID [IRI]'].astype(str).str.split('/').str[-1]

            # Create efficient lookup dictionaries for both directions
            # Handle potential duplicate labels by storing a list of IDs
            self._id_to_label_map = dict(zip(self.cell_ontology_df['ID'], self.cell_ontology_df['LABEL']))

            # For label to ID, group by label to handle one-to-many relationships (multiple IDs for one label)
            self._label_to_id_map = {}
            for index, row in self.cell_ontology_df.iterrows():
                label = row['LABEL']
                cell_id = row['ID']
                if label not in self._label_to_id_map:
                    self._label_to_id_map[label] = []
                self._label_to_id_map[label].append(cell_id)

            print(f"Cell Ontology data loaded successfully from {self.cell_ontology_filepath}")

        except pd.errors.EmptyDataError:
            print(f"Error: Cell Ontology TSV file at {self.cell_ontology_filepath} is empty.")
            self.cell_ontology_df = None
        except FileNotFoundError:
            # This case is already handled by os.path.exists check, but good for robustness
            print(f"Error: Cell Ontology TSV file not found at {self.cell_ontology_filepath}.")
            self.cell_ontology_df = None
        except Exception as e:
            print(f"An unexpected error occurred while loading or processing Cell Ontology data: {e}")
            self.cell_ontology_df = None

    def is_data_loaded(self) -> bool:
        """Checks if the ontology data was successfully loaded."""
        return self.cell_ontology_df is not None

    def get_cell_name_given_cell_id(self, cell_id):
        """
        Retrieves the cell name (LABEL) for a given cell ID (e.g., 'CL_0000000').

        Args:
            cell_id (str): The Cell Ontology ID (e.g., 'CL_0000000').

        Returns:
            Optional[str]: The corresponding cell name, or None if the ID is not found.
        """
        if not self.is_data_loaded():
            print("Warning: Ontology data not loaded. Cannot perform lookup.")
            return None
        return self._id_to_label_map.get(cell_id)

    def get_cell_id_given_cell_name(self, cell_name):
        """
        Retrieves the cell ID(s) for a given cell name (LABEL).
        Note: Cell names may not be unique, so this can return a list of IDs.

        Args:
            cell_name (str): The Cell Ontology label (e.g., 'cell').

        Returns:
            Optional[Union[str, List[str]]]: A string if only one ID is found,
                                              a list of strings if multiple IDs share the same label,
                                              or None if the name is not found.
        """
        if not self.is_data_loaded():
            print("Warning: Ontology data not loaded. Cannot perform lookup.")
            return None

        ids = self._label_to_id_map.get(cell_name)
        if ids:
            return ids[0] if len(ids) == 1 else ids # Return single ID or list of IDs
        return None

    def get_all_cell_ids(self):
        """
        Returns a list of all unique cell IDs present in the loaded ontology.
        """
        if not self.is_data_loaded():
            print("Warning: Ontology data not loaded. Cannot retrieve all IDs.")
            return None
        return self.cell_ontology_df['ID'].drop_duplicates().tolist()

    def get_all_cell_names(self):
        """
        Returns a list of all unique cell names (labels) present in the loaded ontology.
        """
        if not self.is_data_loaded():
            print("Warning: Ontology data not loaded. Cannot retrieve all names.")
            return None
        return self.cell_ontology_df['LABEL'].drop_duplicates().tolist()

    def get_dataframe(self):
        """
        Returns the loaded and processed pandas DataFrame.
        """
        return self.cell_ontology_df
