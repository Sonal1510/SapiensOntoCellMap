#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Parses cellmarker database and uses BaseParser for normalization.
"""
import pandas as pd
import sys
import os
try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.parser.base_parser import BaseParser # Import the new base class
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)


class CellMarkerDBParser:
    """
    A class to parse the Cell Marker database. It prepares the data
    and then uses the BaseParser for ontology mapping and normalization.
    """
    def __init__(self, df: pd.DataFrame):
        """
        Initializes the parser by preparing the data and invoking the BaseParser.

        Args:
            cellmarkerdb_df (pd.DataFrame): The input DataFrame from the CellMarker database.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Prepare the DataFrame with 'db_' prefixed columns ---
        
        # Rename columns to the intermediate schema
        df.rename(columns={
            'tissue_class': 'db_tissue_name',
            'uberonongology_id': 'db_tissue_id',
            'cell_name': 'db_cell_name',
            'cellontology_id': 'db_cell_id',
            'marker': 'gene'
        }, inplace=True)
        
        # --- Step 2: Handle source_info and other metadata ---
        df['source_type'] = df['marker_source']
        df['source_info'] = (
            df['Title'].astype(str) + ";" +
            df['journal'].astype(str) + ";" +
            df['year'].astype(str)
        )
        df["database"] = "cellmarkerdb"
        
        # Keep cancer_type for the BaseParser to handle the suffix
        df['cancer_type'] = df.get('cancer_type', pd.Series(index=df.index, dtype=str))


        # --- Step 3: Use the BaseParser for normalization ---
        normalizer = BaseParser()
        # The normalize_dataframe method now returns both the processed_df and recovery_df
        self.processed_df, self.recovery_df = normalizer.normalize_dataframe(df)