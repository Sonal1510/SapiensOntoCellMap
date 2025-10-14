#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 10/01/2025
Description     : Parses Panglao database and uses BaseParser for normalization.
"""
import pandas as pd
from src.parser.base_parser import BaseParser # Import the new base class

class PanglaoParser:
    """
    A class to parse the Panglao database, preparing the data for the BaseParser.
    """
    def __init__(self, panglaodb_df: pd.DataFrame):
        """
        Initializes the parser by preparing the data and invoking the BaseParser.

        Args:
            panglaodb_df (pd.DataFrame): The input DataFrame from the Panglao database.
        """
        if not isinstance(panglaodb_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Step 1: Filter and Prepare the DataFrame ---
        df = panglaodb_df[panglaodb_df['species'].isin(['Hs', 'Mm Hs'])].copy()

        # Rename columns to the 'db_' intermediate schema
        df.rename(columns={
            'organ': 'db_tissue_name',
            'cell type': 'db_cell_name',
            'official gene symbol': 'gene'
        }, inplace=True)
        
        # PanglaoDB doesn't provide IDs, so these columns will be None
        df['db_tissue_id'] = None
        df['db_cell_id'] = None

        # --- Step 2: Format Source Metadata ---
        df['source_type'] = "Literature"
        df['source_info'] = df.apply(
            lambda row: (
                f"canonical_marker_flag={row.get('canonical marker', '')}"
                f"|human_sensitivity={row.get('sensitivity_human', '')}"
                f"|human_specificity={row.get('specificity_human', '')}"
            ),
            axis=1
        )
        df['database'] = "PanglaoDB"

        # --- Step 3: Use the BaseParser for normalization ---
        normalizer = BaseParser()
        self.processed_df, self.recovery_df = normalizer.normalize_dataframe(df)