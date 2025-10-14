#!/usr.bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 13/10/2025
Description     : Parses WIMMS database and uses BaseParser for normalization.
"""
import pandas as pd
from src.parser.base_parser import BaseParser

class WimmsMelanocyteParser: # Class name corrected for consistency
    """
    A class to parse the WIMMS database for melanocyte markers.
    """
    def __init__(self, wimms_df: pd.DataFrame):
        if not wimms_df.empty:
            df = wimms_df.copy()
            
            if 'Unnamed: 0' in df.columns:
                df = df.drop(columns=['Unnamed: 0'])

            # --- Step 1: Transform from wide to long format ---
            rows = []
            base_cell_name = "Melanocyte"
            
            for col_name in df.columns:
                genes = df[col_name].dropna().tolist()
                subtype = '_'.join(col_name.split("_")[2:])
                
                for gene in genes:
                    rows.append({
                        'db_tissue_name': "Skin",
                        'db_cell_name': f"{base_cell_name}_{subtype}",
                        'gene': gene,
                        'source_info': '_'.join(col_name.split("_")[:2]),
                    })
            
            long_df = pd.DataFrame(rows)
            
            # --- Step 2: Add metadata ---
            long_df["db_tissue_id"] = None
            long_df["db_cell_id"] = None
            long_df["source_type"] = "literature"
            long_df["database"] = "wimms"

            # --- Step 3: Use the BaseParser for normalization ---
            normalizer = BaseParser()
            self.processed_df, self.recovery_df = normalizer.normalize_dataframe(long_df)
        else:
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()