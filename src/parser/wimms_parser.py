#!/usr/bin/python3
"""
Author         : Sonal Rashmi (expert review by Gemini)
Date           : 13/10/2025
Description    : Parses WIMMS database and uses BaseParser for normalization.
"""
import pandas as pd
import sys
import os

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.parser.base_parser import BaseParser
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    print("Please ensure that all dependencies are installed and the script is run from the project's root directory.")
    sys.exit(1)

class WimmsMelanocyteParser:
    """
    A class to parse the WIMMS database for melanocyte markers.
    """
    def __init__(self, df: pd.DataFrame):
        if not df.empty:
            if 'Unnamed: 0' in df.columns:
                df = df.drop(columns=['Unnamed: 0'])

            # --- Mapping from study column to cell state based on the provided image ---
            cell_state_mapping = study_to_category = {
                # AXL
                "Tirosh_2016_AXL": "AXL",
                "Ryu_2011_BRAFV600E_Targets": "AXL",
                "Riesenberg_2015_TNF_Response": "AXL",

                # Neuro
                "Tsoi_2018_Neural_Crest": "Neuro",
                "Rambow_2018_Neuro": "Neuro",
                "Wouters_2020_Intermediate": "Neuro",

                # Invasive
                "Tsoi_2018_Undifferentiated": "Invasive",
                "McNeal_2021_EDN1": "Invasive",
                "Rambow_2018_Invasion": "Invasive",
                "Belote_2021_ADT": "Invasive",
                "Belote_2021_MSC": "Invasive",
                "Widmer_2012_Invasive": "Invasive",
                "Hoek_2006_Invasive": "Invasive",
                "Wouters_2020_Mesenchymal": "Invasive",
                "Andrews_2022_MES": "Invasive",
                "Andrews_2022_NPLAS_Up": "Invasive",
                "Verfaillie_2015_TEAD_Targets": "Invasive",
                "Verfaillie_2015_Invasive": "Invasive",
                "Jeffs_2009_Invasive": "Invasive",

                # Amelanotic
                "Tsoi_2018_Transitory": "Amelanotic",
                "Belote_2021_Vmel": "Amelanotic",

                # Differentiated
                "Belote_2021_Cmel": "Differentiated",
                "Belote_2021_ADT": "Differentiated",
                "Rambow_2018_Pigmentation": "Differentiated",
                "Hoek_2006_Proliferative": "Differentiated",
                "Andrews_2022_MEL": "Differentiated",
                "Verfaillie_2015_Proliferative": "Differentiated",
                "Wouters_2020_Melanocytic": "Differentiated",
                "McNeal_2021_PMA": "Differentiated",
                "Tirosh_2016_MITF": "Differentiated",
                "Rambow_2018_MITF_Targets": "Differentiated",
                "Hoek_2008_MITF_Targets": "Differentiated",
                "Jeffs_2009_Proliferation": "Differentiated",
                "Widmer_2012_Proliferative": "Differentiated",

                # Hypometabolic
                "Rambow_2018_Hypometabolic": "Hypometabolic",

                # Mitotic/MYC
                "Rambow_2018_Mitosis": "Mitotic/MYC",
                "Kauffmann_2008_DNA_Repair": "Mitotic/MYC",
                "Kauffmann_2008_DNA_Replication": "Mitotic/MYC",
            }


            # --- Step 1: Transform from wide to long format ---
            rows = []
            base_cell_name = "melanocyte" # Changed to lowercase for the requested format

            for col_name in df.columns:
                # Look up the cell state from the mapping dictionary
                cell_state = cell_state_mapping.get(col_name)

                # Only process columns that have a corresponding cell state in the mapping
                if cell_state:
                    genes = df[col_name].dropna().tolist()
                    for gene in genes:
                        rows.append({
                            'db_tissue_name': "Skin",
                            'db_cell_name': f"{base_cell_name}_{cell_state}", # Generate name using the mapped cell state
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