#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 21/11/2025
Description     : Parses the Skin Atlas (Excel multi-sheet) database.
                  Splits DEG and Cluster annotation columns, expands abbreviations,
                  and uses BaseParser for ontology normalization.
"""
import pandas as pd
import re
import sys
import os
import numpy as np

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)
    from src.parser.base_parser import BaseParser
except ImportError as e:
    print(f"❌ A critical import error occurred: {e}")
    sys.exit(1)


class SkinAtlasParser:
    """
    Parses multi-sheet Excel files from the Skin Atlas.
    Expects sheets to have a split layout (DEG data | Unnamed Column | Cluster Metadata).
    """

    def __init__(self, df: dict, **kwargs):
        """
        Args:
            df (dict): Dictionary {sheet_name: pd.DataFrame} loaded via pd.read_excel(..., sheet_name=None).
                       (Named 'df' to match the call signature in database_creator.py)
            **kwargs: Capture additional config arguments (e.g., database_name, tissue_name).
        """
        self.processed_df = pd.DataFrame()
        self.recovery_df = pd.DataFrame()
        
        # Extract config with defaults
        self.database_name = kwargs.get("database_name", "SkinAtlas")
        self.tissue_name = kwargs.get("tissue_name", "Skin")

        if not df or not isinstance(df, dict):
            print("❌ Input 'df' must be a dictionary of DataFrames (read with sheet_name=None).")
            return

        # Abbreviation map for cell type expansion
        self.celltype_abbreviation_map = {
            'Adipo': 'Adipocyte', 'Bas': 'Basal', 'Cyc': 'Cycling', 'DC': 'Dendritic',
            'CellDiff': 'Differentiated', 'DP': 'Dermal Papilla', 'DS': 'Dermal Sheath',
            'Ecc': 'Eccrine', 'EC': 'Endothelial Cell', 'Fib': 'Fibroblast',
            'HEC': 'High Endothelial Venules', 'HFE': 'Hair Follicle Epithelia',
            'HS': 'Hair Shaft', 'Imm': 'Immune', 'Inf': 'Infindibulum',
            'IRS': 'Inner Root Sheath', 'LC': 'Langerhans Cell',
            'LEC': 'Lymphatic Endothelial Cell', 'Lym': 'Lymphocyte',
            'Mac': 'Macrophage', 'Melano': 'Melanocyte', 'ORS': 'Outer Root Sheath',
            'Papil': 'Papillary', 'Peri': 'Pericyte', 'Perivasc': 'Perivascular',
            'Retic': 'Reticular', 'Seb': 'Sebaceous', 'SM': 'Smooth Muscle',
            'Spn': 'Spinous', 'Tc': 'T Cytotoxic', 'Th': 'T Helper',
            'Treg': 'T Regulatory', 'VEC': 'Vascular Endothelial Cell',
            'KC': 'Keratinocyte', 'Mono': 'Monocyte'
        }

        processed_dfs_list = []

        for sheet_name, sheet_df in df.items():
            try:
                print(f"\n🔍 Processing sheet: {sheet_name}")
                cleaned_df = self._process_single_sheet(sheet_df, sheet_name)
                if not cleaned_df.empty:
                    processed_dfs_list.append(cleaned_df)
            except Exception as e:
                print(f"⚠️ Error processing sheet '{sheet_name}': {e}")

        if processed_dfs_list:
            combined_df = pd.concat(processed_dfs_list, ignore_index=True)
            
            # --- Use BaseParser for Normalization ---
            print(f"--- Running Ontology Normalization for {self.database_name} ---")
            normalizer = BaseParser()
            self.processed_df, self.recovery_df = normalizer.normalize_dataframe(combined_df)
        else:
            print("❌ No data could be parsed from any sheet.")

    def _process_single_sheet(self, df: pd.DataFrame, sheet_name: str) -> pd.DataFrame:
        """
        Internal method to process a single dataframe/sheet.
        """
        dat_columns = df.columns.tolist()
        unnamed_column = None
        cluster_key_column = None
        cell_type_key_column = None

        # --- 1. Detect key columns ---
        for col in dat_columns:
            col_lower = str(col).lower()
            if unnamed_column is None and "unnamed" in col_lower:
                unnamed_column = col
            if cluster_key_column is None and re.search(r"(leiden|louvain)", col_lower):
                cluster_key_column = col
            if cell_type_key_column is None and re.search(r"(celltype|cell_type)", col_lower):
                cell_type_key_column = col

        # --- 2. Validate essential columns ---
        if not unnamed_column:
            print(f"   Skipping {sheet_name}: No 'Unnamed' split column found.")
            return pd.DataFrame()
        if not cluster_key_column or not cell_type_key_column:
            print(f"   Skipping {sheet_name}: Missing cluster or celltype ID columns.")
            return pd.DataFrame()

        # --- 3. Split DataFrame ---
        split_index = dat_columns.index(unnamed_column)
        deg_columns = dat_columns[:split_index]
        
        # Left side: Gene expression data
        deg_df = df[deg_columns].copy().dropna(how='all')
        
        # Right side: Cluster Metadata
        cluster_df = df[[cluster_key_column, cell_type_key_column]].copy().dropna(how='any')
        
        # Normalize Cluster IDs to string for mapping
        cluster_df[cluster_key_column] = (
            cluster_df[cluster_key_column]
            .astype(float).astype(int).astype(str)
        )

        # --- 4. Expand Cell Type Abbreviations ---
        cluster_df["final_cell_name"] = cluster_df[cell_type_key_column].apply(
            lambda x: self._expand_celltype(str(x))
        )
        
        # Create Map: ID -> Name
        cluster_to_celltype = dict(zip(cluster_df[cluster_key_column], cluster_df["final_cell_name"]))

        # --- 5. Process DEG Data ---
        # Identify the Gene Column (usually the first column or named 'names'/'gene')
        gene_col = None
        for col in deg_df.columns:
            if str(col).lower() in ['names', 'gene', 'symbol']:
                gene_col = col
                break
        if not gene_col and not deg_df.empty:
            # Fallback: Assume first column is gene if it contains strings
            gene_col = deg_df.columns[0]

        # Identify Cluster Column in DEG side
        deg_cluster_col = None
        for col in deg_df.columns:
            if str(col).lower() in ['cluster', 'group']:
                deg_cluster_col = col
                break
        
        if not deg_cluster_col:
            print(f"   Skipping {sheet_name}: No 'cluster' or 'group' column in DEG section.")
            return pd.DataFrame()

        # Filter based on statistical thresholds
        if "auc" in deg_df.columns:
            deg_df = deg_df[(deg_df["auc"] >= 0.75)]
        if "log2FC" in deg_df.columns:
             deg_df = deg_df[(deg_df["log2FC"] >= 0.5)]
        if "padj" in deg_df.columns:
            deg_df = deg_df[(deg_df["padj"] < 0.05)]

        # Map Cell Types
        deg_df[deg_cluster_col] = deg_df[deg_cluster_col].astype(str)
        deg_df['db_cell_name'] = deg_df[deg_cluster_col].map(cluster_to_celltype)

        # Filter out invalid cell types
        deg_df = deg_df[~deg_df['db_cell_name'].isin(['Unknown', 'Doublet', np.nan])]

        # --- 6. Finalize Schema for BaseParser ---
        output_df = pd.DataFrame()
        output_df['gene'] = deg_df[gene_col]
        output_df['db_cell_name'] = deg_df['db_cell_name']
        output_df['db_cell_id'] = None 
        output_df['db_tissue_name'] = self.tissue_name # Use dynamic tissue name from config
        output_df['db_tissue_id'] = None # Let BaseParser find UBERON ID
        output_df['database'] = self.database_name # Use dynamic database name from config
        output_df['source_type'] = "Computational"
        
        # Construct source_info from available stats
        info_parts = []
        if 'padj' in deg_df.columns:
            info_parts.append("padj:" + deg_df['padj'].astype(str))
        if 'log2FC' in deg_df.columns:
            info_parts.append("log2FC:" + deg_df['log2FC'].astype(str))
        if 'auc' in deg_df.columns:
            info_parts.append("auc:" + deg_df['auc'].astype(str))
            
        if info_parts:
            # Vectorized concatenation
            output_df['source_info'] = pd.concat(info_parts, axis=1).apply(lambda x: "; ".join(x), axis=1)
        else:
            output_df['source_info'] = None

        print(f"   ✅ Extracted {len(output_df)} rows from {sheet_name}")
        return output_df

    def _expand_celltype(self, name: str) -> str:
        """
        Expands abbreviations inside complex celltype labels.
        e.g., 'Papil Fib II' -> 'Papillary Fibroblast II'
        """
        if pd.isna(name) or name == 'nan':
            return "Unknown"

        # Split keeping separators: space, /, +
        tokens = re.split(r'([ /+])', name)
        expanded = []

        for tok in tokens:
            # Keep separators as-is
            if tok in [' ', '/', '+', '']:
                expanded.append(tok)
                continue
            
            # Try to match abbreviation at the start
            replaced = False
            for abbr, full in self.celltype_abbreviation_map.items():
                if tok == abbr:
                    expanded.append(full)
                    replaced = True
                    break
                elif tok.startswith(abbr):
                    suffix = tok[len(abbr):]
                    expanded.append(full + suffix)
                    replaced = True
                    break
            
            if not replaced:
                expanded.append(tok)

        return ''.join(expanded)