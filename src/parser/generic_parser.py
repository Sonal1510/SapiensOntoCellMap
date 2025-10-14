#!/usr/bin/python3
"""
Author          : Gemini
Date            : 14/10/2025
Description     : A configurable parser for simple, ad-hoc data files.
"""
import pandas as pd
from typing import Optional
from src.parser.base_parser import BaseParser

class GenericFileParser:
    """
    A generic parser for files that follow a simple structure.
    This class is configured with column names and metadata upon initialization.
    """
    def __init__(self,
                 df: pd.DataFrame,
                 database_name: str,
                 gene_col: str,
                 tissue_name: Optional[str] = None,
                 tissue_name_col: Optional[str] = None,
                 base_cell_name: Optional[str] = None,
                 cell_name_col: Optional[str] = None,
                 cell_subtype_col: Optional[str] = None,
                 source_type: str = "Unspecified",
                 source_info_cols: Optional[list] = None,
                 info_separator: str = "; "):
        """
        Initializes the generic parser and processes the data.

        Args:
            df (pd.DataFrame): The input dataframe to parse.
            database_name (str): The name of the database/source.
            gene_col (str): The name of the column containing gene symbols.
            tissue_name (Optional[str]): A fixed tissue name to apply to all rows.
            tissue_name_col (Optional[str]): The column name containing tissue names.
            base_cell_name (Optional[str]): A fixed base cell name for constructing cell names.
            cell_name_col (Optional[str]): The column name containing full cell names.
            cell_subtype_col (Optional[str]): Column with subtypes to append to 'base_cell_name'.
            source_type (str): The source type (e.g., "Computational").
            source_info_cols (Optional[list]): Columns to concatenate for the 'source_info' field.
            info_separator (str): Separator for joining 'source_info_cols'.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        # --- Validation for clear configuration ---
        if not tissue_name and not tissue_name_col:
            raise ValueError("Must provide either 'tissue_name' (a fixed string) or 'tissue_name_col' (a column name).")
        if tissue_name and tissue_name_col:
            raise ValueError("Cannot provide both 'tissue_name' and 'tissue_name_col'. Please choose one.")
        if not base_cell_name and not cell_name_col:
            raise ValueError("Must provide either 'base_cell_name' (for constructing names) or 'cell_name_col' (a column with full cell names).")
        if cell_name_col and (base_cell_name or cell_subtype_col):
            raise ValueError("If 'cell_name_col' is provided, 'base_cell_name' and 'cell_subtype_col' must be omitted.")

        # --- Step 1: Prepare the DataFrame ---
        df_proc = df.copy()

        # Map required columns to the 'db_' schema
        df_proc.rename(columns={gene_col: 'gene'}, inplace=True)
        df_proc['db_tissue_id'] = None
        df_proc['db_cell_id'] = None

        # Construct tissue name (either from fixed string or a column)
        if tissue_name:
            df_proc['db_tissue_name'] = tissue_name
        else:
            df_proc['db_tissue_name'] = df_proc[tissue_name_col]

        # Construct cell name (either from a column or by combining base + subtype)
        if cell_name_col:
            df_proc['db_cell_name'] = df_proc[cell_name_col]
        else:
            if cell_subtype_col and cell_subtype_col in df_proc.columns:
                df_proc['db_cell_name'] = base_cell_name + "_" + df_proc[cell_subtype_col].astype(str)
            else:
                df_proc['db_cell_name'] = base_cell_name
        
        # --- Step 2: Handle Metadata ---
        df_proc['database'] = database_name
        df_proc['source_type'] = source_type

        if source_info_cols:
            info_df = df_proc[source_info_cols].astype(str)
            df_proc['source_info'] = info_df.apply(
                lambda row: info_separator.join(f"{col}:{val}" for col, val in row.items()),
                axis=1
            )
        else:
            df_proc['source_info'] = None

        # --- Step 3: Use the BaseParser for normalization ---
        normalizer = BaseParser()
        self.processed_df, self.recovery_df = normalizer.normalize_dataframe(df_proc)

# --- USAGE EXAMPLE ---
#
# To use this parser, you would typically do the following in your orchestrator script
# (e.g., database_creator.py), driven by your config file.
#
# import pandas as pd
# from src.parser.generic_parser import GenericFileParser
#
# # --- SCENARIO 1: Hard-coded tissue, constructed cell name (like Epi_Cluster) ---
# data1 = {'GeneSymbol': ['CD4', 'CD8A'], 'ClusterID': ['T_helper', 'T_cyto'], 'Score': [0.9, 0.95]}
# df1 = pd.DataFrame(data1)
#
# parser1 = GenericFileParser(
#     df=df1,
#     database_name="MyFirstDataset",
#     gene_col="GeneSymbol",
#     tissue_name="Blood",  # Fixed tissue name
#     base_cell_name="T-Cell",  # Base for cell name
#     cell_subtype_col="ClusterID",  # Subtype to append
#     source_info_cols=["Score"]
# )
# # This would produce db_cell_names like "T-Cell_T_helper"
# # processed_df1 = parser1.processed_df
#
# # --- SCENARIO 2: Tissue name from a column, cell name is fixed ---
# data2 = {'markers': ['KRT5', 'KRT14'], 'organ': ['Skin', 'Esophagus'], 'pval': [0.01, 0.005]}
# df2 = pd.DataFrame(data2)
#
# parser2 = GenericFileParser(
#     df=df2,
#     database_name="MySecondDataset",
#     gene_col="markers",
#     tissue_name_col="organ",  # Get tissue from the 'organ' column
#     base_cell_name="Basal Epithelial Cell",  # Cell name is the same for all rows
#     source_info_cols=["pval"]
# )
# # This would produce db_tissue_names "Skin", "Esophagus" and a db_cell_name "Basal Epithelial Cell"
# # processed_df2 = parser2.processed_df
#
# # --- SCENARIO 3: Both tissue and cell names are read directly from columns ---
# data3 = {'gene_id': ['EPCAM', 'COL1A1'], 'tissue_type': ['Lung', 'Skin'], 'cell_type_original': ['Epithelial cell', 'Fibroblast']}
# df3 = pd.DataFrame(data3)
#
# parser3 = GenericFileParser(
#     df=df3,
#     database_name="MyThirdDataset",
#     gene_col="gene_id",
#     tissue_name_col="tissue_type",  # Get tissue from a column
#     cell_name_col="cell_type_original"  # Get the full cell name from a column
# )
# # This is the most flexible option for simple tables.
# # processed_df3 = parser3.processed_df