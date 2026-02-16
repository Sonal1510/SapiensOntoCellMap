#!/usr/bin/python3
"""
Author          : Gemini
Date            : 14/10/2025
Description     : A configurable parser for simple, ad-hoc data files.
"""
import pandas as pd
from typing import Optional
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
            base_cell_name (Optional[str]): A fixed base cell name for *querying*.
            cell_name_col (Optional[str]): The column name containing full cell names (for querying).
            cell_subtype_col (Optional[str]): Column with subtypes to append to 'base_cell_name' *after* querying.
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

        # --- (NEW LOGIC) Handle Cell Name Querying vs. Final Name ---
        
        # This variable will hold the *final* specific names, if they are different from the query names.
        final_cell_names: Optional[pd.Series] = None 

        if cell_name_col:
            # Mode 1: A column provides the full cell name.
            # We use this for *both* querying and as the final name.
            df_proc['db_cell_name'] = df_proc[cell_name_col]
        else:
            # Mode 2: We use 'base_cell_name'.
            # The 'db_cell_name' column is *temporarily* set to the base name for querying.
            df_proc['db_cell_name'] = base_cell_name
            
            if cell_subtype_col and cell_subtype_col in df_proc.columns:
                # If a subtype exists, we *store* the final composite name
                # to be used *after* normalization.
                final_cell_names = base_cell_name + "_" + df_proc[cell_subtype_col].astype(str)
            # If no subtype, the query name (base_cell_name) is also the final name,
            # so 'final_cell_names' remains None.
        
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
        # The normalizer will query using the 'db_cell_name' column,
        # which now contains the *general* name (e.g., "T-Cell").
        normalizer = BaseParser()
        processed_df, self.recovery_df = normalizer.normalize_dataframe(df_proc)

        # --- Step 4: (NEW LOGIC) Overwrite with Final Cell Name ---
        # If we stored final composite names, we apply them now to the
        # *processed* dataframe, replacing the general query name.
        if final_cell_names is not None:
            # Use .loc to align the names based on the index,
            # which is crucial if normalization dropped any rows.
            processed_df['db_cell_name'] = final_cell_names.loc[processed_df.index]

        # Store the final processed dataframe
        self.processed_df = processed_df