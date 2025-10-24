#!/usr/bin/python3
"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 21/10/2025
Description     : Parses the Human SCC from Cell 2020 (mmc2.xlsx).
                  Correctly handles potential empty columns and varying column orders within triplets
                  by identifying the gene column first.
                  Uses "Keratinocyte" for lookup while retaining original db_cell_name.
"""
import pandas as pd
import os
import sys
import numpy as np

try:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)

    base_parser_path = os.path.join(project_root, 'src', 'parser')
    if base_parser_path not in sys.path:
        sys.path.append(base_parser_path)
    from base_parser import BaseParser

except ImportError as e:
    print(f"❌ A critical import error occurred in human_scc_cell_2020_parser.py: {e}")
    print(f"   Attempted to import BaseParser from paths including: {base_parser_path}")
    print(f"   Current sys.path: {sys.path}")
    sys.exit(1)

class HumanSccCell2020Parser:
    """
    Parses the wide-format, multi-header file from the Human SCC (Cell 2020).
    It transforms the data into a long format, identifies triplets based on the
    gene column, and passes it to the BaseParser using a temporary lookup name.
    """

    def __init__(self, df: pd.DataFrame, **kwargs):
        """
        Initializes the parser. Assumes 'df' is read from the raw Excel/CSV
        with multi-level headers (like mmc2.xlsx - Sheet 1.csv).
        """
        if df.empty:
            print("❌ Input DataFrame is empty for Human SCC parser.")
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()
            return

        all_markers = []

        # --- Header Processing ---
        if len(df) < 3:
             print("❌ DataFrame has fewer than 3 rows, cannot parse headers for Human SCC.")
             self.processed_df = pd.DataFrame()
             self.recovery_df = pd.DataFrame()
             return

        main_headers = df.iloc[1].ffill()
        sub_headers = df.iloc[2]
        data_df = df.iloc[3:].copy()
        num_cols = len(df.columns)
        # --- End Header Processing ---

        processed_indices = set() # Keep track of columns already processed as part of a triplet
        i = 0
        while i < num_cols:
            # Skip if this index was already part of a previous triplet or has no main header
            if i in processed_indices or pd.isna(main_headers.iloc[i]):
                i += 1
                continue

            # --- Check if current column `i` is a potential Gene Column ---
            sub_header_val = sub_headers.iloc[i]
            # A potential gene column's sub-header is not 'avg_logFC' or 'p_val_adj' (case-insensitive)
            # Also handle NaN sub-headers explicitly - they cannot be logFC or pVal, so could be gene column
            is_potential_gene_col = False
            if pd.isna(sub_header_val):
                 # Consider NaN sub-header as potential gene column *only if* the main header is valid
                 if pd.notna(main_headers.iloc[i]):
                      is_potential_gene_col = True
                      print(f"⚠️ Warning: Column {i} has NaN sub-header but valid main header '{main_headers.iloc[i]}'. Assuming it's a gene column.")
                 else:
                      # If main header is also NaN, skip entirely
                      i += 1
                      continue
            else:
                 sub_header_str = str(sub_header_val).lower()
                 if 'avg_logfc' not in sub_header_str and 'p_val_adj' not in sub_header_str:
                      is_potential_gene_col = True

            if not is_potential_gene_col:
                # If column `i` is definitively logFC or pVal, skip it, it cannot start a block
                processed_indices.add(i) # Mark as processed (as a non-gene column)
                i += 1
                continue
            # --- End Check ---


            # --- Found Potential Gene Column at index `i` ---
            gene_col_idx = i
            gene_sub_header = str(sub_headers.iloc[gene_col_idx]).strip() if pd.notna(sub_headers.iloc[gene_col_idx]) else 'nan' # Use 'nan' string if NaN
            main_header = main_headers.iloc[gene_col_idx]
            base_db_cell_name = str(main_header).strip().replace(' ', '_')

            # Look for partners in the next two columns
            logfc_col_idx, pval_col_idx = -1, -1
            partner_indices = []

            # Check index i+1
            if i + 1 < num_cols and main_headers.iloc[i+1] == main_header: # Ensure partner belongs to the same main header block
                 partner_idx = i + 1
                 partner_sub_header_val = sub_headers.iloc[partner_idx]
                 if pd.notna(partner_sub_header_val):
                      partner_sub_header_str = str(partner_sub_header_val).lower()
                      if 'avg_logfc' in partner_sub_header_str:
                           logfc_col_idx = partner_idx
                           partner_indices.append(partner_idx)
                      elif 'p_val_adj' in partner_sub_header_str:
                           pval_col_idx = partner_idx
                           partner_indices.append(partner_idx)

            # Check index i+2
            if i + 2 < num_cols and main_headers.iloc[i+2] == main_header:
                 partner_idx = i + 2
                 partner_sub_header_val = sub_headers.iloc[partner_idx]
                 if pd.notna(partner_sub_header_val):
                      partner_sub_header_str = str(partner_sub_header_val).lower()
                      # Only assign if not already found at i+1
                      if logfc_col_idx == -1 and 'avg_logfc' in partner_sub_header_str:
                           logfc_col_idx = partner_idx
                           partner_indices.append(partner_idx)
                      elif pval_col_idx == -1 and 'p_val_adj' in partner_sub_header_str:
                           pval_col_idx = partner_idx
                           partner_indices.append(partner_idx)


            # --- Construct Final db_cell_name ---
            final_db_cell_name = base_db_cell_name
            # Check gene_sub_header (which comes from the identified gene column)
            if gene_sub_header and gene_sub_header.lower() != 'nan' and gene_sub_header.lower() != 'gene':
                clean_state = gene_sub_header.replace(' ', '_')
                if clean_state != base_db_cell_name: # Avoid duplication like Normal_Keratinocytes_Normal_Keratinocytes
                    final_db_cell_name = f"{base_db_cell_name}_{clean_state}"
            # --- End Construct Final db_cell_name ---

            # --- Extract Data ---
            genes = data_df.iloc[:, gene_col_idx].dropna()
            if genes.empty:
                 print(f"ℹ️ No gene data found for '{final_db_cell_name}' in column {gene_col_idx}. Skipping block.")
            else:
                chunk_df = pd.DataFrame({'gene': genes})
                chunk_df['db_cell_name'] = final_db_cell_name

                if logfc_col_idx != -1:
                    chunk_df['avg_logFC'] = pd.to_numeric(data_df.iloc[:, logfc_col_idx].reindex(genes.index), errors='coerce')
                else:
                    chunk_df['avg_logFC'] = np.nan
                    print(f"ℹ️ Note: Missing 'avg_logFC' column for '{final_db_cell_name}' (gene col {gene_col_idx}).")

                if pval_col_idx != -1:
                    chunk_df['p_val_adj'] = pd.to_numeric(data_df.iloc[:, pval_col_idx].reindex(genes.index), errors='coerce')
                else:
                    chunk_df['p_val_adj'] = np.nan
                    print(f"ℹ️ Note: Missing 'p_val_adj' column for '{final_db_cell_name}' (gene col {gene_col_idx}).")

                all_markers.append(chunk_df)
            # --- End Extract Data ---

            # Mark processed indices and advance main index `i`
            processed_indices.add(gene_col_idx)
            for partner_idx in partner_indices:
                processed_indices.add(partner_idx)
            # Advance i past the current block (guaranteed to move by at least 1)
            i = max([gene_col_idx] + partner_indices) + 1


        # --- Post-Processing ---
        if not all_markers:
            print("❌ No data parsed from Human SCC file after processing all potential blocks.")
            self.processed_df = pd.DataFrame()
            self.recovery_df = pd.DataFrame()
            return

        combined_df = pd.concat(all_markers, ignore_index=True)

        # Step 2: Add metadata
        combined_df['db_tissue_name'] = "Skin"
        combined_df['db_tissue_id'] = None
        combined_df['db_cell_id'] = None
        combined_df['database'] = "HumanSCC_Cell2020"
        combined_df['source_type'] = "Computational"
        combined_df['source_info'] = (
            "p_val_adj:" + combined_df['p_val_adj'].astype(str).fillna('NA') +
            "; avg_logFC:" + combined_df['avg_logFC'].astype(str).fillna('NA')
        )
        combined_df = combined_df.drop(['p_val_adj', 'avg_logFC'], axis=1)

        # Step 3: Normalize using BaseParser
        original_db_cell_names = combined_df['db_cell_name'].copy()
        combined_df['db_cell_name'] = "Keratinocyte" # Force lookup term

        normalizer = BaseParser()
        processed_df_temp, self.recovery_df = normalizer.normalize_dataframe(combined_df)

        if len(processed_df_temp) == len(original_db_cell_names):
            processed_df_temp['db_cell_name'] = original_db_cell_names.values
        else:
            print(f"⚠️ Warning: Length mismatch after normalization in HumanSCC parser.")
            print(f"   Original length: {len(original_db_cell_names)}, Processed length: {len(processed_df_temp)}")
            print("   Cannot reliably restore original db_cell_name. Check recovery logs.")
            if 'db_cell_name' in processed_df_temp.columns and len(processed_df_temp) > 0:
                 try:
                     processed_df_temp['db_cell_name'] = original_db_cell_names.reindex(processed_df_temp.index).values
                 except Exception as e:
                      print(f"   Error during partial restore: {e}. 'db_cell_name' might be incorrect.")

        self.processed_df = processed_df_temp
        # --- End Post-Processing ---