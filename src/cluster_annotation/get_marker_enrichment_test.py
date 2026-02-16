#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Marker enrichment analysis pipeline.
- Automatically selects 'Conventional' statistics for simple gene lists.
- Automatically selects 'Weighted' statistics for gene-weight dictionaries.
"""

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import os
import logging

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

class MarkerEnrichmentTest:
    def __init__(
        self,
        deg_file,
        marker_db_dict,
        deg_file_type='spatial', 
        p_val_thresh=0.05,
        log2fc_thresh=1.0,
        mean_counts_thresh=0,
        top_genes=None,
    ):
        if not deg_file or not os.path.exists(deg_file):
            raise FileNotFoundError(f"DEG file not found: {deg_file}")

        self.raw_deg_df = pd.read_csv(deg_file)
        self.deg_file_type = deg_file_type
        
        # --- DETECT MODE (Weighted vs Conventional) ---
        # Normalize DB to {cell: {gene: weight}} for internal processing, 
        # but keep a flag to know if we should output weighted columns.
        self.marker_dict_internal = {}
        self.is_weighted_mode = False
        
        all_db_genes = set()
        
        # Check first item to determine mode (assuming consistency)
        first_val = next(iter(marker_db_dict.values())) if marker_db_dict else []
        if isinstance(first_val, dict):
            self.is_weighted_mode = True
            logging.info("Input is dictionary: Activating WEIGHTED Enrichment Statistics.")
        else:
            self.is_weighted_mode = False
            logging.info("Input is list: Activating CONVENTIONAL Enrichment Statistics.")

        for k, v in marker_db_dict.items():
            k_str = str(k).strip()
            if self.is_weighted_mode:
                # {gene: weight}
                clean_v = {str(g).strip(): float(w) for g, w in v.items() if pd.notna(g)}
            else:
                # List of genes -> Convert to {gene: 1.0} for unified math
                clean_v = {str(g).strip(): 1.0 for g in v if pd.notna(g)}
                
            self.marker_dict_internal[k_str] = clean_v
            all_db_genes.update(clean_v.keys())

        self.p_val_thresh = p_val_thresh
        self.log2fc_thresh = log2fc_thresh
        self.mean_counts_thresh = mean_counts_thresh 
        self.top_genes = top_genes
        
        self.deg_df_long = self._normalize_deg_df()
        
        if self.deg_file_type == 'spatial' and self.mean_counts_thresh == 0.0:
            if 'Mean Counts' in self.deg_df_long.columns:
                pos = self.deg_df_long[self.deg_df_long['Mean Counts'] > 0]['Mean Counts']
                if not pos.empty:
                    self.mean_counts_thresh = pos.quantile(0.75)

        self.background_genes = set(self.deg_df_long["Feature Name"].unique())
        self.cluster_markers = {}
        self.results_ = pd.DataFrame()

    def _normalize_deg_df(self):
        # ... (Same normalization logic) ...
        if self.deg_file_type == 'scrna':
            df = self.raw_deg_df.copy()
            rename = {'gene': 'Feature Name', 'cluster': 'Cluster', 'p_val_adj': 'Adjusted p value', 'avg_log2FC': 'Log2 fold change'}
            df.rename(columns=rename, inplace=True)
            if 'Mean Counts' not in df.columns: df['Mean Counts'] = 0
            df['Cluster'] = "Cluster " + df['Cluster'].astype(str)
            return df[['Feature Name', 'Cluster', 'Adjusted p value', 'Log2 fold change', 'Mean Counts']]
        elif self.deg_file_type == 'spatial':
            df = self.raw_deg_df.copy()
            id_vars = [c for c in ['Feature ID', 'Feature Name'] if c in df.columns]
            df_long = pd.melt(df, id_vars=id_vars, var_name='Metric', value_name='Value')
            df_long[['Cluster', 'Metric Type']] = df_long['Metric'].str.extract(r'(Cluster \d+)\s(.*)')
            df_long.dropna(subset=['Cluster', 'Metric Type'], inplace=True)
            df = df_long.pivot_table(index=id_vars+['Cluster'], columns='Metric Type', values='Value', aggfunc='first').reset_index()
            df.rename(columns={'Adjusted p value': 'Adjusted p value', 'Log2 fold change': 'Log2 fold change', 'Mean Counts': 'Mean Counts'}, inplace=True)
            df.fillna({'Adjusted p value': 1.0, 'Log2 fold change': 0.0, 'Mean Counts': 0.0}, inplace=True)
            return df
        else:
            raise ValueError(f"Unknown type: {self.deg_file_type}")

    def fit(self):
        self.cluster_markers = self._filter_and_extract_degs()
        self.results_ = self._run_test()
        return self

    def _filter_and_extract_degs(self):
        # ... (Same logic) ...
        d = {}
        for cid, cdf in self.deg_df_long.groupby('Cluster'):
            mask = (cdf['Adjusted p value'] < self.p_val_thresh) & (cdf['Log2 fold change'] > self.log2fc_thresh)
            if self.mean_counts_thresh > 0 and 'Mean Counts' in cdf.columns:
                mask &= (cdf['Mean Counts'] > self.mean_counts_thresh)
            res = cdf[mask].copy()
            if self.top_genes: res = res.sort_values('Log2 fold change', ascending=False).head(self.top_genes)
            d[cid] = res["Feature Name"].dropna().astype(str).unique().tolist()
        return d

    def _run_test(self):
        N = len(self.background_genes)
        all_results_list = []
        
        for cluster_id, cluster_genes in self.cluster_markers.items():
            n = len(cluster_genes)
            if n == 0: continue

            results = []
            for cell_type, weighted_genes in self.marker_dict_internal.items():
                cell_genes_set = set(weighted_genes.keys()).intersection(self.background_genes)
                K = len(cell_genes_set)
                if K == 0: continue

                overlap_genes = set(cluster_genes).intersection(cell_genes_set)
                k = len(overlap_genes)
                if k == 0: continue

                # 1. P-value (Always Hypergeometric on counts)
                p_val = hypergeom.sf(k - 1, N, K, n)
                
                # 2. Base Metrics
                enrichment_ratio = (k / n) / (K / N) if K > 0 else 0
                
                row = {
                    "Cluster": cluster_id,
                    "Cell_type": cell_type,
                    "p_value": p_val,
                    "Enrichment_ratio": enrichment_ratio,
                    "Overlapping_genes_count": k,
                    "Overlapping_genes": ", ".join(sorted(overlap_genes)),
                }

                # 3. Weighted Metrics (Only calculate/add if in Weighted Mode)
                if self.is_weighted_mode:
                    overlap_w_sum = sum(weighted_genes[g] for g in overlap_genes)
                    ref_w_sum = sum(weighted_genes[g] for g in cell_genes_set)
                    
                    # Weighted Recall
                    row["Weighted_Recall"] = overlap_w_sum / ref_w_sum if ref_w_sum > 0 else 0
                    
                    # Weighted Enrichment
                    # (Sum of Weights in Overlap / n) / (Sum of Weights in Background / N)
                    # NOTE: Background weight sum is needed for true ratio. 
                    # Approximation: (OverlapW / n) * (N / RefW)
                    row["Weighted_Enrichment"] = (overlap_w_sum / n) * (N / ref_w_sum) if ref_w_sum > 0 else 0

                results.append(row)

            if results:
                df = pd.DataFrame(results)
                df["adj_p_value"] = multi.multipletests(df["p_value"], method="fdr_bh")[1]
                all_results_list.append(df)

        if not all_results_list: return pd.DataFrame()
        return pd.concat(all_results_list, ignore_index=True).sort_values(["Cluster", "adj_p_value"])