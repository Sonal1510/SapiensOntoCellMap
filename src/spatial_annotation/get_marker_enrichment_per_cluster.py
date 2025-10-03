#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-08-25
Description :
Marker enrichment analysis pipeline for single-cell DEGs.
"""

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import re
import os
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

class MarkerEnrichmentTest:
    """Pipeline for enrichment analysis of DEGs against a marker database."""

    def __init__(self, deg_file, marker_db_dict,
                 p_val_thresh=0.05, log2fc_thresh=1.0, mean_counts_thresh=0):
        """
        Initialize the pipeline with DEG and marker database files.

        Args:
            deg_file (str): Path to CSV file containing differential expression results
            marker_db_dict (dict): Dict mapping cell_id -> list of marker genes
            p_val_thresh (float): Adjusted p-value cutoff for DEGs
            log2fc_thresh (float): Log2 fold-change cutoff for DEGs
            mean_counts_thresh (float): Mean expression cutoff for DEGs
        """
        if not deg_file or not os.path.exists(deg_file):
            raise FileNotFoundError(f"DEG file not found at: {deg_file}")

        self.deg_file = deg_file
        self.deg_df = pd.read_csv(deg_file)

        # Clean marker dictionary
        self.marker_dict_clean = {
            str(cell_id).strip(): [str(g).strip() for g in set(gene_list) if pd.notna(g)]
            for cell_id, gene_list in marker_db_dict.items()
        }

        # Parameters
        self.p_val_thresh = p_val_thresh
        self.log2fc_thresh = log2fc_thresh
        self.mean_counts_thresh = mean_counts_thresh

        # Internal attributes
        self.background_genes = set(self.deg_df['Feature Name'].unique())
        self.cluster_markers = None
        self.results_ = None

    def fit(self):
        """Run the full pipeline: filter DEGs and compute enrichment tests."""
        logging.info("Filtering DEGs...")
        self.cluster_markers = self._filter_and_extract_degs()

        logging.info("Running hypergeometric enrichment test...")
        self.results_ = self._run_hypergeometric_test()

        if self.results_.empty:
            logging.warning("No significant enrichment results found.")
        else:
            logging.info(f"Enrichment results shape: {self.results_.shape}")

        return self

    def _filter_and_extract_degs(self):
        """Filter DEGs per cluster based on thresholds."""
        cluster_cols = [col for col in self.deg_df.columns if "Cluster" in col and "Adjusted p value" in col]
        cluster_ids = sorted([int(re.search(r"(\d+)", col).group(1)) for col in cluster_cols])

        cluster_markers = {}
        for i in cluster_ids:
            p_val_col = f"Cluster {i} Adjusted p value"
            log2fc_col = f"Cluster {i} Log2 fold change"
            mean_counts_col = f"Cluster {i} Mean Counts"

            sig_genes = self.deg_df[
                (self.deg_df[p_val_col] < self.p_val_thresh) &
                (self.deg_df[log2fc_col] > self.log2fc_thresh) &
                (self.deg_df[mean_counts_col] > self.mean_counts_thresh)
            ]

            ranked = sig_genes.sort_values(log2fc_col, ascending=False)
            genes = ranked['Feature Name'].dropna().astype(str).tolist()
            cluster_markers[f"Cluster {i}"] = list(set(genes))

        logging.info(f"Found significant markers for {len(cluster_markers)} clusters.")
        return cluster_markers

    def _run_hypergeometric_test(self):
        """Perform hypergeometric test for enrichment."""
        N = len(self.background_genes)
        all_results = []

        for cluster_id, cluster_genes in self.cluster_markers.items():
            n = len(cluster_genes)
            if n == 0:
                continue

            results = []
            for cell_type, cell_genes in self.marker_dict_clean.items():
                cell_genes_set = set(cell_genes).intersection(self.background_genes)
                K = len(cell_genes_set)
                if K == 0:
                    continue

                overlap = set(cluster_genes).intersection(cell_genes_set)
                k = len(overlap)
                if k == 0:
                    continue

                p_val = hypergeom.sf(k - 1, N, K, n)
                enrichment = (k / n) / (K / N) if K > 0 and n > 0 else np.nan

                results.append({
                    "Cluster": cluster_id,
                    "Cell Type": cell_type,
                    "p_value": p_val,
                    "Enrichment Ratio": enrichment,
                    "Overlap Size (k)": k,
                    "Cluster Markers (n)": n,
                    "Cell Type Markers (K)": K,
                    "Overlapping Genes": ", ".join(sorted(overlap))
                })

            if results:
                df = pd.DataFrame(results)
                _, pvals_corrected, _, _ = multi.multipletests(df["p_value"], method="fdr_bh")
                df["adj_p_value"] = pvals_corrected
                all_results.append(df)

        if not all_results:
            return pd.DataFrame()

        final_df = pd.concat(all_results, ignore_index=True)
        final_df = final_df.sort_values(["Cluster", "adj_p_value", "Enrichment Ratio"],
                                        ascending=[True, True, False])
        return final_df

