#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-03
Description :
Marker enrichment analysis pipeline for spatial or scRNA-seq DEGs per cluster.
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
        """
        Initializes the MarkerEnrichmentTest class.

        Args:
            deg_file (str): Path to the CSV file containing differentially expressed genes.
            marker_db_dict (dict): Dictionary with cell types as keys and lists of marker genes as values.
            deg_file_type (str): Type of DEG file, either 'spatial' (wide format) or 'scrna' (long format).
            p_val_thresh (float): Adjusted p-value threshold for filtering significant genes.
            log2fc_thresh (float): Log2 fold change threshold for filtering significant genes.
            mean_counts_thresh (float): Mean counts threshold for filtering significant genes.
                                     For 'spatial' data, if this is 0 (default), it will be
                                     re-calculated as the 75th percentile of positive counts.
            top_genes (int, optional): If provided, selects the top N genes based on Log2 fold change per cluster. Defaults to None.
        """
        if not deg_file or not os.path.exists(deg_file):
            raise FileNotFoundError(f"DEG file not found at: {deg_file}")

        self.raw_deg_df = pd.read_csv(deg_file)
        self.deg_file_type = deg_file_type
        self.marker_dict_clean = {
            str(k).strip(): [str(g).strip() for g in set(v) if pd.notna(g)]
            for k, v in marker_db_dict.items()
        }
        self.p_val_thresh = p_val_thresh
        self.log2fc_thresh = log2fc_thresh
        # Store user-provided value. We will modify this *after* normalization.
        self.mean_counts_thresh = mean_counts_thresh 
        self.top_genes = top_genes
        
        # Normalize the DEG dataframe into a consistent long format
        self.deg_df_long = self._normalize_deg_df()
        
        # If spatial data and user left default mean_counts_thresh (0.0),
        # calculate 75th percentile (3rd quartile) from *positive* counts.
        if self.deg_file_type == 'spatial' and self.mean_counts_thresh == 0.0:
            logging.info("Default mean threshold (0.0) provided for spatial data. Calculating 75th percentile...")
            if 'Mean Counts' in self.deg_df_long.columns:
                # Calculate percentile only from positive (expressed) counts
                positive_mean_counts = self.deg_df_long[self.deg_df_long['Mean Counts'] > 0]['Mean Counts']
                
                if not positive_mean_counts.empty:
                    q3_thresh = positive_mean_counts.quantile(0.75)
                    self.mean_counts_thresh = q3_thresh
                    logging.info(f"Using 75th percentile of positive mean counts as threshold: {self.mean_counts_thresh:.4f}")
                else:
                    logging.warning("No positive mean counts found. Defaulting threshold to 0.0")
                    self.mean_counts_thresh = 0.0
            else:
                logging.warning("'Mean Counts' column not found after normalization. Defaulting threshold to 0.0")
                self.mean_counts_thresh = 0.0
        
        elif self.deg_file_type == 'scrna':
            # scRNA-seq doesn't use mean_counts, so ensure threshold is 0
            if self.mean_counts_thresh > 0.0:
                logging.warning(f"Mean counts threshold ({self.mean_counts_thresh}) is not applicable for 'scrna' data. Setting to 0.0.")
            self.mean_counts_thresh = 0.0
        
        # Background genes are all unique genes present in the normalized table
        self.background_genes = set(self.deg_df_long["Feature Name"].unique())
        self.cluster_markers = {}
        self.results_ = pd.DataFrame()

    def _normalize_deg_df(self):
        """
        Normalizes the input DEG dataframe (wide or long) into a consistent long format.
        """
        logging.info(f"Normalizing DEG file with format: '{self.deg_file_type}'")
        if self.deg_file_type == 'scrna':
            # Assumes 'long' format based on the user's image
            # Columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
            df = self.raw_deg_df.copy()
            # Rename columns to a standard internal representation
            rename_map = {
                'gene': 'Feature Name',
                'cluster': 'Cluster',
                'p_val_adj': 'Adjusted p value',
                'avg_log2FC': 'Log2 fold change',
                # Mean Counts is not present in the provided scRNAseq example, so we'll fill with 0
            }
            df.rename(columns=rename_map, inplace=True)
            
            # Ensure required columns exist
            if 'Feature Name' not in df.columns or 'Cluster' not in df.columns:
                raise ValueError("scRNA-seq DEG file must contain 'gene' and 'cluster' columns.")

            # Add a placeholder for mean counts if it doesn't exist
            if 'Mean Counts' not in df.columns:
                 df['Mean Counts'] = 0
                 logging.info("No 'Mean Counts' column found for scrna data. Creating a placeholder column with 0s.")

            # Ensure cluster names are strings
            df['Cluster'] = "Cluster " + df['Cluster'].astype(str)
            
            return df[['Feature Name', 'Cluster', 'Adjusted p value', 'Log2 fold change', 'Mean Counts']]

        elif self.deg_file_type == 'spatial':
            # Assumes 'wide' format and converts to 'long'
            df = self.raw_deg_df.copy()
            id_vars = [c for c in ['Feature ID', 'Feature Name'] if c in df.columns]
            if not id_vars:
                raise ValueError("Spatial DEG file must contain 'Feature ID' or 'Feature Name'.")
            
            # Melt the dataframe
            df_long = pd.melt(df, id_vars=id_vars, var_name='Metric', value_name='Value')
            
            # Extract Cluster and Metric Type
            df_long[['Cluster', 'Metric Type']] = df_long['Metric'].str.extract(r'(Cluster \d+)\s(.*)')
            df_long.dropna(subset=['Cluster', 'Metric Type'], inplace=True)

            # Pivot to get metrics as columns
            df_pivot = df_long.pivot_table(
                index=id_vars + ['Cluster'],
                columns='Metric Type',
                values='Value',
                aggfunc='first'
            ).reset_index()
            df_pivot.columns.name = None
            
            # Rename for consistency
            df_pivot.rename(columns={
                'Log2 fold change': 'Log2 fold change',
                'Adjusted p value': 'Adjusted p value',
                'Mean Counts': 'Mean Counts'
            }, inplace=True)
            
            # Ensure required columns are numeric
            for col in ['Adjusted p value', 'Log2 fold change', 'Mean Counts']:
                 if col in df_pivot.columns:
                    df_pivot[col] = pd.to_numeric(df_pivot[col], errors='coerce')
                 else: # If a column is missing (like Mean Counts), add it
                    logging.warning(f"Column '{col}' not found in spatial DEG file. Creating a placeholder column with 0s.")
                    df_pivot[col] = 0
            
            # Fill NaNs in numeric columns with 0, especially important for Mean Counts
            df_pivot['Adjusted p value'] = df_pivot['Adjusted p value'].fillna(1.0)
            df_pivot['Log2 fold change'] = df_pivot['Log2 fold change'].fillna(0.0)
            df_pivot['Mean Counts'] = df_pivot['Mean Counts'].fillna(0.0)
            
            return df_pivot

        else:
            raise ValueError(f"Unknown deg_file_type: '{self.deg_file_type}'. Choose 'spatial' or 'scrna'.")


    def fit(self):
        """
        Runs the full enrichment analysis pipeline.
        """
        logging.info("Filtering DEGs from normalized dataframe...")
        # Log the final thresholds being used for filtering
        logging.info(f"Filtering with: p_val <= {self.p_val_thresh}, log2fc >= {self.log2fc_thresh}, mean_counts >= {self.mean_counts_thresh}")
        self.cluster_markers = self._filter_and_extract_degs()
        logging.info("Running hypergeometric enrichment test...")
        self.results_ = self._run_hypergeometric_test()
        if self.results_.empty:
            logging.warning("No significant enrichment results found.")
        return self

    def _filter_and_extract_degs(self):
        """
        Filters and extracts significant DEGs for each cluster from the normalized long-format DataFrame.
        """
        cluster_markers_dict = {}
        
        # The dataframe is already in a long, tidy format, so we can group by cluster
        for cluster_id, cluster_df in self.deg_df_long.groupby('Cluster'):
            
            # Build the filter mask
            mask = (cluster_df['Adjusted p value'] < self.p_val_thresh) & \
                   (cluster_df['Log2 fold change'] > self.log2fc_thresh)

            # This condition now uses the dynamically calculated threshold for spatial
            # or 0 for scrna.
            if self.mean_counts_thresh > 0 and 'Mean Counts' in cluster_df.columns:
                mask &= (cluster_df['Mean Counts'] > self.mean_counts_thresh)

            sig_genes_df = cluster_df[mask].copy()

            if self.top_genes and isinstance(self.top_genes, int):
                # Sort by Log2 fold change (descending) to get the "top" genes
                sig_genes_df = sig_genes_df.sort_values(
                    by='Log2 fold change', ascending=False
                ).head(self.top_genes)

            genes = sig_genes_df["Feature Name"].dropna().astype(str).unique().tolist()
            cluster_markers_dict[cluster_id] = genes
            
        return cluster_markers_dict

    def _run_hypergeometric_test(self):
        # This method does not need changes as it operates on the output of _filter_and_extract_degs
        N = len(self.background_genes)
        all_results_list = []
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
                enrichment = (k / n) / (K / N) if K > 0 and n > 0 else 0
                results.append(
                    {
                        "Cluster": cluster_id,
                        "Cell_type": cell_type,
                        "p_value": p_val,
                        "Enrichment_ratio": enrichment,
                        "Overlapping_genes_count": k,
                        "Overlapping_genes": ", ".join(sorted(overlap)),
                    }
                )

            if results:
                df = pd.DataFrame(results)
                df["adj_p_value"] = multi.multipletests(df["p_value"], method="fdr_bh")[1]
                all_results_list.append(df)

        if not all_results_list:
            return pd.DataFrame()
        return pd.concat(all_results_list, ignore_index=True).sort_values(
            ["Cluster", "adj_p_value"]
        )

    def plot_results(
        self, p_val_cutoff=0.05, value_col="Enrichment_ratio", top_n_per_cluster=10
    ):
        # This method does not need changes
        if self.results_.empty:
            logging.warning("No results to plot.")
            return

        plot_df = self.results_[self.results_["adj_p_value"] < p_val_cutoff].copy()
        if plot_df.empty:
            logging.warning(
                f"No significant results to plot with adjusted p-value cutoff {p_val_cutoff}."
            )
            return

        plot_df = (
            plot_df.groupby("Cluster")
            .apply(lambda x: x.nsmallest(top_n_per_cluster, "adj_p_value"))
            .reset_index(drop=True)
        )

        if plot_df.empty:
            logging.warning(f"No results remain after selecting top {top_n_per_cluster} per cluster.")
            return

        plot_df["-log10(adj_p_value)"] = -np.log10(plot_df["adj_p_value"])

        plt.figure(
            figsize=(
                max(8, len(plot_df["Cluster"].unique()) * 0.8),
                max(6, len(plot_df["Cell_type"].unique()) * 0.5),
            )
        )

        scatter = sns.scatterplot(
            data=plot_df,
            x="Cluster",
            y="Cell_type",
            size="Overlapping_genes_count",
            hue=value_col,
            palette="viridis",
            sizes=(50, 500),
            edgecolor="black",
            linewidth=0.5,
        )

        plt.title("Marker Enrichment Analysis Results")
        plt.xlabel("Cluster")
        plt.ylabel("Cell_type")
        plt.xticks(rotation=45, ha="right")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.show()