#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-03
Description :
Marker enrichment analysis pipeline for spatial DEGs per cluster.
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
            p_val_thresh (float): Adjusted p-value threshold for filtering significant genes.
            log2fc_thresh (float): Log2 fold change threshold for filtering significant genes.
            mean_counts_thresh (float): Mean counts threshold for filtering significant genes.
            top_genes (int, optional): If provided, selects the top N genes based on Log2 fold change per cluster. Defaults to None.
        """
        if not deg_file or not os.path.exists(deg_file):
            raise FileNotFoundError(f"DEG file not found at: {deg_file}")

        self.deg_df = pd.read_csv(deg_file)
        self.marker_dict_clean = {
            str(k).strip(): [str(g).strip() for g in set(v) if pd.notna(g)]
            for k, v in marker_db_dict.items()
        }
        self.p_val_thresh = p_val_thresh
        self.log2fc_thresh = log2fc_thresh
        self.mean_counts_thresh = mean_counts_thresh
        self.top_genes = top_genes
        self.background_genes = set(self.deg_df["Feature Name"].unique())
        self.cluster_markers = {}
        self.results_ = pd.DataFrame()

    def fit(self):
        """
        Runs the full enrichment analysis pipeline.
        """
        logging.info("Filtering DEGs...")
        self.cluster_markers = self._filter_and_extract_degs()
        logging.info("Running hypergeometric enrichment test...")
        self.results_ = self._run_hypergeometric_test()
        if self.results_.empty:
            logging.warning("No significant enrichment results found.")
        return self

    def _filter_and_extract_degs(self):
        """
        Filters and extracts significant DEGs for each cluster based on defined thresholds.
        """
        cluster_cols = [
            c for c in self.deg_df.columns if "Cluster" in c and "Adjusted p value" in c
        ]
        if not cluster_cols:
            logging.error(
                "No valid 'Cluster X Adjusted p value' columns found in DEG file."
            )
            return {}

        cluster_ids = sorted(
            [int(re.search(r"(\d+)", c).group(1)) for c in cluster_cols]
        )
        cluster_markers_dict = {}
        for i in cluster_ids:
            p_col = f"Cluster {i} Adjusted p value"
            fc_col = f"Cluster {i} Log2 fold change"
            mc_col = f"Cluster {i} Mean Counts"

            required_cols = [p_col, fc_col]
            # Only require mean counts column if a threshold is set
            if self.mean_counts_thresh > 0:
                required_cols.append(mc_col)

            if not all(col in self.deg_df.columns for col in required_cols):
                logging.warning(
                    f"Skipping Cluster {i} due to missing required columns ({', '.join(required_cols)})."
                )
                continue

            # Build the filter mask
            mask = (self.deg_df[p_col] < self.p_val_thresh) & (
                self.deg_df[fc_col] > self.log2fc_thresh
            )

            if self.mean_counts_thresh > 0:
                mask &= self.deg_df[mc_col] > self.mean_counts_thresh

            sig_genes_df = self.deg_df[mask].copy()

            if self.top_genes and isinstance(self.top_genes, int):
                # Sort by Log2 fold change (descending) to get the "top" genes
                sig_genes_df = sig_genes_df.sort_values(
                    by=fc_col, ascending=False
                ).head(self.top_genes)

            genes = sig_genes_df["Feature Name"].dropna().astype(str).unique().tolist()
            cluster_markers_dict[f"Cluster {i}"] = genes
        return cluster_markers_dict

    def _run_hypergeometric_test(self):
        """
        Performs the hypergeometric test for each cluster against each cell type marker list.
        """
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

                # p-value for enrichment: P(X >= k)
                p_val = hypergeom.sf(k - 1, N, K, n)
                enrichment = (k / n) / (K / N) if K > 0 and n > 0 else 0
                results.append(
                    {
                        "Cluster": cluster_id,
                        "Cell Type": cell_type,
                        "p_value": p_val,
                        "Enrichment Ratio": enrichment,
                        "Overlapping Genes Count": k,
                        "Overlapping Genes": ", ".join(sorted(overlap)),
                    }
                )

            if results:
                df = pd.DataFrame(results)
                df["adj_p_value"] = multi.multipletests(df["p_value"], method="fdr_bh")[
                    1
                ]
                all_results_list.append(df)

        if not all_results_list:
            return pd.DataFrame()
        return pd.concat(all_results_list, ignore_index=True).sort_values(
            ["Cluster", "adj_p_value"]
        )

    def plot_results(
        self, p_val_cutoff=0.05, value_col="Enrichment Ratio", top_n_per_cluster=10
    ):
        """
        Generates a dot plot of the enrichment results.

        Args:
            p_val_cutoff (float): P-value cutoff for results to be plotted.
            value_col (str): Column to use for color scale (e.g., 'Enrichment Ratio').
            top_n_per_cluster (int): Max number of cell types to show per cluster, ranked by p-value.
        """
        if self.results_.empty:
            logging.warning("No results to plot.")
            return

        plot_df = self.results_[self.results_["adj_p_value"] < p_val_cutoff].copy()
        if plot_df.empty:
            logging.warning(
                f"No significant results to plot with adjusted p-value cutoff {p_val_cutoff}."
            )
            return

        # Select top N results per cluster to avoid cluttered plots
        plot_df = (
            plot_df.groupby("Cluster")
            .apply(lambda x: x.nsmallest(top_n_per_cluster, "adj_p_value"))
            .reset_index(drop=True)
        )

        if plot_df.empty:
            logging.warning(f"No results remain after selecting top {top_n_per_cluster} per cluster.")
            return

        # Create -log10(p-value) for better color visualization
        plot_df["-log10(adj_p_value)"] = -np.log10(plot_df["adj_p_value"])

        plt.figure(
            figsize=(
                max(8, len(plot_df["Cluster"].unique()) * 0.8),
                max(6, len(plot_df["Cell Type"].unique()) * 0.5),
            )
        )

        scatter = sns.scatterplot(
            data=plot_df,
            x="Cluster",
            y="Cell Type",
            size="Overlapping Genes Count",
            hue=value_col,
            palette="viridis",
            sizes=(50, 500),
            edgecolor="black",
            linewidth=0.5,
        )

        plt.title("Marker Enrichment Analysis Results")
        plt.xlabel("Cluster")
        plt.ylabel("Cell Type")
        plt.xticks(rotation=45, ha="right")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.show()