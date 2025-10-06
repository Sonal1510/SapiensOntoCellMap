#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Annotate each cluster with cell types using hypergeometric test
and generate an HTML report.
"""

import os
import logging
import argparse
import pandas as pd
from config import PROCESSED_COMBINED_DATA_DIR
from get_marker_enrichment_test import MarkerEnrichmentTest
from get_html_report import (
    plot_interactive_summary_heatmap,
    plot_enrichment_dotplot,
    create_deg_tables_html,
    generate_html_report,
)

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def get_cell_name_marker_reference_map(combined_df):
    combined_df = combined_df.assign(
        database_cell_name=combined_df["database"] + "_" + combined_df["cell_name"]
    )
    return combined_df.groupby("database_cell_name")["gene"].apply(list).to_dict()


def get_spaceranger_differential_cluster_file_path_multi_path_key(path_list, sample_name):
    import os

    if not path_list:
        return {}

    common_prefix = os.path.commonpath(path_list)
    reversed_paths = [p.split(os.sep)[::-1] for p in path_list]
    common_suffix_parts = []
    for parts in zip(*reversed_paths):
        if len(set(parts)) == 1:
            common_suffix_parts.append(parts[0])
        else:
            break
    common_suffix = os.sep.join(common_suffix_parts[::-1]) if common_suffix_parts else ""

    path_dict = {}
    for path in path_list:
        relative_to_prefix = os.path.relpath(path, common_prefix)
        middle_part = relative_to_prefix
        if common_suffix and middle_part.endswith(common_suffix):
            middle_part = middle_part[: -len(common_suffix)].strip(os.sep)
        variable_key_part = middle_part.replace(os.sep, "_") or os.path.basename(os.path.dirname(path))
        final_key = f"{sample_name}_{variable_key_part}"
        path_dict[final_key] = path

    return path_dict


def get_spaceranger_differential_cluster_file_path(spaceranger_out_path, sample_name, spatial_method):
    import os

    path_map = {}

    if os.path.isfile(spaceranger_out_path) and spaceranger_out_path.endswith("differential_expression.csv"):
        path_map[sample_name] = spaceranger_out_path
        return path_map

    if os.path.isdir(spaceranger_out_path) and spaceranger_out_path.rstrip(os.sep).endswith("outs"):
        if spatial_method == "xenium":
            deg_file = os.path.join(
                spaceranger_out_path,
                "analysis/diffexp/gene_expression_graphclust/differential_expression.csv",
            )
            if os.path.exists(deg_file):
                path_map[sample_name] = deg_file
            else:
                logger.warning("Differential expression file not found for xenium at %s", deg_file)

        elif spatial_method == "visium":
            paths = []
            for root, _, files in os.walk(spaceranger_out_path):
                if "differential_expression.csv" in files and "_graphclust" in root:
                    paths.append(os.path.join(root, "differential_expression.csv"))

            if paths:
                path_map = get_spaceranger_differential_cluster_file_path_multi_path_key(paths, sample_name)
            else:
                logger.warning("No '*_graphclust/differential_expression.csv' files found in %s", spaceranger_out_path)
        else:
            logger.error("Invalid spatial_method '%s'", spatial_method)
    else:
        logger.error("Invalid input path: %s", spaceranger_out_path)

    return path_map


def main():
    # --- Basic Argument Parsing ---
    parser = argparse.ArgumentParser(description="Run spatial marker enrichment pipeline")
    parser.add_argument("spaceranger_out_path", help="Path to Space Ranger outs directory or differential_expression.csv")
    parser.add_argument("sample_name", help="Sample name for outputs")
    parser.add_argument("spatial_method", choices=["visium", "xenium"], help="Spatial method")
    parser.add_argument("output_dir", help="Directory to store output reports")

    # Optional thresholds (use defaults if not given)
    parser.add_argument("--pval", type=float, default=0.05, help="Adjusted p-value threshold (default: 0.05)")
    parser.add_argument("--log2fc", type=float, default=1.0, help="Log2 fold-change threshold (default: 1.0)")
    parser.add_argument("--mean", type=float, default=0.0, help="Mean count threshold (default: 0)")
    parser.add_argument("--topgenes", type=int, default=None, help="Top N genes to include (default: None)")

    args = parser.parse_args()

    # --- Assign CLI args directly ---
    spaceranger_out_path = args.spaceranger_out_path
    sample_name = args.sample_name
    spatial_method = args.spatial_method
    spatial_anno_outdir = args.output_dir
    p_val_thresh = args.pval
    log2fc_thresh = args.log2fc
    mean_counts_thresh = args.mean
    top_genes = args.topgenes

    os.makedirs(spatial_anno_outdir, exist_ok=True)

    try:
        combined_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
        if not os.path.exists(combined_path):
            raise FileNotFoundError(f"Missing combined DB at: {combined_path}")

        combined_df = pd.read_csv(combined_path)
        combined_df_dict = get_cell_name_marker_reference_map(combined_df)
        logger.info("Cell marker reference map created successfully.")

        path_dict = get_spaceranger_differential_cluster_file_path(spaceranger_out_path, sample_name, spatial_method)
        if not path_dict:
            logger.error("No differential expression files found.")
            return

        for visium_hd_algo, deg_file in path_dict.items():
            logger.info("Processing %s", visium_hd_algo)
            algo_out_dir = os.path.join(spatial_anno_outdir, visium_hd_algo)
            os.makedirs(algo_out_dir, exist_ok=True)

            pipeline = MarkerEnrichmentTest(
                deg_file=deg_file,
                marker_db_dict=combined_df_dict,
                p_val_thresh=p_val_thresh,
                log2fc_thresh=log2fc_thresh,
                mean_counts_thresh=mean_counts_thresh,
                top_genes=top_genes,
            )
            pipeline.fit()

            all_results = pipeline.results_
            sig_results = all_results[all_results["adj_p_value"] < p_val_thresh].copy()

            all_results.to_csv(os.path.join(algo_out_dir, f"{visium_hd_algo}_all_results.csv"), index=False)
            sig_results.to_csv(os.path.join(algo_out_dir, f"{visium_hd_algo}_sig_results.csv"), index=False)

            plots = {
                "heatmap": plot_dynamic_heatmap_with_bars(sig_results, max_celltypes_per_cluster=1),
                "dotplot": plot_enrichment_dotplot(sig_results, top_n=5),
            }
            deg_tables_html = create_deg_tables_html(pipeline.deg_df, pipeline.cluster_markers)

            analysis_params = {
                "p_val": p_val_thresh,
                "log2fc": log2fc_thresh,
                "mean_counts_thresh": mean_counts_thresh,
                "top_genes": top_genes,
            }

            report_path = os.path.join(algo_out_dir, f"{visium_hd_algo}_report.html")
            generate_html_report(
                sample_name=visium_hd_algo,
                output_path=report_path,
                sig_results_df=sig_results,
                plots_html=plots,
                deg_table_html=deg_tables_html,
                params=analysis_params,
            )

            logger.info("Report saved at: %s", report_path)

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)


if __name__ == "__main__":
    main()
