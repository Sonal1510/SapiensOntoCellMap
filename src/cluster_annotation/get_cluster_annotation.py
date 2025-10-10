#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Annotate each cluster with cell types using hypergeometric test
and generate an HTML report for spatial or scRNA-seq data.
"""

import os
import logging
import argparse
import pandas as pd
# Assuming config.py and other scripts are in the python path
# from config import PROCESSED_COMBINED_DATA_DIR
PROCESSED_COMBINED_DATA_DIR = "../spatial_annotator/SapiensOntoCellMap/data/processed_combined_db/" # Placeholder if config is not available
from get_marker_enrichment_test import MarkerEnrichmentTest
from get_html_report import (
    plot_dynamic_heatmap_with_bars,
    create_deg_tables_html,
    generate_html_report,
    natural_sort_key
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
    """
    Generates a dictionary of paths where the key is the variable part of the path,
    prefixed with the sample name for uniqueness.
    """
    if not path_list:
        return {}

    split_paths = [p.split(os.sep) for p in path_list]

    # Find common prefix using os.path.commonpath for simplicity and reliability
    common_prefix = os.path.commonpath(path_list)

    # Find common suffix
    reversed_paths = [p.split(os.sep)[::-1] for p in path_list]
    common_suffix_parts = []
    for parts in zip(*reversed_paths):
        if len(set(parts)) == 1:
            common_suffix_parts.append(parts[0])
        else:
            break
    common_suffix = os.sep.join(common_suffix_parts[::-1])

    path_dict = {}
    for path in path_list:
        relative_to_prefix = os.path.relpath(path, common_prefix)

        middle_part = relative_to_prefix
        # More robustly remove the suffix
        if middle_part.endswith(common_suffix):
            middle_part = middle_part[:-len(common_suffix)].strip(os.sep)

        variable_key_part = middle_part.replace(os.sep, '_')
        if not variable_key_part:
            variable_key_part = os.path.basename(os.path.dirname(path))

        # Combine the variable part with the sample name for a unique key
        final_key = f"{sample_name}_{variable_key_part}"
        
        path_dict[final_key] = path

    return path_dict

def get_spaceranger_differential_cluster_file_path(spaceranger_out_path, sample_name):
    """
    Finds the path(s) to differential expression CSV files from a Space Ranger output directory.
    """
    path_map = {}

    if spaceranger_out_path.endswith("differential_expression.csv"):
        path_map[sample_name] = spaceranger_out_path
    elif os.path.isdir(spaceranger_out_path):
        paths = []
        for root, dirs, files in os.walk(spaceranger_out_path):
            if 'differential_expression.csv' in files and '_graphclust' in root:
                full_path = os.path.join(root, 'differential_expression.csv')
                paths.append(full_path)
            
        if paths:
            path_map = get_spaceranger_differential_cluster_file_path_multi_path_key(paths, sample_name)
        else:
            print(f"Warning: No '*_graphclust/differential_expression.csv' files found in {spaceranger_out_path}")
    else:
        print(f"Error: Path is not a valid directory ending in 'outs' or a direct CSV file path: {spaceranger_out_path}")

    return path_map

def find_deg_files(input_path, sample_name, deg_type):
    path_map = {}
    if deg_type == 'scrna':
        if os.path.isfile(input_path) and input_path.endswith(".csv"):
            path_map[sample_name] = input_path
        else:
            logger.error(f"For --deg_type scrna, the input path must be a valid .csv file. Got: {input_path}")
    else:
        if os.path.isdir(input_path):
            path_map = get_spaceranger_differential_cluster_file_path(input_path, sample_name)    
    if not path_map:
        logger.warning(f"No spatial differential_expression.csv files found in {input_path}")
        
    return path_map

def main():
    parser = argparse.ArgumentParser(description="Run marker enrichment pipeline for spatial or scRNA-seq data")
    parser.add_argument("input_path", help="Path to DEG CSV file (for scrna) or Space Ranger outs directory (for spatial)")
    parser.add_argument("sample_name", help="Sample name for outputs")
    parser.add_argument("output_dir", help="Directory to store output reports")
    parser.add_argument("--deg_type", choices=["spatial", "scrna"], required=True, help="Type of DEG file: 'spatial' (wide format) or 'scrna' (long format)")
    parser.add_argument("--pval", type=float, default=0.05, help="Adjusted p-value threshold (default: 0.05)")
    parser.add_argument("--log2fc", type=float, default=1.0, help="Log2 fold-change threshold (default: 1.0)")
    parser.add_argument("--mean", type=float, default=0.0, help="Mean count threshold (default: 0)")
    parser.add_argument("--topgenes", type=int, default=None, help="Top N genes to include (default: None)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        combined_path = os.path.join(PROCESSED_COMBINED_DATA_DIR, "master_cell_marker_db.csv")
        if not os.path.exists(combined_path):
             raise FileNotFoundError(f"Missing combined DB at: {combined_path}")
        combined_df = pd.read_csv(combined_path)
        combined_df_dict = get_cell_name_marker_reference_map(combined_df)
        logger.info("Cell marker reference map created successfully.")

        path_dict = find_deg_files(args.input_path, args.sample_name, args.deg_type)
        if not path_dict:
            logger.error("No differential expression files found.")
            return

        for job_name, deg_file in path_dict.items():
            logger.info("Processing %s", job_name)
            algo_out_dir = os.path.join(args.output_dir, job_name)
            os.makedirs(algo_out_dir, exist_ok=True)

            pipeline = MarkerEnrichmentTest(
                deg_file=deg_file,
                marker_db_dict=combined_df_dict,
                deg_file_type=args.deg_type,
                p_val_thresh=args.pval,
                log2fc_thresh=args.log2fc,
                mean_counts_thresh=args.mean,
                top_genes=args.topgenes,
            )
            pipeline.fit()

            all_results = pipeline.results_
            sig_results = all_results[all_results["adj_p_value"] < args.pval].copy()

            all_results.to_csv(os.path.join(algo_out_dir, f"{job_name}_all_results.csv"), index=False)
            sig_results.to_csv(os.path.join(algo_out_dir, f"{job_name}_sig_results.csv"), index=False)

            analysis_params = {
                'p_val': pipeline.p_val_thresh,
                'log2fc': pipeline.log2fc_thresh,
                'mean_counts': pipeline.mean_counts_thresh,
                'top_n_genes': pipeline.top_genes
            }
            
            logger.info("Generating plots and tables for HTML report...")
            plots_html = {'heatmap': {}}
            top_n_options_for_plots = [1, 3, 5, 10]

            for n in top_n_options_for_plots:
                logger.info(f"  - Generating heatmap for Top {n}...")
                plots_html['heatmap'][str(n)] = plot_dynamic_heatmap_with_bars(
                    sig_results_df=sig_results,
                    top_n_celltypes=n
                )
            
            # --- CHANGE: Pass the mean_counts_thresh to the HTML function ---
            deg_tables_html = create_deg_tables_html(
                deg_df=pipeline.raw_deg_df,
                cluster_markers=pipeline.cluster_markers,
                p_val_thresh=pipeline.p_val_thresh,
                log2fc_thresh=pipeline.log2fc_thresh,
                mean_counts_thresh=pipeline.mean_counts_thresh
            )

            report_path = os.path.join(algo_out_dir, f"{job_name}_report.html")
            
            generate_html_report(
                sample_name=job_name,
                output_path=report_path,
                sig_results_df=sig_results,
                plots_html=plots_html,
                deg_table_html=deg_tables_html,
                params=analysis_params
            )
            logger.info(f"Report generated at: {report_path}")

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)


if __name__ == "__main__":
    main()