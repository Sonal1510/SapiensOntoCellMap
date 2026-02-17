#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Annotate each cluster with cell types using:
1. Conventional Hypergeometric Test for 'db_cell_name' (Level 1)
2. Weighted Enrichment Analysis for 'cell_name' (Level 2)

(MODIFIED: Adds a Top Annotation Summary CSV output with logic:
 Selected Tissue (L2) > All Tissue (L2))
"""

import os
import logging
import argparse
import pandas as pd
import re
from get_marker_enrichment_test import MarkerEnrichmentTest
from get_html_report import (
    plot_dynamic_heatmap_with_bars,
    create_deg_tables_html,
    generate_html_report
)

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Helper for natural sorting (cluster 2 before cluster 10)
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

def get_weighted_marker_maps(df, level_col):
    """
    Generates a WEIGHTED marker map (Level 2).
    Returns: { 'CellType': {'geneA': 5, 'geneB': 1} }
    """
    df_clean = df.dropna(subset=[level_col, 'gene']).copy()
    if df_clean.empty: return {}

    # Count frequency of gene per cell type
    counts = df_clean.groupby([level_col, 'gene']).size()
    
    weighted_map = {}
    for (cell_type, gene), count in counts.items():
        if cell_type not in weighted_map: weighted_map[cell_type] = {}
        weighted_map[cell_type][gene] = int(count)
        
    return weighted_map

def get_unweighted_marker_maps(df, level_col, distinct_col=None):
    """
    Generates a CONVENTIONAL (Unweighted) marker map (Level 1).
    Returns: { 'DB_CellType': ['geneA', 'geneB'] }
    """
    df_clean = df.dropna(subset=[level_col, 'gene']).copy()
    if df_clean.empty: return {}

    # Group and just get unique list of genes
    unweighted_map = df_clean.groupby(level_col)['gene'].apply(lambda x: list(set(x))).to_dict()
    return unweighted_map

def get_marker_maps_context(df):
    """
    Generates appropriate maps for Level 1 (Unweighted) and Level 2 (Weighted).
    """
    maps = {}
    
    # --- Level 1: Conventional (Unweighted) ---
    # Key: database_cell_name
    df_l1 = df.dropna(subset=['database', 'db_cell_name']).copy()
    if not df_l1.empty:
        df_l1['unique_db_cell_name'] = df_l1["database"] + "_" + df_l1["db_cell_name"]
        maps['level1'] = get_unweighted_marker_maps(df_l1, 'unique_db_cell_name')
    else:
        maps['level1'] = {}

    # --- Level 2: Weighted ---
    # Key: cell_name
    maps['level2'] = get_weighted_marker_maps(df, 'cell_name')
        
    return maps

def get_spaceranger_differential_cluster_file_path_multi_path_key(path_list, sample_name):
    if not path_list: return {}
    split_paths = [p.split(os.sep) for p in path_list]
    common_prefix = os.path.commonpath(path_list)
    reversed_paths = [p.split(os.sep)[::-1] for p in path_list]
    common_suffix_parts = []
    for parts in zip(*reversed_paths):
        if len(set(parts)) == 1: common_suffix_parts.append(parts[0])
        else: break
    common_suffix = os.sep.join(common_suffix_parts[::-1])
    path_dict = {}
    for path in path_list:
        relative_to_prefix = os.path.relpath(path, common_prefix)
        middle_part = relative_to_prefix
        if middle_part.endswith(common_suffix):
            middle_part = middle_part[:-len(common_suffix)].strip(os.sep)
        variable_key_part = middle_part.replace(os.sep, '_')
        if not variable_key_part: variable_key_part = os.path.basename(os.path.dirname(path))
        final_key = f"{sample_name}_{variable_key_part}"
        path_dict[final_key] = path
    return path_dict

def get_spaceranger_differential_cluster_file_path(spaceranger_out_path, sample_name):
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

def run_enrichment_pipeline(deg_file, marker_db_dict, args, deg_type):
    if not marker_db_dict:
        return (None, pd.DataFrame(), pd.DataFrame(), {'heatmap': {}}, 
                {'p_val': args.pval, 'log2fc': args.log2fc, 'mean_counts': args.mean, 'top_n_genes': args.topgenes})

    pipeline = MarkerEnrichmentTest(
        deg_file=deg_file,
        marker_db_dict=marker_db_dict,
        deg_file_type=deg_type,
        p_val_thresh=args.pval,
        log2fc_thresh=args.log2fc,
        mean_counts_thresh=args.mean,
        top_genes=args.topgenes,
        min_overlap=args.min_overlap,
        background_gene_count=args.background_gene_count,
    )
    pipeline.fit()

    all_results = pipeline.results_
    sig_results = pd.DataFrame()
    if not all_results.empty:
        sig_results = all_results[all_results["adj_p_value"] < args.pval].copy()
    
    analysis_params = {
        'p_val': pipeline.p_val_thresh,
        'log2fc': pipeline.log2fc_thresh,
        'mean_counts': pipeline.mean_counts_thresh,
        'top_n_genes': pipeline.top_genes
    }
            
    plots_html = {'heatmap': {}}
    if not sig_results.empty:
        top_n_options_for_plots = [1, 3, 5, 10]
        for n in top_n_options_for_plots:
            plots_html['heatmap'][str(n)] = plot_dynamic_heatmap_with_bars(
                sig_results_df=sig_results,
                top_n_celltypes=n
            )
    else:
        plots_html['heatmap']['1'] = "<h3>Heatmap</h3><p>No significant enrichments found.</p>"
    
    return pipeline, all_results, sig_results, plots_html, analysis_params

def generate_top_annotation_summary(final_sig_results, output_dir, job_name, tissue_specified):
    """
    Generates a CSV with the top 1 annotation per cluster.
    Priority: Selected Tissue (Level 2) > All Tissue (Level 2).
    """
    summary_data = []
    all_clusters = set()

    # Gather clusters from any context that produced results
    def get_clusters_from_ctx(ctx, lvl):
        if ctx in final_sig_results and lvl in final_sig_results[ctx]:
            df = final_sig_results[ctx][lvl]
            if not df.empty and 'Cluster' in df.columns:
                return set(df['Cluster'].unique())
        return set()

    all_clusters.update(get_clusters_from_ctx('selected_tissue', 'level2'))
    all_clusters.update(get_clusters_from_ctx('all_tissue', 'level2'))
    
    # If we still haven't found clusters (e.g. nothing significant anywhere), 
    # we might want to check the raw DEGs, but sticking to sig results implies 
    # we only list clusters that had *some* result.
    
    sorted_clusters = sorted(list(all_clusters), key=natural_sort_key)
    
    for cluster in sorted_clusters:
        selected_hit = None
        
        # 1. Try Selected Tissue (Level 2)
        if tissue_specified and 'selected_tissue' in final_sig_results:
            df = final_sig_results['selected_tissue'].get('level2')
            if df is not None and not df.empty:
                cluster_res = df[df['Cluster'] == cluster]
                if not cluster_res.empty:
                    # Sort: Adj P (asc), then Weighted Enrichment (desc)
                    sort_cols = ['adj_p_value']
                    asc = [True]
                    if 'Weighted_Enrichment' in cluster_res.columns:
                        sort_cols.append('Weighted_Enrichment')
                        asc.append(False)
                    else:
                        sort_cols.append('Enrichment_ratio')
                        asc.append(False)
                        
                    cluster_res = cluster_res.sort_values(by=sort_cols, ascending=asc)
                    top_row = cluster_res.iloc[0]
                    
                    score = top_row.get('Weighted_Enrichment', top_row.get('Enrichment_ratio'))
                    selected_hit = {
                        'Cluster': cluster,
                        'Top_Cell_Type': top_row['Cell_type'],
                        'Source': 'Selected Tissue',
                        'P_Value': top_row['adj_p_value'],
                        'Score': score,
                        'Genes': top_row.get('Overlapping_genes', '')
                    }

        # 2. If no hit, try All Tissue (Level 2)
        if selected_hit is None:
             if 'all_tissue' in final_sig_results:
                df = final_sig_results['all_tissue'].get('level2')
                if df is not None and not df.empty:
                    cluster_res = df[df['Cluster'] == cluster]
                    if not cluster_res.empty:
                        sort_cols = ['adj_p_value']
                        asc = [True]
                        if 'Weighted_Enrichment' in cluster_res.columns:
                            sort_cols.append('Weighted_Enrichment')
                            asc.append(False)
                        else:
                            sort_cols.append('Enrichment_ratio')
                            asc.append(False)
                            
                        cluster_res = cluster_res.sort_values(by=sort_cols, ascending=asc)
                        top_row = cluster_res.iloc[0]
                        score = top_row.get('Weighted_Enrichment', top_row.get('Enrichment_ratio'))
                        
                        selected_hit = {
                            'Cluster': cluster,
                            'Top_Cell_Type': top_row['Cell_type'],
                            'Source': 'All Tissue',
                            'P_Value': top_row['adj_p_value'],
                            'Score': score,
                            'Genes': top_row.get('Overlapping_genes', '')
                        }

        if selected_hit:
            summary_data.append(selected_hit)
        else:
            summary_data.append({
                'Cluster': cluster,
                'Top_Cell_Type': 'Unannotated',
                'Source': 'None',
                'P_Value': None, 'Score': None, 'Genes': ''
            })

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        out_path = os.path.join(output_dir, f"{job_name}_top_annotation_summary.csv")
        summary_df.to_csv(out_path, index=False)
        logger.info(f"Top annotation summary saved to: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Run hybrid enrichment pipeline")
    parser.add_argument("input_path", help="Path to DEG CSV or Space Ranger dir")
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("output_dir", help="Directory to store output reports")
    parser.add_argument("--deg_type", choices=["spatial", "scrna"], required=True)
    parser.add_argument("--pval", type=float, default=0.05, help="Adj p-value threshold")
    parser.add_argument("--log2fc", type=float, default=1.0, help="Log2FC threshold")
    parser.add_argument("--mean", type=float, default=0.0, help="Mean count threshold")
    parser.add_argument("--topgenes", type=int, default=None, help="Top N genes")
    parser.add_argument("--marker_db", required=True, help="Path to marker DB CSV")
    parser.add_argument("--tissue", type=str, default=None, help="Filter marker DB for tissue")
    parser.add_argument("--min_overlap", type=int, default=2, help="Minimum gene overlap count to report (default: 2)")
    parser.add_argument("--background_gene_count", type=int, default=None, help="Override background gene count N for hypergeometric test")
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # 1. Load DB
        full_df = pd.read_csv(args.marker_db, dtype=str)
        
        # 2. Prepare Contexts
        contexts = {}
        
        logger.info("Preparing maps for 'All Tissue'...")
        contexts['all_tissue'] = get_marker_maps_context(full_df)
        
        tissue_specified = False
        if args.tissue:
            logger.info(f"Preparing maps for Selected Tissue: '{args.tissue}'...")
            subset_df = full_df[full_df['tissue_name'].str.contains(args.tissue, case=False, na=False)]
            if not subset_df.empty:
                contexts['selected_tissue'] = get_marker_maps_context(subset_df)
                tissue_specified = True
            else:
                logger.warning(f"Tissue '{args.tissue}' not found. Skipping.")

        # 3. Find Files
        path_dict = find_deg_files(args.input_path, args.sample_name, args.deg_type)
        if not path_dict: return

        # 4. Process
        for job_name, deg_file in path_dict.items():
            logger.info(f"Processing Job: {job_name}")
            algo_out_dir = os.path.join(args.output_dir, job_name)
            os.makedirs(algo_out_dir, exist_ok=True)

            final_sig = {}
            final_plots = {}
            final_params = {}
            valid_pipe = None

            # Run pipeline for each context/level
            for ctx_name, maps in contexts.items():
                final_sig[ctx_name] = {}
                final_plots[ctx_name] = {}
                final_params[ctx_name] = {}

                for lvl_name, m_map in maps.items():
                    logger.info(f"--- Context: {ctx_name} | Level: {lvl_name} ---")
                    
                    (pipe, all_r, sig_r, plots, params) = run_enrichment_pipeline(
                        deg_file, m_map, args, args.deg_type
                    )

                    if pipe:
                        prefix = f"{job_name}_{ctx_name}_{lvl_name}"
                        all_r.to_csv(os.path.join(algo_out_dir, f"{prefix}_all_results.csv"), index=False)
                        sig_r.to_csv(os.path.join(algo_out_dir, f"{prefix}_sig_results.csv"), index=False)
                        
                        final_sig[ctx_name][lvl_name] = sig_r
                        final_plots[ctx_name][lvl_name] = plots
                        final_params[ctx_name][lvl_name] = params
                        if valid_pipe is None: valid_pipe = pipe

            # --- Generate Top Annotation Summary CSV ---
            # This uses the priority logic: Selected Tissue (Level 2) > All Tissue (Level 2)
            generate_top_annotation_summary(
                final_sig_results=final_sig,
                output_dir=algo_out_dir,
                job_name=job_name,
                tissue_specified=tissue_specified
            )

            # --- DEG Tables ---
            if not valid_pipe:
                deg_html = "<p>Error: No DEGs processed.</p>"
            else:
                deg_html = create_deg_tables_html(
                    deg_df=valid_pipe.raw_deg_df,
                    cluster_markers=valid_pipe.cluster_markers,
                    p_val_thresh=valid_pipe.p_val_thresh,
                    log2fc_thresh=valid_pipe.log2fc_thresh,
                    mean_counts_thresh=valid_pipe.mean_counts_thresh
                )

            # --- HTML Report ---
            report_path = os.path.join(algo_out_dir, f"{job_name}_report.html")
            generate_html_report(
                sample_name=job_name,
                output_path=report_path,
                sig_results_maps=final_sig, 
                plots_html_maps=final_plots,
                deg_table_html=deg_html,
                params_maps=final_params,
                selected_tissue_name=args.tissue if 'selected_tissue' in contexts else None
            )
            logger.info(f"Report: {report_path}")

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)

if __name__ == "__main__":
    main()