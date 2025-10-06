#!/usr/bin/env python3
"""
Author: Sonal Rashmi
Date: 2025-10-06
Description:
A script to perform enrichment analysis and generate a comprehensive, interactive
HTML report. The report includes a clustered heatmap, an interactive dot plot of
significant enrichments, and a browser for differentially expressed genes (DEGs).
"""

# ==============================================================================
# 1. IMPORT LIBRARIES
# ==============================================================================

# --- Standard Libraries ---
import os
import io
import base64
from typing import Dict, List

# --- Third-party Libraries ---
import pandas as pd
import numpy as np
from jinja2 import Template

# --- Plotting Libraries ---
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend for script execution
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import plotly.express as px

# --- Scientific Computing Libraries ---
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, dendrogram


# ==============================================================================
# 2. PLOTTING FUNCTIONS
# ==============================================================================

def plot_dynamic_heatmap_with_bars(
    sig_results_df: pd.DataFrame,
    max_celltypes_per_cluster: int = 1,
    normalize: str = 'row',
    cluster_axes: bool = True
) -> str:
    """
    Generates a publication-quality heatmap with clustering, normalization, and an adjacent bar plot.

    Args:
        sig_results_df (pd.DataFrame): DataFrame with significant enrichment results.
                                       Must contain ['Cluster', 'Cell Type', 'adj_p_value', 'Enrichment Ratio'].
        max_celltypes_per_cluster (int): The number of top enriched cell types to display per cluster,
                                         ranked by adjusted p-value. Defaults to 1.
        normalize (str): Method for normalizing heatmap data ('row', 'col', or None). Defaults to 'row'.
        cluster_axes (bool): If True, hierarchically cluster rows and columns. Defaults to True.

    Returns:
        str: An HTML <img> tag containing the base64-encoded PNG of the plot.
    """
    if sig_results_df.empty:
        return "<h3>Enrichment Heatmap</h3><p>No significant enrichments found.</p>"

    # --- Data Preparation ---
    # Ensure required columns exist
    required_cols = ['Cluster', 'Cell Type', 'adj_p_value', 'Enrichment Ratio']
    if not all(col in sig_results_df.columns for col in required_cols):
        raise ValueError(f"Input DataFrame is missing required columns: {required_cols}")

    # Filter data to top N cell types per cluster based on p-value
    df_top = (
        sig_results_df.copy()
        .sort_values(['Cluster', 'adj_p_value'])
        .groupby('Cluster', group_keys=False)
        .head(max_celltypes_per_cluster)
    )

    # Pivot the table to create a matrix: Clusters (rows) vs. Cell Types (columns)
    heatmap_df = df_top.pivot_table(
        index='Cluster',
        columns='Cell Type',
        values='Enrichment Ratio',
        fill_value=0
    )

    # Attempt to sort clusters numerically for a more intuitive default order
    try:
        sorted_clusters = sorted(heatmap_df.index, key=lambda x: int(x.split(' ')[-1]))
        cat_type = pd.api.types.CategoricalDtype(categories=sorted_clusters, ordered=True)
        heatmap_df.index = heatmap_df.index.astype(cat_type)
        heatmap_df = heatmap_df.sort_index()
    except (ValueError, IndexError):
        heatmap_df = heatmap_df.sort_index()  # Fallback for non-standard cluster names

    # --- Normalization (Optional) ---
    cbar_label = 'Enrichment Ratio'
    if normalize and not heatmap_df.empty:
        if normalize == 'row':
            heatmap_df = pd.DataFrame(zscore(heatmap_df, axis=1), index=heatmap_df.index, columns=heatmap_df.columns)
            cbar_label = 'Row-wise Z-score'
        elif normalize == 'col':
            heatmap_df = pd.DataFrame(zscore(heatmap_df, axis=0), index=heatmap_df.index, columns=heatmap_df.columns)
            cbar_label = 'Column-wise Z-score'
        heatmap_df.fillna(0, inplace=True)  # z-score can produce NaNs if std dev is 0

    # --- Hierarchical Clustering (Optional) ---
    if cluster_axes and not heatmap_df.empty:
        # Cluster columns (cell types) if there's more than one
        if heatmap_df.shape[1] > 1:
            col_linkage = linkage(heatmap_df.T, method='ward', metric='euclidean')
            col_dendrogram = dendrogram(col_linkage, no_plot=True, labels=heatmap_df.columns)
            heatmap_df = heatmap_df[col_dendrogram['ivl']]  # Reorder columns

        # Cluster rows (clusters) if there's more than one
        if heatmap_df.shape[0] > 1:
            row_linkage = linkage(heatmap_df, method='ward', metric='euclidean')
            row_dendrogram = dendrogram(row_linkage, no_plot=True, labels=heatmap_df.index)
            heatmap_df = heatmap_df.reindex(row_dendrogram['ivl'])  # Reorder rows

    # Calculate hit counts per cluster, ensuring order matches the final heatmap
    cluster_freq = sig_results_df.groupby('Cluster')['Cell Type'].nunique().reindex(heatmap_df.index).fillna(0)

    # --- Figure Generation ---
    fig_height = max(8, len(heatmap_df.index) * 0.3)
    fig_width = max(10, len(heatmap_df.columns) * 0.5)
    fig = plt.figure(figsize=(fig_width, fig_height), layout="constrained")

    # Create a grid with space for the bar plot and the heatmap
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 5], wspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharey=ax_bar) # Share Y-axis for alignment

    # Plot the heatmap
    sns.heatmap(
        heatmap_df, ax=ax_heat, cmap='viridis', linewidths=.5, linecolor='lightgray',
        cbar_kws={'label': cbar_label}, yticklabels=False, xticklabels=True
    )
    ax_heat.set_title(
        f'Significant Cluster Annotations (adj p < 0.05)\nTop {max_celltypes_per_cluster} Cell Type(s) per Cluster',
        fontsize=14, weight='bold', pad=20
    )
    ax_heat.set_xlabel('Annotated Cell Type', fontsize=12, labelpad=10)
    ax_heat.set_ylabel('')
    plt.setp(ax_heat.get_xticklabels(), rotation=90, ha='right', rotation_mode='anchor', fontsize=10)

    # Plot the horizontal bar chart for hit counts
    y_pos = np.arange(len(heatmap_df.index))
    bars = ax_bar.barh(y_pos, cluster_freq.values, height=0.8, color="steelblue", edgecolor="black", linewidth=0.5)
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(heatmap_df.index, fontsize=10)
    ax_bar.set_xlabel("Significant Hits", fontsize=11)
    ax_bar.set_ylabel("Spatial Cluster", fontsize=12)
    ax_bar.invert_yaxis()  # Match the heatmap's top-to-bottom order
    ax_bar.grid(axis='x', linestyle='--', alpha=0.6)

    # Add text labels to the bars
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax_bar.text(width * 1.01, bar.get_y() + bar.get_height() / 2., '%d' % int(width),
                        ha='left', va='center', fontsize=8)
    ax_bar.set_xlim(right=cluster_freq.max() * 1.15) # Add space for labels

    # --- Export to Base64 ---
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=200)
    plt.close(fig)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')

    return f'<h3>Enrichment Heatmap</h3><img src="data:image/png;base64,{img_base64}" alt="Enrichment Heatmap" style="width:100%; height:auto;">'

def plot_enrichment_dotplot(sig_results_df: pd.DataFrame, top_n: int = 5) -> str:
    """
    Generates an interactive enrichment dot plot using Plotly.
    Handles missing hover fields defensively to prevent errors.

    Args:
        sig_results_df (pd.DataFrame): DataFrame with significant enrichment results.
        top_n (int): Number of top cell types to show per cluster. Defaults to 5.

    Returns:
        str: An HTML div containing the Plotly graph.
    """
    if sig_results_df.empty:
        return "<h3>Enrichment Dot Plot</h3><p>No significant results to plot.</p>"

    # Select top N results per cluster by adjusted p-value
    df = sig_results_df.sort_values("adj_p_value").groupby("Cluster").head(top_n).copy()
    if df.empty:
        return "<h3>Enrichment Dot Plot</h3><p>No data remains after filtering for top N.</p>"

    # Calculate size metric for dots, adding epsilon for stability
    df['-log10(p_adj)'] = -np.log10(df['adj_p_value'] + 1e-300)

    # Dynamically prepare hover data: only include columns that exist in the dataframe
    hover_data_config = {}
    available_cols = set(df.columns)
    preferred_hover = {
        'Enrichment Ratio': ':.2f',
        'adj_p_value': ':.2e',
        'Overlapping Genes Count': True,
        'Overlapping Genes': True,
        'k': True, 'n': True, 'K': True, 'N': True  # Hypergeometric test parameters
    }
    for col, fmt in preferred_hover.items():
        if col in available_cols:
            hover_data_config[col] = fmt

    # Conditionally set color parameter
    color_col = 'Enrichment Ratio' if 'Enrichment Ratio' in df.columns else None

    # Create the interactive scatter plot
    fig = px.scatter(
        df,
        x='Cluster',
        y='Cell Type',
        size='-log10(p_adj)',
        color=color_col,
        color_continuous_scale='viridis' if color_col else None,
        hover_name='Cell Type',
        hover_data=hover_data_config,
        labels={
            'Cluster': '<b>Cluster</b>',
            'Cell Type': '<b>Annotated Cell Type</b>',
            '-log10(p_adj)': '-log10(Adj. P-Value)'
        },
        title=f'<b>Top {top_n} Enriched Cell Types per Cluster</b>'
    )

    # Adjust layout for a clean, professional look
    fig.update_layout(
        xaxis_type='category',
        font=dict(family="Arial"),
        plot_bgcolor='white',
        height=max(400, df['Cell Type'].nunique() * 25), # Dynamic height
        coloraxis_colorbar=dict(title='Enrichment Ratio')
    )

    # Return partial HTML for embedding in the final report
    html_output = fig.to_html(full_html=False, include_plotlyjs=False)
    return f"<h3>Enrichment Dot Plot</h3>{html_output}"


# ==============================================================================
# 3. HTML COMPONENT GENERATORS
# ==============================================================================

def create_deg_tables_html(deg_df: pd.DataFrame, cluster_markers: Dict[str, List[str]]) -> str:
    """
    Creates interactive HTML tables for browsing DEGs per cluster.

    Args:
        deg_df (pd.DataFrame): DataFrame containing all DEG statistics.
        cluster_markers (Dict[str, List[str]]): Dictionary mapping cluster names to lists of marker genes.

    Returns:
        str: HTML string with a dropdown selector and hidden tables.
    """
    # Start with the dropdown menu
    html_parts = [
        '<label for="cluster_select"><b>Select a Cluster to view its DEGs:</b></label>',
        '<select id="cluster_select" onchange="showTable(this.value)">'
    ]
    cluster_names = sorted(cluster_markers.keys())
    html_parts.append('<option value="">--Select--</option>')
    for i, cluster in enumerate(cluster_names):
        num_genes = len(cluster_markers.get(cluster, []))
        html_parts.append(f'<option value="deg_table_{i}">{cluster} ({num_genes} genes)</option>')
    html_parts.append('</select>')

    # Create a hidden div with a DataTable for each cluster
    for i, cluster in enumerate(cluster_names):
        genes = cluster_markers.get(cluster, [])
        if not genes:
            continue

        # Select and rename columns for the specific cluster
        cols_to_get = ['Feature Name', f"{cluster} Adjusted p value", f"{cluster} Log2 fold change", f"{cluster} Mean Counts"]
        cols_exist = [col for col in cols_to_get if col in deg_df.columns]
        cluster_deg_df = deg_df.loc[deg_df['Feature Name'].isin(genes), cols_exist].copy()
        rename_map = {
            'Feature Name': 'Gene',
            f"{cluster} Adjusted p value": 'Adjusted p-value',
            f"{cluster} Log2 fold change": 'Log2 Fold Change',
            f"{cluster} Mean Counts": 'Mean Counts'
        }
        cluster_deg_df.rename(columns=rename_map, inplace=True)

        # Convert DataFrame to an HTML table with a unique ID
        table_html = cluster_deg_df.to_html(classes="display compact", index=False, table_id=f"deg_table_{i}_data", float_format='{:.2e}'.format)
        html_parts.append(f'<div id="deg_table_{i}" class="deg-table-container" style="display:none;"><h4>DEGs for {cluster}</h4>{table_html}</div>')

    return ''.join(html_parts)

def generate_html_report(sample_name, output_path, sig_results_df, plots_html, deg_table_html, params):
    """
    Generates a standalone HTML report from analysis results and plots using a Jinja2 template.

    Args:
        sample_name (str): The name of the sample being analyzed.
        output_path (str): The file path to save the HTML report.
        sig_results_df (pd.DataFrame): DataFrame of significant enrichment results for the full table tab.
        plots_html (Dict[str, str]): Dictionary containing HTML for the plots.
        deg_table_html (str): HTML string for the DEG browser.
        params (Dict[str, float]): Dictionary of analysis parameters (e.g., p_val, log2fc).
    """
    # Jinja2 template for the HTML report
    template_str = """
    <!DOCTYPE html>
    <html><head><title>{{ sample_name }} Report</title>
        <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"/>
        <script type="text/javascript" src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">
        <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
        <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style> body{font-family:Arial,sans-serif;margin:0;padding:0;background-color:#f8f9fa}.container{width:95%;margin:20px auto;padding:20px;background-color:#fff;box-shadow:0 4px 8px #0000001a;border-radius:8px}header{border-bottom:3px solid #007bff;padding-bottom:15px;margin-bottom:20px}h1,h2,h3{color:#0056b3}h2{border-bottom:2px solid #e9ecef;padding-bottom:8px;margin-top:40px}.tabs{display:flex;border-bottom:1px solid #ccc}.tab-link{padding:10px 20px;cursor:pointer;background:#f1f1f1;border-bottom:none}.tab-link.active{background:#fff;border:1px solid #ccc;border-bottom:1px solid #fff;position:relative;top:1px}.tab-content{display:none;padding:20px;border:1px solid #ddd;border-top:none}.tab-content.active{display:block}.params{background-color:#e9ecef;padding:15px;border-radius:5px;margin-bottom:20px} </style>
    </head><body>
    <div class="container">
        <header><h1>Enrichment Analysis Report</h1><h2>Sample: <strong>{{ sample_name }}</strong></h2></header>
        <div class="params"><strong>Analysis Parameters:</strong> p-value &le; {{ params.p_val }} | log2FC &ge; {{ params.log2fc }}</div>
        <div class="tabs">
            <div class="tab-link active" onclick="openTab(event, 'summary')">Enrichment Plots</div>
            <div class="tab-link" onclick="openTab(event, 'degs')">DEG Browser</div>
            <div class="tab-link" onclick="openTab(event, 'results')">Full Results</div>
        </div>
        <div id="summary" class="tab-content active">{{ plots.heatmap|safe }}<hr>{{ plots.dotplot|safe }}</div>
        <div id="degs" class="tab-content">{{ deg_tables|safe }}</div>
        <div id="results" class="tab-content">{{ sig_table|safe }}</div>
    </div>
    <script>function openTab(e,t){let n,c,l;for(c=document.getElementsByClassName("tab-content"),n=0;n<c.length;n++)c[n].style.display="none";for(l=document.getElementsByClassName("tab-link"),n=0;n<l.length;n++)l[n].className=l[n].className.replace(" active","");document.getElementById(t).style.display="block",e.currentTarget.className+=" active"}function showTable(id){$(".deg-table-container").hide();if(id){$("#"+id).show();var tableId="#"+id+"_data";if(!$.fn.DataTable.isDataTable(tableId)){$(tableId).DataTable({pageLength:10,dom:"Bfrtip",buttons:["copy","csv"]})}}} $(document).ready(function(){$("#results_table").DataTable({pageLength:25,dom:"Bfrtip",buttons:["copy","csv","excel"]});openTab({currentTarget:document.querySelector(".tab-link.active")},"summary")});</script>
    </body></html>
    """
    # Render the template with the provided data
    html_content = Template(template_str).render(
        sample_name=sample_name,
        plots=plots_html,
        sig_table=sig_results_df.to_html(classes="display compact", index=False, table_id="results_table"),
        deg_tables=deg_table_html,
        params=params
    )
    # Write the rendered HTML to the output file
    with open(output_path, "w", encoding='utf-8') as f:
        f.write(html_content)
    print(f"[INFO] Report generated at: {output_path}")
