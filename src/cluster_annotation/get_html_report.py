#!/usr/bin/env python3
"""
Author: Sonal Rashmi
Date: 2025-10-06
Description:
A script to perform enrichment analysis and generate a comprehensive, interactive
HTML report. The report includes a clustered heatmap, an interactive dot plot of
significant enrichments, and a browser for differentially expressed genes (DEGs).
"""

import pandas as pd
import numpy as np
import os
import io
import base64
import re
import logging
import matplotlib
matplotlib.use('Agg') # Use a non-interactive backend for scripts
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from jinja2 import Template
from typing import Dict, List
import plotly.express as px
import plotly.graph_objects as go

# Setup logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

# A small constant to prevent log(0) errors.
_EPSILON = 1e-300

# ==============================================================================
# HELPER FUNCTION FOR SORTING
# ==============================================================================
def natural_sort_key(s: str):
    """
    Creates a key for sorting strings with numbers in a natural way.
    Example: 'Cluster 10' comes after 'Cluster 2'.
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

# ==============================================================================
# CLUSTERED HEATMAP (MODIFIED WITH NEW DYNAMIC LOGIC)
# ==============================================================================
def plot_dynamic_heatmap_with_bars(
    sig_results_df: pd.DataFrame,
    top_n_celltypes: int = 1,
    cluster_axes: bool = True
) -> str:
    """
    Generates a dynamic, publication-quality heatmap based on raw enrichment ratios.
    """
    if sig_results_df.empty:
        return f"<h3>Heatmap (Top {top_n_celltypes})</h3><p>No significant enrichments found.</p>"

    required_cols = ['Cluster', 'Cell Type', 'adj_p_value', 'Enrichment Ratio']
    if not all(col in sig_results_df.columns for col in required_cols):
        raise ValueError(f"Input DataFrame is missing one of the required columns: {required_cols}")

    df = sig_results_df.copy()
    df['Enrichment Ratio'] = pd.to_numeric(df['Enrichment Ratio'], errors='coerce').fillna(0)
    df['adj_p_value'] = pd.to_numeric(df['adj_p_value'], errors='coerce')
    df_sorted = df.sort_values(
        by=['Cluster', 'adj_p_value', 'Enrichment Ratio'],
        ascending=[True, True, False]
    )
    top_hits_per_cluster = df_sorted.groupby('Cluster').head(top_n_celltypes)
    top_cell_types_list = top_hits_per_cluster['Cell Type'].unique().tolist()
    top_sig_results = df[df['Cell Type'].isin(top_cell_types_list)]
    heatmap_df = top_sig_results.pivot_table(
        index='Cluster',
        columns='Cell Type',
        values='Enrichment Ratio',
        fill_value=0
    )
    try:
        all_clusters = df['Cluster'].unique()
        sorted_clusters = sorted(all_clusters, key=natural_sort_key)
        heatmap_df = heatmap_df.reindex(sorted_clusters, fill_value=0)
    except (ValueError, IndexError):
        heatmap_df = heatmap_df.sort_index()

    cbar_label = 'Enrichment Ratio'

    if cluster_axes and not heatmap_df.empty:
        if heatmap_df.shape[1] > 1:
            # --- CHANGE 1: Use correlation metric for better biological clustering ---
            col_linkage = linkage(heatmap_df.T, method='ward', metric='correlation')
            col_dendrogram = dendrogram(col_linkage, no_plot=True, labels=heatmap_df.columns)
            heatmap_df = heatmap_df[col_dendrogram['ivl']]
        if heatmap_df.shape[0] > 1:
            # --- CHANGE 1: Use correlation metric for better biological clustering ---
            row_linkage = linkage(heatmap_df, method='ward', metric='correlation')
            row_dendrogram = dendrogram(row_linkage, no_plot=True, labels=heatmap_df.index)
            heatmap_df = heatmap_df.reindex(row_dendrogram['ivl'])

    cluster_freq = sig_results_df.groupby('Cluster')['Cell Type'].nunique().reindex(heatmap_df.index).fillna(0)
    fig_height = max(8, len(heatmap_df.index) * 0.3)
    fig_width = max(10, len(heatmap_df.columns) * 0.5)
    fig = plt.figure(figsize=(fig_width, fig_height), layout="constrained")
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 5], wspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharey=ax_bar)

    sns.heatmap(
        heatmap_df, ax=ax_heat, cmap='rocket_r', linewidths=.5, linecolor='lightgray',
        cbar_kws={'label': cbar_label}, yticklabels=False, xticklabels=True
    )
    ax_heat.set_title(
        f'Significant Cluster Annotations (adj p < 0.05)\nTop {top_n_celltypes} Cell Type(s) per Cluster',
        fontsize=14, weight='bold', pad=20
    )
    ax_heat.set_xlabel('Annotated Cell Type', fontsize=12, labelpad=10)
    ax_heat.set_ylabel('')
    plt.setp(ax_heat.get_xticklabels(), rotation=90, ha='right', rotation_mode='anchor', fontsize=10)

    y_pos = np.arange(len(heatmap_df.index))
    y_pos_centered = y_pos + 0.5
    bars = ax_bar.barh(y_pos_centered, cluster_freq.values, height=0.8, color="steelblue", edgecolor="black", linewidth=0.5)
    ax_bar.set_yticks(y_pos_centered)
    ax_bar.set_yticklabels(heatmap_df.index, fontsize=10, va='center')
    ax_bar.set_xlabel("Significant Hits", fontsize=11)
    ax_bar.set_ylabel("Cluster", fontsize=12)
    ax_bar.invert_yaxis()
    ax_bar.grid(axis='x', linestyle='--', alpha=0.6)
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax_bar.text(width * 1.01, bar.get_y() + bar.get_height()/2., '%d' % int(width), ha='left', va='center', fontsize=8)
    ax_bar.set_xlim(right=cluster_freq.max() * 1.15 if cluster_freq.max() > 0 else 1)

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=200)
    plt.close(fig)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    return f'<img src="data:image/png;base64,{img_base64}" alt="Enrichment Heatmap" style="width:100%; height:auto;">'

# ==============================================================================
# DEG BROWSER PLOTS & TABLES
# ==============================================================================
def _reshape_deg_df(deg_df: pd.DataFrame) -> pd.DataFrame:
    # This function remains the same as the previous version
    if 'cluster' in deg_df.columns and 'gene' in deg_df.columns:
        logging.info("Input DEG dataframe appears to be in long format. Standardizing columns.")
        df_reshaped = deg_df.copy()
        df_reshaped.rename(columns={'gene': 'Feature Name', 'cluster': 'Cluster', 'p_val_adj': 'adj_p_value', 'avg_log2FC': 'log2FC'}, inplace=True)
        if 'mean_counts' not in df_reshaped.columns: df_reshaped['mean_counts'] = 0
        if 'Feature ID' not in df_reshaped.columns: df_reshaped['Feature ID'] = df_reshaped['Feature Name']
        df_reshaped['Cluster'] = "Cluster " + df_reshaped['Cluster'].astype(str)
        return df_reshaped
    p_val_cols, log2fc_cols, mean_counts_cols = ([c for c in deg_df.columns if s in c] for s in ['Adjusted p value', 'Log2 fold change', 'Mean Counts'])
    if not p_val_cols or not log2fc_cols: raise ValueError("Could not find required columns for reshaping.")
    id_vars = ['Feature ID', 'Feature Name']
    value_vars = p_val_cols + log2fc_cols + mean_counts_cols
    df_long = pd.melt(deg_df, id_vars=id_vars, value_vars=value_vars, var_name='metric', value_name='value')
    df_long[['Cluster', 'Metric Type']] = df_long['metric'].str.extract(r'(Cluster \d+)\s(.*)')
    df_reshaped = df_long.pivot_table(index=id_vars + ['Cluster'], columns='Metric Type', values='value', aggfunc='first').reset_index()
    df_reshaped.rename(columns={'Adjusted p value': 'adj_p_value', 'Log2 fold change': 'log2FC', 'Mean Counts': 'mean_counts'}, inplace=True)
    if 'mean_counts' not in df_reshaped.columns: df_reshaped['mean_counts'] = 0
    for col, dtype in {'adj_p_value': 'float', 'log2FC': 'float', 'mean_counts': 'float'}.items():
        if col in df_reshaped.columns: df_reshaped[col] = pd.to_numeric(df_reshaped[col], errors='coerce')
    df_reshaped['mean_counts'].fillna(0, inplace=True)
    return df_reshaped

def create_deg_violin_plots(deg_df_reshaped: pd.DataFrame, p_val_thresh: float, log2fc_thresh: float, mean_counts_thresh: float) -> str:
    """
    Generate interactive violin plots for adj p-value, log2FC, and mean counts per cluster.
    Log2FC and Mean Counts plots now use a logarithmic y-axis.
    A threshold line is added for Mean Counts.
    """
    filtered_deg_df = deg_df_reshaped[
        (deg_df_reshaped['adj_p_value'] < p_val_thresh) & (abs(deg_df_reshaped['log2FC']) >= log2fc_thresh)
    ].copy()

    if filtered_deg_df.empty:
        return "<h3>DEG Distributions</h3><p>No significant DEGs found with the given thresholds.</p>"

    filtered_deg_df['neg_log10_p_value'] = -np.log10(filtered_deg_df['adj_p_value'] + _EPSILON)
    filtered_deg_df['abs_log2FC'] = filtered_deg_df['log2FC'].abs()
    
    hover_template = '<b>%{customdata[0]}</b><br>Cluster: %{x}<br>Adj p-value: %{customdata[1]:.2e}<br>Log2FC: %{customdata[2]:.2f}<br>Mean Counts: %{customdata[3]:.2f}<extra></extra>'

    sorted_clusters = sorted(filtered_deg_df['Cluster'].unique(), key=natural_sort_key)
    filtered_deg_df['Cluster'] = pd.Categorical(filtered_deg_df['Cluster'], categories=sorted_clusters, ordered=True)
    filtered_deg_df.sort_values('Cluster', inplace=True)
    
    # p-value plot (y-axis is already a log-like scale)
    fig_p_val = px.violin(filtered_deg_df, x='Cluster', y='neg_log10_p_value', box=True, points=False,
                            title='Adjusted p-value Distribution', color='Cluster',
                            custom_data=['Feature Name', 'adj_p_value', 'log2FC', 'mean_counts'])
    fig_p_val.update_traces(hovertemplate=hover_template, box_visible=True, meanline_visible=True)
    fig_p_val.add_hline(y=-np.log10(p_val_thresh), line_dash="dash", line_color="red", annotation_text=f"p-value threshold = {p_val_thresh}", annotation_position="bottom right")
    fig_p_val.update_layout(xaxis_title='', yaxis_title='-log₁₀(Adjusted p-value)')

    # --- CHANGE 2: Add log_y=True for exponential y-axis ---
    fig_log2fc = px.violin(filtered_deg_df, x='Cluster', y='abs_log2FC', box=True, points=False,
                            title='Absolute Log₂ Fold Change Distribution', color='Cluster',
                            custom_data=['Feature Name', 'adj_p_value', 'log2FC', 'mean_counts'],
                            log_y=True)
    fig_log2fc.update_traces(hovertemplate=hover_template, box_visible=True, meanline_visible=True)
    fig_log2fc.add_hline(y=log2fc_thresh, line_dash="dash", line_color="red", annotation_text=f"log2FC threshold = {log2fc_thresh}", annotation_position="bottom right")
    fig_log2fc.update_layout(xaxis_title='', yaxis_title='Absolute Log₂ Fold Change (Log Scale)')
    
    p_val_html = fig_p_val.to_html(full_html=False, include_plotlyjs='cdn')
    log2fc_html = fig_log2fc.to_html(full_html=False, include_plotlyjs=False)
    
    mean_counts_html = ""
    if 'mean_counts' in filtered_deg_df.columns and filtered_deg_df['mean_counts'].sum() > 1e-6:
        # --- CHANGE 2: Add log_y=True for exponential y-axis ---
        fig_mean_counts = px.violin(filtered_deg_df, x='Cluster', y='mean_counts', box=True, points=False,
                                      title='Mean Counts Distribution', color='Cluster',
                                      custom_data=['Feature Name', 'adj_p_value', 'log2FC', 'mean_counts'],
                                      log_y=True)
        fig_mean_counts.update_traces(hovertemplate=hover_template, box_visible=True, meanline_visible=True)
        
        # --- CHANGE 3: Add horizontal line for mean counts threshold ---
        if mean_counts_thresh > 0:
            fig_mean_counts.add_hline(y=mean_counts_thresh, line_dash="dash", line_color="red", annotation_text=f"mean counts threshold = {mean_counts_thresh:.2f}", annotation_position="bottom right")

        fig_mean_counts.update_layout(xaxis_title='Cluster', yaxis_title='Mean Counts (Log Scale)')
        mean_counts_html = fig_mean_counts.to_html(full_html=False, include_plotlyjs=False)

    combined_html = f"""
    <div style="text-align: center;"><h3>Adjusted p-value Distribution</h3>{p_val_html}</div>
    <div style="text-align: center;"><h3>Absolute Log₂ Fold Change Distribution</h3>{log2fc_html}</div>
    """
    if mean_counts_html:
        combined_html += f"""
        <div style="text-align: center;"><h3>Mean Counts Distribution</h3>{mean_counts_html}</div>
        """
    return combined_html


def create_deg_tables_html(deg_df: pd.DataFrame, cluster_markers: Dict[str, List[str]], p_val_thresh: float, log2fc_thresh: float, mean_counts_thresh: float) -> str:
    deg_df_reshaped = _reshape_deg_df(deg_df)
    # Pass the new mean_counts_thresh parameter
    dist_plots_html = create_deg_violin_plots(deg_df_reshaped, p_val_thresh, log2fc_thresh, mean_counts_thresh)
    bar_plot_html = plot_deg_counts_barchart(cluster_markers)
    # ... (rest of the function is unchanged)
    dropdown_parts = ['<label for="cluster_select"><b>Select a Cluster to view its DEGs:</b></label>', '<select id="cluster_select" onchange="showTable(this.value)">']
    cluster_names = sorted(cluster_markers.keys(), key=natural_sort_key)
    dropdown_parts.append('<option value="">--Select--</option>')
    for cluster in cluster_names:
        sanitized_cluster_name = re.sub(r'\s+', '_', str(cluster))
        dropdown_parts.append(f'<option value="deg_table_{sanitized_cluster_name}">{cluster} ({len(cluster_markers.get(cluster,[]))} genes)</option>')
    dropdown_parts.append('</select>')
    table_parts = []
    for cluster in cluster_names:
        genes = cluster_markers.get(cluster, [])
        if not genes: continue
        sanitized_cluster_name = re.sub(r'\s+', '_', str(cluster))
        if 'avg_log2FC' in deg_df.columns:
             cluster_num = str(cluster).split()[-1]
             cluster_deg_df = deg_df.loc[deg_df['gene'].isin(genes) & (deg_df['cluster'].astype(str) == cluster_num)].copy()
             rename_map = {'gene': 'Gene', 'p_val_adj': 'Adjusted p-value', 'avg_log2FC': 'Log2 Fold Change'}
             cluster_deg_df = cluster_deg_df[[k for k in rename_map if k in cluster_deg_df.columns]]
        else:
            cols_to_get = ['Feature Name', f"{cluster} Adjusted p value", f"{cluster} Log2 fold change", f"{cluster} Mean Counts"]
            cols_exist = [col for col in cols_to_get if col in deg_df.columns]
            cluster_deg_df = deg_df.loc[deg_df['Feature Name'].isin(genes), cols_exist].copy()
            rename_map = {'Feature Name': 'Gene', f"{cluster} Adjusted p value": 'Adjusted p-value', f"{cluster} Log2 fold change": 'Log2 Fold Change', f"{cluster} Mean Counts": 'Mean Counts'}
        cluster_deg_df.rename(columns=rename_map, inplace=True)
        table_html = cluster_deg_df.to_html(classes="display compact", index=False, table_id=f"deg_table_{sanitized_cluster_name}_data", float_format='{:.2e}'.format)
        table_parts.append(f'<div id="deg_table_{sanitized_cluster_name}" class="deg-table-container" style="display:none;"><h4>DEGs for {cluster}</h4>{table_html}</div>')
    return dist_plots_html + bar_plot_html + ''.join(dropdown_parts) + ''.join(table_parts)

# The rest of the file (plot_deg_counts_barchart, generate_html_report) remains unchanged.
def plot_deg_counts_barchart(cluster_markers: Dict[str, List[str]]) -> str:
    # ... (unchanged)
    if not cluster_markers: return ""
    deg_counts = pd.Series({cluster: len(genes) for cluster, genes in cluster_markers.items()})
    deg_counts = deg_counts[deg_counts > 0]
    if deg_counts.empty: return "<h3>DEGs per Cluster</h3><p>No differentially expressed genes found for any cluster.</p>"
    sorted_index = sorted(deg_counts.index, key=natural_sort_key)
    deg_counts = deg_counts.reindex(sorted_index)
    fig = px.bar(x=deg_counts.index, y=deg_counts.values, title='Number of Genes per Cluster for Hypergeometric Test', labels={'x': 'Cluster', 'y': 'Number of Genes'}, color=deg_counts.index, text=deg_counts.values)
    fig.update_layout(showlegend=False)
    fig.update_traces(marker_color='steelblue')
    return f'<h3>DEGs per Cluster</h3>{fig.to_html(full_html=False)}<hr style="margin: 25px 0;">'

def generate_html_report(sample_name, output_path, sig_results_df, plots_html, deg_table_html, params):
    # ... (unchanged)
    df_for_html = sig_results_df.copy()
    if 'Cluster' in df_for_html.columns:
        unique_clusters = df_for_html['Cluster'].unique()
        sorted_clusters = sorted(unique_clusters, key=natural_sort_key)
        cat_type = pd.api.types.CategoricalDtype(categories=sorted_clusters, ordered=True)
        df_for_html['Cluster'] = df_for_html['Cluster'].astype(cat_type)
        df_for_html.sort_values(by=['Cluster', 'adj_p_value'], inplace=True)
    has_genes_col = 'Overlapping Genes' in df_for_html.columns
    if has_genes_col:
        gene_limit = 10
        def truncate_for_table(gene_string):
            if not isinstance(gene_string, str) or not gene_string: return ""
            genes = [g.strip() for g in gene_string.replace(';',',').split(',') if g.strip()]
            return f"{', '.join(genes[:gene_limit])}, ..." if len(genes) > gene_limit else ', '.join(genes)
        df_for_html['Full Gene List'] = df_for_html['Overlapping Genes']
        df_for_html['Overlapping Genes'] = df_for_html['Full Gene List'].apply(truncate_for_table)
        df_for_html.insert(0, '', '')
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
        <style> body{font-family:Arial,sans-serif;margin:0;padding:0;background-color:#f8f9fa} .container{width:95%;margin:20px auto;padding:20px;background-color:#fff;box-shadow:0 4px 8px #0000001a;border-radius:8px} header{border-bottom:3px solid #007bff;padding-bottom:15px;margin-bottom:20px} h1,h2,h3{color:#0056b3} h2{border-bottom:2px solid #e9ecef;padding-bottom:8px;margin-top:40px} .tabs{display:flex;border-bottom:1px solid #ccc} .tab-link{padding:10px 20px;cursor:pointer;background:#f1f1f1;border-bottom:none} .tab-link.active{background:#fff;border:1px solid #ccc;border-bottom:1px solid #fff;position:relative;top:1px} .tab-content{display:none;padding:20px;border:1px solid #ddd;border-top:none} .tab-content.active{display:block} .params{background-color:#e9ecef;padding:15px;border-radius:5px;margin-bottom:20px} .plot-controls{margin-bottom:15px;} .plot-controls label{font-weight:bold; margin-right:10px;} td.dt-control { background: url('https://datatables.net/examples/resources/details_open.png') no-repeat center center; cursor: pointer; } tr.dt-hasChild td.dt-control { background: url('https://datatables.net/examples/resources/details_close.png') no-repeat center center; } </style>
    </head><body>
    <div class="container">
        <header><h1>Enrichment Analysis Report</h1><h2>Sample: <strong>{{ sample_name }}</strong></h2></header>
        <div class="params"><strong>DEG Selection Parameters:</strong> Adj. p-value &le; {{ params.p_val }} | Log2FC &ge; {{ params.log2fc }} | Mean Counts &ge; {{ params.mean_counts }}{% if params.top_n_genes and params.top_n_genes > 0 %} | Top {{ params.top_n_genes }} Genes per Cluster{% endif %}</div>
        <div class="tabs"> <div class="tab-link active" onclick="openTab(event, 'degs')">DEG Browser</div> <div class="tab-link" onclick="openTab(event, 'visuals')">Enrichment Visuals</div> <div class="tab-link" onclick="openTab(event, 'results')">Hypergeometric Result</div> </div>
        <div id="degs" class="tab-content active">{{ deg_tables|safe }}</div>
        <div id="visuals" class="tab-content">
            <div class="plot-controls"> <label for="heatmap_top_n_select">Heatmap Top N:</label> <select id="heatmap_top_n_select" onchange="showPlot('heatmap', this.value)"> {% for n_key in plots.heatmap.keys() %}<option value="{{ n_key }}" {% if loop.first %}selected{% endif %}>{{ n_key }}</option>{% endfor %} </select> </div>
            {% for n_key, plot_html in plots.heatmap.items() %}<div id="heatmap_{{ n_key }}" class="plot-container heatmap-plot" style="display: {% if loop.first %}block{% else %}none{% endif %};">{{ plot_html|safe }}</div>{% endfor %}
        </div>
        <div id="results" class="tab-content">
            <div style="margin-bottom: 15px;"> <label for="top_n_select"><b>Show Top N Hits per Cluster:</b></label> <select id="top_n_select" onchange="filterTopNResults()"><option value="1">1</option><option value="3">3</option><option value="5" selected>5</option><option value="10">10</option><option value="all">Show All</option></select> </div>
            {{ sig_table|safe }}
        </div>
    </div>
    <script> function openTab(e,t){let n,c,l;for(c=document.getElementsByClassName("tab-content"),n=0;n<c.length;n++)c[n].style.display="none";for(l=document.getElementsByClassName("tab-link"),n=0;n<l.length;n++)l[n].className=l[n].className.replace(" active","");document.getElementById(t).style.display="block",e.currentTarget.className+=" active"} function showTable(id){$(".deg-table-container").hide();if(id){$("#"+id.replace(/[\\s'"]/g, '_')).show();var tableId="#"+id.replace(/[\\s'"]/g, '_')+"_data";if(!$.fn.DataTable.isDataTable(tableId)){$(tableId).DataTable({pageLength:10,dom:"Bfrtip",buttons:["copy","csv"]})}}} function showPlot(e,t){for(var l=document.getElementsByClassName(e+"-plot"),a=0;a<l.length;a++)l[a].style.display="none";document.getElementById(e+"_"+t).style.display="block"} var resultsTable,allResultsData=[]; function filterTopNResults(){if(!resultsTable)return;var e=$("#top_n_select").val();resultsTable.rows().every(function(){this.child.isShown()&&this.child.hide()});if("all"===e){resultsTable.clear().rows.add(allResultsData).draw();return}var t=parseInt(e,10),l=resultsTable.columns().header().toArray().map(e=>$(e).text()),a=l.indexOf("Cluster"),r=l.indexOf("adj_p_value"),s=l.indexOf("Enrichment Ratio");if(-1===a||-1===r||-1===s)return console.error("Required columns for filtering not found."),void resultsTable.clear().rows.add(allResultsData).draw();var o=allResultsData.reduce((e,t)=>(e[t[a]]=e[t[a]]||[],e[t[a]].push(t),e),{}),n=[];var c=Object.keys(o);c.sort((e,t)=>e.localeCompare(t,void 0,{numeric:!0,sensitivity:"base"}));for(const i of c){var d=o[i];d.sort((e,t)=>{let l=parseFloat(e[r])-parseFloat(t[r]);return 0!==l?l:parseFloat(t[s])-parseFloat(e[s])}),n.push(...d.slice(0,t))}resultsTable.clear().rows.add(n).draw()} $(document).ready(function(){var e={{'true' if has_genes_col else 'false'}};if($("#results_table").length){var t=[],l=-1;if(e){var a=Array.from($("#results_table thead th")).map(e=>$(e).text());(l=a.indexOf("Full Gene List"))>-1&&t.push({targets:l,visible:!1}),t.push({targets:0,className:"dt-control",orderable:!1,data:null,defaultContent:""})}resultsTable=$("#results_table").DataTable({pageLength:25,dom:"Bfrtip",buttons:["copy","csv","excel"],columnDefs:t,order:e?[[1,"asc"]]:[[0,"asc"]]}),allResultsData=resultsTable.rows().data().toArray(),filterTopNResults(),e&&l>-1&&$("#results_table tbody").on("click","td.dt-control",function(){var e=$(this).closest("tr"),t=resultsTable.row(e);if(t.child.isShown())t.child.hide(),e.removeClass("dt-hasChild");else{var a=resultsTable.cell(t.index(),l).data(),s='<div style="padding: 8px 30px; background-color: #f9f9f9; border-left: 3px solid #007bff;"><b>Full List of Overlapping Genes:</b><br><p style="word-break: break-word; white-space: normal; margin-top: 5px;">'+a+"</p></div>";t.child(s).show(),e.addClass("dt-hasChild")}}) }openTab({currentTarget:document.querySelector(".tab-link.active")},"degs")}); </script>
    </body></html>
    """
    html_content = Template(template_str).render(sample_name=sample_name, plots=plots_html, sig_table=df_for_html.to_html(classes="display compact", index=False, table_id="results_table"), deg_tables=deg_table_html, params=params, has_genes_col=has_genes_col)
    with open(output_path, "w", encoding='utf-8') as f: f.write(html_content)