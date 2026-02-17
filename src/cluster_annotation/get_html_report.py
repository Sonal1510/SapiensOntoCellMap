#!/usr/bin/env python3
"""
Author: Sonal Rashmi
Date: 2025-10-06
Description:
A script to perform enrichment analysis and generate a comprehensive, interactive
HTML report.

(MODIFIED: Visuals reverted to 'Old' style, Weighted Logic retained)
"""

import pandas as pd
import numpy as np
import io
import base64
import re
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram
from jinja2 import Template
import plotly.express as px
from typing import Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

def plot_dynamic_heatmap_with_bars(
    sig_results_df: pd.DataFrame,
    top_n_celltypes: int = 1,
    cluster_axes: bool = True
) -> str:
    """
    Generates a heatmap at 150 DPI (screen-optimized) with download button, maintaining the old style/colors.
    """
    if sig_results_df.empty:
        return f"<h3>Heatmap (Top {top_n_celltypes})</h3><p>No significant enrichments found.</p>"

    df = sig_results_df.copy()
    
    # --- DYNAMIC METRIC SELECTION ---
    if 'Weighted_Enrichment' in df.columns:
        value_col = 'Weighted_Enrichment'
        title_metric = "Weighted Enrichment"
    else:
        value_col = 'Enrichment_ratio'
        title_metric = "Enrichment Ratio"
        
    required = ['Cluster', 'Cell_type', 'adj_p_value', value_col]
    if not all(col in df.columns for col in required):
         return f"<p>Error: Missing columns {required}</p>"

    df[value_col] = pd.to_numeric(df[value_col], errors='coerce').fillna(0)
    df['adj_p_value'] = pd.to_numeric(df['adj_p_value'], errors='coerce')
    
    # Sort by P-value first, then Score
    df_sorted = df.sort_values(
        by=['Cluster', 'adj_p_value', value_col],
        ascending=[True, True, False]
    )
    
    top_hits_per_cluster = df_sorted.groupby('Cluster').head(top_n_celltypes)
    top_cell_types_list = top_hits_per_cluster['Cell_type'].unique().tolist()
    top_sig_results = df[df['Cell_type'].isin(top_cell_types_list)]
    
    heatmap_df = top_sig_results.pivot_table(
        index='Cluster', columns='Cell_type', values=value_col, fill_value=0
    )
    
    try:
        sorted_clusters = sorted(df['Cluster'].unique(), key=natural_sort_key)
        # Filter to only existing clusters in heatmap
        existing = [c for c in sorted_clusters if c in heatmap_df.index]
        heatmap_df = heatmap_df.reindex(existing, fill_value=0)
    except:
        heatmap_df = heatmap_df.sort_index()

    if cluster_axes and not heatmap_df.empty:
        rows_var = heatmap_df.index[heatmap_df.std(axis=1) > 0]
        cols_var = heatmap_df.columns[heatmap_df.std(axis=0) > 0]
        clustered = heatmap_df.loc[rows_var, cols_var]
        
        final_row_order = list(heatmap_df.index)
        final_col_order = list(heatmap_df.columns)

        if len(cols_var) > 1:
            try:
                Z = linkage(clustered.T, method='average', metric='correlation')
                D = dendrogram(Z, no_plot=True, labels=clustered.columns)
                final_col_order = D['ivl'] + [c for c in heatmap_df.columns if c not in D['ivl']]
            except: pass

        if len(rows_var) > 1:
            try:
                Z = linkage(clustered, method='average', metric='correlation')
                D = dendrogram(Z, no_plot=True, labels=clustered.index)
                final_row_order = D['ivl'] + [r for r in heatmap_df.index if r not in D['ivl']]
            except: pass
            
        heatmap_df = heatmap_df.reindex(index=final_row_order, columns=final_col_order)

    cluster_freq = sig_results_df.groupby('Cluster')['Cell_type'].nunique().reindex(heatmap_df.index).fillna(0)
    
    # Dynamic dimensions
    fig_h = max(8, len(heatmap_df.index) * 0.3)
    fig_w = max(10, len(heatmap_df.columns) * 0.5)
    
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 5], wspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharey=ax_bar)

    # --- CRITICAL FIX: Wipe index name completely to prevent Seaborn from grabbing it ---
    heatmap_df.index.name = ""
    # -----------------------------------------------------------------------------------

    sns.heatmap(
        heatmap_df, ax=ax_heat, cmap='rocket_r', linewidths=.5, linecolor='lightgray',
        cbar_kws={'label': title_metric}, yticklabels=False, xticklabels=True
    )
    
    # --- [FIX APPLIED HERE] ---
    # Explicitly turn off y-ticks/labels on the heatmap to prevent double labeling
    ax_heat.tick_params(axis='y', which='both', left=False, labelleft=False)
    # --------------------------

    ax_heat.set_title(
        f'Significant Annotations (adj p < 0.05)\nTop {top_n_celltypes} per Cluster',
        fontsize=14, weight='bold', pad=20
    )
    ax_heat.set_xlabel('Cell Type', fontsize=12)
    
    # --- CRITICAL FIX: Force removal of label on heatmap axis AFTER plotting ---
    ax_heat.set_ylabel("") 
    ax_heat.yaxis.set_label_text("")
    # -------------------------------------------------------------------------
    
    plt.setp(ax_heat.get_xticklabels(), rotation=90, ha='right', rotation_mode='anchor', fontsize=10)

    y_pos = np.arange(len(heatmap_df.index)) + 0.5
    bars = ax_bar.barh(y_pos, cluster_freq.values, height=0.8, color="steelblue", edgecolor="black")
    ax_bar.set_yticks(y_pos)
    ax_bar.set_yticklabels(heatmap_df.index)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("# Hits")
    
    # Label is ONLY here on the bar chart
    ax_bar.set_ylabel("Cluster", fontsize=12) 
    
    ax_bar.grid(axis='x', linestyle='--', alpha=0.6)
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax_bar.text(width * 1.01, bar.get_y() + bar.get_height()/2., '%d' % int(width), ha='left', va='center', fontsize=8)
            
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode('utf-8')
    return (
        f'<div style="position:relative;">'
        f'<img src="data:image/png;base64,{b64}" loading="lazy" style="width:100%; height:auto;">'
        f'<br><a download="heatmap.png" href="data:image/png;base64,{b64}" '
        f'style="display:inline-block;margin:8px 0;padding:6px 16px;background:#007bff;'
        f'color:#fff;border-radius:4px;text-decoration:none;font-size:13px;">'
        f'Download Heatmap</a></div>'
    )

def _reshape_deg_df(deg_df: pd.DataFrame) -> pd.DataFrame:
    # Standardize DEG DF to long format
    if 'cluster' in deg_df.columns and 'gene' in deg_df.columns:
        df = deg_df.copy()
        df.rename(columns={'gene':'Feature Name','cluster':'Cluster','p_val_adj':'adj_p_value','avg_log2FC':'log2FC'}, inplace=True)
        if 'mean_counts' not in df.columns: df['mean_counts'] = 0
        if 'Feature ID' not in df.columns: df['Feature ID'] = df['Feature Name']
        df['Cluster'] = "Cluster " + df['Cluster'].astype(str)
        return df
    
    p_cols = [c for c in deg_df.columns if 'Adjusted p value' in c]
    fc_cols = [c for c in deg_df.columns if 'Log2 fold change' in c]
    mc_cols = [c for c in deg_df.columns if 'Mean Counts' in c]
    
    if not p_cols or not fc_cols:
         raise ValueError("Missing P-value or Log2FC columns in spatial DEG file.")

    id_vars = ['Feature ID', 'Feature Name']
    df_long = pd.melt(deg_df, id_vars=id_vars, value_vars=p_cols+fc_cols+mc_cols, var_name='metric', value_name='value')
    df_long[['Cluster', 'Type']] = df_long['metric'].str.extract(r'(Cluster \d+)\s(.*)')
    df_long.dropna(subset=['Cluster', 'Type'], inplace=True)
    
    df = df_long.pivot_table(index=id_vars+['Cluster'], columns='Type', values='value', aggfunc='first').reset_index()
    df.rename(columns={'Adjusted p value':'adj_p_value', 'Log2 fold change':'log2FC', 'Mean Counts':'mean_counts'}, inplace=True)
    df.fillna({'mean_counts':0, 'adj_p_value':1, 'log2FC':0}, inplace=True)
    return df

def create_deg_violin_plots(deg_df_reshaped, p_val_thresh, log2fc_thresh, mean_counts_thresh):
    filtered_deg_df = deg_df_reshaped
    if filtered_deg_df.empty:
        return "<h3>DEG Distributions</h3><p>No significant DEGs found with the given thresholds.</p>"

    # Sort clusters naturally for consistent plot order
    sorted_clusters = sorted(filtered_deg_df['Cluster'].unique(), key=natural_sort_key)
    filtered_deg_df['Cluster'] = pd.Categorical(filtered_deg_df['Cluster'], categories=sorted_clusters, ordered=True)
    filtered_deg_df.sort_values('Cluster', inplace=True)

    custom_data_cols = ['Feature Name', 'adj_p_value', 'log2FC', 'mean_counts']
    hover_template = '<b>%{customdata[0]}</b><br>Cluster: %{x}<br>Adj p-value: %{customdata[1]:.2e}<br>Log2FC: %{customdata[2]:.2f}<br>Mean Counts: %{customdata[3]:.2f}<extra></extra>'

    # Adaptive display: for large datasets, subsample per cluster to keep violin shape
    # accurate while cutting Plotly JSON size by ~90%
    MAX_POINTS_PER_CLUSTER = 1000
    if len(filtered_deg_df) > 5000:
        point_mode = 'outliers'
        sampled_parts = []
        for _, grp in filtered_deg_df.groupby('Cluster', observed=True):
            sampled_parts.append(grp.sample(n=min(len(grp), MAX_POINTS_PER_CLUSTER), random_state=42))
        plot_df = pd.concat(sampled_parts, ignore_index=True)
        # Restore categorical ordering after groupby
        plot_df['Cluster'] = pd.Categorical(plot_df['Cluster'], categories=sorted_clusters, ordered=True)
    else:
        point_mode = 'all'
        plot_df = filtered_deg_df

    # --- Plot 1: Adjusted p-value (Log Scale with Visual Floor) ---
    fig_p_val = px.violin(
        plot_df, x='Cluster', y='adj_p_value', box=True, points=point_mode,
        title='Adjusted p-value Distribution', color='Cluster',
        custom_data=custom_data_cols
    )
    fig_p_val.update_traces(
        hovertemplate=hover_template, box_visible=True, meanline_visible=True,
        points=point_mode, jitter=0.5, pointpos=0, marker_opacity=0.6, marker_size=4
    )
    p_val_floor = 1e-100
    non_zero_pvals = plot_df['adj_p_value'][plot_df['adj_p_value'] > 0]
    yaxis_pval_dict = dict(title_text='Adjusted p-value (Log Scale)', type="log", tickformat=".1e")
    if not non_zero_pvals.empty:
        min_pval_for_display = max(non_zero_pvals.min(), p_val_floor)
        max_pval_range = non_zero_pvals.max() + 10
        yaxis_pval_dict['range'] = [np.log10(min_pval_for_display), np.log10(max_pval_range)]
    fig_p_val.update_layout(xaxis_title='', yaxis=yaxis_pval_dict)
    fig_p_val.add_hline(
        y=p_val_thresh, line_dash="dash", line_color="red",
        annotation_text=f"p-value threshold = {p_val_thresh}", annotation_position="top left"
    )

    # --- Plot 2: Log2 Fold Change (Symlog) ---
    C = 1.0
    symlog_transform = lambda x: np.sign(x) * np.log1p(np.abs(x / C))
    plot_df['log2FC_transformed'] = symlog_transform(plot_df['log2FC'])

    fig_log2fc = px.violin(
        plot_df, x='Cluster', y='log2FC_transformed', box=True, points=point_mode,
        title='Log2 Fold Change Distribution', color='Cluster',
        custom_data=custom_data_cols
    )
    fig_log2fc.update_traces(
        hovertemplate=hover_template, box_visible=True, meanline_visible=True,
        points=point_mode, jitter=0.5, pointpos=0, marker_opacity=0.6, marker_size=4
    )

    max_abs_val = plot_df['log2FC'].abs().max()
    tick_values_orig = [0]
    if max_abs_val > 0 and log2fc_thresh > 0:
        positive_ticks = [10**i for i in range(int(np.floor(np.log10(log2fc_thresh))), int(np.ceil(np.log10(max_abs_val)))+1)]
        tick_values_orig.extend(p for p in positive_ticks if p > 0)
        tick_values_orig.extend([-p for p in positive_ticks if p > 0])
    tick_values_orig = sorted(list(set(tick_values_orig)))
    tick_values_transformed = [symlog_transform(v) for v in tick_values_orig]

    min_fc = plot_df['log2FC'].min()
    max_fc = plot_df['log2FC'].max()
    range_min = symlog_transform(min_fc - 0.1)
    range_max = symlog_transform(max_fc + 1)

    fig_log2fc.update_layout(
        xaxis_title='',
        yaxis=dict(
            title_text='Log₂ Fold Change (Symlog Scale)',
            tickvals=tick_values_transformed,
            ticktext=[f"{v:g}" for v in tick_values_orig],
            range=[range_min, range_max]
        )
    )
    fig_log2fc.add_hline(
        y=symlog_transform(log2fc_thresh), line_dash="dash", line_color="red",
        annotation_text=f"Upregulated Log\u2082FC threshold = {log2fc_thresh}", annotation_position="bottom right"
    )

    # --- Plot 3: Mean Counts (Box Plot, publication quality) ---
    mean_counts_html = ""
    if 'mean_counts' in plot_df.columns and plot_df['mean_counts'].sum() > 1e-6:
        fig_mean_counts = px.box(
            plot_df, x='Cluster', y='mean_counts',
            color='Cluster', points='outliers', custom_data=custom_data_cols
        )
        fig_mean_counts.update_traces(
            hovertemplate=hover_template, marker=dict(opacity=0.7, size=4)
        )
        counts_floor = 1e-4
        non_zero_counts = plot_df['mean_counts'][plot_df['mean_counts'] > 0]
        yaxis_counts_dict = dict(title_text='Mean Counts (Log Scale)', type="log")
        if not non_zero_counts.empty:
            q_low = non_zero_counts.quantile(0.005)
            q_high = non_zero_counts.quantile(0.995)
            min_for_display = max(q_low, counts_floor)
            range_min_val = min_for_display * 0.5
            range_max_val = q_high * 2.0
            if range_max_val > range_min_val:
                yaxis_counts_dict['range'] = [np.log10(range_min_val), np.log10(range_max_val)]
                min_power = np.floor(np.log10(range_min_val))
                max_power = np.ceil(np.log10(range_max_val))
                tick_values = [10**i for i in range(int(min_power), int(max_power) + 1)]
                tick_values = [v for v in tick_values if range_min_val <= v <= range_max_val]
                if tick_values:
                    yaxis_counts_dict['tickvals'] = tick_values
                    yaxis_counts_dict['ticktext'] = [f'{v:g}' for v in tick_values]
        fig_mean_counts.update_layout(
            title_text='Distribution of Mean Gene Counts per Cluster', xaxis_title='Cluster',
            yaxis=yaxis_counts_dict, showlegend=False,
            template='simple_white', font=dict(size=12), title=dict(x=0.5, font=dict(size=16))
        )
        if mean_counts_thresh > 0:
            fig_mean_counts.add_hline(
                y=mean_counts_thresh, line_dash="dash", line_color="red", line_width=2,
                layer="above",
                annotation_text=f"Mean counts threshold = {mean_counts_thresh:.2f}",
                annotation_position="top right"
            )
        mean_counts_html = fig_mean_counts.to_html(full_html=False, include_plotlyjs=False)

    p_val_html = fig_p_val.to_html(full_html=False, include_plotlyjs=False) # JS already loaded in header CDN
    log2fc_html = fig_log2fc.to_html(full_html=False, include_plotlyjs=False)

    combined_html = f"""
    <div style="text-align: center;"><h3>Adjusted p-value Distribution</h3>{p_val_html}</div>
    <div style="text-align: center;"><h3>Log2 Fold Change Distribution</h3>{log2fc_html}</div>
    """
    if mean_counts_html:
        combined_html += f"""
        <div style="text-align: center;"><h3>Mean Counts Distribution</h3>{mean_counts_html}</div>
        """
    return combined_html

def plot_deg_counts_barchart(markers):
    if not markers: return ""
    s = pd.Series({k:len(v) for k,v in markers.items()})
    s = s[s > 0]
    if s.empty: return "<h3>DEGs per Cluster</h3><p>No differentially expressed genes found for any cluster.</p>"
    s = s.reindex(sorted(s.index, key=natural_sort_key))
    fig = px.bar(x=s.index, y=s.values,
                 title='Number of Genes per Cluster for Hypergeometric Test',
                 labels={'x':'Cluster', 'y':'Number of Genes'}, color=s.index, text=s.values)
    fig.update_layout(showlegend=False)
    fig.update_traces(marker_color='steelblue')
    return f'<h3>DEGs per Cluster</h3>{fig.to_html(full_html=False, include_plotlyjs=False)}<hr style="margin: 25px 0;">'

def create_deg_tables_html(deg_df, cluster_markers, p_val_thresh, log2fc_thresh, mean_counts_thresh):
    reshaped = _reshape_deg_df(deg_df)
    
    # Generate Plots
    plots = create_deg_violin_plots(reshaped, p_val_thresh, log2fc_thresh, mean_counts_thresh)
    bar = plot_deg_counts_barchart(cluster_markers)
    
    # Generate Dropdown and Tables
    dropdown = ['<label for="cluster_select"><b>Select a Cluster to view its DEGs:</b></label>', 
                '<select id="cluster_select" onchange="showTable(this.value)">']
    tables = []
    
    for cl in sorted(cluster_markers.keys(), key=natural_sort_key):
        safe_cl = re.sub(r'\s+','_',str(cl))
        dropdown.append(f'<option value="deg_table_{safe_cl}">{cl} ({len(cluster_markers[cl])} genes)</option>')
        
        # Filter
        sub = reshaped[(reshaped['Cluster']==cl) & (reshaped['Feature Name'].isin(cluster_markers[cl]))].copy()
        
        # Format
        sub.rename(columns={'Feature Name':'Gene','adj_p_value':'Adjusted p-value','log2FC':'Log2 Fold Change','mean_counts':'Mean Counts'}, inplace=True)
        cols = ['Gene','Adjusted p-value','Log2 Fold Change']
        if 'Mean Counts' in sub.columns and sub['Mean Counts'].sum() > 0: cols.append('Mean Counts')
        
        tables.append(f'<div id="deg_table_{safe_cl}" class="deg-table-container" style="display:none;"><h4>DEGs for {cl}</h4>{sub[cols].to_html(classes="display compact", index=False, table_id=f"deg_table_{safe_cl}_data", formatters={"Adjusted p-value":"{:.2e}".format, "Log2 Fold Change":"{:.2f}".format, "Mean Counts":"{:.2f}".format})}</div>')
        
    dropdown.append('</select>')
    
    return plots + bar + "".join(dropdown) + "".join(tables)

def _prepare_sig_table(df, tid):
    if df.empty: return "<p>No significant results.</p>", False
    d = df.copy()
    
    if 'Cluster' in d.columns:
        cats = sorted(d['Cluster'].unique(), key=natural_sort_key)
        d['Cluster'] = d['Cluster'].astype(pd.api.types.CategoricalDtype(categories=cats, ordered=True))
        
        sort_cols = ['Cluster', 'adj_p_value']
        asc = [True, True]
        
        if 'Weighted_Enrichment' in d.columns:
            sort_cols.append('Weighted_Enrichment')
            asc.append(False)
        else:
            sort_cols.append('Enrichment_ratio')
            asc.append(False)
            
        d.sort_values(by=sort_cols, ascending=asc, inplace=True)

    has_genes = 'Overlapping_genes' in d.columns
    if has_genes:
        d['Full Gene List'] = d['Overlapping_genes']
        d['Overlapping_genes'] = d['Overlapping_genes'].apply(lambda x: (str(x)[:40] + '...') if len(str(x))>40 else x)
        d.insert(0, '', '')

    fmts = {
        'p_value': '{:.2e}'.format, 
        'adj_p_value': '{:.2e}'.format, 
        'Enrichment_ratio': '{:.2f}'.format, 
        'Weighted_Enrichment': '{:.2f}'.format,
        'Weighted_Recall': '{:.2f}'.format
    }
            
    html = d.to_html(
        classes="display compact", 
        index=False, 
        table_id=tid, 
        formatters=fmts, 
        float_format='{:.3f}'.format
    )
    return html, has_genes

def generate_html_report(sample_name, output_path, sig_results_maps, plots_html_maps, deg_table_html, params_maps, selected_tissue_name=None):
    report_data = {}
    for ctx, levels in sig_results_maps.items():
        report_data[ctx] = {}
        for lvl, res in levels.items():
            tid = f"results_table_{ctx}_{lvl}"
            html, genes = _prepare_sig_table(res, tid)
            report_data[ctx][lvl] = {
                'table': html, 'has_genes': genes, 'id': tid,
                'plots': plots_html_maps[ctx][lvl], 'params': params_maps[ctx][lvl]
            }

    t = """
    <!DOCTYPE html>
    <html><head><title>{{ sample_name }} Report</title>
        <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css"/>
        <script type="text/javascript" src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">
        <script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
        <script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>
        <script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
        <style> 
            body{font-family:Arial,sans-serif;margin:0;padding:0;background-color:#f8f9fa} 
            .container{width:95%;margin:20px auto;padding:20px;background-color:#fff;box-shadow:0 4px 8px #0000001a;border-radius:8px} 
            header{border-bottom:3px solid #007bff;padding-bottom:15px;margin-bottom:20px} 
            h1,h2,h3{color:#0056b3} 
            h2{border-bottom:2px solid #e9ecef;padding-bottom:8px;margin-top:40px} 
            
            /* OLD STYLE TABS */
            .tabs{display:flex;border-bottom:1px solid #ccc} 
            .tab-link{padding:10px 20px;cursor:pointer;background:#f1f1f1;border-bottom:none; margin-right: 2px;} 
            .tab-link.active{background:#fff;border:1px solid #ccc;border-bottom:1px solid #fff;position:relative;top:1px; font-weight: bold;} 
            .tab-content{display:none;padding:20px;border:1px solid #ddd;border-top:none} 
            .tab-content.active{display:block} 
            
            .params{background-color:#e9ecef;padding:15px;border-radius:5px;margin-bottom:20px; display:none;} 
            .plot-controls{margin-bottom:15px;} 
            .plot-controls label, .level-selector label{font-weight:bold; margin-right:10px;} 
            .level-selector{padding: 10px; background-color: #f9f9f9; border-bottom: 1px solid #ddd; margin-bottom: 15px;}
            
            td.dt-control { background: url('https://datatables.net/examples/resources/details_open.png') no-repeat center center; cursor: pointer; } 
            tr.shown td.dt-control { background: url('https://datatables.net/examples/resources/details_close.png') no-repeat center center; } 
            .view-section { display: none; }
        </style>
    </head><body>
    <div class="container">
        <header><h1>Enrichment Analysis Report</h1><h2>Sample: <strong>{{ sample_name }}</strong></h2></header>
        
        {% for context, levels in report_data.items() %}
            {% for level, data in levels.items() %}
            <div id="params_{{ context }}_{{ level }}" class="params">
                <strong>Filters:</strong> Adj. p-value &le; {{ data.params.p_val }} | Log2FC &ge; {{ data.params.log2fc }} | Mean Counts &ge; {{ "%.3f"|format(data.params.mean_counts) }}
            </div>
            {% endfor %}
        {% endfor %}

        <div class="tabs"> 
            <div class="tab-link active" onclick="openTab(event, 'degs')">DEG Browser</div> 
            <div class="tab-link" onclick="openTab(event, 'visuals')">Enrichment Visuals</div> 
            <div class="tab-link" onclick="openTab(event, 'results')">Hypergeometric Result</div> 
        </div>
        
        <div id="degs" class="tab-content active">{{ deg_tables|safe }}</div>
        
        <div id="global_controls" style="display:none;">
             <div class="level-selector">
                <label>1. Tissue Scope:</label>
                <select id="tissue_select" onchange="updateView()">
                    <option value="all_tissue" selected>All Tissue</option>
                    {% if selected_tissue_name %}<option value="selected_tissue">Selected Tissue ({{ selected_tissue_name }})</option>{% endif %}
                </select>
                <span style="margin: 0 15px; border-left: 1px solid #ccc;"></span>
                <label>2. Annotation Level:</label>
                <select id="level_select" onchange="updateView()">
                    <option value="level1" selected>Level 1 (Granular/Conventional)</option>
                    <option value="level2">Level 2 (Broad/Weighted)</option>
                </select>
            </div>
        </div>
        
        <div id="visuals" class="tab-content">
            <div id="visuals_controls_placeholder"></div>
            {% for context, levels in report_data.items() %}
                {% for level, data in levels.items() %}
                <div id="visuals_{{ context }}_{{ level }}" class="view-section visuals-section">
                    <div class="plot-controls"> 
                        <label>Heatmap Top N:</label> 
                        <select onchange="showPlot('{{ context }}_{{ level }}', this.value)"> 
                            {% for n_key in data.plots.heatmap.keys() %}
                            <option value="{{ n_key }}" {% if loop.first %}selected{% endif %}>{{ n_key }}</option>
                            {% endfor %} 
                        </select> 
                    </div>
                    {% for n_key, plot_html in data.plots.heatmap.items() %}
                    <div id="heatmap_{{ context }}_{{ level }}_{{ n_key }}" class="plot-container plot-group-{{ context }}-{{ level }}" style="display: {% if loop.first %}block{% else %}none{% endif %};">
                        {{ plot_html|safe }}
                    </div>
                    {% endfor %}
                </div>
                {% endfor %}
            {% endfor %}
        </div>
        
        <div id="results" class="tab-content">
            <div id="results_controls_placeholder"></div>
            {% for context, levels in report_data.items() %}
                {% for level, data in levels.items() %}
                <div id="results_{{ context }}_{{ level }}" class="view-section results-section" data-has-genes="{{ 'true' if data.has_genes else 'false' }}" data-table-id="{{ data.id }}">
                    {{ data.table|safe }}
                </div>
                {% endfor %}
            {% endfor %}
        </div>
    </div>
    
    <script> 
    function openTab(e, tabName){
        let i, tabcontent, tablinks;
        tabcontent = document.getElementsByClassName("tab-content");
        for(i=0; i<tabcontent.length; i++) tabcontent[i].style.display = "none";
        tablinks = document.getElementsByClassName("tab-link");
        for(i=0; i<tablinks.length; i++) tablinks[i].className = tablinks[i].className.replace(" active", "");
        document.getElementById(tabName).style.display = "block";
        e.currentTarget.className += " active";
        
        var controls = document.getElementById('global_controls');
        if(tabName === 'visuals') {
            controls.style.display = 'block';
            document.getElementById('visuals_controls_placeholder').appendChild(controls);
        } else if (tabName === 'results') {
            controls.style.display = 'block';
            document.getElementById('results_controls_placeholder').appendChild(controls);
        } else {
            controls.style.display = 'none';
        }
        
        // Resize plotly charts when tab becomes visible
        window.dispatchEvent(new Event('resize'));
    } 

    function initResultsTable(tableId, hasGenes) {
        if (!$(tableId).length || $.fn.DataTable.isDataTable(tableId)) return;
        var t = $(tableId).DataTable({
            pageLength: 25, dom: "Bfrtip", buttons: ["copy", "csv", "excel"]
        });
        if (hasGenes) {
            var gIdx = -1;
            $(tableId + ' thead th').each(function(i){ if($(this).text()=='Full Gene List') gIdx=i; });
            if (gIdx > -1) {
                t.column(gIdx).visible(false);
                $(tableId + ' tbody').on('click', 'td.dt-control', function(){
                    var tr = $(this).closest('tr'), row = t.row(tr);
                    if(row.child.isShown()){ row.child.hide(); tr.removeClass('shown'); }
                    else {
                        var g = t.cell(row.index(), gIdx).data();
                        row.child('<div style="background:#f9f9f9;padding:10px;border-left:4px solid #007bff"><b>Genes:</b><br>'+g+'</div>').show();
                        tr.addClass('shown');
                    }
                });
            }
        }
    }

    function updateView() {
        var tissue = document.getElementById('tissue_select').value;
        var level = document.getElementById('level_select').value;
        $('.view-section').hide();
        $('.params').hide();
        $('#visuals_' + tissue + '_' + level).show();
        $('#results_' + tissue + '_' + level).show();
        $('#params_' + tissue + '_' + level).show();

        // Lazy-init DataTable for the now-visible result section
        var wrapper = document.getElementById('results_' + tissue + '_' + level);
        if (wrapper) {
            var tid = wrapper.getAttribute('data-table-id');
            var hasGenes = wrapper.getAttribute('data-has-genes') === 'true';
            if (tid) initResultsTable('#' + tid, hasGenes);
        }
    }

    function showPlot(groupKey, nVal){
        $('.plot-group-' + groupKey.replaceAll('_', '-')).hide();
        $('#heatmap_' + groupKey + '_' + nVal).show();
    }

    function showTable(id){
        $(".deg-table-container").hide();
        if(id){
            var safeId = id.replace(/[\\s'"]/g, '_');
            $("#"+safeId).show();
            var tableId="#"+safeId+"_data";
            if(!$.fn.DataTable.isDataTable(tableId)){
                $(tableId).DataTable({pageLength:10,dom:"Bfrtip",buttons:["copy","csv"]})
            }
        }
    }

    $(document).ready(function(){
        updateView();
        openTab({currentTarget:document.querySelector(".tab-link.active")},"degs");
    });
    </script>
    </body></html>
    """
    
    html_content = Template(t).render(
        sample_name=sample_name, 
        deg_tables=deg_table_html,
        report_data=report_data,
        selected_tissue_name=selected_tissue_name
    )
    
    with open(output_path, "w", encoding='utf-8') as f: 
        f.write(html_content)

if __name__ == "__main__":
    pass