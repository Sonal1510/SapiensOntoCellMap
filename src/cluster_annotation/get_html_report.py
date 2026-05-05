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
import plotly.graph_objects as go
from collections import defaultdict
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
    except (KeyError, TypeError, ValueError):
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
            except (ValueError, np.linalg.LinAlgError):
                pass

        if len(rows_var) > 1:
            try:
                Z = linkage(clustered, method='average', metric='correlation')
                D = dendrogram(Z, no_plot=True, labels=clustered.index)
                final_row_order = D['ivl'] + [r for r in heatmap_df.index if r not in D['ivl']]
            except (ValueError, np.linalg.LinAlgError):
                pass
            
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
    ax_bar.set_xlabel("# Sig. Cell Type Hits")

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

    # Scanpy rank_genes_groups CSV format: names, group, pvals_adj, logfoldchanges
    if 'names' in deg_df.columns and 'group' in deg_df.columns:
        df = deg_df.copy()
        df.rename(columns={
            'names': 'Feature Name',
            'group': 'Cluster',
            'pvals_adj': 'adj_p_value',
            'logfoldchanges': 'log2FC',
        }, inplace=True)
        if 'mean_counts' not in df.columns:
            df['mean_counts'] = 0
        if 'Feature ID' not in df.columns:
            df['Feature ID'] = df['Feature Name']
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
    df.columns.name = None  # Remove "Type" column axis name from pivot
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
    <div style="text-align:center;">
        <h3>Adjusted p-value Distribution (BH-corrected, log scale)</h3>
        <p style="font-size:12px;color:#555;max-width:800px;margin:0 auto 8px;">
        Distribution of BH-adjusted p-values across all DEGs per cluster.
        Red dashed line = significance threshold. Values below the threshold passed filtering.
        </p>
        {p_val_html}
    </div>
    <div style="text-align:center;">
        <h3>Log₂ Fold Change Distribution (symlog scale)</h3>
        <p style="font-size:12px;color:#555;max-width:800px;margin:0 auto 8px;">
        Distribution of log₂ fold changes per cluster (cluster vs. all other clusters).
        Red dashed line = upregulation threshold. Only genes above this line were used for enrichment.
        Symlog scale handles both small and large values without clipping.
        </p>
        {log2fc_html}
    </div>
    """
    if mean_counts_html:
        combined_html += f"""
        <div style="text-align:center;">
            <h3>Mean Gene Counts Distribution (log scale, spatial only)</h3>
            <p style="font-size:12px;color:#555;max-width:800px;margin:0 auto 8px;">
            Average UMI/read count per gene across cells/spots in each cluster.
            Low mean counts indicate lowly-expressed genes that may be filtered out.
            Red dashed line = auto-calibrated mean counts threshold (75th percentile, spatial).
            </p>
            {mean_counts_html}
        </div>
        """
    return combined_html

def plot_deg_counts_barchart(markers):
    if not markers: return ""
    s = pd.Series({k:len(v) for k,v in markers.items()})
    s = s[s > 0]
    if s.empty: return "<h3>DEGs per Cluster</h3><p>No differentially expressed genes found for any cluster.</p>"
    s = s.reindex(sorted(s.index, key=natural_sort_key))
    fig = px.bar(x=s.index, y=s.values,
                 title='DEGs per Cluster used as input to Hypergeometric Enrichment Test<br>'
                       '<sup>Only genes passing adj. p-value, Log₂FC and mean counts filters are counted</sup>',
                 labels={'x': 'Cluster', 'y': 'Number of DEGs'}, color=s.index, text=s.values)
    fig.update_layout(showlegend=False)
    fig.update_traces(marker_color='steelblue', textposition='outside')
    return (
        f'<h3>DEGs per Cluster</h3>'
        f'<p style="font-size:12px;color:#555;">Number of differentially expressed genes per cluster '
        f'that passed all thresholds and were used as input to the hypergeometric enrichment test. '
        f'Higher counts generally yield more robust annotations.</p>'
        f'{fig.to_html(full_html=False, include_plotlyjs=False)}'
        f'<hr style="margin: 25px 0;">'
    )

def create_deg_tables_html(deg_df, cluster_markers, p_val_thresh, log2fc_thresh, mean_counts_thresh):
    reshaped = _reshape_deg_df(deg_df)
    
    # Generate Plots
    plots = create_deg_violin_plots(reshaped, p_val_thresh, log2fc_thresh, mean_counts_thresh)
    bar = plot_deg_counts_barchart(cluster_markers)
    
    # Generate Dropdown and Tables
    dropdown = ['<label for="cluster_select"><b>Select a Cluster to view its DEGs:</b></label>',
                '<select id="cluster_select" onchange="showTable(this.value)">',
                '<option value="">--Select--</option>']
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
        d['Overlapping_genes'] = d['Overlapping_genes'].apply(
            lambda x: (str(x)[:40] + '…') if len(str(x)) > 40 else x)
        d.insert(0, '', '')  # expand arrow column

    # Rename columns to human-readable labels
    col_rename = {
        'Cell_type':            'Cell Type',
        'adj_p_value':          'Adj. P-Value (BH-FDR)',
        'p_value':              'P-Value',
        'Enrichment_ratio':     'Enrichment Ratio (k/n)÷(K/N)',
        'Weighted_Recall':      'Weighted Recall (W_overlap/W_ref)',
        'Weighted_Enrichment':  'Weighted Enrichment',
        'Combined_Score':       'Combined Score',
        'N_Databases':          'N Databases',
        'Overlapping_genes':    'Overlapping Genes ▶',
        'n_overlap':            'k (Overlap)',
        'n_query':              'n (Query)',
        'n_ref':                'K (Reference)',
        'N_background':         'N (Background)',
    }
    d.rename(columns=col_rename, inplace=True)

    # Build formatters using renamed keys
    fmts = {}
    for old, new in col_rename.items():
        if old in ('adj_p_value', 'p_value') and new in d.columns:
            fmts[new] = '{:.2e}'.format
        elif old in ('Enrichment_ratio', 'Weighted_Enrichment',
                     'Weighted_Recall', 'Combined_Score') and new in d.columns:
            fmts[new] = '{:.3f}'.format

    html = d.to_html(
        classes="display compact",
        index=False,
        table_id=tid,
        formatters=fmts,
        float_format='{:.3f}'.format,
    )
    return html, has_genes

def _build_icicle_chart(c_df, cluster_name):
    """Build a Plotly icicle chart for one cluster's hierarchy nodes."""
    rows = c_df.to_dict('records')

    if len(rows) < 2:
        if rows:
            r = rows[0]
            return (f'<p style="margin:10px 0;"><b>Single annotation:</b> '
                    f'{r["Cell_Type"]} ({r["CL_ID"]}) &mdash; '
                    f'Confidence: {r["Confidence"]:.3f}</p>')
        return ""

    # Parse Supporting_Types into sets for subset checking
    for r in rows:
        st = str(r.get('Supporting_Types', ''))
        r['_support_set'] = set(s.strip() for s in st.split(',') if s.strip())

    # Build parent map: for each node, find closest ancestor whose
    # support set is a superset of this node's support set
    parent_map = {}
    for node in rows:
        node_id = node['CL_ID']
        node_depth = node['Depth']
        node_support = node['_support_set']

        best_parent = None
        best_depth = -1
        for candidate in rows:
            if candidate['CL_ID'] == node_id:
                continue
            if candidate['Depth'] >= node_depth:
                continue
            if (node_support and candidate['_support_set']
                    and node_support.issubset(candidate['_support_set'])):
                if candidate['Depth'] > best_depth:
                    best_depth = candidate['Depth']
                    best_parent = candidate['CL_ID']

        parent_map[node_id] = best_parent if best_parent else ""

    # Adjust values for branchvalues="total": parent.value >= sum(children)
    children_of = defaultdict(list)
    for nid, pid in parent_map.items():
        if pid:
            children_of[pid].append(nid)

    node_scores = {r['CL_ID']: max(r.get('Combined_Score', 0), 0.1) for r in rows}

    def get_adjusted_value(nid, visited=None):
        if visited is None:
            visited = set()
        if nid in visited:
            return node_scores.get(nid, 0.1)
        visited.add(nid)
        children = children_of.get(nid, [])
        if not children:
            return node_scores.get(nid, 0.1)
        child_sum = sum(get_adjusted_value(c, visited) for c in children)
        return max(node_scores.get(nid, 0), child_sum + 0.1)

    adjusted = {r['CL_ID']: get_adjusted_value(r['CL_ID']) for r in rows}

    # Build trace data
    ids, labels, parents, values, colors, customdata = [], [], [], [], [], []
    for r in sorted(rows, key=lambda x: x['Depth']):
        ids.append(r['CL_ID'])
        labels.append(r['Cell_Type'])
        parents.append(parent_map.get(r['CL_ID'], ""))
        values.append(adjusted[r['CL_ID']])
        colors.append(r.get('Confidence', 0))
        customdata.append([
            r['Cell_Type'], r['CL_ID'],
            f"{r.get('Confidence', 0):.3f}", r.get('N_Supporting', 0),
            f"{r.get('Combined_Score', 0):.2f}", r.get('Resolution', '')
        ])

    fig = go.Figure(go.Icicle(
        ids=ids, labels=labels, parents=parents, values=values,
        branchvalues="total",
        marker=dict(
            colors=colors,
            colorscale=[
                [0, '#f8d7da'], [0.49, '#f8d7da'],
                [0.5, '#fff3cd'], [0.79, '#fff3cd'],
                [0.8, '#d4edda'], [1.0, '#28a745']
            ],
            cmin=0, cmax=1,
            colorbar=dict(title="Confidence", tickformat=".1f"),
            line=dict(width=1, color='white')
        ),
        customdata=customdata,
        hovertemplate=(
            '<b>%{customdata[0]}</b><br>'
            'CL ID: %{customdata[1]}<br>'
            'Confidence: %{customdata[2]}<br>'
            'N Supporting: %{customdata[3]}<br>'
            'Combined Score: %{customdata[4]}<br>'
            'Resolution: %{customdata[5]}'
            '<extra></extra>'
        ),
        textinfo="label",
        textfont=dict(size=11),
    ))

    height = max(400, len(rows) * 25)
    fig.update_layout(
        title=dict(text=f"Cell Ontology Hierarchy \u2014 {cluster_name}", x=0.5, font=dict(size=14)),
        height=height,
        margin=dict(t=50, l=10, r=10, b=10),
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def _build_hierarchy_html(hierarchical_results):
    """
    Build interactive per-cluster hierarchy view with icicle charts and
    sortable DataTables. Follows the DEG Browser dropdown pattern.

    Args:
        hierarchical_results: dict of {context_name: DataFrame}

    Returns:
        HTML string for the hierarchy section
    """
    if not hierarchical_results:
        return ""

    # Pick the best context (prefer selected_tissue over all_tissue)
    hier_df = None
    for ctx in ['selected_tissue', 'all_tissue']:
        if ctx in hierarchical_results:
            df = hierarchical_results[ctx]
            if isinstance(df, pd.DataFrame) and not df.empty:
                hier_df = df
                break

    if hier_df is None:
        return ""

    clusters = sorted(hier_df['Cluster'].unique(), key=natural_sort_key)
    html_parts = []

    # Dropdown
    html_parts.append(
        '<label for="hierarchy_cluster_select">'
        '<b>Select a Cluster to view its hierarchy:</b></label>')
    html_parts.append(
        '<select id="hierarchy_cluster_select" onchange="showHierarchy(this.value)">')
    html_parts.append('<option value="">--Select a Cluster--</option>')
    for cl in clusters:
        safe_cl = re.sub(r"[\s'\"]+", '_', str(cl))
        n_nodes = len(hier_df[hier_df['Cluster'] == cl])
        html_parts.append(
            f'<option value="hier_{safe_cl}">{cl} ({n_nodes} nodes)</option>')
    html_parts.append('</select>')

    # Legend
    html_parts.append(
        '<div style="margin:10px 0;font-size:12px;">'
        '<span style="background:#d4edda;padding:2px 8px;margin-right:8px;">'
        'High confidence (&ge;0.8)</span>'
        '<span style="background:#fff3cd;padding:2px 8px;margin-right:8px;">'
        'Medium confidence (0.5-0.8)</span>'
        '<span style="background:#f8d7da;padding:2px 8px;">'
        'Low confidence (&lt;0.5)</span>'
        '</div>'
    )

    # Per-cluster containers
    for cl in clusters:
        safe_cl = re.sub(r"[\s'\"]+", '_', str(cl))
        c_df = hier_df[hier_df['Cluster'] == cl].sort_values('Depth')

        html_parts.append(
            f'<div id="hier_{safe_cl}" class="hier-cluster-container"'
            f' style="display:none;">')

        # Icicle chart
        html_parts.append(_build_icicle_chart(c_df, cl))

        # Sortable DataTable
        header_cols = ['Depth', 'CL ID', 'Cell Type', 'N Supporting',
                       'Combined Score', 'Confidence', 'Resolution',
                       'Supporting Types']
        header = '<tr>' + ''.join(f'<th>{c}</th>' for c in header_cols) + '</tr>'

        table_rows = []
        for _, row in c_df.iterrows():
            conf = row.get('Confidence', 0)
            if conf >= 0.8:
                color = '#d4edda'
            elif conf >= 0.5:
                color = '#fff3cd'
            else:
                color = '#f8d7da'

            supporting = str(row.get('Supporting_Types', ''))
            if len(supporting) > 80:
                supporting = supporting[:77] + '...'

            cells = [
                str(row.get('Depth', '')),
                str(row.get('CL_ID', '')),
                f"<b>{row.get('Cell_Type', '')}</b>",
                str(row.get('N_Supporting', '')),
                f"{row.get('Combined_Score', 0):.2f}",
                f"{conf:.3f}",
                str(row.get('Resolution', '')),
                supporting,
            ]
            table_rows.append(
                f'<tr style="background-color:{color};">'
                + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')

        html_parts.append(
            f'<table id="hier_{safe_cl}_table_data" class="display compact"'
            f' style="width:100%;margin-top:15px;">'
            f'<thead>{header}</thead>'
            f'<tbody>{"".join(table_rows)}</tbody>'
            f'</table>')

        html_parts.append('</div>')

    return '\n'.join(html_parts)


def _create_publication_umap(umap_data, sample_name, output_dir):
    """
    Create publication-quality UMAP PNGs matching the BD17 sample summary style.
    Two side-by-side panels:
      Left:  colored by cluster number with centroid labels
      Right: colored by cell type annotation with centroid labels
    Saves standalone PNGs and returns base64-embedded HTML for the report.
    """
    if umap_data is None or umap_data.empty:
        return ""

    logger = logging.getLogger(__name__)

    try:
        import matplotlib.patheffects as pe
        from matplotlib.patches import Patch

        x = umap_data['UMAP-1'].values
        y = umap_data['UMAP-2'].values

        # --- Color palettes ---
        clusters = sorted(umap_data['Cluster'].unique(), key=natural_sort_key)
        cell_types = sorted(umap_data['Top_Cell_Type'].unique(), key=natural_sort_key)

        cluster_palette = sns.color_palette('tab20', n_colors=max(len(clusters), 1))
        cluster_cmap = {c: cluster_palette[i % len(cluster_palette)] for i, c in enumerate(clusters)}

        celltype_palette = sns.color_palette('Set2', n_colors=max(len(cell_types), 1))
        ct_cmap = {ct: celltype_palette[i % len(celltype_palette)] for i, ct in enumerate(cell_types)}
        if 'Unannotated' in ct_cmap:
            ct_cmap['Unannotated'] = (0.8, 0.8, 0.8)

        # --- Compute centroids for labels ---
        cluster_centroids = umap_data.groupby('Cluster')[['UMAP-1', 'UMAP-2']].median()
        ct_centroids = umap_data.groupby('Top_Cell_Type')[['UMAP-1', 'UMAP-2']].median()

        # Text outline for readability
        text_outline = [pe.withStroke(linewidth=2.5, foreground='white')]

        # --- Create figure: two panels side by side ---
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Panel 1: By Cluster
        for cl in clusters:
            mask = umap_data['Cluster'] == cl
            ax1.scatter(x[mask], y[mask], c=[cluster_cmap[cl]], s=1.5,
                        alpha=0.6, edgecolors='none', rasterized=True)
        for cl in clusters:
            if cl in cluster_centroids.index:
                cx, cy = cluster_centroids.loc[cl]
                # Extract cluster number for clean label
                label = cl.replace('Cluster ', '') if cl.startswith('Cluster ') else cl
                ax1.text(cx, cy, label, fontsize=9, fontweight='bold',
                         ha='center', va='center', path_effects=text_outline)

        ax1.set_xlabel('UMAP1', fontsize=11)
        ax1.set_ylabel('UMAP2', fontsize=11, labelpad=0)
        ax1.set_title('Graph-Expression Cluster ID', fontsize=12, fontweight='bold')
        ax1.set_xticks([])
        ax1.set_yticks([])
        for spine in ax1.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color('#333333')

        # Panel 2: By Cell Type
        for ct in cell_types:
            mask = umap_data['Top_Cell_Type'] == ct
            ax2.scatter(x[mask], y[mask], c=[ct_cmap[ct]], s=1.5,
                        alpha=0.6, edgecolors='none', rasterized=True, label=ct)
        for ct in cell_types:
            if ct in ct_centroids.index and ct != 'Unannotated':
                cx, cy = ct_centroids.loc[ct]
                # Truncate long names
                label = ct if len(ct) <= 20 else ct[:18] + '...'
                ax2.text(cx, cy, label, fontsize=8, fontweight='bold',
                         ha='center', va='center', path_effects=text_outline)

        ax2.set_xlabel('UMAP1', fontsize=11)
        ax2.set_ylabel('UMAP2', fontsize=11, labelpad=0)
        ax2.set_title('Cell Type Annotation', fontsize=12, fontweight='bold')
        ax2.set_xticks([])
        ax2.set_yticks([])
        for spine in ax2.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color('#333333')

        fig.suptitle(sample_name, fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()

        # --- Save standalone PNGs ---
        import os
        os.makedirs(output_dir, exist_ok=True)

        # Combined figure
        combined_path = os.path.join(output_dir, f"{sample_name}_umap_celltype.png")
        fig.savefig(combined_path, dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        logger.info(f"UMAP PNG saved: {combined_path}")

        # Also save cell-type-only panel for PDF insertion
        fig_ct, ax_ct = plt.subplots(1, 1, figsize=(7, 6))
        for ct in cell_types:
            mask = umap_data['Top_Cell_Type'] == ct
            ax_ct.scatter(x[mask], y[mask], c=[ct_cmap[ct]], s=1.5,
                          alpha=0.6, edgecolors='none', rasterized=True, label=ct)
        for ct in cell_types:
            if ct in ct_centroids.index and ct != 'Unannotated':
                cx, cy = ct_centroids.loc[ct]
                label = ct if len(ct) <= 20 else ct[:18] + '...'
                ax_ct.text(cx, cy, label, fontsize=9, fontweight='bold',
                           ha='center', va='center', path_effects=text_outline)
        ax_ct.set_xlabel('UMAP1', fontsize=12)
        ax_ct.set_ylabel('UMAP2', fontsize=12, labelpad=0)
        ax_ct.set_xticks([])
        ax_ct.set_yticks([])
        for spine in ax_ct.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.5)
            spine.set_color('#333333')
        plt.tight_layout()

        ct_only_path = os.path.join(output_dir, f"{sample_name}_umap_celltype_only.png")
        fig_ct.savefig(ct_only_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
        plt.close(fig_ct)
        logger.info(f"Cell-type-only UMAP PNG saved: {ct_only_path}")

        # --- Embed combined figure in HTML ---
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        buf.seek(0)
        b64 = base64.b64encode(buf.read()).decode('utf-8')

        return (
            f'<div style="position:relative;">'
            f'<img src="data:image/png;base64,{b64}" loading="lazy" '
            f'style="width:100%; height:auto;">'
            f'<br><a download="{sample_name}_umap_celltype.png" '
            f'href="data:image/png;base64,{b64}" '
            f'style="display:inline-block;margin:8px 0;padding:6px 16px;background:#007bff;'
            f'color:#fff;border-radius:4px;text-decoration:none;font-size:13px;">'
            f'Download UMAP</a>'
            f'<span style="margin-left:10px;font-size:12px;color:#555;">'
            f'300 DPI PNGs also saved to output directory</span></div>'
        )

    except Exception as e:
        logger.warning(f"Publication UMAP failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return ""


def _build_marker_heatmap(top_annotation_df, deg_df):
    """
    Build a marker gene heatmap showing top genes driving each cluster annotation.
    Rows = genes, columns = clusters, values = log2FC.
    """
    if top_annotation_df is None or top_annotation_df.empty:
        return ""
    if deg_df is None or deg_df.empty:
        return ""

    try:
        reshaped = _reshape_deg_df(deg_df)

        # Extract top 5 genes per cluster from the Genes column
        cluster_genes = {}
        for _, row in top_annotation_df.iterrows():
            cluster = row['Cluster']
            genes_str = str(row.get('Genes', ''))
            if genes_str and genes_str != 'nan':
                genes = [g.strip() for g in genes_str.split(',') if g.strip()][:5]
                cluster_genes[cluster] = genes

        if not cluster_genes:
            return ""

        # Pool all unique genes
        all_genes = []
        gene_to_cluster = {}
        for cluster, genes in cluster_genes.items():
            for g in genes:
                if g not in gene_to_cluster:
                    gene_to_cluster[g] = cluster
                all_genes.append(g)
        unique_genes = list(dict.fromkeys(all_genes))  # preserve order, deduplicate

        if not unique_genes:
            return ""

        # Get log2FC values for these genes across all clusters
        reshaped['gene_upper'] = reshaped['Feature Name'].str.upper()
        gene_upper_set = {g.upper() for g in unique_genes}
        gene_data = reshaped[reshaped['gene_upper'].isin(gene_upper_set)].copy()

        if gene_data.empty:
            return ""

        # Pivot: rows=genes, columns=clusters, values=log2FC
        pivot = gene_data.pivot_table(
            index='Feature Name', columns='Cluster',
            values='log2FC', aggfunc='first', fill_value=0
        )
        pivot.columns.name = None

        # Reorder rows by which cluster they support (group genes by annotation)
        sorted_clusters = sorted(cluster_genes.keys(), key=natural_sort_key)
        ordered_genes = []
        for cl in sorted_clusters:
            for g in cluster_genes.get(cl, []):
                matches = [idx for idx in pivot.index if idx.upper() == g.upper()]
                ordered_genes.extend(m for m in matches if m not in ordered_genes)

        # Add any remaining genes not matched
        for g in pivot.index:
            if g not in ordered_genes:
                ordered_genes.append(g)

        pivot = pivot.reindex(index=ordered_genes).fillna(0)

        # Reorder columns naturally
        sorted_cols = sorted(pivot.columns, key=natural_sort_key)
        pivot = pivot.reindex(columns=sorted_cols)

        # Create color annotations for row labels (which cell type each gene supports)
        cluster_to_type = dict(zip(
            top_annotation_df['Cluster'], top_annotation_df['Top_Cell_Type']))
        row_colors_list = []
        unique_types = sorted(set(cluster_to_type.values()), key=natural_sort_key)
        type_palette = sns.color_palette('Set2', n_colors=max(len(unique_types), 1))
        type_color_map = {t: type_palette[i % len(type_palette)] for i, t in enumerate(unique_types)}

        for gene in pivot.index:
            assigned_cluster = gene_to_cluster.get(gene, gene_to_cluster.get(gene.upper(), None))
            if assigned_cluster:
                ct = cluster_to_type.get(assigned_cluster, 'Other')
            else:
                ct = 'Other'
            row_colors_list.append(type_color_map.get(ct, (0.8, 0.8, 0.8)))

        row_colors = pd.Series(row_colors_list, index=pivot.index, name='Cell Type')

        # Create heatmap
        fig_h = max(8, len(pivot.index) * 0.3)
        fig_w = max(10, len(pivot.columns) * 0.5 + 3)

        fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h),
                                  gridspec_kw={'width_ratios': [0.3, 5], 'wspace': 0.02})

        # Color annotation bar
        for i, color in enumerate(row_colors_list):
            axes[0].barh(i, 1, color=color, edgecolor='white', linewidth=0.5)
        axes[0].set_ylim(-0.5, len(pivot.index) - 0.5)
        axes[0].invert_yaxis()
        axes[0].set_xlim(0, 1)
        axes[0].set_xticks([])
        axes[0].set_yticks([])
        axes[0].set_xlabel('Cell Type', fontsize=9)

        # Main heatmap
        vmax = max(abs(pivot.values.min()), abs(pivot.values.max()), 1)
        sns.heatmap(
            pivot, ax=axes[1], cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.5, linecolor='lightgray',
            cbar_kws={'label': 'Log₂ Fold Change', 'shrink': 0.6},
            yticklabels=True, xticklabels=True
        )
        axes[1].set_title('Marker Gene Expression (Top Genes per Annotation)',
                          fontsize=14, weight='bold', pad=15)
        axes[1].set_xlabel('Cluster', fontsize=12)
        axes[1].set_ylabel('')
        plt.setp(axes[1].get_xticklabels(), rotation=45, ha='right', fontsize=10)
        plt.setp(axes[1].get_yticklabels(), fontsize=9)

        # Legend for cell type colors
        from matplotlib.patches import Patch
        legend_patches = [Patch(facecolor=type_color_map[t], label=t) for t in unique_types
                          if t in type_color_map]
        if legend_patches:
            axes[1].legend(handles=legend_patches, title='Annotation',
                          loc='upper left', bbox_to_anchor=(1.15, 1),
                          fontsize=8, title_fontsize=9)

        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        plt.close(fig)
        buf.seek(0)
        b64 = base64.b64encode(buf.read()).decode('utf-8')

        return (
            f'<div style="position:relative;">'
            f'<img src="data:image/png;base64,{b64}" loading="lazy" '
            f'style="width:100%; height:auto;">'
            f'<br><a download="marker_heatmap.png" href="data:image/png;base64,{b64}" '
            f'style="display:inline-block;margin:8px 0;padding:6px 16px;background:#007bff;'
            f'color:#fff;border-radius:4px;text-decoration:none;font-size:13px;">'
            f'Download Marker Heatmap</a></div>'
        )

    except Exception as e:
        logging.getLogger(__name__).warning(f"Marker heatmap failed: {e}")
        return ""


def _build_summary_tab_html(top_annotation_df, umap_html, heatmap_html):
    """
    Assemble the complete Cell Type Summary tab content:
    summary table + UMAP scatter + marker heatmap.
    """
    parts = []

    # Summary DataTable
    if top_annotation_df is not None and not top_annotation_df.empty:
        df = top_annotation_df.copy()

        # --- Proliferative Signature Banner ---
        # Built before column renaming so we can read Proliferative_Flag and
        # Proliferative_Genes from the raw DataFrame.
        prolif_cluster_genes = {}   # {cluster_str: genes_str}
        if 'Proliferative_Flag' in df.columns and 'Cluster' in df.columns:
            for _, pr in df.iterrows():
                if pr.get('Proliferative_Flag') is True or str(pr.get('Proliferative_Flag', '')).lower() == 'true':
                    prolif_cluster_genes[str(pr['Cluster'])] = str(pr.get('Proliferative_Genes', ''))

        if prolif_cluster_genes:
            flagged_list = ', '.join(
                sorted(prolif_cluster_genes.keys(), key=natural_sort_key))
            # Collect all unique proliferative genes seen across flagged clusters
            all_prolif_genes = sorted({
                g.strip() for v in prolif_cluster_genes.values()
                for g in v.split(',') if g.strip()
            })
            gene_examples = ', '.join(all_prolif_genes[:6])
            if len(all_prolif_genes) > 6:
                gene_examples += f' (+{len(all_prolif_genes)-6} more)'
            parts.append(
                '<div style="background:#fff3cd;border:1px solid #ffc107;'
                'border-radius:6px;padding:12px 18px;margin-bottom:14px;'
                'line-height:1.6;">'
                '<strong style="font-size:14px;">&#9888; Proliferative / Potentially Malignant Signature Detected</strong><br>'
                f'Cluster(s) <strong>{flagged_list}</strong> have &ge;2 overlapping genes from the canonical '
                f'cell-cycle gene set (<em>{gene_examples}</em>). '
                'In specimens from <strong>cancer, inflamed, or wound-healing tissue</strong>, '
                'these clusters may represent <strong>cycling or malignant cells</strong> '
                'rather than the annotated resting cell type. '
                'The enrichment test correctly identifies statistical overlap with marker databases, '
                'but cell-cycle genes appear in many cell-type entries — the annotation reflects '
                'proliferative identity, not lineage identity.<br>'
                '<em style="font-size:11px;color:#666;">'
                'Recommended follow-up: Ki-67 (MKI67) IHC, copy-number inference (inferCNV/CopyKAT), '
                'pathologist review, or cell-cycle regression before re-annotation. '
                'References: Tirosh et al., <em>Science</em> 2016 (cell-cycle gene modules); '
                'Whitfield et al., <em>Mol Biol Cell</em> 2002 (MCM complex markers); '
                'clinical Ki-67 scoring guidelines (Dowsett et al., <em>J Clin Oncol</em> 2011).'
                '</em>'
                '</div>'
            )

        # Format columns for display
        fmt_df = df.copy()
        if 'P_Value' in fmt_df.columns:
            fmt_df['P_Value'] = fmt_df['P_Value'].apply(
                lambda x: f'{x:.2e}' if pd.notna(x) else '—')
        if 'Score' in fmt_df.columns:
            fmt_df['Score'] = fmt_df['Score'].apply(
                lambda x: f'{x:.2f}' if pd.notna(x) else '—')
        if 'Confidence' in fmt_df.columns:
            fmt_df['Confidence'] = fmt_df['Confidence'].apply(
                lambda x: f'{x:.3f}' if pd.notna(x) else '—')
        if 'N_Databases' in fmt_df.columns:
            fmt_df['N_Databases'] = fmt_df['N_Databases'].apply(
                lambda x: str(int(x)) if pd.notna(x) else '—')
        if 'Genes' in fmt_df.columns:
            fmt_df['Genes'] = fmt_df['Genes'].apply(
                lambda x: (str(x)[:60] + '...') if len(str(x)) > 60 else str(x))

        # Rename for display (Proliferative_Flag/Genes excluded — shown as badge)
        display_cols = {
            'Cluster': 'Cluster', 'Top_Cell_Type': 'Cell Type',
            'Broad_Type': 'Broad Type', 'Broad_Type_CL_ID': 'Broad Type CL ID',
            'Confidence': 'Confidence',
            'P_Value': 'Adj. P-Value', 'Score': 'Score (Weighted Enrichment)',
            'N_Databases': 'N Databases', 'Source': 'DB Source', 'Genes': 'Top Genes'
        }
        cols_present = [c for c in display_cols if c in fmt_df.columns]
        fmt_df = fmt_df[cols_present].rename(columns=display_cols)

        # Color rows by confidence; inject ⚠ PROLIF badge where flagged
        rows_html = []
        for (_, row), (_, orig_row) in zip(fmt_df.iterrows(), df.iterrows()):
            cluster_key = str(orig_row.get('Cluster', ''))
            conf_str = row.get('Confidence', '—')
            try:
                conf_val = float(conf_str)
                if conf_val >= 0.8:
                    bg = '#d4edda'
                elif conf_val >= 0.5:
                    bg = '#fff3cd'
                else:
                    bg = '#f8d7da'
            except (ValueError, TypeError):
                bg = '#ffffff'

            if cluster_key in prolif_cluster_genes:
                prolif_g = prolif_cluster_genes[cluster_key]
                row = row.copy()
                if 'Cell Type' in row.index:
                    row['Cell Type'] = (
                        f'{row["Cell Type"]} '
                        f'<span title="Proliferative signature genes: {prolif_g}" '
                        f'style="background:#e85d04;color:#fff;border-radius:3px;'
                        f'padding:1px 5px;font-size:10px;cursor:help;white-space:nowrap;">'
                        f'&#9888; PROLIF</span>'
                    )

            cells = ''.join(f'<td>{v}</td>' for v in row.values)
            rows_html.append(f'<tr style="background-color:{bg};">{cells}</tr>')

        header = '<tr>' + ''.join(f'<th>{c}</th>' for c in fmt_df.columns) + '</tr>'
        table_html = (
            f'<table id="summary_annotation_table" class="display compact" '
            f'style="width:100%;">'
            f'<thead>{header}</thead>'
            f'<tbody>{"".join(rows_html)}</tbody>'
            f'</table>'
        )

        parts.append('<h3>Recommended Cell Type per Cluster</h3>')
        parts.append(
            '<p style="font-size:12px;color:#555;margin-bottom:6px;">'
            '<strong>Columns:</strong> '
            '<em>Cell Type</em> — top annotation (Selected Tissue Level&nbsp;2, else All Tissue Level&nbsp;2). '
            'A <span style="background:#e85d04;color:#fff;border-radius:3px;padding:1px 4px;font-size:10px;">&#9888;&nbsp;PROLIF</span> '
            'badge flags clusters where &ge;2 overlapping genes are canonical cell-cycle markers '
            '(see alert above for clinical interpretation). '
            '<em>Broad Type</em> — shallowest meaningful CL ontology ancestor of the annotated cell type (ontology-derived, no hardcoding). '
            '<em>Broad Type CL ID</em> — Cell Ontology identifier for the broad type. '
            '<em>Confidence</em> — fraction of known subtypes supporting this node. '
            '<em>Score</em> — Weighted Enrichment Ratio&nbsp;= (W_overlap/n)&nbsp;÷&nbsp;(W_ref/N). '
            '<em>N&nbsp;Databases</em> — number of independent databases corroborating the overlapping genes. '
            '<em>Top Genes</em> — overlapping marker genes (first 60 chars). '
            'Row color: '
            '<span style="background:#d4edda;padding:1px 6px;border-radius:3px;">green</span> Confidence&nbsp;≥&nbsp;0.8 &nbsp;'
            '<span style="background:#fff3cd;padding:1px 6px;border-radius:3px;">yellow</span> 0.5–0.8 &nbsp;'
            '<span style="background:#f8d7da;padding:1px 6px;border-radius:3px;">red</span> &lt;&nbsp;0.5'
            '</p>'
        )
        parts.append(table_html)

    # UMAP scatter
    if umap_html:
        parts.append('<h3 style="margin-top:30px;">UMAP Visualization</h3>')
        parts.append('<p style="font-size:13px;color:#555;">Left: clusters colored by '
                     'GE cluster ID. Right: clusters colored by recommended cell type '
                     'annotation. 300 DPI PNGs saved to output directory for publication use.</p>')
        parts.append(umap_html)

    # Marker heatmap
    if heatmap_html:
        parts.append('<h3 style="margin-top:30px;">Marker Gene Heatmap</h3>')
        parts.append('<p style="font-size:13px;color:#555;">Log₂ fold change of top '
                     'marker genes driving each cluster annotation. Genes grouped by '
                     'recommended cell type (color bar on left).</p>')
        parts.append(heatmap_html)

    return '\n'.join(parts)


def _build_deconvolution_tab_html(deconv_df):
    """
    Build the Deconvolution tab HTML from a cluster-level proportions DataFrame.

    Shows:
      1. Stacked horizontal bar chart — one bar per cluster, stacked by cell type
         proportion (top 15 cell types by mean proportion across clusters).
      2. Downloadable DataTable of per-cluster proportions.

    Args:
        deconv_df: pd.DataFrame (clusters x cell_types), values = proportions.
                   Index = cluster names (e.g. "Cluster 1").

    Returns:
        HTML string
    """
    if deconv_df is None or deconv_df.empty:
        return ""

    parts = []

    # Select top cell types by mean proportion (skip near-zero columns)
    mean_props = deconv_df.mean(axis=0).sort_values(ascending=False)
    top_ct = mean_props[mean_props > 0.005].head(15).index.tolist()
    if not top_ct:
        top_ct = mean_props.head(10).index.tolist()

    # Sort clusters naturally for consistent ordering
    sorted_clusters = sorted(deconv_df.index.tolist(), key=natural_sort_key)
    plot_df = deconv_df.reindex(sorted_clusters)[top_ct]

    # Compute "Other" = remaining proportion not captured by top_ct
    # This ensures every bar reaches exactly 1.0 (100%)
    other_vals = (1.0 - plot_df.sum(axis=1)).clip(lower=0.0)

    # --- 1. Stacked horizontal bar chart (clusters × cell types) ---
    colors = (
        ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
         '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
         '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5']
    )

    fig = go.Figure()
    for j, ct in enumerate(top_ct):
        ct_short = ct[:50] + ("…" if len(ct) > 50 else "")
        fig.add_trace(go.Bar(
            name=ct_short,
            y=sorted_clusters,
            x=plot_df[ct].values,
            orientation='h',
            marker_color=colors[j % len(colors)],
            hovertemplate=(
                f"{ct_short}<br>%{{y}}<br>Proportion: %{{x:.3f}}<extra></extra>"
            ),
        ))

    # Add "Other" segment so each bar always reaches 1.0
    if other_vals.max() > 0.001:
        fig.add_trace(go.Bar(
            name="Other cell types",
            y=sorted_clusters,
            x=other_vals.values,
            orientation='h',
            marker_color='#d3d3d3',
            hovertemplate="Other cell types<br>%{y}<br>Proportion: %{x:.3f}<extra></extra>",
        ))

    fig.update_layout(
        barmode='stack',
        title=(
            f"Annotation-Derived Cell Type Composition per Cluster"
            f" (Top {len(top_ct)} shown + Other)"
        ),
        xaxis_title="Proportion (bars sum to 1.0)",
        yaxis_title="",
        height=max(350, len(sorted_clusters) * 32 + 120),
        margin=dict(l=110, r=20, t=60, b=50),
        legend=dict(
            orientation='v',
            x=1.01, y=1.0,
            xanchor='left',
            font=dict(size=10),
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=True, gridcolor='#eeeeee',
            showline=True, linecolor='black',
            range=[0, 1],
        ),
        yaxis=dict(
            showgrid=False,
            showline=True, linecolor='black',
            autorange='reversed',
        ),
    )
    stacked_html = fig.to_html(
        full_html=False, include_plotlyjs=False,
        div_id="deconv_stack_bar", default_height="auto",
    )
    parts.append(f'<div style="margin-bottom:30px;">{stacked_html}</div>')

    # --- 2. Per-cluster proportions table (all cell types, not just top) ---
    display_df = deconv_df.round(4).copy()
    # Drop columns that are all-zero across all clusters
    display_df = display_df.loc[:, (display_df > 0).any(axis=0)]
    display_df.index.name = 'Cluster'
    table_html = display_df.to_html(
        table_id="deconv_table",
        classes="display compact",
        border=0,
    )
    parts.append(
        f'<p style="font-size:12px;color:#666;">'
        f'{len(deconv_df)} clusters × {len(display_df.columns)} cell types '
        f'(chart shows top {len(top_ct)} + Other; table shows all non-zero types)</p>'
    )
    parts.append(table_html)

    return "\n".join(parts)


def generate_html_report(sample_name, output_path, sig_results_maps, plots_html_maps, deg_table_html, params_maps, selected_tissue_name=None, hierarchical_results=None, top_annotation_df=None, umap_data=None, deg_df=None, deconv_df=None):
    report_data = {}
    for ctx, levels in sig_results_maps.items():
        report_data[ctx] = {}
        for lvl, res in levels.items():
            if lvl == 'hierarchical':
                continue  # handled separately
            tid = f"results_table_{ctx}_{lvl}"
            html, genes = _prepare_sig_table(res, tid)
            report_data[ctx][lvl] = {
                'table': html, 'has_genes': genes, 'id': tid,
                'plots': plots_html_maps[ctx][lvl], 'params': params_maps[ctx][lvl]
            }

    # Build hierarchy HTML section
    hierarchy_html = _build_hierarchy_html(hierarchical_results or {})

    # Build deconvolution tab
    deconv_html = _build_deconvolution_tab_html(deconv_df)

    # Build summary tab content
    import os
    output_dir = os.path.dirname(output_path)
    summary_html = ""
    if top_annotation_df is not None and not top_annotation_df.empty:
        umap_html = _create_publication_umap(umap_data, sample_name, output_dir) if umap_data is not None else ""
        heatmap_html = _build_marker_heatmap(top_annotation_df, deg_df) if deg_df is not None else ""
        summary_html = _build_summary_tab_html(top_annotation_df, umap_html, heatmap_html)

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
        <script src="https://cdn.plot.ly/plotly-3.1.1.min.js"></script>
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
            .tab-explanation { background-color:#f0f7ff; border-left:4px solid #007bff; padding:12px 16px; margin-bottom:20px; border-radius:0 4px 4px 0; font-size:13px; line-height:1.5; }
            .tab-explanation summary { cursor:pointer; font-size:14px; color:#0056b3; }
            .hier-cluster-container { display:none; }
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
            {% if summary_html %}<div class="tab-link active" onclick="openTab(event, 'summary')">Cell Type Summary</div>{% endif %}
            <div class="tab-link {% if not summary_html %}active{% endif %}" onclick="openTab(event, 'degs')">DEG Browser</div>
            <div class="tab-link" onclick="openTab(event, 'visuals')">Enrichment Visuals</div>
            <div class="tab-link" onclick="openTab(event, 'results')">Hypergeometric Result</div>
            {% if hierarchy_html %}<div class="tab-link" onclick="openTab(event, 'hierarchy')">Hierarchy</div>{% endif %}
            {% if deconv_html %}<div class="tab-link" onclick="openTab(event, 'deconvolution')">Composition</div>{% endif %}
        </div>

        {% if summary_html %}
        <div id="summary" class="tab-content active">
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p><strong>Top annotation per cluster</strong> from Level 2 weighted enrichment
                (Selected Tissue context takes priority over All Tissue). Columns explained in the
                table caption below. Row colors: <span style="background:#d4edda;padding:0 4px;">green</span>
                = high confidence (≥0.8), <span style="background:#fff3cd;padding:0 4px;">yellow</span>
                = moderate (0.5–0.8), <span style="background:#f8d7da;padding:0 4px;">red</span>
                = low (&lt;0.5). UMAP panels (when available) show cluster layout colored by
                cluster ID (left) and cell type (right). Marker heatmap shows log₂FC of top
                genes driving each annotation, grouped by assigned cell type.</p>
            </details>
            {{ summary_html|safe }}
        </div>
        {% endif %}

        <div id="degs" class="tab-content {% if not summary_html %}active{% endif %}">
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p><strong>Differentially expressed genes (DEGs)</strong> for each cluster after
                applying the analysis thresholds (adj. p-value, Log₂FC, mean counts).
                Select a cluster from the dropdown to view its sorted gene table.
                <br>
                <strong>Violin plots:</strong>
                <em>Adj. p-value</em> — BH-corrected p-value on log scale; red dashed line = threshold.
                <em>Log₂FC</em> — symlog-scaled fold change; red dashed line = upregulation threshold.
                <em>Mean Counts</em> — average expression per gene across cells in the cluster (spatial only).
                <br>
                <strong>Bar chart:</strong> number of DEGs per cluster submitted to the enrichment test.</p>
            </details>
            {{ deg_tables|safe }}
        </div>

        <div id="global_controls" style="display:none;">
             <div class="level-selector">
                <label>1. Tissue Scope:</label>
                <select id="tissue_select" onchange="updateView()">
                    <option value="all_tissue" selected>All Tissue (no tissue filter)</option>
                    {% if selected_tissue_name %}<option value="selected_tissue">Selected Tissue: {{ selected_tissue_name }}</option>{% endif %}
                </select>
                <span style="margin: 0 15px; border-left: 1px solid #ccc;"></span>
                <label>2. Annotation Level:</label>
                <select id="level_select" onchange="updateView()">
                    <option value="level1" selected>Level 1 – Granular (per-database, unweighted)</option>
                    <option value="level2">Level 2 – Broad (CL-normalized, evidence-weighted)</option>
                </select>
            </div>
        </div>

        <div id="visuals" class="tab-content">
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p><strong>Enrichment heatmap</strong>: rows = clusters, columns = top-N cell types,
                color = Weighted Enrichment score (higher is more enriched). The bar chart on the
                left shows number of significant cell type hits per cluster.
                Use the <em>Heatmap Top N</em> selector to show 1, 3, 5 or 10 cell types per cluster.
                Use the <em>Tissue Scope</em> and <em>Annotation Level</em> selectors below to switch
                between All-Tissue / Selected-Tissue and granular (Level 1, per-DB) / broad
                (Level 2, CL-normalized) contexts.</p>
            </details>
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
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p><strong>Full hypergeometric enrichment results</strong> — one row per
                (cluster, cell type) pair that passed the significance threshold.
                Table is sortable and exportable (Copy / CSV / Excel buttons).
                Click the <strong>▶</strong> expand arrow on any row to see the full overlapping gene list.
                Use the <em>Tissue Scope</em> and <em>Annotation Level</em> selectors above to switch views.
                <br><strong>Column guide:</strong>
                <em>Adj. P-Value (BH-FDR)</em> — Benjamini–Hochberg corrected p-value across all
                (cluster × cell type) tests; ≤ 0.05 = significant.
                <em>Enrichment Ratio (k/n)÷(K/N)</em> — observed overlap fraction divided by
                expected fraction; &gt;1 = enriched.
                <em>Weighted Recall (W_overlap/W_ref)</em> — evidence-weighted fraction of the
                reference cell type's markers that were found in this cluster's DEGs.
                <em>Weighted Enrichment</em> — weighted analog of Enrichment Ratio
                (W_overlap/n)÷(W_ref/N); primary ranking metric for Level 2.
                <em>Combined Score</em> — Weighted Enrichment × −log₁₀(adj. p-value);
                composite ranking (higher = stronger hit).
                <em>N Databases</em> — number of independent databases corroborating the overlap genes.
                <em>Overlapping Genes ▶</em> — first 40 chars; click ▶ for full list.</p>
            </details>
            <div id="results_controls_placeholder"></div>
            {% for context, levels in report_data.items() %}
                {% for level, data in levels.items() %}
                <div id="results_{{ context }}_{{ level }}" class="view-section results-section" data-has-genes="{{ 'true' if data.has_genes else 'false' }}" data-table-id="{{ data.id }}">
                    {{ data.table|safe }}
                </div>
                {% endfor %}
            {% endfor %}
        </div>

        {% if hierarchy_html %}
        <div id="hierarchy" class="tab-content">
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p>Multi-resolution cell type annotation via Cell Ontology traversal. Significant
                enrichment hits are mapped upward through the ontology tree. Select a cluster
                to see its icicle chart and details table. Confidence = fraction of known
                subtypes supporting the annotation. Green &ge;0.8, yellow 0.5-0.8, red &lt;0.5.</p>
            </details>
            {{ hierarchy_html|safe }}
        </div>
        {% endif %}

        {% if deconv_html %}
        <div id="deconvolution" class="tab-content">
            <details class="tab-explanation" open>
                <summary>About this tab</summary>
                <p>Annotation-derived cell type composition scores per cluster. For each cluster,
                the Combined Score from the hypergeometric enrichment test
                (Weighted_Enrichment &times; &minus;log<sub>10</sub>(adj_p_value)) is normalised
                across all significantly enriched cell types so that scores sum to 1.0 per cluster.
                This approach uses the marker database as intended &mdash; via enrichment testing
                &mdash; and works for both scRNA-seq and spatial data without requiring a full
                expression matrix. The stacked bar chart shows the top cell types per cluster;
                bars reach 1.0 with an &ldquo;Other cell types&rdquo; segment for remaining
                proportions. The table below lists all non-zero scores per cluster and is
                downloadable as CSV.</p>
            </details>
            {{ deconv_html|safe }}
            <script>
            $(document).ready(function(){
                if($("#deconv_table").length && !$.fn.DataTable.isDataTable("#deconv_table")){
                    $("#deconv_table").DataTable({
                        pageLength: 25,
                        dom: "Bfrtip",
                        buttons: ["copy", "csv", "excel"],
                        scrollX: true,
                        order: []
                    });
                }
            });
            </script>
        </div>
        {% endif %}
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

        // Lazy-init summary DataTable
        if(tabName === 'summary') {
            var tid = "#summary_annotation_table";
            if($(tid).length && !$.fn.DataTable.isDataTable(tid)) {
                $(tid).DataTable({
                    pageLength: 50, dom: "Bfrtip", buttons: ["copy", "csv"],
                    order: [[0, "asc"]]
                });
            }
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

    function showHierarchy(id){
        $(".hier-cluster-container").hide();
        if(id){
            var safeId = id.replace(/[\s'"]/g, '_');
            $("#"+safeId).show();
            var tableId = "#"+safeId+"_table_data";
            if(!$.fn.DataTable.isDataTable(tableId)){
                $(tableId).DataTable({
                    pageLength:25, dom:"Bfrtip", buttons:["copy","csv"],
                    order:[[0,"asc"],[5,"desc"]]
                });
            }
            window.dispatchEvent(new Event('resize'));
        }
    }

    $(document).ready(function(){
        updateView();
        var defaultTab = document.getElementById('summary') ? 'summary' : 'degs';
        openTab({currentTarget:document.querySelector(".tab-link.active")}, defaultTab);
    });
    </script>
    </body></html>
    """
    
    html_content = Template(t).render(
        sample_name=sample_name,
        deg_tables=deg_table_html,
        report_data=report_data,
        selected_tissue_name=selected_tissue_name,
        hierarchy_html=hierarchy_html,
        summary_html=summary_html,
        deconv_html=deconv_html,
    )
    
    with open(output_path, "w", encoding='utf-8') as f: 
        f.write(html_content)

if __name__ == "__main__":
    pass