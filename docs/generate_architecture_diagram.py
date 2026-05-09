#!/usr/bin/env python3
"""
Generate publication-quality architecture diagram for SapiensOntoCellMap.
Audience-friendly design: plain-English step labels + technical name as subtitle.
Two-panel layout: Phase 1 (Database Construction) | Phase 2 (Cell Type Annotation)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

# ─── Color palette — modern, warm, accessible ───────────────────────────────
C_PHASE1_BG    = '#faf6ff'   # soft lavender
C_PHASE2_BG    = '#f0f9ff'   # soft sky blue
C_STEP         = '#7c4dff'   # deep violet — processing steps
C_STEP_EDGE    = '#4a148c'
C_STEP_TEXT    = '#ffffff'
C_DATA         = '#fff8e1'   # amber tint — data stores
C_DATA_EDGE    = '#f57f17'
C_DATA_TEXT    = '#4e342e'
C_DB           = '#e8f5e9'   # green tint — outputs
C_DB_EDGE      = '#2e7d32'
C_DB_TEXT      = '#1b5e20'
C_INPUT        = '#fff3e0'   # warm orange — user input
C_INPUT_EDGE   = '#e65100'
C_INPUT_TEXT   = '#bf360c'
C_PARSERS_BG   = '#ede7f6'   # light violet — parser group
C_PARSERS_EDGE = '#7c4dff'
C_VIZ          = '#e3f2fd'   # light blue — visualizer
C_VIZ_EDGE     = '#1565c0'
C_VIZ_TEXT     = '#0d47a1'
C_ARROW        = '#546e7a'
C_ARROW_DB     = '#2e7d32'
C_TITLE        = '#1a237e'
C_STEP_NUM     = '#ff6d00'   # orange step numbers

# ─── Figure ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(22, 14))
ax.set_xlim(0, 22)
ax.set_ylim(0, 14)
ax.axis('off')
fig.patch.set_facecolor('#ffffff')

# ─── Helpers ─────────────────────────────────────────────────────────────────

def rbox(x, y, w, h, fc, ec, title, subtitle='', fontsize=9, bold=True,
         text_color='#ffffff', sub_color=None, radius=0.22, zorder=3,
         step_num=None):
    """Rounded box with a title line and optional subtitle."""
    box = FancyBboxPatch((x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=1.6, zorder=zorder)
    ax.add_patch(box)
    sub_color = sub_color or text_color
    if subtitle:
        # title sits above center, subtitle below
        ax.text(x + w/2, y + h/2 + 0.12, title,
                ha='center', va='center', fontsize=fontsize,
                color=text_color, fontweight='bold' if bold else 'normal',
                zorder=zorder+1, fontfamily='DejaVu Sans')
        ax.text(x + w/2, y + h/2 - 0.16, subtitle,
                ha='center', va='center', fontsize=fontsize - 1.2,
                color=sub_color, fontweight='normal', fontstyle='italic',
                zorder=zorder+1, fontfamily='DejaVu Sans')
    else:
        ax.text(x + w/2, y + h/2, title,
                ha='center', va='center', fontsize=fontsize,
                color=text_color, fontweight='bold' if bold else 'normal',
                zorder=zorder+1, fontfamily='DejaVu Sans')
    if step_num is not None:
        ax.text(x + 0.22, y + h - 0.18, str(step_num),
                ha='center', va='top', fontsize=9, fontweight='bold',
                color=C_STEP_NUM, zorder=zorder+2, fontfamily='DejaVu Sans')


def arr(x0, y0, x1, y1, label='', color=C_ARROW, lw=1.5,
        mutation_scale=14, zorder=5, loffset=(0, 0.13)):
    fa = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='->',
        color=color, linewidth=lw, mutation_scale=mutation_scale, zorder=zorder)
    ax.add_patch(fa)
    if label:
        mx = (x0+x1)/2 + loffset[0]
        my = (y0+y1)/2 + loffset[1]
        ax.text(mx, my, label, ha='center', va='bottom', fontsize=7.5,
                color='#607d8b', fontstyle='italic', zorder=zorder+1,
                fontfamily='DejaVu Sans')


def panel(x, y, w, h, fc, ec, title, title_color='#1a237e'):
    bg = FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=0.35",
        facecolor=fc, edgecolor=ec, linewidth=2.0, zorder=1)
    ax.add_patch(bg)
    ax.text(x + 0.28, y + h - 0.3, title,
            ha='left', va='top', fontsize=11.5, fontweight='bold',
            color=title_color, zorder=2, fontfamily='DejaVu Sans')


# ═══════════════════════════════════════════════════════════════════════════════
# TITLE
# ═══════════════════════════════════════════════════════════════════════════════
ax.text(11, 13.72, 'SapiensOntoCellMap', ha='center', va='top',
        fontsize=18, fontweight='bold', color=C_TITLE,
        fontfamily='DejaVu Sans')
ax.text(11, 13.27, 'Ontology-Aware Cell Type Annotation for Single-Cell & Spatial Transcriptomics',
        ha='center', va='top', fontsize=10, color='#546e7a',
        fontfamily='DejaVu Sans')

# thin rule under title
ax.axhline(13.0, xmin=0.02, xmax=0.98, color='#cfd8dc', linewidth=0.8)

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 PANEL  (left)
# ═══════════════════════════════════════════════════════════════════════════════
panel(0.3, 0.28, 10.1, 12.55, C_PHASE1_BG, '#9575cd',
      'Phase 1 — Build the Marker Database', '#4527a0')
ax.text(0.55, 12.55, '(one-time setup — run once, reuse forever)',
        ha='left', va='top', fontsize=8, color='#7e57c2', fontstyle='italic',
        fontfamily='DejaVu Sans')

CX1 = 5.35   # centre-x of Phase 1 column

# Step A: Configuration
rbox(1.1, 11.25, 4.5, 0.72, C_STEP, C_STEP_EDGE,
     'Configure data sources',
     subtitle='config/config.py  ·  14 public databases defined',
     fontsize=8.8, text_color='#ffffff', sub_color='#e0cfff', step_num='A')

# Step B: Download raw data
rbox(1.1, 10.28, 4.5, 0.72, C_STEP, C_STEP_EDGE,
     'Download raw marker data',
     subtitle='bio_database_downloader.py',
     fontsize=8.8, text_color='#ffffff', sub_color='#e0cfff', step_num='B')

arr(CX1, 11.25, CX1, 11.0, label='reads config')
arr(CX1, 10.28, CX1, 10.02, label='fetches files')

# Data/raw box
rbox(2.7, 9.48, 3.4, 0.52, C_DATA, C_DATA_EDGE,
     'data/raw/   (raw database files)',
     fontsize=8, text_color=C_DATA_TEXT, subtitle='', bold=False)

arr(CX1, 10.28, CX1, 10.0)

# Step C: Orchestrate parsing
rbox(1.1, 8.6, 4.5, 0.7, C_STEP, C_STEP_EDGE,
     'Orchestrate database construction',
     subtitle='database_creator.py',
     fontsize=8.8, text_color='#ffffff', sub_color='#e0cfff', step_num='C')

arr(CX1, 9.48, CX1, 9.3, label='reads')
arr(CX1, 8.6, CX1, 8.35, label='calls')

# Parsers group
parsers_bg = FancyBboxPatch((0.55, 5.52), 9.15, 2.72,
    boxstyle="round,pad=0.0,rounding_size=0.2",
    facecolor=C_PARSERS_BG, edgecolor=C_PARSERS_EDGE,
    linewidth=1.2, linestyle='--', zorder=2)
ax.add_patch(parsers_bg)
ax.text(CX1, 8.14, '14 Source-Specific Parsers',
        ha='center', va='top', fontsize=8.5, fontweight='bold',
        color='#4527a0', zorder=3, fontfamily='DejaVu Sans')

parsers = [
    ('PanglaoDB', 'CellMarkerDB'),
    ('HuBMap HRA', 'CellxGene'),
    ('WIMMS Melanocyte', 'Human SCC (Ji 2020)'),
    ('Skin Fibroblast Atlas', 'GenericFileParser'),
]
for row_i, (left, right) in enumerate(parsers):
    py = 7.62 - row_i * 0.54
    rbox(0.72, py, 3.95, 0.44, '#ede7f6', '#9575cd', left,
         fontsize=7.6, text_color='#311b92', bold=False, radius=0.12, zorder=4)
    rbox(4.88, py, 4.65, 0.44, '#ede7f6', '#9575cd', right,
         fontsize=7.6, text_color='#311b92', bold=False, radius=0.12, zorder=4)

# Step D: Normalise names & ontology IDs
rbox(0.72, 4.66, 4.3, 0.72, C_STEP, C_STEP_EDGE,
     'Normalise gene names & ontology IDs',
     subtitle='base_parser.py  ·  CL / UBERON IDs',
     fontsize=8.2, text_color='#ffffff', sub_color='#e0cfff', step_num='D')

rbox(5.3, 4.66, 3.9, 0.72, C_STEP, C_STEP_EDGE,
     'Build Cell Ontology graph',
     subtitle='ontology_utils.py  ·  obonet CL DAG',
     fontsize=8.2, text_color='#ffffff', sub_color='#e0cfff')

arr(2.5, 5.52, 2.5, 5.38, label='uses')
arr(7.0, 5.52, 7.0, 5.38, label='uses')

# Intermediate outputs
rbox(0.55, 3.7, 4.1, 0.7, C_DATA, C_DATA_EDGE,
     'Per-parser processed CSVs',
     subtitle='data/processed_db_dfs/',
     fontsize=7.8, text_color=C_DATA_TEXT, sub_color='#795548', bold=False)
rbox(4.85, 3.7, 4.8, 0.7, C_DATA, C_DATA_EDGE,
     'HGNC alias-resolved gene IDs',
     subtitle='data/recovered_ids_dfs/',
     fontsize=7.8, text_color=C_DATA_TEXT, sub_color='#795548', bold=False)

arr(2.5, 4.66, 2.5, 4.4, label='outputs')
arr(7.05, 4.66, 7.05, 4.4, label='outputs')

# Master DB — centred prominent box
rbox(1.0, 2.6, 8.3, 0.82, C_DB, C_DB_EDGE,
     '442,000 marker–cell type associations  ·  14 sources  ·  1,208 cell types',
     subtitle='master_cell_marker_db.csv   (the reference database)',
     fontsize=8.5, text_color=C_DB_TEXT, sub_color='#33691e', bold=True, radius=0.22, zorder=4)

arr(2.5, 3.7, 2.5, 3.42)
arr(7.05, 3.7, 7.05, 3.42)

# Visualizer
rbox(1.0, 1.55, 8.3, 0.78, C_VIZ, C_VIZ_EDGE,
     'Interactive Database Visualizer',
     subtitle='sapiens_visualizer.html  ·  Dashboard · Ontology Browser · Cell Explorer · Gene Explorer',
     fontsize=8.5, text_color=C_VIZ_TEXT, sub_color='#1565c0', bold=True, radius=0.22, zorder=4)

arr(5.15, 2.6, 5.15, 2.33, label='generates', color=C_ARROW_DB,
    loffset=(0.7, 0.09))

ax.text(5.15, 1.42, '→ Used as reference in Phase 2',
        ha='center', va='top', fontsize=7.8, color='#2e7d32',
        fontstyle='italic', fontfamily='DejaVu Sans')

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 PANEL  (right)
# ═══════════════════════════════════════════════════════════════════════════════
panel(10.6, 0.28, 11.1, 12.55, C_PHASE2_BG, '#29b6f6',
      'Phase 2 — Annotate Your Data', '#01579b')
ax.text(10.85, 12.55, '(run per experiment — minutes to complete)',
        ha='left', va='top', fontsize=8, color='#0288d1', fontstyle='italic',
        fontfamily='DejaVu Sans')

CX2 = 16.15   # centre-x of Phase 2 column

def p2arr(y0, y1, label='', loffset=(0.55, 0.09)):
    arr(CX2, y0, CX2, y1, label=label, loffset=loffset)

# User input
rbox(CX2 - 3.5, 11.25, 7.0, 0.82, C_INPUT, C_INPUT_EDGE,
     'Your differential expression results (DEGs)',
     subtitle='Seurat · Scanpy · SpaceRanger/CellRanger · generic CSV',
     fontsize=9, text_color=C_INPUT_TEXT, sub_color='#bf360c', bold=True, radius=0.22)

p2arr(11.25, 10.9, label='input')

# Step 1: Parse DEGs
rbox(CX2 - 3.5, 10.25, 7.0, 0.62, C_STEP, C_STEP_EDGE,
     'Parse & standardise DEGs',
     subtitle='auto-detects format: Seurat | Scanpy | generic | SpaceRanger',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='1')

p2arr(10.25, 9.92, label='gene list')

# Step 2: Gene alias resolution
rbox(CX2 - 3.5, 9.27, 7.0, 0.62, C_STEP, C_STEP_EDGE,
     'Resolve gene aliases',
     subtitle='HGNC alias map · maps historical / synonym symbols → approved names',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='2')

p2arr(9.27, 8.94, label='standardised genes')

# Step 3: Enrichment testing
rbox(CX2 - 3.5, 8.05, 7.0, 0.86, C_STEP, C_STEP_EDGE,
     'Test for enriched cell types',
     subtitle='L1: Hypergeometric test (per-database)  ·  L2: Weighted Enrichment + BH FDR\n'
              'Combined Score = Weighted Enrichment × −log₁₀(adj p-value)',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='3')

p2arr(8.05, 7.7, label='significant hits')

# Step 4: Hierarchical annotation
rbox(CX2 - 3.5, 6.92, 7.0, 0.72, C_STEP, C_STEP_EDGE,
     'Traverse the Cell Ontology hierarchy',
     subtitle='CL DAG walk · assigns confidence at every level (broad → specific subtype)',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='4')

p2arr(6.92, 6.56, label='ranked annotations')

# Step 5: Select top annotation
rbox(CX2 - 3.5, 5.82, 7.0, 0.68, C_STEP, C_STEP_EDGE,
     'Select best cell type per cluster',
     subtitle='tissue-specific result preferred · proliferation flag · lineage conflict check',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='5')

p2arr(5.82, 5.45, label='annotation scores')

# Step 6: Composition scoring
rbox(CX2 - 3.5, 4.67, 7.0, 0.68, C_STEP, C_STEP_EDGE,
     'Compute cell type composition per cluster',
     subtitle='CL ancestor pruning · normalise Combined_Score to sum = 1.0',
     fontsize=8.5, text_color='#fff', sub_color='#e0cfff', step_num='6')

p2arr(4.67, 4.3, label='results')

# Step 7: Generate report
rbox(CX2 - 3.3, 3.6, 6.6, 0.65, C_STEP, C_STEP_EDGE,
     'Generate interactive HTML report',
     subtitle='get_html_report.py  ·  6 tabs: Summary · DEG Browser · Enrichment · Hierarchy · Composition',
     fontsize=8.3, text_color='#fff', sub_color='#e0cfff', step_num='7')

# Outputs — fan out below the report
OUT_Y = 0.48
out_items = [
    ('Top annotation\nper cluster', '*.csv'),
    ('All enrichment\nresults', '*.csv'),
    ('Ontology\nhierarchy', '*.csv'),
    ('Composition\nscores', '*.csv'),
    ('Interactive\nHTML report', '*.html'),
]
n_out = len(out_items)
out_w = 1.72
out_gap = 0.14
total_out_w = n_out * out_w + (n_out - 1) * out_gap
out_start_x = CX2 - total_out_w / 2

for i, (title_o, sub_o) in enumerate(out_items):
    ox = out_start_x + i * (out_w + out_gap)
    rbox(ox, OUT_Y, out_w, 0.86, C_DB, C_DB_EDGE,
         title_o, subtitle=sub_o,
         fontsize=6.8, text_color=C_DB_TEXT, sub_color='#558b2f',
         bold=True, radius=0.15, zorder=4)
    arr(CX2, 3.6, ox + out_w/2, OUT_Y + 0.86,
        color='#78909c', lw=0.9, mutation_scale=10)

# ═══════════════════════════════════════════════════════════════════════════════
# CROSS-PANEL CONNECTOR: master DB → enrichment test (Step 3)
# ═══════════════════════════════════════════════════════════════════════════════
DB_MID_Y = 2.6 + 0.41
ax.annotate('',
    xy=(CX2 - 3.5, 8.05 + 0.43),
    xytext=(9.3, DB_MID_Y),
    arrowprops=dict(
        arrowstyle='->', color=C_ARROW_DB, lw=2.0,
        connectionstyle='arc3,rad=-0.22'
    ), zorder=7)
ax.text(10.4, 6.5, 'marker\ndatabase',
        ha='center', va='center', fontsize=8, color=C_ARROW_DB,
        fontstyle='italic', fontweight='bold',
        fontfamily='DejaVu Sans', zorder=8)

# ═══════════════════════════════════════════════════════════════════════════════
# LEGEND (bottom centre, between panels)
# ═══════════════════════════════════════════════════════════════════════════════
LX, LY = 4.3, 1.15
legend_items = [
    (C_STEP,  C_STEP_EDGE,  '#fff',       'Processing step'),
    (C_INPUT, C_INPUT_EDGE, C_INPUT_TEXT, 'User input'),
    (C_DATA,  C_DATA_EDGE,  C_DATA_TEXT,  'Intermediate data'),
    (C_DB,    C_DB_EDGE,    C_DB_TEXT,    'Output / result'),
    (C_VIZ,   C_VIZ_EDGE,   C_VIZ_TEXT,   'Interactive tool'),
]
ax.text(LX, LY + 0.1, 'Legend:',
        ha='left', va='bottom', fontsize=8, fontweight='bold',
        color='#37474f', fontfamily='DejaVu Sans')
for i, (fc, ec, tc, lbl) in enumerate(legend_items):
    lx_i = LX + i * 1.9
    rbox(lx_i, LY - 0.35, 1.7, 0.32, fc, ec, lbl,
         fontsize=7, text_color=tc, bold=False, radius=0.1, zorder=5)

# ═══════════════════════════════════════════════════════════════════════════════
# SAVE
# ═══════════════════════════════════════════════════════════════════════════════
out_dir = '/Users/srashmi/Desktop/Personal/UCSF_projects/SapiensOntoCellMap/docs'
fig.savefig(f'{out_dir}/architecture_diagram.png', dpi=200,
            bbox_inches='tight', facecolor='white')
fig.savefig(f'{out_dir}/architecture_diagram.pdf',
            bbox_inches='tight', facecolor='white')
print(f'Saved: {out_dir}/architecture_diagram.png')
print(f'Saved: {out_dir}/architecture_diagram.pdf')
plt.close(fig)
