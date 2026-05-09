#!/usr/bin/env python3
"""
Generate publication-quality architecture diagram for SapiensOntoCellMap.
Nature Methods figure quality: clean, minimal, two-panel layout.
Phase 1 (Database Construction) | Phase 2 (Cell Type Annotation)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# ─── Color palette ─────────────────────────────────────────────────────────────
C_PHASE1_BG    = '#faf6ff'   # soft lavender
C_PHASE2_BG    = '#f0f9ff'   # soft sky blue
C_STEP         = '#7c4dff'   # deep violet — processing steps
C_STEP_EDGE    = '#4a148c'
C_STEP_TEXT    = '#ffffff'
C_DATA         = '#fff8e1'   # amber tint — data stores
C_DATA_EDGE    = '#f57f17'
C_DATA_TEXT    = '#4e342e'
C_DB           = '#e8f5e9'   # green tint — outputs / database
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

# ─── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(20, 17))
ax.set_xlim(0, 20)
ax.set_ylim(0, 17)
ax.axis('off')
fig.patch.set_facecolor('#ffffff')

# ─── Helpers ───────────────────────────────────────────────────────────────────

def rbox(x, y, w, h, fc, ec, title, subtitle='', fontsize=10, bold=True,
         text_color='#ffffff', sub_color=None, radius=0.22, zorder=3):
    """Rounded box with a title and optional subtitle."""
    box = FancyBboxPatch((x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=1.6, zorder=zorder)
    ax.add_patch(box)
    sub_color = sub_color or text_color
    if subtitle:
        ax.text(x + w / 2, y + h / 2 + 0.13, title,
                ha='center', va='center', fontsize=fontsize,
                color=text_color, fontweight='bold' if bold else 'normal',
                zorder=zorder + 1, fontfamily='DejaVu Sans')
        ax.text(x + w / 2, y + h / 2 - 0.16, subtitle,
                ha='center', va='center', fontsize=fontsize - 1.5,
                color=sub_color, fontweight='normal',
                zorder=zorder + 1, fontfamily='DejaVu Sans')
    else:
        ax.text(x + w / 2, y + h / 2, title,
                ha='center', va='center', fontsize=fontsize,
                color=text_color, fontweight='bold' if bold else 'normal',
                zorder=zorder + 1, fontfamily='DejaVu Sans')


def arr(x0, y0, x1, y1, color=C_ARROW, lw=1.4, mutation_scale=13, zorder=5):
    """Clean arrow — no inline label."""
    fa = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='->',
        color=color, linewidth=lw, mutation_scale=mutation_scale, zorder=zorder)
    ax.add_patch(fa)


def panel(x, y, w, h, fc, ec, title, title_color='#1a237e'):
    """Panel background with a clearly demarcated title band at top."""
    # Main background
    bg = FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=0.35",
        facecolor=fc, edgecolor=ec, linewidth=2.0, zorder=1)
    ax.add_patch(bg)

    # Title band: a slightly darker strip at the top of the panel
    title_band_h = 0.70
    title_band = FancyBboxPatch((x, y + h - title_band_h), w, title_band_h,
        boxstyle="round,pad=0.0,rounding_size=0.30",
        facecolor=ec, edgecolor=ec, linewidth=0, alpha=0.18, zorder=2)
    ax.add_patch(title_band)

    # Separator line below title band
    ax.plot([x + 0.1, x + w - 0.1], [y + h - title_band_h, y + h - title_band_h],
            color=ec, linewidth=1.2, alpha=0.6, zorder=3)

    # Title text centred in the band
    ax.text(x + w / 2, y + h - title_band_h / 2, title,
            ha='center', va='center', fontsize=12, fontweight='bold',
            color=title_color, zorder=4, fontfamily='DejaVu Sans')


# ═══════════════════════════════════════════════════════════════════════════════
# TITLE  (top of figure)
# ═══════════════════════════════════════════════════════════════════════════════
ax.text(10.0, 16.60, 'SapiensOntoCellMap', ha='center', va='top',
        fontsize=17, fontweight='bold', color=C_TITLE,
        fontfamily='DejaVu Sans')
ax.text(10.0, 16.18, 'Ontology-Aware Cell Type Annotation for Single-Cell & Spatial Transcriptomics',
        ha='center', va='top', fontsize=9.5, color='#546e7a',
        fontfamily='DejaVu Sans')

ax.axhline(15.90, xmin=0.02, xmax=0.98, color='#cfd8dc', linewidth=0.8)

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 PANEL  (left)   y: 0.25 → 15.65  (h=15.40, title band top at 14.95)
# Content boxes start at 13.80 — well below the title band bottom at 14.95
# ═══════════════════════════════════════════════════════════════════════════════
PANEL_Y  = 0.25
PANEL_H  = 15.40
panel(0.28, PANEL_Y, 9.3, PANEL_H, C_PHASE1_BG, '#9575cd',
      'Phase 1 — Database Construction', '#4527a0')

CX1 = 4.93   # centre-x of Phase 1 column

# ── config.py  (top box — starts at y=13.80, title band bottom at 14.95 → 1.15 gap) ──
rbox(1.0, 13.60, 7.8, 0.72, C_STEP, C_STEP_EDGE,
     'config.py',
     subtitle='14 source databases',
     fontsize=10, text_color='#ffffff', sub_color='#e0cfff')

arr(CX1, 13.60, CX1, 13.33)

# ── Database Downloader ──
rbox(1.0, 12.66, 7.8, 0.65, C_STEP, C_STEP_EDGE,
     'Database Downloader',
     fontsize=10, text_color='#ffffff')

arr(CX1, 12.66, CX1, 12.38)

# ── data/raw/ ──
rbox(2.4, 11.86, 5.0, 0.50, C_DATA, C_DATA_EDGE,
     'data/raw/',
     fontsize=9, text_color=C_DATA_TEXT, bold=False)

arr(CX1, 11.86, CX1, 11.60)

# ── Database Creator ──
rbox(1.0, 10.95, 7.8, 0.63, C_STEP, C_STEP_EDGE,
     'Database Creator',
     fontsize=10, text_color='#ffffff')

arr(CX1, 10.95, CX1, 10.70)

# ── Parser group box ──
parsers_bg = FancyBboxPatch((0.45, 8.00), 8.68, 2.65,
    boxstyle="round,pad=0.0,rounding_size=0.2",
    facecolor=C_PARSERS_BG, edgecolor=C_PARSERS_EDGE,
    linewidth=1.2, linestyle='--', zorder=2)
ax.add_patch(parsers_bg)

# Group label
ax.text(CX1, 10.52, 'Source Parsers (×8)',
        ha='center', va='top', fontsize=9, fontweight='bold',
        color='#4527a0', zorder=3, fontfamily='DejaVu Sans')

# Parser names as a 2-column text grid
parser_names_left  = ['PanglaoDB', 'HuBMap HRA', 'WIMMS Melanocyte', 'Skin Fibroblast Atlas']
parser_names_right = ['CellMarkerDB', 'CellxGene', 'Human SCC (Ji 2020)', 'GenericFileParser']
for row_i, (left, right) in enumerate(zip(parser_names_left, parser_names_right)):
    py = 10.06 - row_i * 0.52
    ax.text(2.55, py, left,  ha='center', va='center', fontsize=7.5,
            color='#311b92', zorder=4, fontfamily='DejaVu Sans')
    ax.text(7.30, py, right, ha='center', va='center', fontsize=7.5,
            color='#311b92', zorder=4, fontfamily='DejaVu Sans')

# ── Base Parser + Ontology Utils (side by side) ──
rbox(0.45, 7.10, 4.1, 0.72, C_STEP, C_STEP_EDGE,
     'Base Parser',
     subtitle='CL / UBERON normalisation',
     fontsize=9, text_color='#ffffff', sub_color='#e0cfff')

rbox(4.75, 7.10, 3.9, 0.72, C_STEP, C_STEP_EDGE,
     'Ontology Utils',
     subtitle='CL DAG (obonet)',
     fontsize=9, text_color='#ffffff', sub_color='#e0cfff')

arr(2.3, 8.00, 2.3, 7.82)
arr(6.55, 8.00, 6.55, 7.82)

# ── Intermediate data outputs ──
rbox(0.45, 6.20, 3.85, 0.68, C_DATA, C_DATA_EDGE,
     'Per-parser CSVs',
     subtitle='data/processed_db_dfs/',
     fontsize=8.5, text_color=C_DATA_TEXT, sub_color='#795548', bold=False)

rbox(4.55, 6.20, 4.1, 0.68, C_DATA, C_DATA_EDGE,
     'HGNC alias map',
     subtitle='data/recovered_ids_dfs/',
     fontsize=8.5, text_color=C_DATA_TEXT, sub_color='#795548', bold=False)

arr(2.3, 7.10, 2.3, 6.88)
arr(6.55, 7.10, 6.55, 6.88)

# ── Master DB ── (prominent anchor)
rbox(0.55, 5.16, 8.28, 0.82, C_DB, C_DB_EDGE,
     'master_cell_marker_db.csv',
     subtitle='442K associations  ·  14 sources  ·  1,208 cell types',
     fontsize=10, text_color=C_DB_TEXT, sub_color='#33691e', bold=True, radius=0.22, zorder=4)

arr(2.3, 6.20, 2.3, 5.98)
arr(6.55, 6.20, 6.55, 5.98)

# ── Dynamic DB Update — selling-point callout ─────────────────────────────────
# master_db bottom = 5.16;  place this box y=3.94–5.04 (h=1.10)
# Interactive Visualizer shifted down to y=2.85
DYN_X = 0.55
DYN_Y = 3.76
DYN_W = 8.28
DYN_H = 1.28

# Outer frame — slightly more opaque fill to make it pop
dyn_frame = FancyBboxPatch((DYN_X, DYN_Y), DYN_W, DYN_H,
    boxstyle="round,pad=0.0,rounding_size=0.22",
    facecolor='#fff8f0', edgecolor=C_INPUT_EDGE, linewidth=2.2, zorder=4)
ax.add_patch(dyn_frame)

# Accent bar on left edge
ax.add_patch(FancyBboxPatch((DYN_X, DYN_Y), 0.18, DYN_H,
    boxstyle="round,pad=0.0,rounding_size=0.10",
    facecolor=C_INPUT_EDGE, edgecolor='none', linewidth=0, zorder=5))

# Title
ax.text(DYN_X + 0.36, DYN_Y + DYN_H - 0.22,
        '★  Flexible Database Updates',
        ha='left', va='center', fontsize=10, fontweight='bold',
        color=C_INPUT_EDGE, fontfamily='DejaVu Sans', zorder=6)

# Separator
ax.plot([DYN_X + 0.28, DYN_X + DYN_W - 0.14],
        [DYN_Y + DYN_H - 0.36, DYN_Y + DYN_H - 0.36],
        color=C_INPUT_EDGE, linewidth=0.9, alpha=0.45, zorder=6)

# Bullet points — 3 key selling points
bullets = [
    '+  Add any new source:  drop in a CSV + one config entry',
    '♻  Re-run DatabaseCreate  →  updated master DB in minutes',
    '■  Existing annotations unaffected; new cell types auto-integrated',
]
for bi, btxt in enumerate(bullets):
    ax.text(DYN_X + 0.40, DYN_Y + DYN_H - 0.62 - bi * 0.26,
            btxt, ha='left', va='center', fontsize=8.2,
            color='#6d3a00', fontfamily='DejaVu Sans', zorder=6)

# master_db bottom → Dynamic Update top
arr(4.69, 5.16, 4.69, DYN_Y + DYN_H, color=C_INPUT_EDGE, lw=1.3, mutation_scale=11)
# Dynamic Update bottom → Interactive Visualizer top
arr(4.69, DYN_Y, 4.69, 2.68 + 0.84, color=C_ARROW_DB)

# ── Interactive Visualizer ──
rbox(0.55, 2.68, 8.28, 0.84, C_VIZ, C_VIZ_EDGE,
     'Interactive Visualizer',
     subtitle='sapiens_visualizer.html  ·  4 tabs',
     fontsize=10, text_color=C_VIZ_TEXT, sub_color='#1565c0', bold=True, radius=0.22, zorder=4)

# ─── Legend — framed box at bottom of Phase 1 ─────────────────────────────────
LEG_X, LEG_Y = 0.40, 0.35
LEG_W, LEG_H = 9.10, 1.40

# Legend frame
legend_frame = FancyBboxPatch((LEG_X, LEG_Y), LEG_W, LEG_H,
    boxstyle="round,pad=0.0,rounding_size=0.18",
    facecolor='#f5f0ff', edgecolor='#9575cd', linewidth=1.4, zorder=2)
ax.add_patch(legend_frame)

# Legend title bar
legend_title_band = FancyBboxPatch((LEG_X, LEG_Y + LEG_H - 0.34), LEG_W, 0.34,
    boxstyle="round,pad=0.0,rounding_size=0.15",
    facecolor='#9575cd', edgecolor='#9575cd', linewidth=0, alpha=0.25, zorder=3)
ax.add_patch(legend_title_band)

ax.text(LEG_X + LEG_W / 2, LEG_Y + LEG_H - 0.17, 'LEGEND',
        ha='center', va='center', fontsize=8.5, fontweight='bold',
        color='#4527a0', fontfamily='DejaVu Sans', zorder=4)

# Separator line
ax.plot([LEG_X + 0.1, LEG_X + LEG_W - 0.1],
        [LEG_Y + LEG_H - 0.34, LEG_Y + LEG_H - 0.34],
        color='#9575cd', linewidth=1.0, alpha=0.5, zorder=4)

legend_items = [
    (C_STEP,  C_STEP_EDGE,  '#fff',       'Processing step'),
    (C_INPUT, C_INPUT_EDGE, C_INPUT_TEXT, 'User input'),
    (C_DATA,  C_DATA_EDGE,  C_DATA_TEXT,  'Intermediate data'),
    (C_DB,    C_DB_EDGE,    C_DB_TEXT,    'Output / result'),
    (C_VIZ,   C_VIZ_EDGE,   C_VIZ_TEXT,   'Interactive tool'),
]
item_w = 1.68
item_gap = 0.10
total_items_w = len(legend_items) * item_w + (len(legend_items) - 1) * item_gap
item_start_x = LEG_X + (LEG_W - total_items_w) / 2
for i, (fc, ec, tc, lbl) in enumerate(legend_items):
    lx_i = item_start_x + i * (item_w + item_gap)
    rbox(lx_i, LEG_Y + 0.12, item_w, 0.36, fc, ec, lbl,
         fontsize=7, text_color=tc, bold=False, radius=0.09, zorder=5)

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 PANEL  (right)
# ═══════════════════════════════════════════════════════════════════════════════
panel(9.82, PANEL_Y, 9.9, PANEL_H, C_PHASE2_BG, '#29b6f6',
      'Phase 2 — Cell Type Annotation', '#01579b')

CX2 = 14.77   # centre-x of Phase 2 column

def p2arr(y0, y1):
    arr(CX2, y0, CX2, y1)

# ── Differential Expression File  (starts at 13.60, title band bottom at 14.95 → 1.35 gap) ──
rbox(CX2 - 3.3, 13.60, 6.6, 0.78, C_INPUT, C_INPUT_EDGE,
     'Differential Expression File',
     subtitle='Seurat  ·  Scanpy  ·  SpaceRanger  ·  generic CSV',
     fontsize=10, text_color=C_INPUT_TEXT, sub_color='#bf360c', bold=True, radius=0.22)

p2arr(13.60, 13.34)

# ── DEG Parser ──
rbox(CX2 - 3.3, 12.68, 6.6, 0.63, C_STEP, C_STEP_EDGE,
     'DEG Parser',
     fontsize=10, text_color='#ffffff')

p2arr(12.68, 12.41)

# ── Gene Alias Resolution ──
rbox(CX2 - 3.3, 11.76, 6.6, 0.63, C_STEP, C_STEP_EDGE,
     'Gene Alias Resolution',
     fontsize=10, text_color='#ffffff')

p2arr(11.76, 11.48)

# ── Marker Enrichment Test ──
rbox(CX2 - 3.3, 10.58, 6.6, 0.82, C_STEP, C_STEP_EDGE,
     'Marker Enrichment Test',
     subtitle='L1: Hypergeometric  ·  L2: Weighted + BH FDR',
     fontsize=10, text_color='#ffffff', sub_color='#e0cfff')

p2arr(10.58, 10.32)

# ── Hierarchical Annotator ──
rbox(CX2 - 3.3, 9.46, 6.6, 0.78, C_STEP, C_STEP_EDGE,
     'Hierarchical Annotator',
     subtitle='CL DAG traversal',
     fontsize=10, text_color='#ffffff', sub_color='#e0cfff')

p2arr(9.46, 9.18)

# ── Top Annotation Summary ──
rbox(CX2 - 3.3, 8.40, 6.6, 0.68, C_STEP, C_STEP_EDGE,
     'Top Annotation Summary',
     fontsize=10, text_color='#ffffff')

p2arr(8.40, 8.13)

# ── Composition Scoring ──
rbox(CX2 - 3.3, 7.35, 6.6, 0.68, C_STEP, C_STEP_EDGE,
     'Composition Scoring',
     fontsize=10, text_color='#ffffff')

p2arr(7.35, 7.08)

# ── HTML Report Generator ──
rbox(CX2 - 3.1, 6.32, 6.2, 0.65, C_STEP, C_STEP_EDGE,
     'HTML Report Generator',
     fontsize=10, text_color='#ffffff')

# ── Output fan-out ──
OUT_Y = 4.62
out_items = [
    'Annotation\nSummary',
    'Enrichment\nResults',
    'CL Hierarchy',
    'Composition\nScores',
    'Interactive\nHTML Report',
]
n_out = len(out_items)
out_w = 1.52
out_gap = 0.16
total_out_w = n_out * out_w + (n_out - 1) * out_gap
out_start_x = CX2 - total_out_w / 2

for i, title_o in enumerate(out_items):
    ox = out_start_x + i * (out_w + out_gap)
    rbox(ox, OUT_Y, out_w, 0.86, C_DB, C_DB_EDGE,
         title_o,
         fontsize=7.5, text_color=C_DB_TEXT,
         bold=True, radius=0.14, zorder=4)
    arr(CX2, 6.32, ox + out_w / 2, OUT_Y + 0.86,
        color='#78909c', lw=0.9, mutation_scale=10)

# ─── Outputs legend — framed box at bottom of Phase 2 ─────────────────────────
OUT_LEG_X, OUT_LEG_Y = 9.94, 0.35
OUT_LEG_W, OUT_LEG_H = 9.60, 1.40

out_legend_frame = FancyBboxPatch((OUT_LEG_X, OUT_LEG_Y), OUT_LEG_W, OUT_LEG_H,
    boxstyle="round,pad=0.0,rounding_size=0.18",
    facecolor='#e8f5e9', edgecolor='#29b6f6', linewidth=1.4, zorder=2)
ax.add_patch(out_legend_frame)

# Title band
out_legend_title_band = FancyBboxPatch(
    (OUT_LEG_X, OUT_LEG_Y + OUT_LEG_H - 0.34), OUT_LEG_W, 0.34,
    boxstyle="round,pad=0.0,rounding_size=0.15",
    facecolor='#29b6f6', edgecolor='#29b6f6', linewidth=0, alpha=0.25, zorder=3)
ax.add_patch(out_legend_title_band)

ax.text(OUT_LEG_X + OUT_LEG_W / 2, OUT_LEG_Y + OUT_LEG_H - 0.17,
        'PIPELINE STEPS',
        ha='center', va='center', fontsize=8.5, fontweight='bold',
        color='#01579b', fontfamily='DejaVu Sans', zorder=4)

ax.plot([OUT_LEG_X + 0.1, OUT_LEG_X + OUT_LEG_W - 0.1],
        [OUT_LEG_Y + OUT_LEG_H - 0.34, OUT_LEG_Y + OUT_LEG_H - 0.34],
        color='#29b6f6', linewidth=1.0, alpha=0.5, zorder=4)

callout = (
    '① DEG Parser   ② Gene Alias Resolution   ③ Marker Enrichment Test\n'
    '④ Hierarchical Annotation   ⑤ Top Annotation Summary   ⑥ Composition Scoring   ⑦ HTML Report'
)
ax.text(OUT_LEG_X + OUT_LEG_W / 2, OUT_LEG_Y + 0.62, callout,
        ha='center', va='center', fontsize=7.5, color='#37474f',
        fontfamily='DejaVu Sans', zorder=6, linespacing=1.6)

# ═══════════════════════════════════════════════════════════════════════════════
# CROSS-PANEL CONNECTOR: master_db → MarkerEnrichmentTest
# ═══════════════════════════════════════════════════════════════════════════════
CONN_Y = 10.99   # mid-y of Marker Enrichment Test box (10.58 + 0.82/2)

X_DB_RIGHT   = 0.55 + 8.28       # 8.83 — right edge of master_db
X_ME_LEFT    = CX2 - 3.3         # 11.47 — left edge of MarkerEnrichmentTest
DB_MID_X     = 0.55 + 8.28 / 2   # 4.69
DB_MID_Y     = 5.16 + 0.82 / 2   # 5.57 — mid-y of master_db box

# Vertical segment: right edge of master_db up to connector y-level
ax.plot([X_DB_RIGHT, X_DB_RIGHT], [DB_MID_Y, CONN_Y],
        color=C_ARROW_DB, lw=1.8, zorder=6, solid_capstyle='round')

# Horizontal arrow across to Phase 2
fa_cross = FancyArrowPatch(
    (X_DB_RIGHT, CONN_Y), (X_ME_LEFT, CONN_Y),
    arrowstyle='->', color=C_ARROW_DB, linewidth=2.0,
    mutation_scale=14, zorder=7
)
ax.add_patch(fa_cross)

# Junction dot
ax.plot(X_DB_RIGHT, DB_MID_Y, 'o', color=C_ARROW_DB, ms=4.5, zorder=8)

# Label above midpoint of horizontal arrow
MID_CONN_X = (X_DB_RIGHT + X_ME_LEFT) / 2
ax.text(MID_CONN_X, CONN_Y + 0.12, 'reference marker database',
        ha='center', va='bottom', fontsize=8, color=C_ARROW_DB,
        fontweight='normal', fontfamily='DejaVu Sans', zorder=8)

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
