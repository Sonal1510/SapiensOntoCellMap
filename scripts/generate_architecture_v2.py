#!/usr/bin/env python3
"""
SapiensOntoCellMap — Architecture v1 vs v2 comparison diagram
Output: SapiensOntoCellMap_architecture_v2.png  (22×16 in, 180 dpi)

Layout (left → right, top → bottom):
  Column A  [0.02–0.31]  V1 Baseline (what existed)
  Column B  [0.38–0.62]  Shared / Core components (with change callouts)
  Column C  [0.69–0.98]  V2 New additions

A vertical arrow spine runs through column B.
Horizontal arrows connect A ↔ B and B ↔ C.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import os

# ── Palette ────────────────────────────────────────────────────────────────────
BG        = "#F8F9FB"
V1C       = "#D6E8F8"   # blue tint   — v1 box fill
V1E       = "#3A78AA"   # blue edge
V2C       = "#FEF0CC"   # amber tint  — v2 new box fill
V2E       = "#CC8800"   # amber edge
V2STAR    = "#B85C00"
SHAREC    = "#EAF3E6"   # green tint  — shared / improved
SHAREE    = "#4A8A58"
RETIREDC  = "#FCE8E8"   # pink        — retired in v2
RETIREDE  = "#BB3333"
OUTC      = "#EDE5F8"   # lavender    — output
OUTE      = "#6644AA"
BGCOL     = "#EDEFF4"   # section bg
HDRCOL    = "#1C2D40"
ARROWC    = "#44556A"

def rbox(ax, x, y, w, h, title, sub=None,
         fc=SHAREC, ec=SHAREE, tsz=8.5, bold=False,
         tag=None, tag_fc="#BB3333", radius=0.012):
    """Rounded rectangle with optional subtitle and corner tag."""
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle=f"round,pad=0.004,rounding_size={radius}",
                           facecolor=fc, edgecolor=ec, linewidth=1.3, zorder=3)
    ax.add_patch(patch)
    n_lines = 1 + (title.count('\n')) + (1 if sub else 0)
    cy = y + h / 2 + (0.010 if sub else 0)
    ax.text(x + w/2, cy, title, ha='center', va='center',
            fontsize=tsz, fontweight='bold' if bold else 'normal',
            color='#111122', zorder=4, linespacing=1.3)
    if sub:
        ax.text(x + w/2, y + h/2 - 0.013, sub,
                ha='center', va='center', fontsize=6.8,
                color='#445566', style='italic', zorder=4, linespacing=1.25)
    if tag:
        tw = 0.060;  th = 0.018
        tx = x + w - tw - 0.003;  ty = y + h - th - 0.003
        t = FancyBboxPatch((tx, ty), tw, th,
                           boxstyle="round,pad=0.002,rounding_size=0.003",
                           facecolor=tag_fc, edgecolor='none', zorder=5)
        ax.add_patch(t)
        ax.text(tx + tw/2, ty + th/2, tag, ha='center', va='center',
                fontsize=5.8, color='white', fontweight='bold', zorder=6)

def harrow(ax, x0, x1, y, col=ARROWC, lw=1.3):
    ax.annotate('', xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle='->', color=col, lw=lw), zorder=2)

def varrow(ax, x, y0, y1, col=ARROWC, lw=1.4):
    ax.annotate('', xy=(x, y1), xytext=(x, y0),
                arrowprops=dict(arrowstyle='->', color=col, lw=lw), zorder=2)

def sect(ax, x, y, w, h, title):
    bg = FancyBboxPatch((x, y), w, h,
                        boxstyle="round,pad=0.006,rounding_size=0.018",
                        facecolor=BGCOL, edgecolor='#C8CAD0',
                        linewidth=0.7, zorder=1)
    ax.add_patch(bg)
    ax.text(x + w/2, y + h + 0.010, title, ha='center', va='bottom',
            fontsize=10.5, fontweight='bold', color=HDRCOL, zorder=4)

def divider(ax, y, label=None, col='#AABBCC'):
    ax.axhline(y, xmin=0.01, xmax=0.99, color=col, lw=0.6, ls='--', zorder=1)
    if label:
        ax.text(0.5, y + 0.005, label, ha='center', va='bottom',
                fontsize=7.5, color=col, style='italic')

# ── Canvas ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(22, 16))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_aspect('equal'); ax.axis('off')
fig.patch.set_facecolor(BG); ax.set_facecolor(BG)

# ── Title ──────────────────────────────────────────────────────────────────────
ax.text(0.50, 0.982, 'SapiensOntoCellMap — Architecture: v1 vs v2',
        ha='center', va='top', fontsize=17, fontweight='bold', color=HDRCOL)
ax.text(0.50, 0.963, 'Sonal Rashmi  ·  Shain Lab  ·  UCSF Dermatology  ·  2026',
        ha='center', va='top', fontsize=9.5, color='#556677')

# ── Column headers ─────────────────────────────────────────────────────────────
ax.text(0.165, 0.944, 'V1  BASELINE', ha='center', fontsize=12,
        fontweight='bold', color=V1E,
        bbox=dict(fc=V1C, ec=V1E, boxstyle='round,pad=0.3', lw=1.2))
ax.text(0.500, 0.944, 'CORE  ENGINE', ha='center', fontsize=12,
        fontweight='bold', color=SHAREE,
        bbox=dict(fc=SHAREC, ec=SHAREE, boxstyle='round,pad=0.3', lw=1.2))
ax.text(0.835, 0.944, 'V2  NEW  ADDITIONS', ha='center', fontsize=12,
        fontweight='bold', color=V2STAR,
        bbox=dict(fc=V2C, ec=V2E, boxstyle='round,pad=0.3', lw=1.2))

# ═══════════════════════════════════════════════════════════════════════════════
# COLUMN GEOMETRY
# ═══════════════════════════════════════════════════════════════════════════════
A_x, A_w = 0.020, 0.290   # V1 baseline
B_x, B_w = 0.380, 0.240   # Core engine
C_x, C_w = 0.690, 0.300   # V2 new

# Section backgrounds
sect(ax, A_x - 0.005, 0.025, A_w + 0.010, 0.895, 'V1 Baseline')
sect(ax, B_x - 0.005, 0.025, B_w + 0.010, 0.895, 'Shared / Core')
sect(ax, C_x - 0.005, 0.025, C_w + 0.010, 0.895, 'V2 Additions')

# ── Central spine arrow (data flows downward through B) ────────────────────────
spine_x = B_x + B_w / 2
varrow(ax, spine_x, 0.900, 0.035, col=SHAREE, lw=2.0)

# ═══════════════════════════════════════════════════════════════════════════════
# ROW DEFINITIONS  (top → bottom)
# Rows (y centres): 0.875, 0.790, 0.700, 0.605, 0.510, 0.415, 0.310, 0.190, 0.090
# ═══════════════════════════════════════════════════════════════════════════════
BH = 0.062   # standard box height
BH_L = 0.075  # tall box
BH_XL = 0.090 # extra tall

# ── ROW 1: Databases ──────────────────────────────────────────────────────────
R1 = 0.845
# V1 — 8 databases
rbox(ax, A_x, R1, A_w, BH_XL,
     'Marker Databases  (8 sources)',
     sub='CellMarker 2.0 · PanglaoDB · CellxGene\nHuBMAP HRA · WIMMS · SCC 2020\nFibroblast Atlas · Epi Clusters',
     fc=V1C, ec=V1E, tsz=8.5, bold=True)
# Core — BaseParser (shared)
rbox(ax, B_x, R1, B_w, BH_XL,
     'BaseParser',
     sub='CL + UBERON ontology normalisation\nobsolete term replacement\nCellxGeneOntologyParser (OBO)',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 — 6 new skin-specific databases
rbox(ax, C_x, R1, C_w, BH_XL,
     '+6 Skin-Specific Databases',
     sub='ScarCellMarker × 4 (GSE 130973/163973/138669/156326)\nSkin Atlas MERFISH + scRNA-seq\n→ Total: 14 sources,  442,375 marker associations',
     fc=V2C, ec=V2E, tsz=8.5, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R1 + BH_XL/2)
harrow(ax, B_x + B_w, C_x, R1 + BH_XL/2)

divider(ax, R1 - 0.010)

# ── ROW 2: Database Validator / Master DB ────────────────────────────────────
R2 = 0.745
# V1 — no validator, basic concat
rbox(ax, A_x, R2, A_w, BH_L,
     'DB Assembly (no validation)',
     sub='concat parsers → master_cell_marker_db.csv\nNo schema enforcement\nDuplicates not removed (~32K extra rows)',
     fc=V1C, ec=V1E, tsz=8)
# Core — master DB
rbox(ax, B_x, R2, B_w, BH_L,
     'master_cell_marker_db.csv',
     sub='14 cols · UBERON / CL IDs throughout\n22,469 genes · 1,208 cell types · 266 tissues',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 — DatabaseValidator
rbox(ax, C_x, R2, C_w, BH_L,
     'DatabaseValidator',
     sub='UBERON/CL prefix rules · controlled vocab\ndeduplication (–32,720 duplicate rows removed)\nquarantine_log.csv for rejected rows',
     fc=V2C, ec=V2E, tsz=8.5, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R2 + BH_L/2)
harrow(ax, B_x + B_w, C_x, R2 + BH_L/2)

divider(ax, R2 - 0.010)

# ── ROW 3: DEG Input / Pre-processing ─────────────────────────────────────────
R3 = 0.650
# V1
rbox(ax, A_x, R3, A_w, BH,
     'DEG Input  (Seurat only)',
     sub='Hard-coded column names: gene, p_val_adj, avg_log2FC, cluster\nNo format detection · No alias resolution',
     fc=V1C, ec=V1E, tsz=8)
# Core
rbox(ax, B_x, R3, B_w, BH,
     'DEG Pre-Processing',
     sub='Filter: adj_p ≤ 0.05 · log2FC ≥ 1.0\nmean_counts auto-calibration (spatial 75th pctile)',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 additions
rbox(ax, C_x, R3, C_w, BH,
     'Format Auto-Detect  +  HGNC Aliases',
     sub='Seurat / Scanpy / generic schemas (column-name match)\n58K HGNC aliases: alias_symbol + prev_symbol\nCase normalisation → upstream + downstream',
     fc=V2C, ec=V2E, tsz=8.5, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R3 + BH/2)
harrow(ax, B_x + B_w, C_x, R3 + BH/2)

divider(ax, R3 - 0.010)

# ── ROW 4: Enrichment Test ─────────────────────────────────────────────────────
R4 = 0.548
# V1
rbox(ax, A_x, R4, A_w, BH_L,
     'MarkerEnrichmentTest  (v1)',
     sub='Hypergeometric  P(X≥k | N,K,n)\n✗  Per-cluster BH-FDR  (anti-conservative)\n✗  N = DEG count only  (underestimated)\n✗  No min-overlap filter  (k=1 allowed)',
     fc=V1C, ec=V1E, tsz=8)
# Core
rbox(ax, B_x, R4, B_w, BH_L,
     'MarkerEnrichmentTest',
     sub='scipy.stats.hypergeom.sf(k-1, N, K, n)\nWeighted_Enrichment  ·  Weighted_Recall\nCombined_Score = WE × −log₁₀(adj_p)',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 fixes
rbox(ax, C_x, R4, C_w, BH_L,
     'Statistical Corrections  (v2)',
     sub='★  Global BH-FDR across ALL (cluster×celltype) pairs\n★  min_overlap k≥2  (removes single-gene coincidences)\n★  Robust N detection: warn if N<15K + >80% sig\n★  Source weights: Exp=4 · SC=3 · Lit/Rev=2 · Comp=0.5\n★  Cross-DB multiplier: 1 + log₂(D databases)',
     fc=V2C, ec=V2E, tsz=8, bold=True, tag='FIXED')
harrow(ax, A_x + A_w, B_x, R4 + BH_L/2)
harrow(ax, B_x + B_w, C_x, R4 + BH_L/2)

divider(ax, R4 - 0.010)

# ── ROW 5: Hierarchical Annotation ────────────────────────────────────────────
R5 = 0.435
# V1 — flat only
rbox(ax, A_x, R5, A_w, BH_L,
     'Flat Enrichment Only  (v1)',
     sub='Level 1: per-database raw cell names\nLevel 2: CL-normalised cell names\nNo hierarchy traversal · No confidence scoring\nNo icicle chart',
     fc=V1C, ec=V1E, tsz=8)
# Core
rbox(ax, B_x, R5, B_w, BH_L,
     'Level 1 + Level 2\nAnnotation',
     sub='Level 1: database × cell_name (per-DB hits)\nLevel 2: CL cell_name (cross-DB consensus)\nSelected-tissue  +  all-tissue contexts',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 — HierarchicalAnnotator
rbox(ax, C_x, R5, C_w, BH_L,
     'HierarchicalAnnotator  (new class)',
     sub='CL ontology DAG traversal  (ancestor → leaf)\nConfidence(v) = k_v / K_v  per CL node\nResolution labels: broad / intermediate / fine / uncertain\nCombined_Score(v) = Σ −log₁₀(p_c) supporting descendants\nIcicle chart  per cluster in HTML report',
     fc=V2C, ec=V2E, tsz=8, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R5 + BH_L/2)
harrow(ax, B_x + B_w, C_x, R5 + BH_L/2)

divider(ax, R5 - 0.010)

# ── ROW 6: Composition ────────────────────────────────────────────────────────
R6 = 0.315
# V1 — NNLS
rbox(ax, A_x, R6, A_w, BH_L,
     'Composition:  NNLS  (retired)',
     sub='Marker DB treated as sparse expression matrix\n✗  L1-norm amplifies narrow cell types (1 marker)\n✗  Cell types with 42 markers penalised to 0.024/gene\n✗  Biologically implausible dominance of rare types',
     fc=RETIREDC, ec=RETIREDE, tsz=8, tag='RETIRED', tag_fc=RETIREDE)
# Core — shared composition concept
rbox(ax, B_x, R6, B_w, BH_L,
     'Cell Type\nComposition',
     sub='composition(c) = raw_score(c) / Σ raw_score(c)\nSums to 1.0 per cluster\nNo expression matrix required',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2 — Annotation-derived
rbox(ax, C_x, R6, C_w, BH_L,
     'Annotation-Derived Composer  (new)',
     sub='Uses enrichment scores, not marker matrix\nCL ancestor pruning via N_Supporting column\nTop-1 safety net: re-adds pruned dominant type\nFilters: obsolete terms  ·  depth < 2\n14/16 expected types confirmed on test data',
     fc=V2C, ec=V2E, tsz=8, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R6 + BH_L/2, col=RETIREDE)   # retired → replaced
harrow(ax, B_x + B_w, C_x, R6 + BH_L/2)

divider(ax, R6 - 0.010)

# ── ROW 7: Proliferative Flag + SignatureMatrix ────────────────────────────────
R7 = 0.205
# V1 — nothing
rbox(ax, A_x, R7, A_w, BH,
     'Not implemented in v1',
     sub='No cycling-cluster detection\nNo deconvolution reference output',
     fc=V1C, ec=V1E, tsz=8, tag='V1 GAP', tag_fc='#8899AA')
# Core — two outputs
rbox(ax, B_x, R7, B_w, BH,
     'Additional Outputs',
     sub='Proliferative flag  ·  SignatureMatrix\nfor downstream RCTD/PRISM integration',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)
# V2
rbox(ax, C_x, R7, C_w, BH,
     'Proliferative Flag  +  SignatureMatrix',
     sub='G2M + E2F Hallmark sets (Tirosh 2016) · ⚠ PROLIF badge in report\nKi-67 IHC alert · prevents false cycling-cell annotations\nSignatureMatrix: gene × cell_type, evidence-weighted\nfor PRISM Pillar 3 RCTD deconvolution',
     fc=V2C, ec=V2E, tsz=8, bold=True, tag='NEW')
harrow(ax, A_x + A_w, B_x, R7 + BH/2)
harrow(ax, B_x + B_w, C_x, R7 + BH/2)

divider(ax, R7 - 0.010)

# ── ROW 8: Outputs / HTML Report ──────────────────────────────────────────────
R8 = 0.058
BH8 = 0.115
# V1 HTML report
rbox(ax, A_x, R8, A_w, BH8,
     'HTML Report  (3 tabs)',
     sub='Tab 1: DEG Browser\nTab 2: Enrichment Visuals\nTab 3: Hypergeometric Result\n\n+ Top Annotation Summary CSV\n+ All / Sig Results CSVs\n(15–19 MB scRNA,  20–24 MB spatial)\nPlotly 2.26.0',
     fc=V1C, ec=V1E, tsz=8)

# Core — shared report infra
rbox(ax, B_x, R8, B_w, BH8,
     'Report Infrastructure',
     sub='Jinja2 templating\nDataTables (lazy init on tab switch)\nPlotly CDN (single load, include_plotlyjs=False)\nBase64-encoded PNG heatmaps',
     fc=SHAREC, ec=SHAREE, tsz=8.5, bold=True)

# V2 HTML report
rbox(ax, C_x, R8, C_w, BH8,
     'HTML Report  (6 tabs)',
     sub='★ Tab 1: Cell Type Summary  [UMAP auto-detect + marker heatmap]\n   Tab 2: DEG Browser\n   Tab 3: Enrichment Visuals\n   Tab 4: Hypergeometric Results\n★ Tab 5: Hierarchy  [icicle chart per cluster]\n★ Tab 6: Composition  [stacked bar + table]\n★ –75% report size for spatial data\n★ Additional CSVs: hierarchical CSV · composition_scores.csv\n★ Pub-quality matplotlib UMAP panels\nPlotly 3.1.1',
     fc=V2C, ec=V2E, tsz=7.8, bold=False, tag='NEW')
harrow(ax, A_x + A_w, B_x, R8 + BH8/2, col=OUTE)
harrow(ax, B_x + B_w, C_x, R8 + BH8/2, col=OUTE)

# ── DB Explorer (standalone, far right) ───────────────────────────────────────
# Small callout overlaid at bottom right of C column
ex = C_x + C_w + 0.005
# (no room — put a note inside C column about explorer being a separate tool)

# ═══════════════════════════════════════════════════════════════════════════════
# STATS SUMMARY BAR (bottom)
# ═══════════════════════════════════════════════════════════════════════════════
stats = [
    ('Databases',       '8 → 14',   '+75%'),
    ('Marker rows',     '~280K → 442K', '+58%'),
    ('HTML report tabs','3 → 6',    '+100%'),
    ('FDR correction',  'per-cluster', '→ global ★'),
    ('Gene matching',   'exact string', '→ HGNC aliases ★'),
    ('Composition',     'NNLS (matrix)', '→ enrichment-derived ★'),
    ('Hierarchy',       'none',     '→ CL ontology DAG ★'),
    ('Prolif. flag',    'none',     '→ G2M/E2F sets ★'),
    ('Benchmarking',    'none',     '→ PBMC3k done ★'),
]
n = len(stats)
bar_y = 0.028;  bar_h = 0.020
cell_w = 0.96 / n
for i, (lbl, v1_val, v2_val) in enumerate(stats):
    bx3 = 0.020 + i * cell_w
    is_new = '★' in v2_val
    bg_fc = V2C if is_new else '#F0F0F0'
    bg_ec = V2E if is_new else '#BBBBBB'
    patch = FancyBboxPatch((bx3, bar_y), cell_w - 0.003, bar_h,
                           boxstyle='round,pad=0.002,rounding_size=0.003',
                           facecolor=bg_fc, edgecolor=bg_ec, lw=0.8, zorder=3)
    ax.add_patch(patch)
    ax.text(bx3 + (cell_w-0.003)/2, bar_y + bar_h - 0.005, lbl,
            ha='center', va='top', fontsize=6.2, fontweight='bold',
            color='#1A2A3A', zorder=4)
    ax.text(bx3 + (cell_w-0.003)/2, bar_y + 0.010, f'{v1_val}  →  {v2_val.replace(" ★","")}',
            ha='center', va='center', fontsize=5.8,
            color=V2STAR if is_new else '#445566', zorder=4,
            fontweight='bold' if is_new else 'normal')

# ═══════════════════════════════════════════════════════════════════════════════
# LEGEND
# ═══════════════════════════════════════════════════════════════════════════════
legend_items = [
    (V1C, V1E, 'V1 baseline'),
    (SHAREC, SHAREE, 'Shared / improved'),
    (V2C, V2E, 'New in v2  ★'),
    (RETIREDC, RETIREDE, 'Retired / replaced in v2'),
    (OUTC, OUTE, 'Output artefacts'),
]
for i, (fc, ec, lbl) in enumerate(legend_items):
    lbx = 0.020 + i * 0.190
    lby = 0.932
    p = FancyBboxPatch((lbx, lby), 0.018, 0.014,
                       boxstyle='round,pad=0.002,rounding_size=0.003',
                       facecolor=fc, edgecolor=ec, lw=1.0, zorder=5)
    ax.add_patch(p)
    ax.text(lbx + 0.021, lby + 0.007, lbl, fontsize=8, va='center',
            color='#333344', zorder=5)

# ── Save ──────────────────────────────────────────────────────────────────────
out_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', 'SapiensOntoCellMap_architecture_v2.png'))
plt.tight_layout(pad=0.1)
plt.savefig(out_path, dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print(f'Saved: {out_path}')
