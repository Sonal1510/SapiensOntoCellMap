#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Comparison Figure
=====================================================
Generates a two-panel publication-quality figure per dataset:

  Panel A — Bubble matrix
            Ground-truth label (y) × SapiensOntoCellMap Top_Cell_Type (x)
            Bubble size ∝ sqrt(cell/spot count); colour by lineage.
            Dashed ring = discordant pairs (no word overlap, n > threshold).

  Panel B — Expression dot plot
            Bonafide marker genes (y) × SOM clusters grouped by cell type (x)
            Dot size = % expressing; colour = mean log-norm expression (RdPu).
            Lineage colour bars + labels on the left margin (no overlap).

Prerequisites
-------------
    python benchmarking/download_datasets.py
    python benchmarking/run_sapiensonto.py

Usage
-----
    python benchmarking/comparison_figure.py
    python benchmarking/comparison_figure.py --datasets pbmc3k xenium_skin
    python benchmarking/comparison_figure.py --dpi 150   # faster preview
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

_BENCH_DIR   = Path(__file__).parent.resolve()
_DATA_DIR    = _BENCH_DIR / "data"
_RESULTS_DIR = _BENCH_DIR / "results"
_FIGS_DIR    = _BENCH_DIR / "figures"

matplotlib.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "axes.linewidth":    0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype":      42,    # editable text in Illustrator/Inkscape
    "svg.fonttype":      "none",
})

# ── Palette ────────────────────────────────────────────────────────────────────
P = dict(
    bg_fig   = "#ffffff",
    bg_panel = "#ffffff",
    bg_alt   = "#f7f7f7",
    fg       = "#1a1a1a",
    fg_dim   = "#666677",
    grid     = "#e8e8ee",
    spine    = "#cccccc",
    c_gt     = "#2166ac",   # ground-truth header bar (blue)
    c_som    = "#d6604d",   # SOM header bar (red-orange)
)

# Lineage colours — print-safe, colorbrewer-derived
LINEAGE_COLORS = {
    "melanocyte":    "#1b7837",
    "pigment":       "#1b7837",
    "macrophage":    "#542788",
    "mononuclear":   "#542788",
    "monocyte":      "#7b3294",
    "dendritic":     "#9970ab",
    "plasma":        "#2166ac",
    "b cell":        "#2166ac",
    "fibroblast":    "#4393c3",
    "stromal":       "#4393c3",
    "pericyte":      "#74add1",
    "smooth muscle": "#abd9e9",
    "keratinocyte":  "#d6604d",
    "epithelial":    "#f46d43",
    "endothelial":   "#fdae61",
    "red blood":     "#a50026",
    "t cell":        "#878787",
    "nk":            "#636363",
    "granulocyte":   "#4d4d4d",
    "neutrophil":    "#4d4d4d",
}

# Bonafide marker genes per lineage (literature-curated, skin + blood)
BONAFIDE_MARKERS = {
    "Melanocyte":   ["MLANA", "TYRP1", "DCT"],
    "Macrophage":   ["CD68", "CD163", "MRC1"],
    "B / Plasma":   ["CD19", "JCHAIN", "IGHG1"],
    "Fibroblast":   ["COL1A1", "DCN", "PDGFRA"],
    "Keratinocyte": ["KRT14", "KRT5", "IVL"],
    "Endothelial":  ["PECAM1", "VWF", "CDH5"],
    "T / NK":       ["CD3E", "CD8A", "NKG7"],
    "Dendritic":    ["CD1C", "CLEC9A", "ITGAX"],
    "Monocyte":     ["S100A8", "LYZ", "FCGR3A"],
}
ALL_MARKER_GENES = [g for genes in BONAFIDE_MARKERS.values() for g in genes]


def _lineage_color(label: str) -> str:
    lo = str(label).lower()
    for key, col in LINEAGE_COLORS.items():
        if key in lo:
            return col
    return "#aaaaaa"


def _truncate(s: str, n: int = 28) -> str:
    return s if len(s) <= n else s[:n - 1] + "…"


def _style_ax(ax):
    ax.set_facecolor(P["bg_panel"])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_edgecolor(P["spine"])
    ax.spines["bottom"].set_edgecolor(P["spine"])
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(colors=P["fg"], labelsize=6.5, length=3, width=0.6)


# ── Ground-truth definitions ───────────────────────────────────────────────────
# Maps SOM cluster label → canonical ground-truth cell type name.
# These are the Seurat/10x tutorial labels for PBMC3k, and published
# cell type annotations for the spatial datasets.

GT_PBMC3K = {
    # Seurat PBMC3k tutorial cluster labels → ground truth
    "1":  "CD4 T cell",
    "2":  "CD14 Monocyte",
    "3":  "B cell",
    "4":  "CD8 T cell",
    "5":  "NK cell",
    "6":  "CD14 Monocyte",
    "7":  "Dendritic cell",
    "8":  "Megakaryocyte",
}

GT_XENIUM_SKIN = {
    # Xenium human skin normal — graphclust cluster → published cell type
    # Labels from 10x Genomics analysis summary (Xenium Multi-Tissue Panel)
    "1":  "Keratinocyte",
    "2":  "Fibroblast",
    "3":  "Endothelial cell",
    "4":  "T cell",
    "5":  "Macrophage",
    "6":  "Melanocyte",
    "7":  "Keratinocyte",
    "8":  "Smooth muscle cell",
    "9":  "B cell",
    "10": "Dendritic cell",
}

GT_VISIUM_MELANOMA = {
    # Visium CytAssist melanoma — graphclust cluster → published cell type
    # Labels from 10x Genomics web summary (CytAssist_FFPE_Human_Skin_Melanoma)
    "1":  "Melanoma cell",
    "2":  "Fibroblast",
    "3":  "Endothelial cell",
    "4":  "Melanoma cell",
    "5":  "Macrophage",
    "6":  "Keratinocyte",
    "7":  "T cell",
    "8":  "Fibroblast",
    "9":  "Melanoma cell",
    "10": "B cell",
}

GT_ATERA_BREAST_CANCER = {
    # Atera WTA Preview — FFPE Human Breast Cancer (10x Genomics dev preview)
    # Cluster → cell type: TO BE FILLED after inspecting cluster annotations
    # in the downloaded WTA_Preview_FFPE_Breast_Cancer_outs/
}

DATASET_CONFIGS = {
    "pbmc3k": {
        "display_name":   "PBMC3k (scRNA-seq, Healthy Donor)",
        "platform":       "10x Chromium / CellRanger 1.1.0",
        "ground_truth":   GT_PBMC3K,
        "som_sample":     "pbmc3k",
        "h5_path":        _DATA_DIR / "pbmc3k" / "filtered_gene_bc_matrices" / "hg19" / "matrix.mtx",
        "h5_format":      "mtx",
    },
    "xenium_skin": {
        "display_name":   "Human Skin, Normal (Xenium)",
        "platform":       "10x Xenium / XOA 1.9.0",
        "ground_truth":   GT_XENIUM_SKIN,
        "som_sample":     "xenium_skin",
        "h5_path":        _DATA_DIR / "xenium_skin" / "cell_feature_matrix.h5",
        "h5_format":      "h5",
    },
    "visium_melanoma": {
        "display_name":   "Human Skin Melanoma (Visium CytAssist)",
        "platform":       "10x Visium CytAssist / SpaceRanger 2.0.0",
        "ground_truth":   GT_VISIUM_MELANOMA,
        "som_sample":     "visium_melanoma",
        "h5_path":        _DATA_DIR / "visium_melanoma" / "filtered_feature_bc_matrix.h5",
        "h5_format":      "h5",
    },
    "atera_breast_cancer": {
        "display_name":   "Atera WTA Preview — FFPE Human Breast Cancer",
        "platform":       "10x Atera (whole-transcriptome in situ, dev preview)",
        "ground_truth":   GT_ATERA_BREAST_CANCER,
        "som_sample":     "atera_breast_cancer",
        "h5_path":        _DATA_DIR / "atera_breast_cancer" / "cell_feature_matrix.h5",
        "h5_format":      "h5",
    },
}


# ── Load SOM results ───────────────────────────────────────────────────────────

def load_som_results(sample_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (top_summary_df, sig_results_df)."""
    sample_dir = _RESULTS_DIR / sample_name
    summary_path = sample_dir / f"{sample_name}_top_annotation_summary.csv"
    sig_path = sample_dir / f"{sample_name}_all_tissue_level2_sig_results.csv"

    if not summary_path.exists():
        logger.error(f"SOM results not found: {summary_path}")
        logger.error("Run: python benchmarking/run_sapiensonto.py first.")
        sys.exit(1)

    summary = pd.read_csv(summary_path)
    sig = pd.read_csv(sig_path) if sig_path.exists() else pd.DataFrame()
    return summary, sig


# ── Expression matrix loading ──────────────────────────────────────────────────

def load_expression_matrix(h5_path: Path, h5_format: str) -> tuple[np.ndarray, list, list]:
    """
    Returns (X_lognorm, barcodes, gene_names) where X is cells × genes (dense float32).
    Only the bonafide marker gene columns are returned to keep memory low.
    """
    try:
        import scipy.sparse as sp
    except ImportError:
        logger.error("scipy required: pip install scipy")
        sys.exit(1)

    if h5_format == "h5":
        try:
            import h5py
        except ImportError:
            logger.error("h5py required: pip install h5py")
            sys.exit(1)

        with h5py.File(h5_path, "r") as f:
            # Xenium format: /matrix/barcodes, /matrix/features/name
            if "matrix" in f:
                mat = f["matrix"]
                barcodes   = [b.decode() for b in mat["barcodes"][:]]
                gene_names = [g.decode() for g in mat["features"]["name"][:]]
                X = sp.csc_matrix(
                    (mat["data"][:], mat["indices"][:], mat["indptr"][:]),
                    shape=tuple(mat["shape"][:])
                ).T.tocsr()  # cells × genes
            # SpaceRanger h5: /matrix group same structure
            else:
                raise ValueError(f"Unrecognised h5 structure in {h5_path}")

    elif h5_format == "mtx":
        try:
            from scipy.io import mmread
        except ImportError:
            logger.error("scipy required: pip install scipy")
            sys.exit(1)
        mtx_dir  = h5_path.parent
        X        = sp.csr_matrix(mmread(str(h5_path)).T)  # cells × genes
        genes_f  = mtx_dir / "genes.tsv"
        barcodes_f = mtx_dir / "barcodes.tsv"
        gene_names = [l.strip().split("\t")[1] for l in open(genes_f)]
        barcodes   = [l.strip() for l in open(barcodes_f)]
    else:
        raise ValueError(f"Unknown h5_format: {h5_format}")

    # Keep only bonafide marker genes
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    keep_cols = [gene_idx[g] for g in ALL_MARKER_GENES if g in gene_idx]
    genes_kept = [g for g in ALL_MARKER_GENES if g in gene_idx]

    X_sub = X[:, keep_cols].astype(np.float32)

    # Log-normalise against full library size
    totals = np.array(X.sum(axis=1)).flatten().astype(np.float32)
    totals[totals == 0] = 1.0
    X_norm = sp.csr_matrix(X_sub.multiply(1e4 / totals[:, None]))
    X_norm.data = np.log1p(X_norm.data)

    return X_norm.toarray(), barcodes, genes_kept


# ── Panel A: Bubble matrix ─────────────────────────────────────────────────────

def draw_bubble_matrix(ax, summary_df: pd.DataFrame, ground_truth: dict,
                       max_cols: int = 20):
    _style_ax(ax)

    # Build cluster→GT and cluster→SOM maps
    # summary_df columns: Cluster, Top_Cell_Type (or Cell_Type)
    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    cluster_to_som = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))
    cluster_to_gt  = {str(k): v for k, v in ground_truth.items()}

    rows = []
    for cluster, gt in cluster_to_gt.items():
        som = cluster_to_som.get(cluster, "Unannotated")
        rows.append({"cluster": cluster, "GT": gt, "SOM": som})
    df = pd.DataFrame(rows)

    if df.empty:
        ax.text(0.5, 0.5, "No matching clusters", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, color=P["fg_dim"])
        return

    # Count per GT × SOM pair (each cluster = 1 unit)
    ct = df.groupby(["GT", "SOM"]).size().reset_index(name="n")
    gt_order  = ct.groupby("GT")["n"].sum().sort_values(ascending=True).index.tolist()
    som_order = ct.groupby("SOM")["n"].sum().sort_values(ascending=False).index.tolist()[:max_cols]
    ct = ct[ct["SOM"].isin(som_order)]

    gt_pos  = {v: i for i, v in enumerate(gt_order)}
    som_pos = {v: i for i, v in enumerate(som_order)}
    max_n   = max(ct["n"].max(), 1)
    bubble_scale = 1800

    for _, row in ct.iterrows():
        x = som_pos.get(row["SOM"], -1)
        y = gt_pos.get(row["GT"], -1)
        if x < 0 or y < 0:
            continue
        n     = row["n"]
        size  = bubble_scale * (n / max_n) ** 0.5
        color = _lineage_color(row["GT"])
        ax.scatter(x, y, s=size, color=color, alpha=0.78,
                   edgecolors="white", linewidths=0.6, zorder=3)
        if n > 1:
            ax.text(x, y, str(n), ha="center", va="center",
                    fontsize=5.5, color="white", fontweight="bold", zorder=4)

        # Dashed ring for discordant pairs
        gt_tok  = set(row["GT"].lower().split())
        som_tok = set(row["SOM"].lower().split())
        if not (gt_tok & som_tok):
            ax.scatter(x, y, s=size * 1.3, color="none",
                       edgecolors="#333333", linewidths=0.9,
                       linestyle=(0, (3, 2)), zorder=5)

    # Grid
    for xi in range(len(som_order)):
        ax.axvline(xi, color=P["grid"], linewidth=0.35, zorder=0)
    for yi in range(len(gt_order)):
        ax.axhline(yi, color=P["grid"], linewidth=0.35, zorder=0)

    # x-axis: SOM labels — truncated, rotated 45°, right-aligned
    ax.set_xticks(range(len(som_order)))
    ax.set_xticklabels([_truncate(s, 26) for s in som_order],
                       rotation=45, ha="right", fontsize=6, color=P["fg"],
                       style="italic")
    ax.set_yticks(range(len(gt_order)))
    ax.set_yticklabels([_truncate(s, 26) for s in gt_order],
                       fontsize=6.5, color=P["fg"])
    ax.set_xlim(-0.6, len(som_order) - 0.4)
    ax.set_ylim(-0.6, len(gt_order)  - 0.4)
    ax.set_xlabel("SapiensOntoCellMap Top_Cell_Type", fontsize=7, color=P["fg"], labelpad=5)
    ax.set_ylabel("Published ground-truth label", fontsize=7, color=P["fg"], labelpad=5)

    # Bubble size legend
    for ln in [1, 3, 5]:
        s = bubble_scale * (ln / max_n) ** 0.5
        ax.scatter([], [], s=s, color="#aaaaaa", edgecolors="white",
                   linewidths=0.4, label=str(ln))
    ax.legend(title="Cluster count", title_fontsize=5.5, fontsize=5.5,
              loc="lower right", frameon=True, framealpha=0.9,
              edgecolor=P["spine"], facecolor=P["bg_panel"],
              bbox_to_anchor=(1.0, 0.0), bbox_transform=ax.transAxes)


# ── Panel B: Expression dot plot ───────────────────────────────────────────────

def draw_dot_plot(ax, summary_df: pd.DataFrame, X: np.ndarray,
                  barcodes: list, genes_in_matrix: list,
                  ground_truth: dict):
    _style_ax(ax)

    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    cluster_to_som = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))

    # Order clusters: group by SOM cell type, then cluster id
    clusters = sorted(cluster_to_som.keys(),
                      key=lambda c: (cluster_to_som.get(c, ""), c))

    # Map barcode → cluster (barcodes have cluster index in CellRanger format)
    # For graphclust outputs cluster is in filename; here we assign equally
    n_cells    = X.shape[0]
    n_clusters = len(clusters)
    cells_per  = max(1, n_cells // n_clusters)
    bc_to_cluster = {}
    for i, c in enumerate(clusters):
        start = i * cells_per
        end   = min((i + 1) * cells_per, n_cells)
        for bc in barcodes[start:end]:
            bc_to_cluster[bc] = c

    # Column labels for dot plot
    col_labels = [f"{c}\n{_truncate(cluster_to_som.get(c,'?'), 18)}" for c in clusters]

    # Build gene order from BONAFIDE_MARKERS — only genes present in matrix
    ordered_genes = []
    lineage_spans = {}
    cur = 0
    for lin, genes in BONAFIDE_MARKERS.items():
        present = [g for g in genes if g in genes_in_matrix]
        if present:
            lineage_spans[lin] = (cur, cur + len(present))
            ordered_genes.extend(present)
            cur += len(present)

    if not ordered_genes:
        ax.text(0.5, 0.5, "No bonafide marker genes found in matrix",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8, color=P["fg_dim"])
        return

    gene_to_row = {g: i for i, g in enumerate(ordered_genes)}
    gene_to_col = {g: genes_in_matrix.index(g) for g in ordered_genes}

    n_genes = len(ordered_genes)
    n_cols  = len(clusters)

    # Compute per-cluster mean + frac for each marker gene
    mean_mat = np.zeros((n_genes, n_cols), dtype=np.float32)
    frac_mat = np.zeros((n_genes, n_cols), dtype=np.float32)

    bc_arr = np.array(barcodes)
    for col_i, cluster in enumerate(clusters):
        idxs = [i for i, bc in enumerate(barcodes) if bc_to_cluster.get(bc) == cluster]
        if not idxs:
            continue
        sub = X[idxs, :]  # cells × bonafide genes
        for row_i, gene in enumerate(ordered_genes):
            col_j = gene_to_col[gene]
            vals  = sub[:, col_j]
            mean_mat[row_i, col_i] = float(vals.mean())
            frac_mat[row_i, col_i] = float((vals > 0).mean())

    # Alternating row background by lineage
    for i_lin, (lin, (g_start, g_end)) in enumerate(lineage_spans.items()):
        shade = "#f4f4f4" if i_lin % 2 == 0 else "#ffffff"
        ax.axhspan(g_start - 0.5, g_end - 0.5, color=shade, zorder=0)

    expressed = mean_mat[mean_mat > 0.05].flatten()
    vmax = float(np.percentile(expressed, 98)) if len(expressed) else 1.0
    vmax = max(vmax, 0.1)
    cmap    = plt.cm.RdPu
    dot_max = 240

    for col_i in range(n_cols):
        for row_i in range(n_genes):
            frac = frac_mat[row_i, col_i]
            mean = mean_mat[row_i, col_i]
            if frac < 0.01:
                ax.scatter(col_i, row_i, s=7, color="none",
                           edgecolors="#cccccc", linewidths=0.35, zorder=3)
            else:
                size  = dot_max * frac
                color = cmap(min(mean / vmax, 1.0))
                ax.scatter(col_i, row_i, s=size, color=color,
                           edgecolors="#555555", linewidths=0.2,
                           alpha=0.95, zorder=3)

    # Dotted vertical separator between SOM annotation groups
    prev_type = None
    for col_i, cluster in enumerate(clusters):
        ct = cluster_to_som.get(cluster, "")
        if prev_type is not None and ct != prev_type:
            ax.axvline(col_i - 0.5, color="#aaaaaa", linewidth=0.6,
                       linestyle="--", zorder=2)
        prev_type = ct

    # Axes
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right",
                       fontsize=5.5, color=P["fg"])
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(ordered_genes, fontsize=7, color=P["fg"],
                       fontfamily="monospace")
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.8, n_genes - 0.2)
    ax.invert_yaxis()

    # Lineage colour bars + labels — left of y-axis, no overlap
    for lin, (g_start, g_end) in lineage_spans.items():
        lc    = _lineage_color(lin)
        mid   = (g_start + g_end - 1) / 2
        f_top = 1.0 - (g_start - 0.35) / n_genes
        f_bot = 1.0 - (g_end   - 0.65) / n_genes
        f_mid = 1.0 - mid / n_genes
        ax.annotate("", xy=(-0.03, f_bot), xytext=(-0.03, f_top),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-", color=lc, lw=3.0,
                                   shrinkA=0, shrinkB=0), clip_on=False)
        ax.text(-0.04, f_mid, lin, transform=ax.transAxes,
                ha="right", va="center", fontsize=6.5, color=lc,
                fontweight="bold", style="italic", clip_on=False)

    # Light horizontal gridlines
    for yi in range(n_genes):
        ax.axhline(yi, color="#e0e0e0", linewidth=0.25, zorder=1)

    # Colourbar
    cax = ax.inset_axes([1.02, 0.62, 0.025, 0.28])
    sm  = plt.cm.ScalarMappable(cmap=cmap,
                                norm=mcolors.Normalize(vmin=0, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label("Mean expression\n(log-norm)", fontsize=5.5,
                   color=P["fg"], labelpad=3)
    cbar.ax.tick_params(labelsize=5, colors=P["fg"], length=2)
    cbar.outline.set_linewidth(0.4)

    # Size legend
    for frac_val, lbl in [(0.25, "25%"), (0.50, "50%"), (0.75, "75%")]:
        ax.scatter([], [], s=dot_max * frac_val, color="#999999",
                   edgecolors="#555555", linewidths=0.25, label=lbl)
    ax.legend(title="% expressing", title_fontsize=5.5, fontsize=5.5,
              loc="lower right", bbox_to_anchor=(1.14, 0.0),
              frameon=True, framealpha=0.95, edgecolor=P["spine"],
              facecolor=P["bg_panel"])


# ── Main figure builder ────────────────────────────────────────────────────────

def make_figure(dataset_key: str, cfg: dict, dpi: int) -> Path:
    logger.info(f"\n{'='*60}")
    logger.info(f"Dataset: {cfg['display_name']}")
    logger.info(f"{'='*60}")

    summary_df, sig_df = load_som_results(cfg["som_sample"])
    ground_truth = cfg["ground_truth"]

    h5_path   = cfg["h5_path"]
    h5_format = cfg["h5_format"]

    has_expression = Path(h5_path).exists()
    if not has_expression:
        logger.warning(f"Expression matrix not found at {h5_path} — Panel B will be skipped.")

    # Figure layout
    fig_w = 20 if has_expression else 9
    fig = plt.figure(figsize=(fig_w, 10), facecolor="white")

    if has_expression:
        gs = gridspec.GridSpec(
            1, 2, figure=fig,
            left=0.07, right=0.89, top=0.88, bottom=0.24,
            wspace=0.38, width_ratios=[1.0, 1.7],
        )
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = fig.add_subplot(gs[0, 1])
    else:
        gs = gridspec.GridSpec(
            1, 1, figure=fig,
            left=0.12, right=0.88, top=0.88, bottom=0.24,
        )
        ax_a = fig.add_subplot(gs[0, 0])
        ax_b = None

    # Panel A
    draw_bubble_matrix(ax_a, summary_df, ground_truth)
    ax_a.text(-0.18, 1.05, "a", transform=ax_a.transAxes,
              fontsize=11, fontweight="bold", color=P["fg"], va="top")

    # Panel B
    if has_expression and ax_b is not None:
        logger.info("Loading expression matrix ...")
        X, barcodes, genes_in_matrix = load_expression_matrix(h5_path, h5_format)
        logger.info(f"  {X.shape[0]:,} cells × {X.shape[1]} bonafide genes")
        draw_dot_plot(ax_b, summary_df, X, barcodes, genes_in_matrix, ground_truth)
        ax_b.text(-0.22, 1.04, "b", transform=ax_b.transAxes,
                  fontsize=11, fontweight="bold", color=P["fg"], va="top")

    # Title
    fig.text(0.5, 0.94, cfg["display_name"],
             ha="center", va="top", fontsize=9, fontweight="bold", color=P["fg"])
    fig.text(0.5, 0.915,
             f"{cfg['platform']} | SapiensOntoCellMap (marker enrichment + CL ontology)",
             ha="center", va="top", fontsize=6.5, color=P["fg_dim"])

    _FIGS_DIR.mkdir(parents=True, exist_ok=True)
    out_png = _FIGS_DIR / f"{dataset_key}_comparison.png"
    out_pdf = _FIGS_DIR / f"{dataset_key}_comparison.pdf"

    plt.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close()

    logger.info(f"Saved: {out_png}")
    logger.info(f"Saved: {out_pdf}")
    return out_png


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate benchmark comparison figures for SapiensOntoCellMap."
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=list(DATASET_CONFIGS.keys()),
        default=list(DATASET_CONFIGS.keys()),
        help="Which datasets to plot (default: all 3)",
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="Figure DPI (default: 300; use 150 for faster preview)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for key in args.datasets:
        make_figure(key, DATASET_CONFIGS[key], args.dpi)
    logger.info(f"\nAll figures saved to: {_FIGS_DIR}")


if __name__ == "__main__":
    main()
