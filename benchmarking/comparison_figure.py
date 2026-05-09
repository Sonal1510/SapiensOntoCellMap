#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking — Comparison Figure
=====================================================
Generates a publication-quality figure per dataset:

  Bubble matrix
    Ground-truth label (y) x SapiensOntoCellMap Top_Cell_Type (x, shown as S1/S2/S3 codes)
    Bubble size proportional to sqrt(cell/spot count); colour by GT lineage.
    Dashed ring = discordant pairs (no word overlap, n > threshold).
    Reference table below matrix: Code -> full SOM label.
    Lineage colour legend to the right of the matrix.
    Accuracy badge in top-right corner of the matrix.

  Expression dot plot (Panel B -- only drawn when expression h5 is found)
    Bonafide marker genes (y) x SOM clusters grouped by cell type (x)
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
import matplotlib.patches as mpatches
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
    "font.size":         8,
    "axes.linewidth":    0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype":      42,    # editable text in Illustrator/Inkscape
    "svg.fonttype":      "none",
})

# -- Palette -------------------------------------------------------------------
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

# Lineage colours -- print-safe, colorbrewer-derived
LINEAGE_COLORS = {
    # Tumour / epithelial
    "tumor":         "#d73027",
    "invasive":      "#d73027",
    "dcis":          "#f46d43",
    "luminal":       "#fdae61",
    "myoepithelial": "#fee090",
    "apocrine":      "#ffffbf",
    "keratinocyte":  "#d6604d",
    "epithelial":    "#f46d43",
    # Stroma
    "fibroblast":    "#4393c3",
    "caf":           "#2166ac",
    "stromal":       "#4393c3",
    "pericyte":      "#74add1",
    "smooth muscle": "#abd9e9",
    "endothelial":   "#fdae61",
    # Immune
    "macrophage":    "#542788",
    "mononuclear":   "#542788",
    "monocyte":      "#7b3294",
    "myeloid":       "#9970ab",
    "dendritic":     "#c2a5cf",
    "mast":          "#1b7837",
    "plasma":        "#2166ac",
    "b cell":        "#4575b4",
    "t cell":        "#878787",
    "lymphocyte":    "#636363",
    "nk":            "#636363",
    "neutrophil":    "#4d4d4d",
    "granulocyte":   "#4d4d4d",
    # Other
    "melanocyte":    "#1b7837",
    "red blood":     "#a50026",
}

# Grouping of lineage keys into legend categories
LINEAGE_CATEGORIES = {
    "Tumor / Epithelial": ["tumor", "invasive", "dcis", "luminal", "myoepithelial",
                           "apocrine", "keratinocyte", "epithelial"],
    "Stroma":             ["fibroblast", "caf", "stromal", "pericyte",
                           "smooth muscle", "endothelial"],
    "Immune":             ["macrophage", "mononuclear", "monocyte", "myeloid",
                           "dendritic", "mast", "plasma", "b cell", "t cell",
                           "lymphocyte", "nk", "neutrophil", "granulocyte"],
    "Other":              ["melanocyte", "red blood"],
}

# Bonafide marker genes per lineage (literature-curated, breast cancer + blood)
BONAFIDE_MARKERS = {
    "Tumor":        ["ERBB2", "ESR1", "MKI67"],
    "Myoepithelial":["ACTA2", "TP63", "KRT17"],
    "Luminal":      ["KRT8", "KRT18", "FOXA1"],
    "CAF":          ["COL1A1", "FAP", "PDGFRA"],
    "Endothelial":  ["PECAM1", "VWF", "CDH5"],
    "Pericyte":     ["RGS5", "PDGFRB", "MCAM"],
    "T cell":       ["CD3E", "CD8A", "FOXP3"],
    "B / Plasma":   ["CD19", "JCHAIN", "IGHG1"],
    "Macrophage":   ["CD68", "CD163", "MRC1"],
    "Dendritic":    ["CD1C", "CLEC9A", "ITGAX"],
    "Mast":         ["TPSAB1", "CPA3", "KIT"],
}
ALL_MARKER_GENES = [g for genes in BONAFIDE_MARKERS.values() for g in genes]


def _lineage_color(label: str) -> str:
    lo = str(label).lower()
    for key, col in LINEAGE_COLORS.items():
        if key in lo:
            return col
    return "#aaaaaa"


def _truncate(s: str, n: int = 28) -> str:
    return s if len(s) <= n else s[:n - 1] + "..."


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


# -- Ground-truth definitions --------------------------------------------------
# Maps SOM cluster label -> canonical ground-truth cell type name.
# These are the Seurat/10x tutorial labels for PBMC3k, and published
# cell type annotations for the spatial datasets.

GT_PBMC3K = {
    # CellRanger 1.x kmeans/8_clusters ordering -- verified against mean counts in
    # benchmarking/results/pbmc3k/pbmc3k_deg_converted.csv.
    # Keys match summary_df["Cluster"] format ("Cluster N").
    #
    # Cluster sizes: C1=267, C2=470, C3=348, C4=1387, C5=1, C6=10, C7=1, C8=216
    # C5/C6/C7 are tiny platelet/megakaryocyte micro-clusters (kmeans k=8 artefact).
    #
    # C4 label: large mixed T cell blob (naive CD4, memory CD4, CD8 all merged).
    #   CD4 is NOT significantly upregulated vs other clusters (LFC=-92, p=1.0).
    #   Key upregulated genes: CD3D, CD3E, IL7R, CCR7, LEF1, TCF7, CD8A, CD8B.
    #   -> "alpha-beta T cell" is the broadest correct label at this resolution.
    #
    # C5 (1 cell): GP9=16, GP1BA=2, ITGA2B=3, SDPR=44 -> megakaryocyte (GP9+ distinguishes)
    # C6 (10 cells): CLU=8.7, PPBP=46, GP9=3.9 -> platelet (CLU dominant, GP9 low)
    # C7 (1 cell): GP9=0, GP1BA=0, PPBP=22, TMSB4X=10, GNG11=9 -> platelet (no megakaryocyte markers)
    "Cluster 1":  "natural killer cell",   # NKG7, GNLY, GZMB, PRF1, FGFBP2
    "Cluster 2":  "CD14-positive monocyte",  # S100A9, S100A8, LYZ, CD14, FCN1
    "Cluster 3":  "B cell",                # CD79A, MS4A1, CD79B, HLA-DQA1, TCL1A
    "Cluster 4":  "alpha-beta T cell",     # CD3D, CD3E, IL7R, CCR7, LEF1, TCF7, CD8A, CD8B
    "Cluster 5":  "megakaryocyte",         # GP9=16, GP1BA=2, ITGA2B=3, SDPR=44, PPBP=36
    "Cluster 6":  "platelet",              # PPBP=46, CLU=8.7, GNG11=10.6, GP9=3.9
    "Cluster 7":  "platelet",              # PPBP=22, TMSB4X=10, GNG11=9 -- GP9=0, GP1BA=0
    "Cluster 8":  "non-classical monocyte",  # FCGR3A=6.4, LST1=15.6, MS4A7, AIF1, CDKN1C
}

# Atera GT: maps cell_groups.csv "group" label -> canonical cell type name for accuracy matching.
# "Unassigned" excluded from accuracy calculation (mapped to None).
GT_ATERA_BREAST_CANCER = {
    "11q13 Invasive Tumor Cells":           "invasive tumor cell",
    "11q13 Invasive Tumor Cells (Mitotic)": "invasive tumor cell",
    "11q13 Invasive Tumor Cells (G1/S)":    "invasive tumor cell",
    "CAFs, DCIS Associated":                "cancer associated fibroblast",
    "CAFs, Invasive Associated":            "cancer associated fibroblast",
    "CXCL14+ Fibroblasts":                  "fibroblast",
    "Luminal-like Amorphous DCIS Cells":    "luminal epithelial cell of mammary gland",
    "Basal-like Structured DCIS Cells":     "basal cell",
    "Myoepithelial Cells":                  "myoepithelial cell",
    "Apocrine Cells":                       "apocrine cell",
    "Endothelial Cells":                    "endothelial cell",
    "Pericytes":                            "pericyte",
    "T Lymphocytes":                        "T cell",
    "B Cells":                              "B cell",
    "Plasma Cells":                         "plasmablast",
    "Macrophages":                          "macrophage",
    "Myeloid Cells":                        "myeloid cell",
    "Dendritic Cells":                      "dendritic cell",
    "Mast Cells":                           "mast cell",
    "Unassigned":                           None,
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
    "atera_breast_cancer": {
        "display_name":   "Atera WTA Preview -- FFPE Human Breast Cancer",
        "platform":       "10x Atera (whole-transcriptome in situ, dev preview)",
        "ground_truth":   GT_ATERA_BREAST_CANCER,
        "gt_type":        "barcode",
        "gt_csv":         _RESULTS_DIR / "atera_breast_cancer" / "atera_cell_groups.csv",
        "gt_barcode_col": "cell_id",
        "gt_label_col":   "group",
        "clusters_csv":   _RESULTS_DIR / "atera_breast_cancer" / "atera_clusters.csv",
        "som_sample":     "atera_breast_cancer",
        "h5_path":        Path("/Volumes/shainlab/Sonal/sapiensontocellmap_atera/cell_feature_matrix.h5"),
        "h5_format":      "h5",
    },
}


# -- Load SOM results ----------------------------------------------------------

def load_som_results(sample_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (top_summary_df, sig_results_df).

    Handles two possible output directory layouts:
      - results/<sample>/<sample>_top_annotation_summary.csv  (flat)
      - results/<sample>/<sample>/<sample>_top_annotation_summary.csv (nested, from run_annotation)
    """
    sample_dir_flat   = _RESULTS_DIR / sample_name
    sample_dir_nested = _RESULTS_DIR / sample_name / sample_name

    summary_path = None
    for candidate_dir in [sample_dir_flat, sample_dir_nested]:
        p = candidate_dir / f"{sample_name}_top_annotation_summary.csv"
        if p.exists():
            summary_path = p
            sample_dir   = candidate_dir
            break

    if summary_path is None:
        logger.error(f"SOM results not found: {sample_dir_flat / (sample_name + '_top_annotation_summary.csv')}")
        logger.error("Run: python benchmarking/run_sapiensonto.py first.")
        sys.exit(1)

    sig_path = sample_dir / f"{sample_name}_all_tissue_level2_sig_results.csv"
    summary = pd.read_csv(summary_path)
    sig = pd.read_csv(sig_path) if sig_path.exists() else pd.DataFrame()
    return summary, sig


# -- Expression matrix loading -------------------------------------------------

def load_expression_matrix(h5_path: Path, h5_format: str) -> tuple[np.ndarray, list, list]:
    """
    Returns (X_lognorm, barcodes, gene_names) where X is cells x genes (dense float32).
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
                ).T.tocsr()  # cells x genes
            # SpaceRanger h5: /matrix group same structure
            else:
                raise ValueError(f"Unrecognised h5 structure in {h5_path}")

    elif h5_format in ("mtx", "mtx_gz"):
        try:
            from scipy.io import mmread
        except ImportError:
            logger.error("scipy required: pip install scipy")
            sys.exit(1)
        import gzip
        mtx_dir  = h5_path.parent
        if h5_format == "mtx_gz":
            with gzip.open(str(h5_path), "rb") as fh:
                X = sp.csr_matrix(mmread(fh).T)  # cells x genes
            # barcodes and features are also .gz in 10x MTX gz bundles
            barcodes_f = mtx_dir / "barcodes.tsv.gz"
            features_f = mtx_dir / "features.tsv.gz"
            if barcodes_f.exists():
                barcodes = [l.strip() for l in gzip.open(str(barcodes_f), "rt")]
            else:
                barcodes = [l.strip() for l in open(mtx_dir / "barcodes.tsv")]
            if features_f.exists():
                gene_names = [l.strip().split("\t")[1] for l in gzip.open(str(features_f), "rt")]
            else:
                gene_names = [l.strip().split("\t")[1] for l in open(mtx_dir / "features.tsv")]
        else:
            X        = sp.csr_matrix(mmread(str(h5_path)).T)  # cells x genes
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


# -- Panel A: Bubble matrix ----------------------------------------------------

def draw_bubble_matrix(ax, summary_df: pd.DataFrame, ground_truth: dict,
                       max_cols: int = 20, accuracy: float = None):
    """Draw bubble matrix with short x-axis codes (S1, S2, ...).

    Returns:
        som_code_map: dict mapping code (e.g. 'S1') -> full SOM label string
        gt_labels_used: list of GT labels actually plotted (for legend filtering)
    """
    _style_ax(ax)

    # Build cluster->GT and cluster->SOM maps
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
        return {}, []

    # Count per GT x SOM pair (each cluster = 1 unit)
    ct = df.groupby(["GT", "SOM"]).size().reset_index(name="n")
    gt_order  = ct.groupby("GT")["n"].sum().sort_values(ascending=True).index.tolist()
    som_order = ct.groupby("SOM")["n"].sum().sort_values(ascending=False).index.tolist()[:max_cols]
    ct = ct[ct["SOM"].isin(som_order)]

    # Build short code map: S1, S2, ... in order of som_order
    som_code_map = {f"S{i+1}": label for i, label in enumerate(som_order)}
    som_to_code  = {label: f"S{i+1}" for i, label in enumerate(som_order)}

    gt_pos  = {v: i for i, v in enumerate(gt_order)}
    som_pos = {v: i for i, v in enumerate(som_order)}
    max_n   = max(ct["n"].max(), 1)
    bubble_scale = 1800

    gt_labels_used = []
    for _, row in ct.iterrows():
        x = som_pos.get(row["SOM"], -1)
        y = gt_pos.get(row["GT"], -1)
        if x < 0 or y < 0:
            continue
        if row["GT"] not in gt_labels_used:
            gt_labels_used.append(row["GT"])
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

    # x-axis: short codes only (S1, S2, ...) -- no rotation needed
    ax.set_xticks(range(len(som_order)))
    ax.set_xticklabels([som_to_code[s] for s in som_order],
                       fontsize=7.5, color=P["fg"])
    ax.set_yticks(range(len(gt_order)))
    ax.set_yticklabels(gt_order, fontsize=7.5, color=P["fg"])
    ax.set_xlim(-0.6, len(som_order) - 0.4)
    ax.set_ylim(-0.6, len(gt_order)  - 0.4)
    ax.set_xlabel("SapiensOntoCellMap Annotation (code)", fontsize=8.5,
                  color=P["fg"], labelpad=5)
    ax.set_ylabel("Published ground-truth label", fontsize=8.5,
                  color=P["fg"], labelpad=5)

    # Bubble size legend -- top-left corner of the matrix (avoids overlapping bubbles)
    for ln in [1, 3, 5]:
        s = bubble_scale * (ln / max_n) ** 0.5
        ax.scatter([], [], s=s, color="#aaaaaa", edgecolors="white",
                   linewidths=0.4, label=str(ln))
    ax.legend(title="Cluster count", title_fontsize=7, fontsize=7,
              loc="upper left", frameon=True, framealpha=0.9,
              edgecolor=P["spine"], facecolor=P["bg_panel"],
              bbox_to_anchor=(0.0, 1.0), bbox_transform=ax.transAxes)

    # Accuracy badge -- top-right corner of matrix axes
    if accuracy is not None:
        badge_txt = f"Accuracy: {accuracy:.0%}"
        ax.text(0.98, 0.98, badge_txt,
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8.5, fontweight="bold", color=P["fg"],
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#333333", linewidth=0.8))

    return som_code_map, gt_labels_used


def draw_reference_table(fig, ax_matrix, som_code_map: dict,
                         table_top: float, table_height: float,
                         left: float, right: float):
    """Draw a code->label reference table in figure coordinates below the matrix.

    Args:
        fig: matplotlib Figure
        ax_matrix: the bubble matrix Axes (used to compute left/right bounds)
        som_code_map: {code: full_label} ordered dict
        table_top: figure y coordinate for top of the table (0-1)
        table_height: figure height to allocate for table (0-1)
        left: figure x coordinate for left edge of table
        right: figure x coordinate for right edge of table
    """
    codes  = list(som_code_map.keys())
    labels = list(som_code_map.values())
    n_rows = len(codes)
    if n_rows == 0:
        return

    # Create a hidden axes for the table
    ax_tbl = fig.add_axes([left, table_top - table_height, right - left, table_height])
    ax_tbl.set_axis_off()

    # Build table data
    table_data = [[code, label] for code, label in zip(codes, labels)]
    col_labels = ["Code", "SapiensOntoCellMap Annotation"]

    tbl = ax_tbl.table(
        cellText=table_data,
        colLabels=col_labels,
        loc="upper left",
        cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.0)

    # Style header row
    for col_i in range(2):
        cell = tbl[0, col_i]
        cell.set_facecolor("#dddddd")
        cell.set_text_props(fontweight="bold", color="#1a1a1a", fontsize=7.0)
        cell.set_edgecolor("#bbbbbb")
        cell.set_linewidth(0.5)

    # Style data rows with alternating shading
    for row_i in range(1, n_rows + 1):
        bg = "#f5f5f5" if row_i % 2 == 1 else "#ffffff"
        for col_i in range(2):
            cell = tbl[row_i, col_i]
            cell.set_facecolor(bg)
            cell.set_edgecolor("#dddddd")
            cell.set_linewidth(0.4)
            cell.set_text_props(fontsize=7.0, color="#1a1a1a")

    # Column widths: code column narrow, label column wide
    tbl.auto_set_column_width([0, 1])
    for row_i in range(n_rows + 1):
        tbl[row_i, 0].set_width(0.07)
        tbl[row_i, 1].set_width(0.93)

    # Section header text above the table
    ax_tbl.text(0.0, 1.02,
                "Reference: SapiensOntoCellMap Annotation Codes",
                transform=ax_tbl.transAxes,
                fontsize=7.0, fontweight="bold", color=P["fg_dim"],
                va="bottom",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="none",
                          edgecolor="none"))


def draw_lineage_legend(fig, gt_labels_used: list,
                        legend_left: float, legend_bottom: float,
                        legend_width: float, legend_height: float):
    """Draw a color-swatch lineage legend panel to the right of the matrix.

    Only shows lineages that appear in gt_labels_used.
    """
    ax_leg = fig.add_axes([legend_left, legend_bottom, legend_width, legend_height])
    ax_leg.set_axis_off()

    # Determine which lineage keys appear in gt_labels_used
    gt_labels_lower = [l.lower() for l in gt_labels_used]

    def label_matches_key(key):
        return any(key in gl for gl in gt_labels_lower)

    # Build legend entries: category header + matching lineages
    entries = []  # list of (label, color, is_header)
    for cat_name, keys in LINEAGE_CATEGORIES.items():
        cat_entries = []
        for key in keys:
            if label_matches_key(key):
                cat_entries.append((key.title(), LINEAGE_COLORS[key]))
        if cat_entries:
            entries.append((cat_name, None, True))
            entries.extend([(lbl, col, False) for lbl, col in cat_entries])

    if not entries:
        # Fallback: show all used colors
        seen = set()
        for gl in gt_labels_used:
            col = _lineage_color(gl)
            if col not in seen:
                seen.add(col)
                entries.append((gl.title(), col, False))

    n_entries = len(entries)
    if n_entries == 0:
        return

    row_h = 1.0 / max(n_entries + 2, 1)
    y     = 1.0 - row_h * 0.5

    ax_leg.text(0.0, y, "Bubble color = GT lineage",
                transform=ax_leg.transAxes,
                fontsize=7, fontweight="bold", color=P["fg"],
                va="center")
    y -= row_h

    for label, color, is_header in entries:
        if is_header:
            ax_leg.text(0.0, y, label,
                        transform=ax_leg.transAxes,
                        fontsize=7, fontweight="bold", color=P["fg_dim"],
                        va="center", style="italic")
        else:
            patch = mpatches.Rectangle((0.0, y - row_h * 0.35),
                                       0.08, row_h * 0.7,
                                       transform=ax_leg.transAxes,
                                       clip_on=False,
                                       facecolor=color, edgecolor="white",
                                       linewidth=0.4)
            ax_leg.add_patch(patch)
            ax_leg.text(0.12, y, label,
                        transform=ax_leg.transAxes,
                        fontsize=7, color=P["fg"], va="center")
        y -= row_h


# -- Panel B: Expression dot plot ----------------------------------------------

def draw_dot_plot(ax, summary_df: pd.DataFrame, X: np.ndarray,
                  barcodes: list, genes_in_matrix: list,
                  ground_truth: dict):
    _style_ax(ax)

    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    cluster_to_som = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))

    # Order clusters: group by SOM cell type, then cluster id
    clusters = sorted(cluster_to_som.keys(),
                      key=lambda c: (cluster_to_som.get(c, ""), c))

    # Map barcode -> cluster (barcodes have cluster index in CellRanger format)
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

    # Build gene order from BONAFIDE_MARKERS -- only genes present in matrix
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
        sub = X[idxs, :]  # cells x bonafide genes
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
                       fontsize=7, color=P["fg"])
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(ordered_genes, fontsize=8, color=P["fg"],
                       fontfamily="monospace")
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-0.8, n_genes - 0.2)
    ax.invert_yaxis()

    # Lineage colour bars + labels -- left of y-axis, no overlap
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
                ha="right", va="center", fontsize=8, color=lc,
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


# -- Barcode-level GT join -----------------------------------------------------

def build_cluster_gt_from_barcodes(gt_csv: Path, clusters_csv: Path,
                                   gt_label_map: dict,
                                   gt_barcode_col: str = "barcode",
                                   gt_label_col: str = "cell_type",
                                   cluster_barcode_col: str = "Barcode",
                                   cluster_id_col: str = "Cluster") -> dict:
    """
    Join barcode-level ground truth with cluster assignments.

    For each cluster, assign the majority GT label (mapped through gt_label_map).
    Labels that map to None in gt_label_map are excluded from majority voting.
    Returns a dict: "Cluster N" -> canonical_gt_label_str.
    """
    gt_df = pd.read_csv(gt_csv)
    cc_df = pd.read_csv(clusters_csv)

    # Rename to common internal column names
    gt_df = gt_df.rename(columns={gt_barcode_col: "barcode", gt_label_col: "raw_label"})
    cc_df = cc_df.rename(columns={cluster_barcode_col: "barcode", cluster_id_col: "cluster"})

    merged = cc_df.merge(gt_df[["barcode", "raw_label"]], on="barcode", how="left")
    merged["raw_label"] = merged["raw_label"].fillna("Unassigned")

    cluster_to_gt = {}
    for cluster, grp in merged.groupby("cluster"):
        # Exclude labels that map to None
        labelled = grp[grp["raw_label"].map(lambda x: gt_label_map.get(x) is not None)]
        if labelled.empty:
            continue
        majority_raw   = labelled["raw_label"].mode()[0]
        majority_frac  = (labelled["raw_label"] == majority_raw).sum() / len(grp)
        canonical      = gt_label_map.get(majority_raw)
        if canonical is not None and majority_frac > 0.1:
            cluster_to_gt[f"Cluster {cluster}"] = canonical
    logger.info(f"Barcode GT join: {len(cluster_to_gt)} clusters assigned a majority GT label")
    return cluster_to_gt


# -- Accuracy computation ------------------------------------------------------

def compute_accuracy(summary_df: pd.DataFrame, ground_truth: dict) -> float:
    """
    Compute per-cluster accuracy using word-overlap between GT label and SOM Top_Cell_Type.
    Returns fraction of clusters where GT and SOM share at least one word token.
    Also logs per-cluster breakdown.
    """
    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    broad_col = "Broad_Type" if "Broad_Type" in summary_df.columns else None
    cluster_to_som = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))
    cluster_to_broad = {}
    if broad_col:
        cluster_to_broad = dict(zip(summary_df["Cluster"].astype(str), summary_df[broad_col]))

    cluster_to_gt = {str(k): v for k, v in ground_truth.items() if v is not None}

    n_match = 0
    total   = 0
    for cluster, gt in sorted(cluster_to_gt.items()):
        som   = cluster_to_som.get(cluster, "Unannotated")
        broad = cluster_to_broad.get(cluster, "")
        gt_tok  = set(gt.lower().replace("-", " ").replace("+", " ").split())
        som_tok = set(som.lower().replace("-", " ").replace("+", " ").split())
        broad_tok = set(str(broad).lower().replace("-", " ").replace("+", " ").split())
        match = bool(gt_tok & som_tok) or bool(gt_tok & broad_tok)
        n_match += int(match)
        total   += 1
        logger.info(
            f"  {cluster:12s} | GT: {gt:40s} | SOM: {_truncate(som, 45):45s} | "
            f"Broad: {_truncate(str(broad), 35):35s} | {'MATCH' if match else 'MISMATCH'}"
        )

    accuracy = n_match / total if total > 0 else 0.0
    logger.info(f"  Accuracy: {n_match}/{total} = {accuracy:.1%}")
    return accuracy


# -- Main figure builder -------------------------------------------------------

def make_figure(dataset_key: str, cfg: dict, dpi: int) -> Path:
    logger.info(f"\n{'='*60}")
    logger.info(f"Dataset: {cfg['display_name']}")
    logger.info(f"{'='*60}")

    summary_df, sig_df = load_som_results(cfg["som_sample"])
    ground_truth = cfg["ground_truth"]

    # For barcode-level GT datasets, build cluster->GT map via majority-vote join
    if cfg.get("gt_type") == "barcode":
        gt_csv       = cfg.get("gt_csv")
        clusters_csv = cfg.get("clusters_csv",
                          _RESULTS_DIR / cfg["som_sample"] / f"{cfg['som_sample']}_cell_clusters.csv")
        gt_barcode_col = cfg.get("gt_barcode_col", "barcode")
        gt_label_col   = cfg.get("gt_label_col", "cell_type")
        if gt_csv and Path(gt_csv).exists() and Path(clusters_csv).exists():
            ground_truth = build_cluster_gt_from_barcodes(
                Path(gt_csv), Path(clusters_csv), ground_truth,
                gt_barcode_col=gt_barcode_col,
                gt_label_col=gt_label_col,
            )
        else:
            logger.warning(
                f"Barcode GT join skipped -- missing file(s): "
                f"gt_csv={gt_csv}, clusters_csv={clusters_csv}"
            )

    logger.info("Computing per-cluster accuracy ...")
    accuracy = compute_accuracy(summary_df, ground_truth)
    logger.info(f"Overall accuracy for {dataset_key}: {accuracy:.1%}")

    h5_path   = cfg["h5_path"]
    h5_format = cfg["h5_format"]

    has_expression = Path(h5_path).exists()
    if not has_expression:
        logger.warning(f"Expression matrix not found at {h5_path} -- Panel B will be skipped.")

    # -- Dynamic figure height ------------------------------------------------
    # Estimate GT row count for bubble matrix sizing
    n_gt_rows = len({v for v in ground_truth.values() if v is not None})
    # Estimate SOM columns (up to 20)
    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    n_som_cols = min(20, summary_df[ct_col].nunique())
    # Estimate table rows (= number of unique SOM labels shown)
    n_table_rows = n_som_cols + 1  # +1 for header

    # Height breakdown (inches):
    #   title block: 0.6
    #   matrix: max(3.5, n_gt_rows * 0.35)
    #   gap: 0.3
    #   table: n_table_rows * 0.18 + 0.4 (header + padding)
    #   bottom margin: 0.5
    matrix_h    = max(3.5, n_gt_rows * 0.35)
    table_h_in  = n_table_rows * 0.18 + 0.4
    fig_h       = max(8.0, 0.6 + matrix_h + 0.3 + table_h_in + 0.5)
    fig_w       = 10.0  # fixed single-panel width

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    # -- Figure coordinate layout ---------------------------------------------
    # All values in figure fraction [0, 1].
    # Legend strip on the right: 14% of width
    # Matrix: left margin 13%, right edge at 80% (leaving 3% gap + ~17% legend)
    fig_left   = 0.13   # left margin for matrix (room for GT y-labels)
    fig_right  = 0.80   # right edge of matrix
    fig_top    = 1.0 - (0.6 / fig_h)   # below title block
    fig_top    = min(fig_top, 0.93)

    # Matrix height fraction
    matrix_frac  = matrix_h / fig_h
    table_frac   = table_h_in / fig_h
    bottom_frac  = 0.5 / fig_h
    gap_frac     = 0.3 / fig_h

    matrix_bottom = fig_top - matrix_frac
    matrix_bottom = max(matrix_bottom, table_frac + bottom_frac + gap_frac)

    # Legend strip: right of matrix
    legend_left   = fig_right + 0.03
    legend_right  = 0.99
    legend_width  = legend_right - legend_left
    legend_bottom = matrix_bottom
    legend_height = matrix_frac

    # Table: below matrix
    table_top_frac    = matrix_bottom - gap_frac
    table_bottom_frac = table_top_frac - table_frac
    table_bottom_frac = max(table_bottom_frac, bottom_frac)

    # Place main axes for bubble matrix
    ax_a = fig.add_axes([fig_left, matrix_bottom, fig_right - fig_left, matrix_frac])

    # Panel A: Bubble matrix
    som_code_map, gt_labels_used = draw_bubble_matrix(
        ax_a, summary_df, ground_truth, accuracy=accuracy
    )
    # Panel letter
    ax_a.text(-0.12, 1.05, "a", transform=ax_a.transAxes,
              fontsize=12, fontweight="bold", color=P["fg"], va="top")

    # Reference table below the matrix
    draw_reference_table(
        fig, ax_a, som_code_map,
        table_top=table_top_frac,
        table_height=table_frac,
        left=fig_left,
        right=fig_right,
    )

    # Lineage colour legend to the right of the matrix
    draw_lineage_legend(
        fig, gt_labels_used,
        legend_left=legend_left,
        legend_bottom=legend_bottom,
        legend_width=legend_width,
        legend_height=legend_height,
    )

    # Panel B -- expression dot plot (only when h5 found)
    if has_expression:
        logger.info("Loading expression matrix ...")
        X, barcodes, genes_in_matrix = load_expression_matrix(h5_path, h5_format)
        logger.info(f"  {X.shape[0]:,} cells x {X.shape[1]} bonafide genes")
        # Panel B shares the same figure -- placed below table for completeness
        dot_bottom = max(0.02, table_bottom_frac - 0.02 - matrix_frac * 0.9)
        ax_b = fig.add_axes([fig_left, dot_bottom, fig_right - fig_left, matrix_frac * 0.9])
        draw_dot_plot(ax_b, summary_df, X, barcodes, genes_in_matrix, ground_truth)
        ax_b.text(-0.12, 1.04, "b", transform=ax_b.transAxes,
                  fontsize=12, fontweight="bold", color=P["fg"], va="top")

    # -- Title block ----------------------------------------------------------
    title_y = min(fig_top + 0.04, 0.97)
    fig.text(0.5, title_y, cfg["display_name"],
             ha="center", va="top", fontsize=9, fontweight="bold", color=P["fg"])
    fig.text(0.5, title_y - 0.035,
             f"{cfg['platform']} | SapiensOntoCellMap (marker enrichment + CL ontology)",
             ha="center", va="top", fontsize=6.5, color=P["fg_dim"])

    # -- Save -----------------------------------------------------------------
    _FIGS_DIR.mkdir(parents=True, exist_ok=True)
    out_png = _FIGS_DIR / f"{dataset_key}_comparison.png"
    out_pdf = _FIGS_DIR / f"{dataset_key}_comparison.pdf"

    plt.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close()

    logger.info(f"Saved: {out_png}")
    logger.info(f"Saved: {out_pdf}")
    logger.info(f"Figure dimensions: {fig_w:.1f} x {fig_h:.1f} inches @ {dpi} DPI")
    return out_png


# -- CLI -----------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate benchmark comparison figures for SapiensOntoCellMap."
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=list(DATASET_CONFIGS.keys()),
        default=list(DATASET_CONFIGS.keys()),
        help="Which datasets to plot (default: all)",
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
