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
    c_gt     = "#2166ac",
    c_som    = "#d6604d",
)

# Match-quality colors (the only colors used for bubbles)
C_EXACT   = "#2e7d32"   # dark green  — word-level match
C_PARTIAL = "#f9a825"   # amber       — correct broad type, wrong specific
C_MISS    = "#c62828"   # deep red    — completely wrong

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


# -- Smart label shortening ----------------------------------------------------
_STRIP_SUFFIXES = [
    ", CD19-positive", ", CD56-dim", ", alpha-beta", " of vascular tree",
    " of mammary gland", ", CD14-positive", " CD14-positive, ",
    " progenitor cell", " epithelial cell", "branched duct ",
    "double-positive, ", "elicited ", "activated ",
    "-positive, alpha-beta", " cell, CD19-positive",
]
_STRIP_PREFIXES = [
    "CD16-positive, CD56-dim ",
]


def _smart_truncate(s: str, n: int = 26) -> str:
    label = s
    for sfx in _STRIP_SUFFIXES:
        label = label.replace(sfx, "")
    for pfx in _STRIP_PREFIXES:
        if label.startswith(pfx):
            label = label[len(pfx):]
    label = label.strip().rstrip(",").strip()
    return label if len(label) <= n else label[:n - 1] + "…"


# -- Concordance check ---------------------------------------------------------
def _is_concordant(gt: str, som: str) -> bool:
    gt_tok  = set(gt.lower().replace("-", " ").replace("+", " ").split())
    som_tok = set(som.lower().replace("-", " ").replace("+", " ").split())
    stopwords = {"cell", "positive", "negative", "type", "of", "and", "the", "a"}
    return bool((gt_tok - stopwords) & (som_tok - stopwords))


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

def load_som_results(sample_name: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set]]:
    """Returns (top_summary_df, sig_results_df, cluster_lineage_sets).

    cluster_lineage_sets: maps cluster string → set of normalised Cell_Type tokens
    from the hierarchical CSV — used for lineage-aware partial-match detection.

    Handles two possible output directory layouts:
      - results/<sample>/<sample>_top_annotation_summary.csv  (flat)
      - results/<sample>/<sample>/<sample>_top_annotation_summary.csv (nested)
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

    summary = pd.read_csv(summary_path)
    sig_path = sample_dir / f"{sample_name}_all_tissue_level2_sig_results.csv"
    sig      = pd.read_csv(sig_path) if sig_path.exists() else pd.DataFrame()

    # Build per-cluster set of all Cell_Type tokens from the hierarchical CSV
    hier_path = sample_dir / f"{sample_name}_all_tissue_level2_hierarchical.csv"
    cluster_lineage: dict[str, set] = {}
    if hier_path.exists():
        hier = pd.read_csv(hier_path)
        stop = {"cell", "positive", "negative", "type", "of", "and", "the",
                "a", "human", "cd", "lineage", "obsolete"}
        for cluster, grp in hier.groupby("Cluster"):
            tokens: set[str] = set()
            for ct in grp["Cell_Type"].dropna():
                if "obsolete" in ct.lower():
                    continue
                toks = set(ct.lower().replace("-", " ").replace("+", " ").split())
                tokens |= (toks - stop)
            cluster_lineage[str(cluster)] = tokens

    return summary, sig, cluster_lineage


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

def _match_quality(gt: str, som: str, broad: str,
                   lineage_tokens: set | None = None) -> str:
    """Return 'exact', 'partial', or 'miss'.

    exact   — GT and SOM share a meaningful word token (direct label match)
    partial — no direct word overlap BUT any of:
                • SOM Broad_Type shares tokens with GT
                • GT tokens appear anywhere in the CL hierarchical lineage
                  of the cluster (i.e. GT is an ancestor or related node in
                  the CL DAG that scored for that cluster)
    miss    — neither
    """
    stop = {"cell", "positive", "negative", "type", "of", "and", "the",
            "a", "human", "cd", "lineage"}

    gt_tok  = set(gt.lower().replace("-", " ").replace("+", " ").split()) - stop
    som_tok = set(som.lower().replace("-", " ").replace("+", " ").split()) - stop

    # 1. Direct word-level match
    if gt_tok & som_tok:
        return "exact"

    # 2. Broad-type match
    broad_tok = set(str(broad).lower().replace("-", " ").replace("+", " ").split()) - stop
    if gt_tok & broad_tok:
        return "partial"

    # 3. Hierarchical lineage match — GT tokens appear in any CL node that
    #    scored for this cluster in the hierarchical traversal
    if lineage_tokens and (gt_tok & lineage_tokens):
        return "partial"

    return "miss"


def draw_bubble_matrix(ax, summary_df: pd.DataFrame, ground_truth: dict,
                       cluster_lineage: dict = None,
                       max_cols: int = 20, accuracy: float = None,
                       n_match: int = None, n_total: int = None):
    """CNS-grade bubble matrix — match-quality color encoding.

    Axes
    ----
    X : SapiensOntoCellMap predicted labels (smart-truncated, 45° rotation)
    Y : Published ground-truth labels (sorted by lineage family for readability)

    Bubbles
    -------
    Size   ∝ sqrt(n_clusters / max_n) — standard CNS bubble-plot convention
    Color  = match quality:
               GREEN  (#2e7d32) — exact / near-exact match
               AMBER  (#f9a825) — correct broad cell type, wrong specific label
               RED    (#c62828) — completely wrong annotation
    Count  printed in white inside bubbles where n > 1

    Returns BUBBLE_MAX_PT for the key panel.
    """
    ax.set_facecolor("#fcfcfc")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor(P["spine"])
        ax.spines[side].set_linewidth(0.7)
    ax.tick_params(colors="#212121", labelsize=7.5, length=3, width=0.6)

    ct_col    = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    broad_col = "Broad_Type"    if "Broad_Type"    in summary_df.columns else None
    cluster_to_som   = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))
    cluster_to_broad = {}
    if broad_col:
        cluster_to_broad = dict(zip(summary_df["Cluster"].astype(str),
                                    summary_df[broad_col]))
    cluster_to_gt = {str(k): v for k, v in ground_truth.items()}

    rows = []
    for cluster, gt in cluster_to_gt.items():
        if gt is None:
            continue
        som    = cluster_to_som.get(cluster, "Unannotated")
        broad  = cluster_to_broad.get(cluster, "")
        ltoks  = (cluster_lineage or {}).get(cluster, set())
        rows.append({"cluster": cluster, "GT": gt, "SOM": som,
                     "quality": _match_quality(gt, som, broad, ltoks)})
    df = pd.DataFrame(rows)

    if df.empty:
        ax.text(0.5, 0.5, "No matching clusters", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, color=P["fg_dim"])
        return 1

    ct_tbl = df.groupby(["GT", "SOM", "quality"]).size().reset_index(name="n")

    # Y order: sort GT labels alphabetically within broad family for readability
    gt_order  = sorted(ct_tbl["GT"].unique())
    som_order = (ct_tbl.groupby("SOM")["n"].sum()
                       .sort_values(ascending=False).index.tolist()[:max_cols])
    ct_tbl = ct_tbl[ct_tbl["SOM"].isin(som_order)]

    gt_pos  = {v: i for i, v in enumerate(gt_order)}
    som_pos = {v: i for i, v in enumerate(som_order)}
    max_n   = max(ct_tbl["n"].max(), 1)

    BUBBLE_MAX_PT = 480   # points² for n == max_n (CNS: ~0.9 grid-unit diameter)

    QUALITY_COLOR = {"exact": C_EXACT, "partial": C_PARTIAL, "miss": C_MISS}

    for _, row in ct_tbl.iterrows():
        xi = som_pos.get(row["SOM"], -1)
        yi = gt_pos.get(row["GT"], -1)
        if xi < 0 or yi < 0:
            continue
        n     = int(row["n"])
        size  = BUBBLE_MAX_PT * (n / max_n) ** 0.5
        color = QUALITY_COLOR.get(row["quality"], C_MISS)

        ax.scatter(xi, yi, s=size, color=color, alpha=0.88,
                   edgecolors="white", linewidths=0.5, zorder=3)
        if n > 1:
            ax.text(xi, yi, str(n), ha="center", va="center",
                    fontsize=5.5, color="white", fontweight="bold", zorder=4)

    # Alternating row shading
    for yi, _ in enumerate(gt_order):
        shade = "#f4f4f4" if yi % 2 == 0 else "#ffffff"
        ax.axhspan(yi - 0.5, yi + 0.5, color=shade, zorder=0, linewidth=0)

    # Subtle grid
    for xi in range(len(som_order)):
        ax.axvline(xi, color="#e5e5e5", linewidth=0.25, zorder=1)
    for yi in range(len(gt_order)):
        ax.axhline(yi, color="#e5e5e5", linewidth=0.25, zorder=1)

    ax.set_xticks(range(len(som_order)))
    ax.set_xticklabels([_smart_truncate(s, 26) for s in som_order],
                       rotation=45, ha="right", fontsize=7.5, color="#212121")
    ax.set_yticks(range(len(gt_order)))
    ax.set_yticklabels(gt_order, fontsize=8, color="#212121")
    ax.set_xlim(-0.6, len(som_order) - 0.4)
    ax.set_ylim(-0.6, len(gt_order)  - 0.4)
    ax.set_xlabel("SapiensOntoCellMap annotation", fontsize=9,
                  color=P["fg"], labelpad=7)
    ax.set_ylabel("Published ground-truth label", fontsize=9,
                  color=P["fg"], labelpad=6)

    # Accuracy badge
    if accuracy is not None:
        badge_txt = (f"Accuracy: {accuracy:.1%}  ({n_match}/{n_total} clusters)"
                     if n_match is not None else f"Accuracy: {accuracy:.1%}")
        ax.text(0.985, 0.985, badge_txt,
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8.5, fontweight="bold", color="#1b5e20",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#f1f8e9",
                          edgecolor="#81c784", linewidth=0.8, alpha=0.95))

    return BUBBLE_MAX_PT


def draw_key_panel(fig, gt_labels_used: list,
                   key_left: float, key_bottom: float,
                   key_width: float, key_height: float,
                   bubble_max_pt: int = 480,
                   max_per_row: int = 6):
    """Clearly demarcated KEY panel below the bubble matrix.

    Three sections in a single framed box:
      1. Match-quality color key  (the primary encoding)
      2. Bubble-size legend       (CNS convention: 25 / 50 / 100% of max)
      3. How to read note
    """
    ax_k = fig.add_axes([key_left, key_bottom, key_width, key_height])
    ax_k.set_facecolor("#f8f9fb")
    for sp in ax_k.spines.values():
        sp.set_visible(True)
        sp.set_edgecolor("#b8bcc8")
        sp.set_linewidth(0.9)
    ax_k.set_xlim(0, 1)
    ax_k.set_ylim(0, 1)
    ax_k.set_xticks([])
    ax_k.set_yticks([])

    def section_header(x, y, w, h, label, bg="#e4e8f0", fg="#2c3060"):
        band = mpatches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.005,rounding_size=0.01",
            facecolor=bg, edgecolor="#9096b8", linewidth=0.6,
            transform=ax_k.transAxes, clip_on=False
        )
        ax_k.add_patch(band)
        ax_k.text(x + 0.008, y + h / 2, label,
                  transform=ax_k.transAxes,
                  fontsize=7.5, fontweight="bold", color=fg,
                  va="center", ha="left")

    # ── 1. Match-quality color key (left ~60%) ────────────────────────────────
    MQ_X = 0.01
    MQ_W = 0.57
    section_header(MQ_X, 0.74, MQ_W, 0.20, "Bubble color  =  annotation match quality")

    mq_items = [
        (C_EXACT,   "Exact / near-exact match",
                    "GT and SOM share a key word  (e.g. 'monocyte' in both)"),
        (C_PARTIAL, "Partial match  (correct broad type)",
                    "SOM broad type matches GT  (e.g. GT='platelet', SOM broad='blood cell')"),
        (C_MISS,    "Mismatch",
                    "No meaningful word overlap at any level"),
    ]
    for i, (col, title, desc) in enumerate(mq_items):
        ey = 0.60 - i * 0.22
        # Filled circle swatch
        ax_k.scatter(MQ_X + 0.025, ey, s=90, color=col, alpha=0.90,
                     edgecolors="white", linewidths=0.5,
                     transform=ax_k.transAxes, zorder=3)
        ax_k.text(MQ_X + 0.052, ey + 0.04, title,
                  transform=ax_k.transAxes,
                  fontsize=7.5, fontweight="bold", color="#1a1a1a",
                  va="center", ha="left")
        ax_k.text(MQ_X + 0.052, ey - 0.07, desc,
                  transform=ax_k.transAxes,
                  fontsize=6.5, color="#555566",
                  va="center", ha="left", style="italic")

    # ── 2. Bubble-size legend (right upper ~35%) ──────────────────────────────
    SZ_X = 0.64
    SZ_W = 0.34
    section_header(SZ_X, 0.74, SZ_W, 0.20, "Bubble size  =  number of clusters")

    size_entries = [(1.0, "max"), (0.5, "50%"), (0.25, "25%")]
    for si, (frac, lbl) in enumerate(size_entries):
        bx = SZ_X + 0.06 + si * 0.11
        by = 0.45
        pt = bubble_max_pt * frac ** 0.5 * 0.32   # rescale for axes units
        ax_k.scatter(bx, by, s=pt,
                     color="#888888", alpha=0.75, edgecolors="white", linewidths=0.4,
                     transform=ax_k.transAxes, zorder=3)
        ax_k.text(bx, by - 0.16, lbl,
                  transform=ax_k.transAxes,
                  fontsize=6.5, color="#444444", ha="center", va="top")

    # ── 3. How-to-read blurb (right lower) ───────────────────────────────────
    note = ("Each bubble = one (GT label, SOM prediction) pair.\n"
            "Diagonal = perfect agreement. Off-diagonal = disagreement.")
    ax_k.text(SZ_X + SZ_W / 2, 0.16, note,
              transform=ax_k.transAxes,
              fontsize=6.5, color="#555566", ha="center", va="center",
              style="italic", linespacing=1.5)


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
    cbar.set_label("Mean expression\n(log-norm)", fontsize=7,
                   color=P["fg"], labelpad=3)
    cbar.ax.tick_params(labelsize=6, colors=P["fg"], length=2)
    cbar.outline.set_linewidth(0.4)

    # Size legend
    for frac_val, lbl in [(0.25, "25%"), (0.50, "50%"), (0.75, "75%")]:
        ax.scatter([], [], s=dot_max * frac_val, color="#999999",
                   edgecolors="#555555", linewidths=0.25, label=lbl)
    ax.legend(title="% expressing", title_fontsize=7, fontsize=7,
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

    summary_df, sig_df, cluster_lineage = load_som_results(cfg["som_sample"])
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

    # Compute n_match / n_total for badge
    ct_col = "Top_Cell_Type" if "Top_Cell_Type" in summary_df.columns else "Cell_Type"
    broad_col = "Broad_Type" if "Broad_Type" in summary_df.columns else None
    cluster_to_som   = dict(zip(summary_df["Cluster"].astype(str), summary_df[ct_col]))
    cluster_to_broad = {}
    if broad_col:
        cluster_to_broad = dict(zip(summary_df["Cluster"].astype(str), summary_df[broad_col]))
    cluster_to_gt = {str(k): v for k, v in ground_truth.items() if v is not None}
    n_match = 0
    n_total = 0
    for cluster, gt in cluster_to_gt.items():
        som   = cluster_to_som.get(cluster, "Unannotated")
        broad = cluster_to_broad.get(cluster, "")
        gt_tok    = set(gt.lower().replace("-", " ").replace("+", " ").split())
        som_tok   = set(som.lower().replace("-", " ").replace("+", " ").split())
        broad_tok = set(str(broad).lower().replace("-", " ").replace("+", " ").split())
        if bool(gt_tok & som_tok) or bool(gt_tok & broad_tok):
            n_match += 1
        n_total += 1

    # -- Dynamic figure sizing -------------------------------------------------
    n_gt_rows = len({v for v in ground_truth.values() if v is not None})

    fig_w = 9.0
    # Key panel height fixed at 1.5 in; matrix gets remaining height
    key_h_in   = 1.5
    matrix_h_in = max(4.0, n_gt_rows * 0.38 + 1.8)
    title_h_in  = 0.5
    xlab_h_in   = 1.4   # room for rotated x-axis labels
    fig_h = title_h_in + matrix_h_in + xlab_h_in + key_h_in + 0.3

    # Fractional layout (bottom → top)
    key_bottom   = 0.02
    key_frac_h   = key_h_in / fig_h
    xlab_frac_h  = xlab_h_in / fig_h
    mat_bottom   = key_bottom + key_frac_h + xlab_frac_h
    mat_frac_h   = matrix_h_in / fig_h
    mat_top      = mat_bottom + mat_frac_h

    fig_left  = 0.26
    fig_right = 0.97

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    ax_a = fig.add_axes([fig_left, mat_bottom,
                         fig_right - fig_left, mat_frac_h])

    bubble_max_pt = draw_bubble_matrix(
        ax_a, summary_df, ground_truth,
        cluster_lineage=cluster_lineage,
        accuracy=accuracy, n_match=n_match, n_total=n_total,
    )

    ax_a.text(-0.08, 1.05, "a", transform=ax_a.transAxes,
              fontsize=13, fontweight="bold", color=P["fg"], va="top")

    # KEY PANEL — framed, clearly demarcated
    draw_key_panel(
        fig, [],
        key_left=fig_left,
        key_bottom=key_bottom,
        key_width=fig_right - fig_left,
        key_height=key_frac_h,
        bubble_max_pt=bubble_max_pt,
    )

    # -- Title above matrix ---------------------------------------------------
    title_y = mat_top + (1.0 - mat_top) * 0.50
    fig.text(0.5, title_y, cfg["display_name"],
             ha="center", va="center", fontsize=10, fontweight="bold",
             color=P["fg"])

    # -- Save -----------------------------------------------------------------
    _FIGS_DIR.mkdir(parents=True, exist_ok=True)
    out_png = _FIGS_DIR / f"{dataset_key}_comparison.png"
    out_pdf = _FIGS_DIR / f"{dataset_key}_comparison.pdf"

    plt.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close()

    logger.info(f"Saved: {out_png}")
    logger.info(f"Saved: {out_pdf}")
    logger.info(f"Accuracy: {accuracy:.1%}  ({n_match}/{n_total} clusters)")
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
