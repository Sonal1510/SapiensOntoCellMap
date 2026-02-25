#!/usr/bin/env python3
"""
SapiensOntoCellMap Benchmarking Suite — PBMC3k
===============================================
Author : Sonal Rashmi
Date   : 2026-02-24

Benchmarks SapiensOntoCellMap against CellTypist on the Zheng et al. 2017
PBMC3k dataset (2,638 PBMCs, 8 cell types). This is the standard scRNA-seq
annotation benchmark used by every published cell type annotation paper.

Comparators
-----------
  • SapiensOntoCellMap  — hypergeometric enrichment, 14+ databases, CL hierarchy
  • CellTypist          — logistic regression, Human Cell Atlas reference
                          (Dominguez Conde et al. Science 2022)
  • scType              — marker-score based, curated DB (Ianevski et al. Nat Commun 2022)
                          NOTE: scType requires R; results must be pre-loaded manually
                          or via rpy2 (see MANUAL_COMPARATORS dict below)

Metrics
-------
  • Top-1 accuracy        — top predicted cell type matches ground truth
  • CL-level exact match  — predicted CL ID == ground truth CL ID
  • Hierarchical accuracy — predicted CL is ancestor/descendant of ground truth CL
  • Broad-type accuracy   — Broad_Type matches top-level lineage
  • Concordance           — SapiensOntoCellMap ↔ CellTypist agreement rate

Usage
-----
  # Full pipeline (requires: pip install scanpy celltypist)
  python benchmarking/benchmark_pbmc3k.py

  # Use pre-computed Seurat FindAllMarkers CSV
  python benchmarking/benchmark_pbmc3k.py --deg_csv path/to/markers.csv --deg_format seurat

  # Skip CellTypist
  python benchmarking/benchmark_pbmc3k.py --no_celltypist

  # Skip scanpy DEG step (DEG CSV already produced)
  python benchmarking/benchmark_pbmc3k.py --deg_csv benchmarking/results/pbmc3k_degs.csv

  # Use different background N (default 20,000 for whole-transcriptome scRNA-seq)
  python benchmarking/benchmark_pbmc3k.py --background_n 20000
"""

import argparse
import logging
import os
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths (auto-resolved relative to this file)
# ---------------------------------------------------------------------------
BENCH_DIR  = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(BENCH_DIR)
RESULTS_DIR = os.path.join(BENCH_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Config paths
sys.path.insert(0, PROJECT_DIR)
sys.path.insert(0, os.path.join(PROJECT_DIR, "config"))
try:
    from config.config import (
        PROCESSED_COMBINED_DATABASE_FILE,
        HGNC_COMPLETE_SET_FILE,
    )
    MARKER_DB = PROCESSED_COMBINED_DATABASE_FILE
    HGNC_FILE = HGNC_COMPLETE_SET_FILE
except ImportError:
    MARKER_DB = os.path.join(PROJECT_DIR, "data", "processed_combined_db", "master_cell_marker_db.csv")
    HGNC_FILE = os.path.join(PROJECT_DIR, "data", "reference", "hgnc_complete_set.txt")

ANNOTATION_SCRIPT = os.path.join(
    PROJECT_DIR, "src", "cluster_annotation", "get_cluster_annotation.py"
)

PBMC_H5AD_URL = (
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/"
    "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)

# ---------------------------------------------------------------------------
# Ground truth: PBMC3k Seurat tutorial cell type labels
# These ARE the louvain cluster labels returned by sc.datasets.pbmc3k_processed().
# Keys are exact louvain label strings; values are (display_name, CL_ID, broad_lineage).
# Source: Seurat PBMC3k tutorial (Hao et al. 2021 Nat Biotech)
# ---------------------------------------------------------------------------
PBMC3K_GROUND_TRUTH = {
    "CD4 T cells":         ("CD4 T cells",      "CL:0000624", "T cell"),
    "CD14+ Monocytes":     ("CD14+ Monocytes",  "CL:0000860", "Monocyte"),
    "B cells":             ("B cells",           "CL:0000236", "B cell"),
    "CD8 T cells":         ("CD8 T cells",       "CL:0000625", "T cell"),
    "NK cells":            ("NK cells",          "CL:0000623", "NK cell"),
    "FCGR3A+ Monocytes":   ("FCGR3A+ Monocytes","CL:0000875", "Monocyte"),
    "Dendritic cells":     ("Dendritic cells",   "CL:0000451", "Dendritic cell"),
    "Megakaryocytes":      ("Megakaryocytes",    "CL:0000556", "Megakaryocyte"),
}

# Acceptable fuzzy matches for Top-1 string accuracy
# Keys are ground truth cell type names (lowercase); values are acceptable substrings
GT_ALIASES = {
    # Top-1 cluster label → acceptable substrings in the predicted cell type name
    "cd4 t cells":         ["cd4", "helper t", "helper cd4"],       # "CD4-positive" → matches "cd4"
    "naive cd4 t cells":   ["cd4", "t cell", "helper", "naive"],
    "memory cd4 t cells":  ["cd4", "t cell", "helper", "memory"],
    "cd14+ monocytes":     ["monocyte", "cd14", "classical monocyte"],
    "b cells":             ["b cell", "b-cell", "b-2 b"],
    "cd8 t cells":         ["cd8", "cytotoxic", "killer"],
    "fcgr3a+ monocytes":   ["monocyte", "non-classical", "fcgr3a", "cd16"],
    "nk cells":            ["natural killer", "nk cell"],
    "dendritic cells":     ["dendritic", "dc"],
    "megakaryocytes":      ["megakaryocyte", "platelet"],
    # Lineage aliases used by hierarchical_match via _string_match(predicted, gt_lineage)
    "megakaryocyte":       ["megakaryocyte", "platelet"],           # platelet is a megakaryocyte product
}

# CL ancestor terms for broad-lineage matching
CL_LINEAGE_MAP = {
    "T cell":          {"CL:0000084", "CL:0000624", "CL:0000625", "CL:0000798",
                        "CL:0000909", "CL:0000910"},
    "B cell":          {"CL:0000236", "CL:0000785", "CL:0000946"},
    "NK cell":         {"CL:0000623", "CL:0000825"},
    "Monocyte":        {"CL:0000576", "CL:0000860", "CL:0000875"},
    "Dendritic cell":  {"CL:0000451", "CL:0000990"},
    "Megakaryocyte":   {"CL:0000556"},
}

# Optional: manually fill in scType results here (cluster → cell_type_string)
# to add scType as a comparator without running R
MANUAL_COMPARATORS = {
    # "scType": {
    #     "0": "CD4+ T cells",
    #     "1": "CD14+ Monocytes",
    #     ...
    # }
}


# ===========================================================================
# STEP 1: Download and preprocess PBMC3k
# ===========================================================================

def get_pbmc3k_adata(output_dir: str):
    """
    Load PBMC3k processed dataset via scanpy.datasets.pbmc3k_processed().

    This returns the dataset as processed in the Seurat tutorial with pre-annotated
    louvain cluster labels (e.g. 'CD4 T cells', 'CD14+ Monocytes', etc.).
    Using these cluster labels as both grouping key and ground truth avoids the
    cluster-ID mapping problem that arises when running leiden fresh.

    Requires: pip install scanpy
    """
    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy not installed. Run: pip install scanpy")
        logger.error("OR provide pre-computed DEG CSV via --deg_csv")
        sys.exit(1)

    h5ad_path = os.path.join(output_dir, "pbmc3k_processed.h5ad")
    if os.path.exists(h5ad_path):
        logger.info(f"Loading cached PBMC3k from {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        return adata

    logger.info("Downloading PBMC3k processed (Seurat tutorial) via scanpy...")
    # pbmc3k_processed() returns the fully preprocessed PBMC3k with Seurat tutorial
    # louvain labels in adata.obs['louvain'] as cell type strings.
    adata = sc.datasets.pbmc3k_processed()

    # adata.X is scaled data; we need log-normalized for DEG fold changes + CellTypist.
    # Re-derive from raw counts: fetch the raw unprocessed data and preprocess properly.
    logger.info("Fetching raw PBMC3k to rebuild log-normalized layer...")
    adata_raw = sc.datasets.pbmc3k()
    sc.pp.filter_cells(adata_raw, min_genes=200)
    sc.pp.filter_genes(adata_raw, min_cells=3)
    adata_raw.var["mt"] = adata_raw.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata_raw = adata_raw[adata_raw.obs.pct_counts_mt < 5].copy()
    sc.pp.normalize_total(adata_raw, target_sum=1e4)
    sc.pp.log1p(adata_raw)

    # Match cells between processed (has louvain) and raw-normalized (has all genes)
    common_cells = adata.obs_names.intersection(adata_raw.obs_names)
    adata_matched = adata_raw[common_cells].copy()
    # Transfer louvain labels
    adata_matched.obs["louvain"] = adata.obs.loc[common_cells, "louvain"]
    # Save log-normalized full-gene-space as raw for DEG + CellTypist
    adata_matched.raw = adata_matched.copy()

    # Also run HVG + scale + PCA + UMAP for structure (optional, for viz)
    sc.pp.highly_variable_genes(adata_matched, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_hvg = adata_matched[:, adata_matched.var.highly_variable].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack")
    sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_hvg)

    # Transfer UMAP back to adata_matched
    adata_matched.obsm["X_umap"] = adata_hvg.obsm["X_umap"]

    adata_matched.write_h5ad(h5ad_path)
    logger.info(f"Saved preprocessed PBMC3k to {h5ad_path}")
    logger.info(f"  {adata_matched.n_obs} cells, {adata_matched.n_vars} genes, "
                f"{adata_matched.obs['louvain'].nunique()} cell type clusters")
    return adata_matched


# ===========================================================================
# STEP 2: Compute DEGs and export as Scanpy-format CSV
# ===========================================================================

def compute_and_export_degs(adata, output_dir: str) -> str:
    """
    Run Wilcoxon rank-sum test per louvain cluster and export to
    Scanpy rank_genes_groups CSV format (names, group, pvals_adj, logfoldchanges).

    Returns path to the saved DEG CSV.
    """
    try:
        import scanpy as sc
    except ImportError:
        logger.error("scanpy required for DEG computation")
        sys.exit(1)

    deg_path = os.path.join(output_dir, "pbmc3k_degs.csv")
    if os.path.exists(deg_path):
        logger.info(f"Using cached DEG file: {deg_path}")
        return deg_path

    logger.info("Computing Wilcoxon DEGs per cluster (one-vs-rest)...")
    # use_raw=True: computes fold changes on the log-normalized full gene space
    # (saved in adata.raw) — avoids NaN logFC from scaled data
    sc.tl.rank_genes_groups(
        adata,
        groupby="louvain",
        method="wilcoxon",
        n_genes=adata.raw.var.shape[0] if adata.raw is not None else adata.n_vars,
        use_raw=True,
    )

    # Export in Scanpy CSV format (columns: names, group, pvals_adj, logfoldchanges)
    records = []
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    for grp in groups:
        names      = adata.uns["rank_genes_groups"]["names"][grp]
        pvals_adj  = adata.uns["rank_genes_groups"]["pvals_adj"][grp]
        logfc      = adata.uns["rank_genes_groups"]["logfoldchanges"][grp]
        for g, p, lfc in zip(names, pvals_adj, logfc):
            records.append({
                "names":          g,
                "group":          grp,
                "pvals_adj":      float(p),
                "logfoldchanges": float(lfc),
            })

    deg_df = pd.DataFrame(records)
    deg_df.to_csv(deg_path, index=False)
    logger.info(f"DEGs exported: {len(deg_df):,} rows → {deg_path}")
    return deg_path


# ===========================================================================
# STEP 3: Run SapiensOntoCellMap via subprocess
# ===========================================================================

def run_sapiensontocellmap(deg_csv: str, output_dir: str,
                           deg_format: str = "scanpy",
                           background_n: int = 20000,
                           tissue: str = None) -> str:
    """
    Call get_cluster_annotation.py via subprocess.
    Returns path to the top_annotation_summary.csv produced.
    """
    sample_name = "pbmc3k_benchmark"
    job_out_dir = os.path.join(output_dir, "sapiensonto_out")
    os.makedirs(job_out_dir, exist_ok=True)

    # Resolve to absolute path — subprocess CWD differs from the caller's CWD
    deg_csv_abs = os.path.abspath(deg_csv)

    cmd = [
        sys.executable, ANNOTATION_SCRIPT,
        deg_csv_abs, sample_name, job_out_dir,
        "--deg_type", "scrna",
        "--marker_db", MARKER_DB,
        "--deg_format", deg_format,
        "--background_gene_count", str(background_n),
        "--log2fc", "0.25",   # lower threshold for PBMC (small panel)
        "--pval", "0.05",
        "--min_overlap", "2",
    ]
    if HGNC_FILE and os.path.exists(HGNC_FILE):
        cmd += ["--hgnc_map", HGNC_FILE]
    if tissue:
        cmd += ["--tissue", tissue]

    logger.info(f"Running SapiensOntoCellMap: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=os.path.join(PROJECT_DIR, "src", "cluster_annotation"))
    if result.returncode != 0:
        logger.warning("SapiensOntoCellMap exited non-zero (may be HTML-only error):")
        logger.warning(result.stderr[-2000:])

    # Find top_annotation_summary.csv
    summary_csv = None
    for root, dirs, files in os.walk(job_out_dir):
        for f in files:
            if f.endswith("_top_annotation_summary.csv"):
                summary_csv = os.path.join(root, f)
                break

    if not summary_csv or not os.path.exists(summary_csv):
        logger.error(f"top_annotation_summary.csv not found under {job_out_dir}")
        sys.exit(1)

    logger.info(f"SapiensOntoCellMap results: {summary_csv}")
    return summary_csv


# ===========================================================================
# STEP 4: Run CellTypist (optional)
# ===========================================================================

def run_celltypist(adata, output_dir: str) -> pd.DataFrame:
    """
    Run CellTypist on the AnnData object.
    Requires: pip install celltypist

    Returns DataFrame with columns [cluster, celltypist_label].
    """
    try:
        import celltypist
        from celltypist import models as ct_models
    except ImportError:
        logger.warning("celltypist not installed. Skipping. Run: pip install celltypist")
        return pd.DataFrame()

    logger.info("Running CellTypist (Immune_All_High model)...")
    ct_models.download_models(force_update=False)

    # CellTypist requires log1p-normalised counts (NOT scaled).
    # adata.raw was saved before HVG subsetting and scaling, so it has
    # the full gene space in log-normalised form — exactly what CellTypist needs.
    import scanpy as sc
    import anndata
    if adata.raw is not None:
        # Convert raw back to AnnData with full gene space
        ct_adata = anndata.AnnData(
            X=adata.raw.X,
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )
    else:
        ct_adata = adata.copy()
        sc.pp.normalize_total(ct_adata, target_sum=1e4)
        sc.pp.log1p(ct_adata)

    predictions = celltypist.annotate(
        ct_adata,
        model="Immune_All_High.pkl",
        majority_voting=True,
        over_clustering="louvain",   # adata.obs["louvain"] set from leiden above
    )
    ct_adata = predictions.to_adata()

    # Aggregate to cluster level by majority vote
    # louvain contains cell type strings (e.g. "CD4 T cells") — use them directly as cluster key
    cluster_labels = (
        ct_adata.obs.groupby("louvain")["majority_voting"]
        .agg(lambda x: x.value_counts().index[0])
        .reset_index()
        .rename(columns={"louvain": "cluster", "majority_voting": "celltypist_label"})
    )
    cluster_labels["cluster"] = cluster_labels["cluster"].astype(str).str.strip()

    out_path = os.path.join(output_dir, "celltypist_results.csv")
    cluster_labels.to_csv(out_path, index=False)
    logger.info(f"CellTypist results: {out_path}")
    return cluster_labels


# ===========================================================================
# STEP 5: Compute metrics
# ===========================================================================

def _string_match(predicted: str, ground_truth_key: str) -> bool:
    """Fuzzy string match: any ground truth alias substring in predicted (case-insensitive)."""
    if not isinstance(predicted, str) or not predicted:
        return False
    pred_lower = predicted.lower()
    aliases = GT_ALIASES.get(ground_truth_key.lower(), [ground_truth_key.lower()])
    return any(alias in pred_lower for alias in aliases)


def _cl_exact_match(predicted_cl: str, gt_cl: str) -> bool:
    if not isinstance(predicted_cl, str):
        return False
    return predicted_cl.strip() == gt_cl.strip()


def _cl_hierarchical_match(predicted_cell_type: str, gt_lineage: str,
                            cell_name_to_cl: dict, cl_lineage_map: dict) -> bool:
    """
    Check if the predicted cell type belongs to the correct broad lineage
    by checking if its CL ID (looked up via cell_name_to_cl) is in the
    known CL ID set for that lineage.
    """
    if not isinstance(predicted_cell_type, str):
        return False
    cl_id = cell_name_to_cl.get(predicted_cell_type.upper())
    if not cl_id:
        return False
    lineage_cls = cl_lineage_map.get(gt_lineage, set())
    return cl_id in lineage_cls


def _broad_type_match(broad_type: str, gt_lineage: str) -> bool:
    """Check if Broad_Type column (from hierarchical annotation) contains the correct lineage."""
    if not isinstance(broad_type, str) or not broad_type:
        return False
    gt_lower = gt_lineage.lower()
    return any(tok in broad_type.lower() for tok in gt_lower.split())


def compute_metrics(
    sapiensonto_csv: str,
    gt_dict: dict,
    cell_name_to_cl: dict = None,
    celltypist_df: pd.DataFrame = None,
    manual_comparators: dict = None,
) -> dict:
    """
    Compute all benchmark metrics.

    Parameters
    ----------
    sapiensonto_csv : path to top_annotation_summary.csv
    gt_dict         : {cluster_str: (gt_cell_type, gt_cl_id, gt_lineage)}
    cell_name_to_cl : {cell_type_name_upper: cl_id} from master DB
    celltypist_df   : DataFrame with [cluster, celltypist_label] or empty
    manual_comparators : dict of {tool_name: {cluster: cell_type}}

    Returns
    -------
    dict with per-cluster results and aggregate metrics
    """
    df = pd.read_csv(sapiensonto_csv, dtype=str)
    # SapiensOntoCellMap prefixes cluster names with "Cluster " (e.g. "Cluster B cells")
    df["_cluster_key"] = df["Cluster"].str.replace("^Cluster ", "", regex=True).str.strip()

    results = []
    for cluster_id, (gt_name, gt_cl, gt_lineage) in gt_dict.items():
        # Match by cluster label after stripping the "Cluster " prefix
        row = df[df["_cluster_key"] == str(cluster_id)]
        if row.empty:
            logger.warning(f"Cluster {cluster_id} not found in SapiensOntoCellMap output")
            results.append({
                "cluster": cluster_id,
                "gt_cell_type": gt_name,
                "gt_cl_id": gt_cl,
                "gt_lineage": gt_lineage,
                "predicted": None,
                "predicted_cl": None,
                "broad_type": None,
                "confidence": None,
                "top1_match": False,
                "cl_exact_match": False,
                "hierarchical_match": False,
                "broad_type_match": False,
            })
            continue

        row = row.iloc[0]
        predicted     = str(row.get("Top_Cell_Type", "")).strip()
        predicted_cl  = str(row.get("Source", "")).strip()   # Source holds CL ID in some formats
        broad_type    = str(row.get("Broad_Type", "")).strip()
        confidence    = row.get("Confidence", None)

        # Top-1 string match
        top1 = _string_match(predicted, gt_name)

        # CL exact match (Source may be "cell_name [CL:...]" — extract CL ID if present)
        cl_match = False
        import re
        cl_hits = re.findall(r"CL:\d+", predicted_cl)
        if cl_hits:
            cl_match = _cl_exact_match(cl_hits[0], gt_cl)

        # Hierarchical lineage match
        hier_match = _broad_type_match(broad_type, gt_lineage) or \
                     _string_match(predicted, gt_lineage)

        # Broad type match
        bt_match = _broad_type_match(broad_type, gt_lineage)

        results.append({
            "cluster":            cluster_id,
            "gt_cell_type":       gt_name,
            "gt_cl_id":           gt_cl,
            "gt_lineage":         gt_lineage,
            "predicted":          predicted,
            "predicted_cl":       predicted_cl,
            "broad_type":         broad_type,
            "confidence":         confidence,
            "top1_match":         top1,
            "cl_exact_match":     cl_match,
            "hierarchical_match": hier_match,
            "broad_type_match":   bt_match,
        })

    results_df = pd.DataFrame(results)

    n = len(results_df)
    aggregate = {
        "n_clusters":           n,
        "top1_accuracy":        results_df["top1_match"].sum() / n,
        "cl_exact_accuracy":    results_df["cl_exact_match"].sum() / n,
        "hierarchical_accuracy":results_df["hierarchical_match"].sum() / n,
        "broad_type_accuracy":  results_df["broad_type_match"].sum() / n,
    }

    # CellTypist concordance
    if celltypist_df is not None and not celltypist_df.empty:
        ct_map = dict(zip(celltypist_df["cluster"].astype(str), celltypist_df["celltypist_label"]))
        results_df["celltypist_label"] = results_df["cluster"].astype(str).map(ct_map)
        results_df["celltypist_top1"]  = results_df.apply(
            lambda r: _string_match(str(r["celltypist_label"]), r["gt_cell_type"])
            if pd.notna(r["celltypist_label"]) else False,
            axis=1
        )
        results_df["concordance"] = (
            results_df["top1_match"] == results_df["celltypist_top1"]
        )
        aggregate["celltypist_top1_accuracy"] = results_df["celltypist_top1"].sum() / n
        aggregate["sapiensonto_celltypist_concordance"] = results_df["concordance"].mean()

    # Manual comparators (e.g., scType)
    for tool_name, pred_map in (manual_comparators or {}).items():
        results_df[f"{tool_name}_label"] = results_df["cluster"].astype(str).map(pred_map)
        results_df[f"{tool_name}_top1"]  = results_df.apply(
            lambda r: _string_match(str(r[f"{tool_name}_label"]), r["gt_cell_type"])
            if pd.notna(r[f"{tool_name}_label"]) else False,
            axis=1
        )
        aggregate[f"{tool_name}_top1_accuracy"] = results_df[f"{tool_name}_top1"].sum() / n

    return {"per_cluster": results_df, "aggregate": aggregate}


# ===========================================================================
# STEP 6: Generate comparison figures
# ===========================================================================

def generate_figures(metrics: dict, output_dir: str):
    """
    Produce benchmark figures:
      1. Per-cluster accuracy heatmap (all metrics × all clusters)
      2. Tool comparison bar chart (top-1 accuracy per tool)
      3. Confidence vs. accuracy scatter
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    per_df = metrics["per_cluster"].copy()
    agg    = metrics["aggregate"]

    # ---- Figure 1: Per-cluster metric heatmap ----------------------------------
    metric_cols = ["top1_match", "cl_exact_match", "hierarchical_match", "broad_type_match"]
    labels = ["Top-1 String Match", "CL Exact Match", "Hierarchical Match", "Broad Type Match"]
    heatmap_data = per_df.set_index("gt_cell_type")[metric_cols].astype(int)
    heatmap_data.columns = labels

    fig, ax = plt.subplots(figsize=(10, max(4, len(per_df) * 0.6)))
    sns.heatmap(
        heatmap_data, annot=True, fmt="d", cmap="RdYlGn",
        linewidths=0.5, cbar=False, ax=ax,
        vmin=0, vmax=1,
    )
    ax.set_title("SapiensOntoCellMap — Per-Cluster Benchmark (PBMC3k)", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("Ground Truth Cell Type", fontsize=11)
    plt.tight_layout()
    fig1_path = os.path.join(output_dir, "benchmark_per_cluster_heatmap.png")
    fig.savefig(fig1_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 1 saved: {fig1_path}")

    # ---- Figure 2: Tool comparison bar chart -----------------------------------
    tool_metrics = {
        "SapiensOntoCellMap\n(Top-1)":       agg["top1_accuracy"],
        "SapiensOntoCellMap\n(Hierarchical)": agg["hierarchical_accuracy"],
        "SapiensOntoCellMap\n(Broad Type)":   agg["broad_type_accuracy"],
    }
    if "celltypist_top1_accuracy" in agg:
        tool_metrics["CellTypist\n(Top-1)"] = agg["celltypist_top1_accuracy"]
    for key in agg:
        if key.endswith("_top1_accuracy") and not key.startswith("celltypist"):
            tool_name = key.replace("_top1_accuracy", "")
            tool_metrics[f"{tool_name}\n(Top-1)"] = agg[key]

    fig, ax = plt.subplots(figsize=(max(6, len(tool_metrics) * 1.4), 5))
    colors = ["#2c7bb6" if "SapiensOnto" in k else "#d7191c" if "CellTypist" in k
              else "#fdae61" for k in tool_metrics]
    bars = ax.bar(list(tool_metrics.keys()), list(tool_metrics.values()),
                  color=colors, edgecolor="white", width=0.6)
    ax.bar_label(bars, fmt="{:.1%}", padding=3, fontsize=10)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Accuracy", fontsize=11)
    ax.set_title("Tool Comparison — PBMC3k (n=9 clusters)", fontsize=13)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig2_path = os.path.join(output_dir, "benchmark_tool_comparison.png")
    fig.savefig(fig2_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 2 saved: {fig2_path}")

    # ---- Figure 3: Confidence vs. accuracy scatter ---------------------------
    per_df["accuracy"] = per_df["top1_match"].astype(int)
    conf_vals = pd.to_numeric(per_df["confidence"], errors="coerce")
    if conf_vals.notna().sum() >= 3:
        fig, ax = plt.subplots(figsize=(6, 5))
        scatter_df = per_df.dropna(subset=["confidence"])
        scatter_df = scatter_df.copy()
        scatter_df["confidence"] = pd.to_numeric(scatter_df["confidence"])
        colors_pt = scatter_df["accuracy"].map({1: "#2c7bb6", 0: "#d7191c"})
        ax.scatter(scatter_df["confidence"], scatter_df["accuracy"] + 0.05 * np.random.randn(len(scatter_df)),
                   c=colors_pt, s=90, alpha=0.8, edgecolors="white")
        rho, pval = spearmanr(scatter_df["confidence"], scatter_df["accuracy"])
        ax.set_xlabel("Confidence Score (HierarchicalAnnotator)", fontsize=11)
        ax.set_ylabel("Correct Annotation (1=Yes, 0=No)", fontsize=11)
        ax.set_title(f"Confidence vs. Accuracy (ρ={rho:.2f}, p={pval:.3f})", fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig3_path = os.path.join(output_dir, "benchmark_confidence_vs_accuracy.png")
        fig.savefig(fig3_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Figure 3 saved: {fig3_path}")


# ===========================================================================
# STEP 7: Print and save summary
# ===========================================================================

def print_summary(metrics: dict, output_dir: str):
    agg = metrics["aggregate"]
    per = metrics["per_cluster"]

    print("\n" + "=" * 70)
    print("  SapiensOntoCellMap — PBMC3k Benchmark Summary")
    print("=" * 70)
    print(f"  Clusters evaluated : {agg['n_clusters']}")
    print(f"  Top-1 accuracy     : {agg['top1_accuracy']:.1%}")
    print(f"  CL exact accuracy  : {agg['cl_exact_accuracy']:.1%}")
    print(f"  Hierarchical acc.  : {agg['hierarchical_accuracy']:.1%}")
    print(f"  Broad-type acc.    : {agg['broad_type_accuracy']:.1%}")
    if "celltypist_top1_accuracy" in agg:
        print(f"  CellTypist Top-1   : {agg['celltypist_top1_accuracy']:.1%}")
        print(f"  Concordance (SOCM↔CT): {agg['sapiensonto_celltypist_concordance']:.1%}")
    for key, val in agg.items():
        if key.endswith("_top1_accuracy") and "celltypist" not in key:
            tool = key.replace("_top1_accuracy", "")
            print(f"  {tool} Top-1      : {val:.1%}")
    print("=" * 70)

    print("\nPer-cluster results:")
    display_cols = ["cluster", "gt_cell_type", "predicted", "broad_type",
                    "top1_match", "hierarchical_match", "confidence"]
    available = [c for c in display_cols if c in per.columns]
    print(per[available].to_string(index=False))

    # Save full results
    results_path = os.path.join(output_dir, "benchmark_results.csv")
    per.to_csv(results_path, index=False)

    summary_path = os.path.join(output_dir, "benchmark_summary.txt")
    with open(summary_path, "w") as f:
        f.write("SapiensOntoCellMap — PBMC3k Benchmark\n")
        f.write(f"Date: 2026-02-24\n\n")
        for k, v in agg.items():
            f.write(f"{k}: {v:.4f}\n" if isinstance(v, float) else f"{k}: {v}\n")

    logger.info(f"Full results: {results_path}")
    logger.info(f"Summary:      {summary_path}")


# ===========================================================================
# Main
# ===========================================================================

def parse_args():
    p = argparse.ArgumentParser(description="SapiensOntoCellMap PBMC3k Benchmark")
    p.add_argument("--deg_csv", type=str, default=None,
                   help="Pre-computed DEG CSV (skip scanpy step)")
    p.add_argument("--deg_format", choices=["scanpy", "seurat", "generic"],
                   default="scanpy",
                   help="DEG input format (default: scanpy)")
    p.add_argument("--background_n", type=int, default=20000,
                   help="Background gene count N for hypergeometric test (default: 20000)")
    p.add_argument("--tissue", type=str, default=None,
                   help="Filter marker DB to this tissue (e.g. 'blood'). "
                        "Default: all tissues (recommended for PBMC — broad DB)")
    p.add_argument("--no_celltypist", action="store_true",
                   help="Skip CellTypist comparison")
    p.add_argument("--results_dir", type=str, default=RESULTS_DIR,
                   help=f"Output directory (default: {RESULTS_DIR})")
    return p.parse_args()


def main():
    args = parse_args()
    output_dir = args.results_dir
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(MARKER_DB):
        logger.error(f"Marker DB not found: {MARKER_DB}")
        logger.error("Run test/test_classes.py to build the database first.")
        sys.exit(1)

    # ---- Step 1-2: Data prep -------------------------------------------------
    adata = None
    deg_csv = args.deg_csv

    if deg_csv and os.path.exists(deg_csv):
        logger.info(f"Using pre-computed DEG CSV: {deg_csv}")
    else:
        logger.info("Step 1: Downloading/loading PBMC3k...")
        adata = get_pbmc3k_adata(output_dir)
        logger.info(f"  PBMC3k: {adata.n_obs} cells, {adata.n_vars} genes, "
                    f"{adata.obs['louvain'].nunique()} cell type clusters: "
                    f"{sorted(adata.obs['louvain'].unique().tolist())}")

        logger.info("Step 2: Computing DEGs...")
        deg_csv = compute_and_export_degs(adata, output_dir)

    # ---- Step 3: SapiensOntoCellMap ------------------------------------------
    logger.info("Step 3: Running SapiensOntoCellMap...")
    summary_csv = run_sapiensontocellmap(
        deg_csv=deg_csv,
        output_dir=output_dir,
        deg_format=args.deg_format,
        background_n=args.background_n,
        tissue=args.tissue,
    )

    # ---- Step 4: CellTypist (optional) ---------------------------------------
    celltypist_df = pd.DataFrame()
    if not args.no_celltypist:
        if adata is None:
            logger.warning("AnnData not loaded (--deg_csv provided). "
                           "Re-loading PBMC3k for CellTypist...")
            try:
                adata = get_pbmc3k_adata(output_dir)
            except SystemExit:
                logger.warning("Skipping CellTypist (scanpy not available)")

        if adata is not None:
            logger.info("Step 4: Running CellTypist...")
            celltypist_df = run_celltypist(adata, output_dir)

    # ---- Step 5: Build cell_name → CL map from master DB --------------------
    logger.info("Step 5: Computing metrics...")
    try:
        db_df = pd.read_csv(MARKER_DB, usecols=["cell_name", "cell_id"], dtype=str)
        db_df = db_df.dropna().drop_duplicates()
        cell_name_to_cl = {
            str(r["cell_name"]).strip().upper(): str(r["cell_id"]).strip()
            for _, r in db_df.iterrows()
            if str(r["cell_id"]).startswith("CL:")
        }
    except Exception as e:
        logger.warning(f"Could not load cell_name→CL map: {e}")
        cell_name_to_cl = {}

    metrics = compute_metrics(
        sapiensonto_csv=summary_csv,
        gt_dict=PBMC3K_GROUND_TRUTH,
        cell_name_to_cl=cell_name_to_cl,
        celltypist_df=celltypist_df if not celltypist_df.empty else None,
        manual_comparators=MANUAL_COMPARATORS if MANUAL_COMPARATORS else None,
    )

    # ---- Step 6: Figures -----------------------------------------------------
    logger.info("Step 6: Generating figures...")
    generate_figures(metrics, output_dir)

    # ---- Step 7: Summary -----------------------------------------------------
    print_summary(metrics, output_dir)


if __name__ == "__main__":
    main()
