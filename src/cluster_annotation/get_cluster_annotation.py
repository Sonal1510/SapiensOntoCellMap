#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Annotate each cluster with cell types using:
1. Conventional Hypergeometric Test for 'db_cell_name' (Level 1)
2. Weighted Enrichment Analysis for 'cell_name' (Level 2)

(MODIFIED: Adds a Top Annotation Summary CSV output with logic:
 Selected Tissue (L2) > All Tissue (L2))
"""

import os
import sys
import logging
import argparse
import pandas as pd
import re
from get_marker_enrichment_test import MarkerEnrichmentTest
from get_html_report import (
    plot_dynamic_heatmap_with_bars,
    create_deg_tables_html,
    generate_html_report
)
from hierarchical_annotation import HierarchicalAnnotator

# --- Setup Logging ---
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# --- Load path defaults from config (project root two levels up from this file) ---
_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
try:
    from config.config import (
        PROCESSED_COMBINED_DATABASE_FILE as _DEFAULT_MARKER_DB,
        HGNC_COMPLETE_SET_FILE as _DEFAULT_HGNC_MAP,
        MSIGDB_G2M_FILE as _DEFAULT_G2M_FILE,
        MSIGDB_E2F_FILE as _DEFAULT_E2F_FILE,
    )
except ImportError:
    _DEFAULT_MARKER_DB = None
    _DEFAULT_HGNC_MAP = None
    _DEFAULT_G2M_FILE = None
    _DEFAULT_E2F_FILE = None

# Source-type quality weights for cross-database agreement scoring
SOURCE_TYPE_WEIGHTS = {
    'Experiment': 1.0,
    'Single-Cell Sequencing': 0.9,
    'Company': 0.8,
    'Literature': 0.7,
    'Review': 0.6,
    'Computational': 0.5,
}

# ---------------------------------------------------------------------------
# Proliferative / Cell-cycle gene signature
# ---------------------------------------------------------------------------
# Source: MSigDB Hallmark gene sets downloaded at pipeline build time via
#   download_reference_data() in src/download/bio_database_downloader.py.
#   • HALLMARK_G2M_CHECKPOINT — 200 genes (G2/M transition)
#   • HALLMARK_E2F_TARGETS    — 200 genes (E2F targets; canonical S-phase / cell-cycle
#     entry gene set in the MSigDB Hallmark collection — the Hallmark collection does
#     not contain a separate G1/S checkpoint gene set)
#   Liberzon et al., Cell Systems 2015 (doi:10.1016/j.cels.2015.12.004)
#   Subramanian et al., PNAS 2005 (doi:10.1073/pnas.0506580102)
#
# If the .grp files are absent, PROLIFERATIVE_GENES is an empty frozenset and
# Proliferative_Flag is set to None (not computed) in the output CSV.
#
# Biological rationale:
#   When a cluster's overlapping marker genes are dominated by cell-cycle genes
#   the enrichment test may match a specific cell type (e.g. "large pre-B-II
#   cell") only because that cell type's database entry contains genes shared
#   with the cell cycle programme.  In cancer or inflamed tissue, the true
#   identity may be a CYCLING or MALIGNANT population.  Flagging this allows
#   the user to correlate with Ki-67 IHC, pathology review, or copy-number
#   inference before accepting the annotation.
#
# Thresholds (_PROLIF_FLAG_*):
#   Flag if ≥3 cycle genes in the overlap  (absolute — clear proliferating cluster)
#   OR  ≥2 cycle genes AND they represent ≥33 % of the total overlap
#       (catches small Xenium panels where a 4-gene overlap of 2 MKI67+TOP2A = 50 %)
# ---------------------------------------------------------------------------

import re as _re
_GENE_SYMBOL_RE = _re.compile(r'^[A-Z][A-Z0-9\-\.]{0,31}$')


def _load_proliferative_genes(g2m_path=None, e2f_path=None):
    """
    Load the proliferative gene set from MSigDB Hallmark .grp files.

    Combines:
      • HALLMARK_G2M_CHECKPOINT (G2/M transition genes)
      • HALLMARK_E2F_TARGETS    (E2F target genes; canonical S-phase set in Hallmark)

    .grp format: plain text, one HGNC gene symbol per line; comment/header
    lines start with '#' or 'HALLMARK_'. Lines that do not match a valid HGNC
    gene symbol pattern (uppercase, 1-32 chars, letters/digits/hyphen/dot) are
    discarded with a warning — this guards against accidentally parsing HTML
    content that some MSigDB endpoints return for invalid gene set names.

    Returns an empty frozenset if files are missing or unreadable.
    Proliferative_Flag is set to None (not computed) in that case.

    Args:
        g2m_path: Path to HALLMARK_G2M_CHECKPOINT.grp (or None)
        e2f_path: Path to HALLMARK_E2F_TARGETS.grp (or None)

    Returns:
        frozenset of uppercase HGNC gene symbols, or frozenset() if unavailable
    """
    genes = set()
    for label, path in [('G2M (HALLMARK_G2M_CHECKPOINT)', g2m_path),
                        ('E2F (HALLMARK_E2F_TARGETS)', e2f_path)]:
        if path and os.path.exists(path):
            try:
                with open(path) as fh:
                    raw_lines = [line.strip() for line in fh]

                valid_genes = set()
                n_rejected = 0
                for s in raw_lines:
                    if not s or s.startswith('#') or s.startswith('HALLMARK_'):
                        continue
                    sym = s.upper()
                    if _GENE_SYMBOL_RE.match(sym):
                        valid_genes.add(sym)
                    else:
                        n_rejected += 1

                if n_rejected > 0:
                    logger.warning(
                        f"MSigDB {label}: discarded {n_rejected} lines that did not "
                        "match HGNC gene symbol format (possible HTML/redirect response)"
                    )

                if valid_genes:
                    genes.update(valid_genes)
                    logger.info(
                        f"Loaded {len(valid_genes)} genes from MSigDB {label} "
                        f"({os.path.basename(path)})"
                    )
                else:
                    logger.warning(
                        f"MSigDB {label} file yielded no valid gene symbols: {path}"
                    )
            except OSError as exc:
                logger.warning(f"MSigDB {label} file could not be read ({path}): {exc}")
        else:
            logger.warning(
                f"MSigDB {label} file not found: {path}. "
                "Run download_reference_data() to fetch it."
            )

    if not genes:
        logger.warning(
            "Proliferative gene sets not loaded — Proliferative_Flag will not be computed. "
            "Run the downloader (src/download/bio_database_downloader.py) to fetch "
            "HALLMARK_G2M_CHECKPOINT.grp and HALLMARK_E2F_TARGETS.grp from MSigDB."
        )

    return frozenset(genes)


PROLIFERATIVE_GENES = _load_proliferative_genes(_DEFAULT_G2M_FILE, _DEFAULT_E2F_FILE)
if PROLIFERATIVE_GENES:
    logger.info(f"Proliferative gene set loaded: {len(PROLIFERATIVE_GENES)} genes")

_PROLIF_FLAG_MIN_GENES = 3    # absolute: ≥3 cycle genes in overlap → flag
_PROLIF_FLAG_MIN_GENES2 = 2   # conditional: ≥2 AND high fraction → flag
_PROLIF_FLAG_MIN_FRAC = 0.33  # fraction threshold for conditional path

# Helper for natural sorting (cluster 2 before cluster 10)
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

def get_weighted_marker_maps(df, level_col):
    """
    Generates a WEIGHTED marker map (Level 2) using cross-database agreement
    and source-type quality weighting.

    Deduplicates to one entry per (cell_type, gene, database) to remove
    intra-database tissue redundancy. Each unique database contributes its
    source_type quality weight to the gene's total score.

    Returns:
        weighted_map:   { 'CellType': {'geneA': 2.3, 'geneB': 0.7} }
        n_databases_map: { 'CellType': {'geneA': 3, 'geneB': 1} }  — per-gene count
        db_names_map:   { 'CellType': {'geneA': {'DB1', 'DB2'}} }  — per-gene db sets
                         Used downstream to compute N_Databases as the union of databases
                         across all overlapping genes (i.e., how many independent studies
                         contributed at least one overlapping marker for this annotation).
    """
    df_clean = df.dropna(subset=[level_col, 'gene']).copy()
    if df_clean.empty:
        return {}, {}, {}

    # Deduplicate: one entry per (cell_type, gene, database)
    deduped = df_clean.drop_duplicates(subset=[level_col, 'gene', 'database']).copy()

    # Map source_type to quality weight
    deduped['_src_w'] = deduped['source_type'].map(SOURCE_TYPE_WEIGHTS).fillna(0.5)

    # Sum weights per (cell_type, gene) — each unique database contributes once
    weights = deduped.groupby([level_col, 'gene'])['_src_w'].sum()

    # Count unique databases per (cell_type, gene)
    db_counts = deduped.groupby([level_col, 'gene'])['database'].nunique()

    # Collect the actual set of database names per (cell_type, gene)
    # This is used downstream to compute N_Databases as the union of databases
    # across all overlapping genes — the correct measure of cross-study corroboration.
    db_sets = deduped.groupby([level_col, 'gene'])['database'].apply(set)

    weighted_map = {}
    n_databases_map = {}
    db_names_map = {}
    for (cell_type, gene), weight in weights.items():
        if cell_type not in weighted_map:
            weighted_map[cell_type] = {}
            n_databases_map[cell_type] = {}
            db_names_map[cell_type] = {}
        weighted_map[cell_type][gene] = round(float(weight), 2)
        n_databases_map[cell_type][gene] = int(db_counts.get((cell_type, gene), 1))
        db_names_map[cell_type][gene] = db_sets.get((cell_type, gene), set())

    return weighted_map, n_databases_map, db_names_map

def get_unweighted_marker_maps(df, level_col, distinct_col=None):
    """
    Generates a CONVENTIONAL (Unweighted) marker map (Level 1).
    Returns: { 'DB_CellType': ['geneA', 'geneB'] }
    """
    df_clean = df.dropna(subset=[level_col, 'gene']).copy()
    if df_clean.empty: return {}

    # Group and just get unique list of genes
    unweighted_map = df_clean.groupby(level_col)['gene'].apply(lambda x: list(set(x))).to_dict()
    return unweighted_map

def get_marker_maps_context(df):
    """
    Generates appropriate maps for Level 1 (Unweighted) and Level 2 (Weighted).
    Returns:
        maps: {'level1': {...}, 'level2': {...}}
        n_databases_maps: {'level2': {...}}  — per-gene database counts
        db_names_maps: {'level2': {...}}     — per-gene database name sets (for N_Databases)
    """
    maps = {}
    n_databases_maps = {}
    db_names_maps = {}

    # --- Level 1: Conventional (Unweighted) ---
    # Key: database_cell_name
    df_l1 = df.dropna(subset=['database', 'db_cell_name']).copy()
    if not df_l1.empty:
        df_l1['unique_db_cell_name'] = df_l1["database"] + "_" + df_l1["db_cell_name"]
        maps['level1'] = get_unweighted_marker_maps(df_l1, 'unique_db_cell_name')
    else:
        maps['level1'] = {}

    # --- Level 2: Weighted ---
    # Key: cell_name
    maps['level2'], n_databases_maps['level2'], db_names_maps['level2'] = \
        get_weighted_marker_maps(df, 'cell_name')

    return maps, n_databases_maps, db_names_maps

def get_spaceranger_differential_cluster_file_path_multi_path_key(path_list, sample_name):
    if not path_list: return {}
    split_paths = [p.split(os.sep) for p in path_list]
    common_prefix = os.path.commonpath(path_list)
    reversed_paths = [p.split(os.sep)[::-1] for p in path_list]
    common_suffix_parts = []
    for parts in zip(*reversed_paths):
        if len(set(parts)) == 1: common_suffix_parts.append(parts[0])
        else: break
    common_suffix = os.sep.join(common_suffix_parts[::-1])
    path_dict = {}
    for path in path_list:
        relative_to_prefix = os.path.relpath(path, common_prefix)
        middle_part = relative_to_prefix
        if middle_part.endswith(common_suffix):
            middle_part = middle_part[:-len(common_suffix)].strip(os.sep)
        variable_key_part = middle_part.replace(os.sep, '_')
        if not variable_key_part: variable_key_part = os.path.basename(os.path.dirname(path))
        final_key = f"{sample_name}_{variable_key_part}"
        path_dict[final_key] = path
    return path_dict

def get_spaceranger_differential_cluster_file_path(spaceranger_out_path, sample_name):
    path_map = {}
    if spaceranger_out_path.endswith("differential_expression.csv"):
        path_map[sample_name] = spaceranger_out_path
    elif os.path.isdir(spaceranger_out_path):
        paths = []
        for root, dirs, files in os.walk(spaceranger_out_path):
            if 'differential_expression.csv' in files and '_graphclust' in root:
                full_path = os.path.join(root, 'differential_expression.csv')
                paths.append(full_path)
        if paths:
            path_map = get_spaceranger_differential_cluster_file_path_multi_path_key(paths, sample_name)
        else:
            print(f"Warning: No '*_graphclust/differential_expression.csv' files found in {spaceranger_out_path}")
    else:
        print(f"Error: Path is not a valid directory ending in 'outs' or a direct CSV file path: {spaceranger_out_path}")
    return path_map

def find_deg_files(input_path, sample_name, deg_type):
    path_map = {}
    if deg_type == 'scrna':
        if os.path.isfile(input_path) and input_path.endswith(".csv"):
            path_map[sample_name] = input_path
        else:
            logger.error(f"For --deg_type scrna, the input path must be a valid .csv file. Got: {input_path}")
    else:
        if os.path.isdir(input_path):
            path_map = get_spaceranger_differential_cluster_file_path(input_path, sample_name)    
    if not path_map:
        logger.warning(f"No spatial differential_expression.csv files found in {input_path}")
    return path_map

def run_enrichment_pipeline(deg_file, marker_db_dict, args, deg_type, n_databases_map=None, db_names_map=None):
    if not marker_db_dict:
        return (None, pd.DataFrame(), pd.DataFrame(), {'heatmap': {}},
                {'p_val': args.pval, 'log2fc': args.log2fc, 'mean_counts': args.mean, 'top_n_genes': args.topgenes})

    pipeline = MarkerEnrichmentTest(
        deg_file=deg_file,
        marker_db_dict=marker_db_dict,
        deg_file_type=deg_type,
        p_val_thresh=args.pval,
        log2fc_thresh=args.log2fc,
        mean_counts_thresh=args.mean,
        top_genes=args.topgenes,
        min_overlap=args.min_overlap,
        min_db_markers=getattr(args, 'min_db_markers', 0),
        min_cluster_degs=getattr(args, 'min_cluster_degs', 0),
        background_gene_count=args.background_gene_count,
        hgnc_map=getattr(args, 'hgnc_map', None),
        deg_format=getattr(args, 'deg_format', None),
        auto_spatial_mean_filter=not getattr(args, 'no_auto_spatial_filter', False),
    )
    pipeline.fit()

    all_results = pipeline.results_
    sig_results = pd.DataFrame()
    if not all_results.empty:
        sig_results = all_results[all_results["adj_p_value"] < args.pval].copy()

    def _compute_n_databases(df, cell_type_col, genes_col, n_databases_map, db_names_map):
        """
        Compute N_Databases per row = number of independent databases/studies that
        contributed at least one of the overlapping marker genes for this annotation.

        Uses the union of database name sets across all overlapping genes, which is
        the biologically correct measure of cross-study corroboration. A cell type
        supported by overlapping genes from 3 different scRNA-seq studies scores 3,
        regardless of how many genes each study contributed.

        Falls back to per-gene max if db_names_map is unavailable.
        """
        n_db_values = []
        for _, row in df.iterrows():
            cell_type = row[cell_type_col]
            genes_str = row.get(genes_col, '')
            if genes_str and cell_type in n_databases_map:
                genes = [g.strip().upper() for g in str(genes_str).split(',')]
                if db_names_map and cell_type in db_names_map:
                    # Correct: count unique databases across all overlapping genes
                    overlap_dbs = set()
                    ct_db = db_names_map[cell_type]
                    for g in genes:
                        overlap_dbs.update(ct_db.get(g, ct_db.get(g.upper(), set())))
                    n_db_values.append(len(overlap_dbs))
                else:
                    # Fallback: max of per-gene counts (better than min)
                    db_counts = [n_databases_map[cell_type].get(g, n_databases_map[cell_type].get(g.upper(), 1)) for g in genes]
                    n_db_values.append(max(db_counts) if db_counts else 0)
            else:
                n_db_values.append(0)
        return n_db_values

    # Add N_Databases column if n_databases_map is available
    if n_databases_map and not sig_results.empty and 'Cell_type' in sig_results.columns and 'Overlapping_genes' in sig_results.columns:
        sig_results['N_Databases'] = _compute_n_databases(
            sig_results, 'Cell_type', 'Overlapping_genes', n_databases_map, db_names_map)

    if n_databases_map and not all_results.empty and 'Cell_type' in all_results.columns and 'Overlapping_genes' in all_results.columns:
        all_results['N_Databases'] = _compute_n_databases(
            all_results, 'Cell_type', 'Overlapping_genes', n_databases_map, db_names_map)

    analysis_params = {
        'p_val': pipeline.p_val_thresh,
        'log2fc': pipeline.log2fc_thresh,
        'mean_counts': pipeline.mean_counts_thresh,
        'top_n_genes': pipeline.top_genes
    }

    plots_html = {'heatmap': {}}
    if not sig_results.empty:
        top_n_options_for_plots = [1, 3, 5, 10]
        for n in top_n_options_for_plots:
            plots_html['heatmap'][str(n)] = plot_dynamic_heatmap_with_bars(
                sig_results_df=sig_results,
                top_n_celltypes=n
            )
    else:
        plots_html['heatmap']['1'] = "<h3>Heatmap</h3><p>No significant enrichments found.</p>"

    return pipeline, all_results, sig_results, plots_html, analysis_params

def generate_top_annotation_summary(final_sig_results, output_dir, job_name, tissue_specified,
                                    hierarchical_annotator=None, tissue_priority_ratio=0.0):
    """
    Generates a CSV with the top 1 annotation per cluster.

    Priority logic: Selected Tissue (Level 2) > All Tissue (Level 2), with
    score-gated fallback — if selected_tissue Combined_Score <
    tissue_priority_ratio * all_tissue Combined_Score, all_tissue wins.

    New output columns vs. prior version:
    - Runner_Up_Cell_Type / Score_Gap : annotation uncertainty signal
    - Broad_Type_Consensus            : True if runner-up agrees on Broad_Type lineage
    - Proliferating_Lineage           : Broad_Type + '_proliferating' when flag is True
    - Lineage_Conflict                : True if Top_Cell_Type is not a CL descendant of Broad_Type
    """
    summary_data = []
    all_clusters = set()

    def _get_clusters_from_ctx(ctx, lvl):
        if ctx in final_sig_results and lvl in final_sig_results[ctx]:
            df = final_sig_results[ctx][lvl]
            if isinstance(df, pd.DataFrame) and not df.empty and 'Cluster' in df.columns:
                return set(df['Cluster'].unique())
        return set()

    all_clusters.update(_get_clusters_from_ctx('selected_tissue', 'level2'))
    all_clusters.update(_get_clusters_from_ctx('all_tissue', 'level2'))

    sorted_clusters = sorted(list(all_clusters), key=natural_sort_key)

    hier_results = {}
    for ctx_name in ['selected_tissue', 'all_tissue']:
        if ctx_name in final_sig_results:
            h = final_sig_results[ctx_name].get('hierarchical')
            if isinstance(h, pd.DataFrame) and not h.empty:
                hier_results[ctx_name] = h

    def _sort_cluster_results(cluster_res):
        if cluster_res.empty:
            return cluster_res
        if 'Combined_Score' in cluster_res.columns and cluster_res['Combined_Score'].notna().any():
            return cluster_res.sort_values('Combined_Score', ascending=False, na_position='last')
        elif 'Weighted_Enrichment' in cluster_res.columns:
            return cluster_res.sort_values(
                ['adj_p_value', 'Weighted_Enrichment'], ascending=[True, False]
            )
        return cluster_res.sort_values(
            ['adj_p_value', 'Enrichment_ratio'], ascending=[True, False]
        )

    # Suppression window: only consider the top-N candidates by Combined_Score.
    # Suppression across the full significant list (300+ terms) is over-aggressive —
    # every broad term has some distant descendant with a better p-value somewhere.
    # Restricting to the top window ensures suppression only resolves genuine
    # head-to-head competition between closely-ranked terms.
    _SUPPRESSION_WINDOW = 10

    def _apply_descendant_suppression(cluster_res):
        """
        Within the top-N candidates (by Combined_Score), suppress any CL ancestor
        term when a more-specific descendant in that same window has an equal-or-better
        adjusted p-value.

        Conditional suppression: ancestor is only removed when descendant adj_p ≤
        ancestor adj_p — so a parent genuinely more significant than all its children
        (e.g. a mixed-lineage cluster) is preserved.

        Operates on the cell_name_to_id map from hierarchical_annotator; if the
        annotator is unavailable the function is a no-op.
        """
        if hierarchical_annotator is None or cluster_res.empty:
            return cluster_res

        name_to_id = hierarchical_annotator.cell_name_to_id
        p_col = 'adj_p_value'

        # Restrict competition window to top-N by Combined_Score
        window = cluster_res.head(_SUPPRESSION_WINDOW)

        # Build {cl_id: adj_p} for terms in the window with a resolvable CL ID
        cl_id_to_p = {}
        for _, row in window.iterrows():
            cid = name_to_id.get(str(row['Cell_type']).strip().upper())
            if cid:
                cl_id_to_p[cid] = float(row[p_col]) if pd.notna(row[p_col]) else 1.0

        if len(cl_id_to_p) < 2:
            return cluster_res

        # Mark an ancestor for suppression only when a window-peer descendant
        # has adj_p ≤ ancestor adj_p
        suppress_ids = set()
        cl_ids = list(cl_id_to_p.keys())
        for anc_id in cl_ids:
            anc_p = cl_id_to_p[anc_id]
            for desc_id in cl_ids:
                if desc_id == anc_id:
                    continue
                desc_ancestors = hierarchical_annotator._get_ancestors(desc_id)
                if anc_id in desc_ancestors and cl_id_to_p[desc_id] <= anc_p:
                    suppress_ids.add(anc_id)
                    break

        if not suppress_ids:
            return cluster_res

        suppress_names = set()
        id_to_name = {v: k for k, v in name_to_id.items()}
        for sid in suppress_ids:
            name = id_to_name.get(sid)
            if name:
                suppress_names.add(name.upper())

        filtered = cluster_res[
            ~cluster_res['Cell_type'].str.strip().str.upper().isin(suppress_names)
        ]

        if filtered.empty:
            return cluster_res

        n_suppressed = len(cluster_res) - len(filtered)
        if n_suppressed > 0:
            suppressed_labels = [id_to_name.get(s, s) for s in suppress_ids]
            logger.info(
                f"Descendant suppression (window={_SUPPRESSION_WINDOW}): "
                f"suppressed {n_suppressed} ancestor(s): {suppressed_labels}"
            )
        return filtered

    def _extract_top_rows(ctx_name, cluster):
        """Return (top_row, runner_up_row) from a context's level2 df, or (None, None)."""
        df = final_sig_results.get(ctx_name, {}).get('level2')
        if not isinstance(df, pd.DataFrame) or df.empty:
            return None, None
        cluster_res = df[df['Cluster'] == cluster]
        if cluster_res.empty:
            return None, None
        cluster_res = _sort_cluster_results(cluster_res)
        cluster_res = _apply_descendant_suppression(cluster_res)
        top = cluster_res.iloc[0]
        runner = cluster_res.iloc[1] if len(cluster_res) > 1 else None
        return top, runner

    def _row_score(row):
        if row is None:
            return None
        cs = row.get('Combined_Score')
        if cs is not None and not pd.isna(cs):
            return float(cs)
        we = row.get('Weighted_Enrichment')
        if we is not None and not pd.isna(we):
            return float(we)
        return float(row.get('Enrichment_ratio', 0) or 0)

    def _build_hit(cluster, top_row, runner_up_row, source_label):
        score = _row_score(top_row)
        runner_score = _row_score(runner_up_row)
        score_gap = (
            round(score - runner_score, 4)
            if score is not None and runner_score is not None
            else None
        )
        n_db = (
            int(top_row['N_Databases'])
            if 'N_Databases' in top_row.index and pd.notna(top_row.get('N_Databases'))
            else None
        )
        return {
            'Cluster': cluster,
            'Top_Cell_Type': top_row['Cell_type'],
            'Runner_Up_Cell_Type': runner_up_row['Cell_type'] if runner_up_row is not None else None,
            'Score_Gap': score_gap,
            'Source': source_label,
            'P_Value': top_row['adj_p_value'],
            'Score': score,
            'Genes': top_row.get('Overlapping_genes', ''),
            'N_Databases': n_db,
        }

    for cluster in sorted_clusters:
        st_top, st_runner = (
            _extract_top_rows('selected_tissue', cluster)
            if tissue_specified else (None, None)
        )
        at_top, at_runner = _extract_top_rows('all_tissue', cluster)

        # Score-gated tissue priority (Problem 5): fall back to all_tissue when
        # selected_tissue score is disproportionately weak.
        use_all_tissue = False
        if st_top is not None and at_top is not None and tissue_priority_ratio > 0:
            st_score = _row_score(st_top) or 0
            at_score = _row_score(at_top) or 0
            if st_score < tissue_priority_ratio * at_score:
                use_all_tissue = True
                logger.info(
                    f"{cluster}: selected_tissue score {st_score:.2f} < "
                    f"{tissue_priority_ratio} × all_tissue score {at_score:.2f} — using all_tissue"
                )

        if st_top is not None and not use_all_tissue:
            selected_hit = _build_hit(cluster, st_top, st_runner, 'Selected Tissue')
            source_ctx = 'selected_tissue'
        elif at_top is not None:
            selected_hit = _build_hit(cluster, at_top, at_runner, 'All Tissue')
            source_ctx = 'all_tissue'
        else:
            summary_data.append({
                'Cluster': cluster,
                'Top_Cell_Type': 'Unannotated',
                'Runner_Up_Cell_Type': None, 'Score_Gap': None,
                'Source': 'None', 'P_Value': None, 'Score': None, 'Genes': '',
                'Confidence': None, 'Broad_Type': None, 'Broad_Type_CL_ID': None,
                'Broad_Type_Consensus': None,
                'N_Databases': None, 'Proliferative_Flag': False,
                'Proliferative_Genes': '', 'Proliferating_Lineage': '',
                'Lineage_Conflict': False,
            })
            continue

        # Hierarchical context: confidence, broad type, consensus, lineage check
        confidence = None
        broad_type = None
        broad_type_cl_id = None
        broad_type_consensus = None
        lineage_conflict = False

        if source_ctx and source_ctx in hier_results and hierarchical_annotator:
            hier_df = hier_results[source_ctx]
            confidence_info = hierarchical_annotator.get_best_resolution(hier_df, cluster)
            if confidence_info:
                confidence = confidence_info.get('Confidence')

            broad_result = hierarchical_annotator.get_broad_type_with_consensus(
                hier_df, cluster,
                top_cell_type=selected_hit.get('Top_Cell_Type'),
                runner_up_cell_type=selected_hit.get('Runner_Up_Cell_Type'),
            )
            broad_type_consensus = broad_result['Broad_Type_Consensus']

            # Replace enrichment-derived Broad_Type with a pure CL DAG walk:
            # look up the Top_Cell_Type's CL ID, then find the shallowest
            # non-scaffold ancestor in the ontology — no hardcoding of cell biology.
            top_cl_id_for_broad = hierarchical_annotator.cell_name_to_id.get(
                str(selected_hit.get('Top_Cell_Type', '')).strip().upper()
            )
            if top_cl_id_for_broad:
                dag_name, dag_cl_id = hierarchical_annotator.get_broad_type_from_dag(
                    top_cl_id_for_broad
                )
                if dag_name:
                    broad_type = dag_name
                    broad_type_cl_id = dag_cl_id

            # Lineage conflict: Top_Cell_Type must be a CL descendant of Broad_Type.
            if broad_type_cl_id and top_cl_id_for_broad and top_cl_id_for_broad != broad_type_cl_id:
                top_ancestors = hierarchical_annotator._get_ancestors(top_cl_id_for_broad)
                lineage_conflict = broad_type_cl_id not in top_ancestors

        selected_hit['Confidence'] = confidence
        selected_hit['Broad_Type'] = broad_type
        selected_hit['Broad_Type_CL_ID'] = broad_type_cl_id
        selected_hit['Broad_Type_Consensus'] = broad_type_consensus
        selected_hit['Lineage_Conflict'] = lineage_conflict

        # --- Proliferative Signature Detection ---
        # Checks whether the overlapping marker genes are dominated by canonical
        # cell cycle genes (MSigDB HALLMARK_G2M_CHECKPOINT + HALLMARK_G1S_CHECKPOINT;
        # Liberzon et al., Cell Systems 2015).  In cancer/inflamed tissue, such clusters
        # may represent cycling or malignant cells rather than the annotated cell type.
        # Two thresholds are applied:
        #   (a) ≥3 cycle genes in overlap → unconditional flag
        #   (b) ≥2 cycle genes AND they represent ≥33 % of total overlap → flag
        # Requires MSigDB .grp files downloaded via download_reference_data().
        # If the gene set was not loaded, flag is set to None (not computed).
        if not PROLIFERATIVE_GENES:
            selected_hit['Proliferative_Flag'] = None
            selected_hit['Proliferative_Genes'] = ''
            selected_hit['Proliferating_Lineage'] = ''
        else:
            genes_str = selected_hit.get('Genes', '')
            if genes_str:
                overlap_genes = [g.strip().upper() for g in str(genes_str).split(',') if g.strip()]
                prolif_hits = [g for g in overlap_genes if g in PROLIFERATIVE_GENES]
                n_prolif = len(prolif_hits)
                n_total = len(overlap_genes)
                frac_prolif = n_prolif / n_total if n_total > 0 else 0.0
                is_prolif = (n_prolif >= _PROLIF_FLAG_MIN_GENES) or \
                            (n_prolif >= _PROLIF_FLAG_MIN_GENES2 and frac_prolif >= _PROLIF_FLAG_MIN_FRAC)
                selected_hit['Proliferative_Flag'] = is_prolif
                selected_hit['Proliferative_Genes'] = ', '.join(prolif_hits) if is_prolif else ''
                selected_hit['Proliferating_Lineage'] = (
                    f"{broad_type}_proliferating" if is_prolif and broad_type else ''
                )
            else:
                selected_hit['Proliferative_Flag'] = False
                selected_hit['Proliferative_Genes'] = ''
                selected_hit['Proliferating_Lineage'] = ''

        summary_data.append(selected_hit)

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        col_order = [
            'Cluster', 'Top_Cell_Type', 'Runner_Up_Cell_Type', 'Score_Gap',
            'Source', 'P_Value', 'Score', 'Genes', 'Confidence',
            'Broad_Type', 'Broad_Type_CL_ID', 'Broad_Type_Consensus', 'N_Databases',
            'Proliferative_Flag', 'Proliferative_Genes',
            'Proliferating_Lineage', 'Lineage_Conflict',
        ]
        summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]
        out_path = os.path.join(output_dir, f"{job_name}_top_annotation_summary.csv")
        summary_df.to_csv(out_path, index=False)
        logger.info(f"Top annotation summary saved to: {out_path}")
        return summary_df

    return pd.DataFrame()



def _compute_annotation_composition(final_sig, tissue_specified):
    """
    Derive per-cluster cell type composition scores from annotation enrichment results.

    For each cluster, collects all Level 2 significant enrichment results (selecting
    the tissue-specific context when available, falling back to all-tissue), then
    applies CL ontology ancestor pruning before normalising Combined_Score so that
    only the most specific (leaf) cell types contribute to the final proportions.

    Ancestor pruning:
        The Level 2 enrichment test returns all significantly enriched cell types,
        including generic CL ontological ancestors (e.g., "epithelial cell",
        "barrier cell", "phagocyte") whose marker gene sets are supersets of more
        specific descendants. Without pruning, these ancestors dilute the composition
        of the specific cell type (e.g., "keratinocyte" gets 12% instead of 50%).

        Pruning uses the hierarchical annotation output (N_Supporting column):
        - N_Supporting == 1 → "leaf" node: this type has no more specific significant
          descendant → keep for composition.
        - N_Supporting > 1 → "ancestor" node: this type has significant descendants
          that are more specific → exclude from composition.

        The Supporting_Types column of leaf rows (N_Supporting == 1) contains the
        database-facing cell type name (the sig_results Cell_type value). These are
        used as the "keep set" per cluster. Cell types absent from hierarchical
        (unmapped CL IDs) are kept by default.

    Args:
        final_sig      : dict {ctx_name: {lvl_name: DataFrame}} from run_enrichment_pipeline
        tissue_specified: bool — whether a tissue-specific context was run

    Returns:
        pd.DataFrame (n_clusters × n_cell_types), values = normalised composition scores,
        index = cluster names (e.g. "Cluster 1"), or None if no enrichment data available.
    """
    import numpy as np

    # Priority: tissue-specific Level 2 > all-tissue Level 2
    # Track which context we chose so we can load its matching hierarchical data
    sig_df = None
    chosen_ctx = None
    for ctx in (['selected_tissue', 'all_tissue'] if tissue_specified else ['all_tissue']):
        candidate = final_sig.get(ctx, {}).get('level2')
        if isinstance(candidate, pd.DataFrame) and not candidate.empty:
            sig_df = candidate
            chosen_ctx = ctx
            break

    if sig_df is None or sig_df.empty:
        return None
    if 'Cluster' not in sig_df.columns or 'Cell_type' not in sig_df.columns:
        return None

    # Use Combined_Score directly if present; otherwise compute it
    if 'Combined_Score' in sig_df.columns:
        score_col = 'Combined_Score'
    elif 'Weighted_Enrichment' in sig_df.columns and 'adj_p_value' in sig_df.columns:
        sig_df = sig_df.copy()
        sig_df['_comp_score'] = (
            sig_df['Weighted_Enrichment'].astype(float)
            * -np.log10(sig_df['adj_p_value'].astype(float).clip(lower=1e-300))
        )
        score_col = '_comp_score'
    else:
        sig_df = sig_df.copy()
        sig_df['_comp_score'] = sig_df.get('Enrichment_ratio', pd.Series(1.0, index=sig_df.index))
        score_col = '_comp_score'

    # --- CL Ancestor Pruning ---
    # Build per-cluster "leaf_names" from the hierarchical annotation:
    #
    # A leaf row (N_Supporting == 1) means this CL node has only one significant
    # contributor — itself. Its Supporting_Types field holds the database-facing
    # cell type name from the enrichment results. Collecting these names gives the
    # set of "most specific" types per cluster; any sig type NOT in this set is an
    # ancestor that has been superseded by a more specific descendant.
    #
    # Two additional filters prevent garbage passing the N_Supporting==1 check:
    #  • OBSOLETE terms: CL retires old classes (e.g. "barrier cell" → CL:0000215
    #    marked obsolete). The hierarchical Cell_Type starts with "obsolete" (case-
    #    insensitive); their database names would pollute the leaf set.
    #  • Depth < 2: root-adjacent terms (eukaryotic cell, animal cell, etc.) are not
    #    connected to the main CL sub-graph and should never appear in composition.
    #
    # Safety net — top-1 re-inclusion:
    #  If the highest-scoring sig type for a cluster was pruned as an ancestor (e.g.
    #  "macrophage" gets pruned because "inflammatory macrophage" is also significant),
    #  it is re-added to the pool. This ensures the dominant biology is always visible
    #  even when more specific descendants are also enriched.
    per_cluster_leaf_names = {}  # cluster → set of UPPERCASE cell type names (leaves)
    hier_df = final_sig.get(chosen_ctx, {}).get('hierarchical') if chosen_ctx else None
    if isinstance(hier_df, pd.DataFrame) and not hier_df.empty \
            and 'N_Supporting' in hier_df.columns \
            and 'Supporting_Types' in hier_df.columns \
            and 'Cluster' in hier_df.columns:
        for cluster, grp in hier_df.groupby('Cluster'):
            leaf_names = set()
            for _, row in grp.iterrows():
                if int(row.get('N_Supporting', 0)) != 1:
                    continue
                # Skip obsolete CL terms (Cell_Type starts with "obsolete")
                cl_name = str(row.get('Cell_Type', '') or '')
                if cl_name.lower().startswith('obsolete'):
                    continue
                # Skip root-adjacent terms (Depth < 2 means not integrated into main ontology)
                depth = int(row.get('Depth', 0) or 0)
                if depth < 2:
                    continue
                raw = str(row.get('Supporting_Types', '') or '')
                for name in raw.split(','):
                    name = name.strip()
                    if name:
                        leaf_names.add(name.upper())
            per_cluster_leaf_names[str(cluster)] = leaf_names

    # Build composition per cluster
    all_compositions = {}
    total_pruned = 0
    for cluster, grp in sig_df.groupby('Cluster'):
        leaf_names = per_cluster_leaf_names.get(str(cluster), set())

        if leaf_names:
            grp_pruned = grp[grp['Cell_type'].str.upper().isin(leaf_names)].copy()
            n_pruned = len(grp) - len(grp_pruned)

            # Safety net: if the top-scoring type was pruned as an ancestor, re-include it.
            # This keeps the dominant biological signal visible while still removing generic
            # ancestors further down the Combined_Score ranking.
            if not grp.empty:
                top_row = grp.nlargest(1, score_col)
                top_name = top_row['Cell_type'].values[0]
                if top_name.upper() not in leaf_names:
                    grp_pruned = pd.concat(
                        [grp_pruned, top_row], ignore_index=True
                    ).drop_duplicates(subset='Cell_type')
                    n_pruned = max(n_pruned - 1, 0)

            if grp_pruned.empty:
                grp_pruned = grp  # full fallback if pruning removed everything
                n_pruned = 0
            total_pruned += n_pruned
            grp = grp_pruned

        scores = grp.set_index('Cell_type')[score_col].astype(float).clip(lower=0)
        total = scores.sum()
        if total > 0:
            all_compositions[cluster] = (scores / total).to_dict()

    if not all_compositions:
        logger.warning("Annotation composition: no clusters with positive enrichment scores.")
        return None

    if total_pruned > 0:
        logger.info(
            f"Ancestor pruning: removed {total_pruned} generic CL ancestor types "
            f"from composition pool (top-scoring type per cluster always retained)"
        )

    comp_df = pd.DataFrame.from_dict(all_compositions, orient='index').fillna(0.0)
    comp_df.index.name = 'Cluster'
    comp_df = comp_df.reindex(sorted(comp_df.index, key=natural_sort_key))

    logger.info(
        f"Annotation composition: {comp_df.shape[0]} clusters × "
        f"{comp_df.shape[1]} cell types"
    )
    return comp_df


def _auto_detect_umap_path(deg_file_path):
    """
    Given a spaceranger DEG file path, reconstruct UMAP and cluster assignment paths.
    Returns (umap_path, clusters_path) if both exist, else (None, None).
    """
    from pathlib import Path
    try:
        analysis_dir = Path(deg_file_path).parent.parent.parent  # .../analysis/
        umap_path = analysis_dir / 'umap' / 'gene_expression_2_components' / 'projection.csv'
        clusters_path = analysis_dir / 'clustering' / 'gene_expression_graphclust' / 'clusters.csv'
        if umap_path.exists() and clusters_path.exists():
            logger.info(f"Auto-detected UMAP: {umap_path}")
            return str(umap_path), str(clusters_path)
    except Exception as e:
        logger.debug(f"UMAP auto-detection failed: {e}")
    return None, None


def _load_umap_data(umap_csv, clusters_csv, top_annotation_df, max_points=20000):
    """
    Load UMAP coordinates and cluster assignments, merge with top annotation data.
    Returns merged DataFrame or None on failure.
    """
    try:
        umap_df = pd.read_csv(umap_csv)
        clusters_df = pd.read_csv(clusters_csv)

        # Normalize column names
        umap_cols = {c: c for c in umap_df.columns}
        for c in umap_df.columns:
            if c.lower() == 'barcode':
                umap_cols[c] = 'Barcode'
            elif 'umap' in c.lower() and '1' in c:
                umap_cols[c] = 'UMAP-1'
            elif 'umap' in c.lower() and '2' in c:
                umap_cols[c] = 'UMAP-2'
        umap_df.rename(columns=umap_cols, inplace=True)

        cluster_cols = {c: c for c in clusters_df.columns}
        for c in clusters_df.columns:
            if c.lower() == 'barcode':
                cluster_cols[c] = 'Barcode'
            elif c.lower() == 'cluster':
                cluster_cols[c] = 'Cluster'
        clusters_df.rename(columns=cluster_cols, inplace=True)

        if 'Barcode' not in umap_df.columns or 'Barcode' not in clusters_df.columns:
            logger.warning("UMAP or clusters CSV missing 'Barcode' column")
            return None

        merged = umap_df[['Barcode', 'UMAP-1', 'UMAP-2']].merge(
            clusters_df[['Barcode', 'Cluster']], on='Barcode', how='inner')

        # Format cluster names to match enrichment output (e.g., "Cluster 1")
        merged['Cluster'] = 'Cluster ' + merged['Cluster'].astype(str)

        # Map cluster -> Top_Cell_Type from summary
        if top_annotation_df is not None and not top_annotation_df.empty:
            cluster_to_type = dict(zip(
                top_annotation_df['Cluster'], top_annotation_df['Top_Cell_Type']))
            merged['Top_Cell_Type'] = merged['Cluster'].map(cluster_to_type).fillna('Unannotated')
        else:
            merged['Top_Cell_Type'] = 'Unannotated'

        # Stratified subsample if too large
        if len(merged) > max_points:
            logger.info(f"Subsampling UMAP from {len(merged)} to {max_points} points")
            sampled_parts = []
            for _, grp in merged.groupby('Cluster'):
                n = max(50, int(len(grp) / len(merged) * max_points))
                sampled_parts.append(grp.sample(n=min(len(grp), n), random_state=42))
            merged = pd.concat(sampled_parts, ignore_index=True)

        logger.info(f"UMAP data loaded: {len(merged)} cells, {merged['Cluster'].nunique()} clusters")
        return merged

    except Exception as e:
        logger.warning(f"Failed to load UMAP data: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Run hybrid enrichment pipeline")
    parser.add_argument("input_path", help="Path to DEG CSV or Space Ranger dir")
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("output_dir", help="Directory to store output reports")
    parser.add_argument("--deg_type", choices=["spatial", "scrna"], required=True)
    parser.add_argument("--pval", type=float, default=0.05, help="Adj p-value threshold")
    parser.add_argument("--log2fc", type=float, default=1.0,
        help="Log2FC threshold (default: 1.0). For spatial data at 8 µm resolution, "
             "consider 0.5 — immune markers are attenuated when bins contain mixed cell types.")
    parser.add_argument("--mean", type=float, default=0.0, help="Mean count threshold")
    parser.add_argument("--topgenes", type=int, default=None, help="Top N genes")
    parser.add_argument("--marker_db", type=str, default=_DEFAULT_MARKER_DB,
        required=(_DEFAULT_MARKER_DB is None),
        help="Path to master_cell_marker_db.csv "
             f"(default: {_DEFAULT_MARKER_DB or 'REQUIRED — config path not found'})")
    parser.add_argument("--tissue", type=str, default=None, help="Filter marker DB for tissue")
    parser.add_argument("--min_overlap", type=int, default=2, help="Minimum gene overlap count to report (default: 2)")
    parser.add_argument("--min_db_markers", type=int, default=5,
        help="Minimum number of database marker genes a cell type must have (after background "
             "intersection) to be tested. Prevents rare cell types with tiny reference sets from "
             "inflating Weighted_Enrichment and winning Combined_Score. "
             "(default: 5; recommended: 5 for spatial, 0 to disable)")
    parser.add_argument("--min_cluster_degs", type=int, default=0,
        help="Minimum cluster DEG count required to compute Weighted Enrichment. "
             "Clusters with fewer DEGs fall back to unweighted scoring, avoiding N/n inflation. "
             "(default: 0 = no filter; recommended: 10 for spatial data)")
    parser.add_argument("--tissue_priority_ratio", type=float, default=0.0,
        help="Score-gated tissue priority: if selected_tissue Combined_Score < ratio × "
             "all_tissue Combined_Score, fall back to all_tissue annotation. "
             "(default: 0.0 = hard priority; recommended: 0.3)")
    parser.add_argument("--background_gene_count", type=int, default=None, help="Override background gene count N for hypergeometric test")
    parser.add_argument("--hgnc_map", type=str, default=_DEFAULT_HGNC_MAP,
        help="Path to HGNC complete set file for gene alias resolution "
             f"(default: {_DEFAULT_HGNC_MAP or 'None'})")
    parser.add_argument("--deg_format", type=str, default=None, choices=['seurat', 'scanpy', 'generic'], help="Force a specific DEG input format (default: auto-detect)")
    parser.add_argument("--no_hierarchy", action="store_true", help="Skip hierarchical annotation")
    parser.add_argument("--umap_csv", type=str, default=None,
        help="CSV with UMAP coordinates (columns: Barcode, UMAP-1, UMAP-2)")
    parser.add_argument("--cell_cluster_csv", type=str, default=None,
        help="CSV with cell-cluster mapping (columns: Barcode, Cluster)")
    parser.add_argument("--no_deconvolution", action="store_true",
        help="Disable the Cell Type Composition tab in the HTML report.")
    parser.add_argument("--no_auto_spatial_filter", action="store_true",
        help="Disable automatic mean_counts_thresh calibration for spatial data.")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # 1. Load DB
        full_df = pd.read_csv(args.marker_db, dtype=str)
        
        # 2. Prepare Contexts
        contexts = {}
        n_databases_contexts = {}

        logger.info("Preparing maps for 'All Tissue'...")
        db_names_contexts = {}
        contexts['all_tissue'], n_databases_contexts['all_tissue'], db_names_contexts['all_tissue'] = \
            get_marker_maps_context(full_df)

        tissue_specified = False
        if args.tissue:
            logger.info(f"Preparing maps for Selected Tissue: '{args.tissue}'...")
            subset_df = full_df[full_df['tissue_name'].str.contains(args.tissue, case=False, na=False)]
            if not subset_df.empty:
                contexts['selected_tissue'], n_databases_contexts['selected_tissue'], db_names_contexts['selected_tissue'] = \
                    get_marker_maps_context(subset_df)
                tissue_specified = True
            else:
                logger.warning(f"Tissue '{args.tissue}' not found. Skipping.")

        # 2b. Initialize Hierarchical Annotator (if ontology available)
        hierarchical_annotator = None
        if not args.no_hierarchy:
            try:
                import sys
                sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'parser'))
                from ontology_utils import CellxGeneOntologyParser
                logger.info("Initializing ontology parser for hierarchical annotation...")
                ontology_parser = CellxGeneOntologyParser()
                hierarchical_annotator = HierarchicalAnnotator(
                    ontology_parser=ontology_parser,
                    master_db_df=full_df,
                )
                logger.info("Hierarchical annotator initialized successfully.")
            except Exception as e:
                logger.warning(f"Could not initialize hierarchical annotator: {e}. Skipping hierarchy.")

        # 3. Find Files
        path_dict = find_deg_files(args.input_path, args.sample_name, args.deg_type)
        if not path_dict: return

        # 4. Process
        for job_name, deg_file in path_dict.items():
            logger.info(f"Processing Job: {job_name}")
            algo_out_dir = os.path.join(args.output_dir, job_name)
            os.makedirs(algo_out_dir, exist_ok=True)

            final_sig = {}
            final_plots = {}
            final_params = {}
            valid_pipe = None

            # Run pipeline for each context/level
            for ctx_name, maps in contexts.items():
                final_sig[ctx_name] = {}
                final_plots[ctx_name] = {}
                final_params[ctx_name] = {}

                for lvl_name, m_map in maps.items():
                    logger.info(f"--- Context: {ctx_name} | Level: {lvl_name} ---")

                    # Pass n_databases_map and db_names_map for N_Databases computation
                    ndb_map = n_databases_contexts.get(ctx_name, {}).get(lvl_name, None)
                    ndb_names = db_names_contexts.get(ctx_name, {}).get(lvl_name, None)

                    (pipe, all_r, sig_r, plots, params) = run_enrichment_pipeline(
                        deg_file, m_map, args, args.deg_type,
                        n_databases_map=ndb_map, db_names_map=ndb_names,
                    )

                    if pipe:
                        prefix = f"{job_name}_{ctx_name}_{lvl_name}"
                        all_r.to_csv(os.path.join(algo_out_dir, f"{prefix}_all_results.csv"), index=False)
                        sig_r.to_csv(os.path.join(algo_out_dir, f"{prefix}_sig_results.csv"), index=False)

                        final_sig[ctx_name][lvl_name] = sig_r
                        final_plots[ctx_name][lvl_name] = plots
                        final_params[ctx_name][lvl_name] = params
                        if valid_pipe is None: valid_pipe = pipe

                        # Run hierarchical annotation on Level 2 significant results
                        if hierarchical_annotator and lvl_name == 'level2' and not sig_r.empty:
                            try:
                                hier_results = hierarchical_annotator.annotate_all_clusters(sig_r)
                                if not hier_results.empty:
                                    hier_path = os.path.join(algo_out_dir, f"{prefix}_hierarchical.csv")
                                    hier_results.to_csv(hier_path, index=False)
                                    logger.info(f"Hierarchical annotation saved: {hier_path}")
                                    # Store for use in top annotation summary
                                    if 'hierarchical' not in final_sig[ctx_name]:
                                        final_sig[ctx_name]['hierarchical'] = {}
                                    final_sig[ctx_name]['hierarchical'] = hier_results
                            except Exception as e:
                                logger.warning(f"Hierarchical annotation failed for {ctx_name}/{lvl_name}: {e}")

            # --- Generate Top Annotation Summary CSV ---
            # This uses the priority logic: Selected Tissue (Level 2) > All Tissue (Level 2)
            top_annotation_df = generate_top_annotation_summary(
                final_sig_results=final_sig,
                output_dir=algo_out_dir,
                job_name=job_name,
                tissue_specified=tissue_specified,
                hierarchical_annotator=hierarchical_annotator,
                tissue_priority_ratio=getattr(args, 'tissue_priority_ratio', 0.0),
            )

            # --- Load UMAP data (auto-detect for spatial, CLI args for scRNA) ---
            umap_data = None
            if args.deg_type == 'spatial':
                umap_path, clusters_path = _auto_detect_umap_path(deg_file)
                if umap_path:
                    umap_data = _load_umap_data(umap_path, clusters_path, top_annotation_df)
            elif args.umap_csv and args.cell_cluster_csv:
                umap_data = _load_umap_data(args.umap_csv, args.cell_cluster_csv, top_annotation_df)

            # --- DEG Tables (spatial) ---
            if not valid_pipe:
                deg_html = "<p>Error: No DEGs processed.</p>"
            else:
                deg_html = create_deg_tables_html(
                    deg_df=valid_pipe.raw_deg_df,
                    cluster_markers=valid_pipe.cluster_markers,
                    p_val_thresh=valid_pipe.p_val_thresh,
                    log2fc_thresh=valid_pipe.log2fc_thresh,
                    mean_counts_thresh=valid_pipe.mean_counts_thresh
                )

            # --- Cell Type Composition (annotation-derived, works for scRNA-seq + spatial) ---
            # Derives per-cluster cell type proportions from enrichment scores:
            #   score(cluster, cell_type) = Weighted_Enrichment × −log10(adj_p_value)
            #   normalised across all significant cell types per cluster → sums to 1.0
            # This is scientifically appropriate for marker-list databases: uses them
            # via enrichment testing (as designed), not as expression profile surrogates.
            deconv_df = None
            if not getattr(args, 'no_deconvolution', False):
                try:
                    deconv_df = _compute_annotation_composition(
                        final_sig=final_sig,
                        tissue_specified=tissue_specified,
                    )
                    if deconv_df is not None:
                        comp_out = os.path.join(
                            algo_out_dir, f"{job_name}_composition_scores.csv"
                        )
                        deconv_df.to_csv(comp_out)
                        logger.info(f"Composition scores saved: {comp_out}")
                except Exception as e:
                    logger.warning(f"Composition scoring failed (skipping tab): {e}")
                    deconv_df = None

            # --- HTML Report ---
            # Collect hierarchical results for HTML report
            hier_for_report = {}
            for ctx_name in final_sig:
                h = final_sig[ctx_name].get('hierarchical')
                if isinstance(h, pd.DataFrame) and not h.empty:
                    hier_for_report[ctx_name] = h

            report_path = os.path.join(algo_out_dir, f"{job_name}_report.html")
            generate_html_report(
                sample_name=job_name,
                output_path=report_path,
                sig_results_maps=final_sig,
                plots_html_maps=final_plots,
                deg_table_html=deg_html,
                params_maps=final_params,
                selected_tissue_name=args.tissue if 'selected_tissue' in contexts else None,
                hierarchical_results=hier_for_report,
                top_annotation_df=top_annotation_df,
                umap_data=umap_data,
                deg_df=valid_pipe.raw_deg_df if valid_pipe else None,
                deconv_df=deconv_df,
            )
            logger.info(f"Report: {report_path}")

    except Exception as e:
        logger.exception("Pipeline failed: %s", e)

if __name__ == "__main__":
    main()