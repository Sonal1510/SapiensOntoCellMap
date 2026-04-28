#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2025-10-06
Description :
Marker enrichment analysis pipeline.
- Automatically selects 'Conventional' statistics for simple gene lists.
- Automatically selects 'Weighted' statistics for gene-weight dictionaries.
"""

import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import statsmodels.stats.multitest as multi
import os
import logging

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

# Known column schemas for auto-detecting DEG input formats
SCRNA_COLUMN_SCHEMAS = {
    'seurat': {
        'Feature Name': ['gene'],
        'Cluster': ['cluster'],
        'Adjusted p value': ['p_val_adj'],
        'Log2 fold change': ['avg_log2FC'],
    },
    'scanpy': {
        'Feature Name': ['names'],
        'Cluster': ['group'],
        'Adjusted p value': ['pvals_adj'],
        'Log2 fold change': ['logfoldchanges'],
    },
    'generic': {
        'Feature Name': ['Gene', 'gene_name', 'feature', 'symbol', 'Gene.name'],
        'Cluster': ['cluster_id', 'ident', 'celltype', 'group_id'],
        'Adjusted p value': ['padj', 'FDR', 'q_value', 'adj.P.Val', 'p_val_adj'],
        'Log2 fold change': ['log2FoldChange', 'logFC', 'avg_logFC', 'log2fc'],
    },
}


def _detect_deg_format(columns, forced_format=None):
    """
    Auto-detect DEG input format by matching column names against known schemas.

    Args:
        columns: list of column names from the input DataFrame
        forced_format: if set, skip auto-detection and use this format

    Returns:
        (format_name, rename_dict) where rename_dict maps original → standard names
    """
    col_set = set(columns)

    if forced_format and forced_format in SCRNA_COLUMN_SCHEMAS:
        schema = SCRNA_COLUMN_SCHEMAS[forced_format]
        rename = {}
        for std_name, possible_names in schema.items():
            for pn in possible_names:
                if pn in col_set:
                    rename[pn] = std_name
                    break
        return forced_format, rename

    # Try each schema in priority order: seurat, scanpy, then generic
    for fmt_name in ['seurat', 'scanpy', 'generic']:
        schema = SCRNA_COLUMN_SCHEMAS[fmt_name]
        rename = {}
        matched_all = True
        for std_name, possible_names in schema.items():
            found = False
            for pn in possible_names:
                if pn in col_set:
                    rename[pn] = std_name
                    found = True
                    break
            if not found:
                matched_all = False
                break
        if matched_all:
            return fmt_name, rename

    # No match found — return None
    return None, {}

def _build_hgnc_alias_map(hgnc_path):
    """
    Build a gene alias → approved symbol mapping from the HGNC complete set file.

    Parses the 'symbol', 'alias_symbol', and 'prev_symbol' columns to create
    a dictionary mapping all known aliases (uppercased) to their approved
    HGNC symbol (uppercased).

    Args:
        hgnc_path: Path to hgnc_complete_set.txt (tab-separated)

    Returns:
        dict: {alias_upper: approved_symbol_upper}
    """
    try:
        df = pd.read_csv(hgnc_path, sep='\t', usecols=['symbol', 'alias_symbol', 'prev_symbol'],
                         dtype=str, na_values=[''])
    except (ValueError, FileNotFoundError) as e:
        logging.warning(f"Could not load HGNC file: {e}")
        return {}

    alias_map = {}
    for _, row in df.iterrows():
        approved = str(row['symbol']).strip().upper() if pd.notna(row['symbol']) else None
        if not approved:
            continue
        for alias_col in ['alias_symbol', 'prev_symbol']:
            raw = row[alias_col]
            if pd.isna(raw):
                continue
            # HGNC uses '|' as delimiter within alias columns
            aliases = str(raw).split('|')
            for alias in aliases:
                alias_clean = alias.strip().upper()
                if alias_clean and alias_clean != approved:
                    # First mapping wins — approved symbols shouldn't be overridden
                    if alias_clean not in alias_map:
                        alias_map[alias_clean] = approved
    return alias_map


class MarkerEnrichmentTest:
    def __init__(
        self,
        deg_file,
        marker_db_dict,
        deg_file_type='spatial',
        p_val_thresh=0.05,
        log2fc_thresh=1.0,
        mean_counts_thresh=0,
        top_genes=None,
        min_overlap=2,
        min_db_markers=0,
        min_cluster_degs=0,
        background_gene_count=None,
        hgnc_map=None,
        deg_format=None,
        auto_spatial_mean_filter=True,
    ):
        if not deg_file or not os.path.exists(deg_file):
            raise FileNotFoundError(f"DEG file not found: {deg_file}")

        self.raw_deg_df = pd.read_csv(deg_file)
        self.deg_file_type = deg_file_type
        
        # --- DETECT MODE (Weighted vs Conventional) ---
        # Normalize DB to {cell: {gene: weight}} for internal processing, 
        # but keep a flag to know if we should output weighted columns.
        self.marker_dict_internal = {}
        self.is_weighted_mode = False
        
        all_db_genes = set()
        
        # Check first item to determine mode (assuming consistency)
        first_val = next(iter(marker_db_dict.values())) if marker_db_dict else []
        if isinstance(first_val, dict):
            self.is_weighted_mode = True
            logging.info("Input is dictionary: Activating WEIGHTED Enrichment Statistics.")
        else:
            self.is_weighted_mode = False
            logging.info("Input is list: Activating CONVENTIONAL Enrichment Statistics.")

        self._gene_display_name = {}  # UPPER -> original case

        for k, v in marker_db_dict.items():
            k_str = str(k).strip()
            if self.is_weighted_mode:
                clean_v = {}
                for g, w in v.items():
                    if pd.notna(g):
                        original = str(g).strip()
                        normalized = original.upper()
                        clean_v[normalized] = float(w)
                        self._gene_display_name[normalized] = original
            else:
                clean_v = {}
                for g in v:
                    if pd.notna(g):
                        original = str(g).strip()
                        normalized = original.upper()
                        clean_v[normalized] = 1.0
                        self._gene_display_name[normalized] = original

            self.marker_dict_internal[k_str] = clean_v
            all_db_genes.update(clean_v.keys())

        self.p_val_thresh = p_val_thresh
        self.log2fc_thresh = log2fc_thresh
        self.mean_counts_thresh = mean_counts_thresh
        self.top_genes = top_genes
        self.min_overlap = min_overlap
        self.min_db_markers = min_db_markers
        self.min_cluster_degs = min_cluster_degs
        self.N_override = background_gene_count
        self.deg_format = deg_format
        self.auto_spatial_mean_filter = auto_spatial_mean_filter

        # --- HGNC Gene Alias Resolution ---
        self.alias_map = {}
        hgnc_path = hgnc_map
        if hgnc_path is None:
            # Try config path first, then fall back to relative path
            try:
                import sys
                _proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
                if _proj_root not in sys.path:
                    sys.path.insert(0, _proj_root)
                from config.config import HGNC_COMPLETE_SET_FILE
                if os.path.exists(HGNC_COMPLETE_SET_FILE):
                    hgnc_path = HGNC_COMPLETE_SET_FILE
            except ImportError:
                pass
            if hgnc_path is None:
                default_hgnc = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'reference', 'hgnc_complete_set.txt')
                if os.path.exists(default_hgnc):
                    hgnc_path = default_hgnc
        if hgnc_path and os.path.exists(hgnc_path):
            self.alias_map = _build_hgnc_alias_map(hgnc_path)
            logging.info(f"Loaded {len(self.alias_map)} gene aliases from HGNC map")
        elif hgnc_map:
            logging.warning(f"HGNC map file not found: {hgnc_map}. Skipping alias resolution.")

        # Resolve marker DB genes to approved HGNC symbols
        if self.alias_map:
            resolved_marker_count = 0
            for cell_type in list(self.marker_dict_internal.keys()):
                genes = self.marker_dict_internal[cell_type]
                resolved = {}
                for gene, weight in genes.items():
                    canonical = self.alias_map.get(gene, gene)
                    if canonical != gene:
                        resolved_marker_count += 1
                    if canonical in resolved:
                        resolved[canonical] = max(resolved[canonical], weight)
                    else:
                        resolved[canonical] = weight
                    # Keep display name for canonical
                    if canonical not in self._gene_display_name:
                        self._gene_display_name[canonical] = self._gene_display_name.get(gene, gene)
                self.marker_dict_internal[cell_type] = resolved
            # Update all_db_genes
            all_db_genes = set()
            for v in self.marker_dict_internal.values():
                all_db_genes.update(v.keys())
            logging.info(f"Gene alias resolution: {resolved_marker_count} marker genes resolved to approved HGNC symbols")

        self.deg_df_long = self._normalize_deg_df()

        if self.auto_spatial_mean_filter and self.deg_file_type == 'spatial' and self.mean_counts_thresh == 0.0:
            if 'Mean Counts' in self.deg_df_long.columns:
                pos = self.deg_df_long[self.deg_df_long['Mean Counts'] > 0]['Mean Counts']
                if not pos.empty:
                    self.mean_counts_thresh = pos.quantile(0.75)
                    logging.info(
                        f"Auto-calibrated spatial mean_counts_thresh to 75th percentile: "
                        f"{self.mean_counts_thresh:.2f}. Pass auto_spatial_mean_filter=False to disable."
                    )

        self.background_genes = set(self.deg_df_long["Feature Name"].unique())

        # Heuristic warning for scRNA-seq with likely underestimated N
        if self.N_override is None and self.deg_file_type == 'scrna':
            n_genes = len(self.background_genes)
            if n_genes < 15000 and 'Adjusted p value' in self.deg_df_long.columns:
                pvals = self.deg_df_long['Adjusted p value'].dropna()
                if len(pvals) > 0 and (pvals < 0.05).mean() > 0.80:
                    logging.warning(
                        f"Detected N={n_genes} genes. For whole-transcriptome scRNA-seq, "
                        f"N should be ~20,000. Over 80% of adjusted p-values are < 0.05, "
                        f"suggesting the DEG file contains only significant genes. "
                        f"Use --background_gene_count to override."
                    )

        if self.N_override is not None:
            logging.info(f"Using user-specified background gene count N={self.N_override}")

        self.cluster_markers = {}
        self.results_ = pd.DataFrame()

    def _normalize_deg_df(self):
        if self.deg_file_type == 'scrna':
            df = self.raw_deg_df.copy()

            # Auto-detect or force input format
            fmt_name, rename = _detect_deg_format(df.columns.tolist(), forced_format=self.deg_format)

            if fmt_name is None:
                expected = []
                for schema_name, schema in SCRNA_COLUMN_SCHEMAS.items():
                    cols = [aliases[0] for aliases in schema.values()]
                    expected.append(f"  {schema_name}: {cols}")
                raise ValueError(
                    f"Could not auto-detect DEG input format. Columns found: {df.columns.tolist()}\n"
                    f"Expected one of:\n" + "\n".join(expected) + "\n"
                    f"Use --deg_format to force a specific format."
                )

            logging.info(f"Auto-detected DEG format: {fmt_name} (columns: {list(rename.keys())})")
            df.rename(columns=rename, inplace=True)

            if 'Mean Counts' not in df.columns:
                df['Mean Counts'] = 0
            df['Cluster'] = "Cluster " + df['Cluster'].astype(str)
            df['Feature Name'] = df['Feature Name'].astype(str).str.strip().str.upper()

            # Resolve DEG gene aliases
            if self.alias_map:
                resolved_count = 0
                resolved_names = []
                for gene in df['Feature Name']:
                    canonical = self.alias_map.get(gene, gene)
                    if canonical != gene:
                        resolved_count += 1
                    resolved_names.append(canonical)
                df['Feature Name'] = resolved_names
                if resolved_count > 0:
                    logging.info(f"Gene alias resolution: {resolved_count} DEG genes resolved to approved HGNC symbols")

            return df[['Feature Name', 'Cluster', 'Adjusted p value', 'Log2 fold change', 'Mean Counts']]
        elif self.deg_file_type == 'spatial':
            df = self.raw_deg_df.copy()
            id_vars = [c for c in ['Feature ID', 'Feature Name'] if c in df.columns]
            df_long = pd.melt(df, id_vars=id_vars, var_name='Metric', value_name='Value')
            df_long[['Cluster', 'Metric Type']] = df_long['Metric'].str.extract(r'(Cluster \d+)\s(.*)')
            df_long.dropna(subset=['Cluster', 'Metric Type'], inplace=True)
            df = df_long.pivot_table(index=id_vars+['Cluster'], columns='Metric Type', values='Value', aggfunc='first').reset_index()
            df.rename(columns={'Adjusted p value': 'Adjusted p value', 'Log2 fold change': 'Log2 fold change', 'Mean Counts': 'Mean Counts'}, inplace=True)
            df.fillna({'Adjusted p value': 1.0, 'Log2 fold change': 0.0, 'Mean Counts': 0.0}, inplace=True)
            df['Feature Name'] = df['Feature Name'].astype(str).str.strip().str.upper()

            # Resolve spatial DEG gene aliases
            if self.alias_map:
                resolved_count = 0
                resolved_names = []
                for gene in df['Feature Name']:
                    canonical = self.alias_map.get(gene, gene)
                    if canonical != gene:
                        resolved_count += 1
                    resolved_names.append(canonical)
                df['Feature Name'] = resolved_names
                if resolved_count > 0:
                    logging.info(f"Gene alias resolution: {resolved_count} spatial DEG genes resolved to approved HGNC symbols")

            return df
        else:
            raise ValueError(f"Unknown type: {self.deg_file_type}")

    def fit(self):
        self.cluster_markers = self._filter_and_extract_degs()
        self.results_ = self._run_test()
        return self

    def _filter_and_extract_degs(self):
        # ... (Same logic) ...
        d = {}
        for cid, cdf in self.deg_df_long.groupby('Cluster'):
            mask = (cdf['Adjusted p value'] < self.p_val_thresh) & (cdf['Log2 fold change'] > self.log2fc_thresh)
            if self.mean_counts_thresh > 0 and 'Mean Counts' in cdf.columns:
                mask &= (cdf['Mean Counts'] > self.mean_counts_thresh)
            res = cdf[mask].copy()
            if self.top_genes: res = res.sort_values('Log2 fold change', ascending=False).head(self.top_genes)
            d[cid] = res["Feature Name"].dropna().astype(str).unique().tolist()
        return d

    def _run_test(self):
        N = self.N_override if self.N_override else len(self.background_genes)
        all_results_list = []

        for cluster_id, cluster_genes in self.cluster_markers.items():
            n = len(cluster_genes)
            if n == 0: continue

            # WE inflates when n is small (N/n term dominates). Skip weighted
            # scoring below threshold so these clusters rank by adj_p only.
            use_weighted = self.is_weighted_mode and (
                self.min_cluster_degs == 0 or n >= self.min_cluster_degs
            )
            if self.is_weighted_mode and not use_weighted:
                logging.warning(
                    f"{cluster_id}: {n} DEGs < min_cluster_degs={self.min_cluster_degs}; "
                    f"WE skipped — using unweighted enrichment for this cluster."
                )

            results = []
            for cell_type, weighted_genes in self.marker_dict_internal.items():
                cell_genes_set = set(weighted_genes.keys()).intersection(self.background_genes)
                K = len(cell_genes_set)
                if K < self.min_db_markers: continue

                overlap_genes = set(cluster_genes).intersection(cell_genes_set)
                k = len(overlap_genes)
                if k < self.min_overlap: continue

                # 1. P-value (Always Hypergeometric on counts)
                p_val = hypergeom.sf(k - 1, N, K, n)

                # 2. Base Metrics
                enrichment_ratio = (k / n) / (K / N) if K > 0 else 0

                # Display original-case gene names in output
                display_genes = [self._gene_display_name.get(g, g) for g in sorted(overlap_genes)]

                row = {
                    "Cluster": cluster_id,
                    "Cell_type": cell_type,
                    "p_value": p_val,
                    "Enrichment_ratio": enrichment_ratio,
                    "Overlapping_genes_count": k,
                    "Overlapping_genes": ", ".join(display_genes),
                }

                # 3. Weighted Metrics (Only if weighted mode AND cluster n is sufficient)
                if use_weighted:
                    overlap_w_sum = sum(weighted_genes[g] for g in overlap_genes)
                    ref_w_sum = sum(weighted_genes[g] for g in cell_genes_set)

                    # Weighted Recall: fraction of total reference evidence weight captured.
                    # = overlap_w_sum / ref_w_sum
                    # Values in [0,1]; 1 means all evidence-weighted markers were found.
                    row["Weighted_Recall"] = round(
                        overlap_w_sum / ref_w_sum if ref_w_sum > 0 else 0, 4
                    )

                    # Weighted Enrichment Ratio: observed evidence density / expected evidence density.
                    # = (overlap_w_sum / n) / (ref_w_sum / N)
                    # Equivalent to Weighted_Recall / (n/N) — the weighted analog of the
                    # standard enrichment ratio (k/n)/(K/N). Values > 1 indicate enrichment.
                    row["Weighted_Enrichment"] = round(
                        (overlap_w_sum / n) / (ref_w_sum / N)
                        if ref_w_sum > 0 and N > 0 else 0, 4
                    )

                results.append(row)

            if results:
                all_results_list.append(pd.DataFrame(results))

        if not all_results_list: return pd.DataFrame()

        # Global FDR correction across ALL clusters and cell types.
        # Per-cluster FDR was tested and reverted: it introduced false positives
        # (e.g. Sertoli cell in granular epidermis) because the looser per-cluster
        # threshold passes borderline hits that global BH correctly rejects.
        combined = pd.concat(all_results_list, ignore_index=True)
        combined["adj_p_value"] = multi.multipletests(combined["p_value"], method="fdr_bh")[1]

        # Combined_Score: integrates statistical significance with evidence-weighted effect size.
        # = Weighted_Enrichment * -log10(adj_p_value)
        # Provides a single ranking metric that rewards both low p-values and high evidence weight.
        # Only computed in weighted mode.
        if self.is_weighted_mode and "Weighted_Enrichment" in combined.columns:
            adj_p_clipped = combined["adj_p_value"].clip(lower=1e-300)
            combined["Combined_Score"] = (
                combined["Weighted_Enrichment"] * (-np.log10(adj_p_clipped))
            ).round(4)

        return combined.sort_values(["Cluster", "adj_p_value"])