#!/usr/bin/env python3

"""
SapiensOntoCellMap Database Visualizer Generator (v9.0)

Generates:
1. Interactive HTML database explorer (Dashboard, Ontology Browser,
   Cell Explorer, Gene Explorer) — self-contained single-file HTML.
2. Static publication-quality summary plot (PNG).

Data pipeline:
  CSV → clean → dashboard stats → cell lineages → ontology sunburst
      → tissue-cell matrix → cell marker evidence → gene summaries
      → JSON embed into HTML template.
"""

import pandas as pd
import numpy as np
import json
import logging
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from typing import Dict, Any, List
from collections import defaultdict, Counter

# --- Path Management & Imports ---
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if project_root not in sys.path:
        sys.path.append(project_root)

    from src.parser.ontology_utils import CellxGeneOntologyParser

    from config.config import (
        PROCESSED_COMBINED_DATABASE_FILE,
        PROCESSED_COMBINED_DATABASE_FILE_HTML
    )
except ImportError as e:
    print(f"CRITICAL ERROR: Failed to import modules. {e}")
    print("Ensure project structure contains: src/parser/ontology_utils.py and config/config.py")
    exit(1)
except FileNotFoundError as e:
    print(f"CRITICAL ERROR: {e}")
    exit(1)

# --- Configuration ---
PLOT_OUTPUT_FILENAME = "database_summary_2x2_publishable_v4.png"

# Gene-centric data: genes with >= this many DBs get full source details
GENE_DETAIL_THRESHOLD = 3
GENE_TOP_N_SOURCES = 15
GENE_TOP_N_BARS = 8


class SapiensMapGenerator:
    """
    Generates the SapiensOntoCellMap interactive HTML visualizer
    and static publication summary plot from the master marker database.
    """

    LOG_LEVEL = logging.INFO

    def __init__(self, input_csv: str, template_html: str, output_html: str):
        self.input_csv_path = input_csv
        self.template_html_path = template_html
        self.output_html_path = output_html

        logging.basicConfig(level=self.LOG_LEVEL,
                            format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

        # Core data
        self.df: pd.DataFrame = pd.DataFrame()
        self.all_genes: List[str] = []
        self.all_tissues_list: List[Dict] = []
        self.all_cells_list: List[Dict] = []

        # Computed data containers
        self.dashboard_stats: Dict[str, Any] = {}
        self.db_contribution: List[Dict] = []
        self.source_type_dist: List[Dict] = []
        self.cross_db_agreement: Dict[str, Any] = {}
        self.cell_lineage_map: Dict[str, Any] = {}
        self.sunburst_data: Dict[str, Any] = {}
        self.tissue_cell_matrix: Dict[str, Any] = {}
        self.cell_markers: Dict[str, List] = {}
        self.cell_tissue_summary: Dict[str, List] = {}
        self.gene_to_summary: Dict[str, Any] = {}
        self.hierarchy_tree: Dict[str, Any] = {}
        self.dot_matrix: Dict[str, Any] = {}
        self.nlp_quality_map: Dict[str, str] = {}
        self.source_type_per_tissue: Dict[str, Any] = {}
        self.final_data: Dict[str, Any] = {}

        # Initialize Ontology Parser
        try:
            self.ontology = CellxGeneOntologyParser()
            self.logger.info("CellxGeneOntologyParser initialized successfully.")
        except Exception as e:
            self.logger.error(f"Failed to initialize CellxGeneOntologyParser: {e}",
                              exc_info=True)
            raise

    # =================================================================
    # 1. DATA LOADING
    # =================================================================

    def _load_and_clean_data(self) -> bool:
        """Load raw CSV, clean IDs, populate helper lists."""
        self.logger.info(f"Loading raw data from {self.input_csv_path}...")
        try:
            df = pd.read_csv(self.input_csv_path, low_memory=False)
        except FileNotFoundError:
            self.logger.error(f"Input file not found: {self.input_csv_path}")
            return False

        self.logger.info(f"Loaded {len(df)} total rows.")

        df.dropna(subset=['tissue_id', 'cell_id'], inplace=True)

        # Standardize separators (CL_123 -> CL:123)
        df['cell_id'] = df['cell_id'].astype(str).str.replace('_', ':')
        df['tissue_id'] = df['tissue_id'].astype(str).str.replace('_', ':')

        # Extract valid ontology patterns only
        df['cell_id'] = df['cell_id'].str.extract(r'(CL:\d+)', expand=False)
        df['tissue_id'] = df['tissue_id'].str.extract(r'(UBERON:\d+)', expand=False)
        df.dropna(subset=['tissue_id', 'cell_id'], inplace=True)

        # Fill metadata NaNs
        df['tissue_name'] = df['tissue_name'].fillna(
            df['db_tissue_name']).fillna('Unknown Tissue')
        df['cell_name'] = df['cell_name'].fillna(
            df['db_cell_name']).fillna('Unknown Cell')
        df['gene'] = df['gene'].fillna('Unknown Gene').astype(str)
        df['database'] = df['database'].fillna('N/A').astype(str)
        df['source_type'] = df['source_type'].fillna('N/A').astype(str)
        df['source_info'] = df['source_info'].fillna('N/A').astype(str)

        # Global gene list — filter out non-gene entries (dates, multi-word strings, empty)
        # Date strings from source DBs (e.g. "2019-09-02 00:00:00") contain spaces;
        # valid gene symbols do not.
        all_genes_set = set(df['gene'].unique()) - {'Unknown Gene'}
        all_genes_set = {g for g in all_genes_set if g and ' ' not in g and len(g) <= 60}
        self.all_genes = sorted(list(all_genes_set))

        # Global tissue list
        tissue_df = df[['tissue_id', 'tissue_name']].drop_duplicates().sort_values(
            by='tissue_name')
        self.all_tissues_list = tissue_df.to_dict('records')

        # Global cell list (canonical name per cell_id)
        cell_name_map = df.groupby('cell_id')['cell_name'].apply(
            lambda x: x.mode().iloc[0]).to_dict()
        self.all_cells_list = sorted(
            [{'cell_id': cid, 'cell_name': cname}
             for cid, cname in cell_name_map.items()],
            key=lambda x: x['cell_name']
        )

        self.logger.info(f"Cleaned data: {len(df)} rows, "
                         f"{len(self.all_genes)} genes, "
                         f"{len(self.all_cells_list)} cell types, "
                         f"{len(self.all_tissues_list)} tissues.")
        self.df = df
        return True

    # =================================================================
    # 2. DASHBOARD DATA
    # =================================================================

    def _compute_dashboard_data(self):
        """Compute summary statistics and chart data for Dashboard tab."""
        self.logger.info("Computing dashboard data...")

        # Summary stats
        self.dashboard_stats = {
            "total_markers": int(len(self.df)),
            "unique_genes": int(self.df['gene'].nunique()),
            "unique_cell_types": int(self.df['cell_id'].nunique()),
            "unique_tissues": int(self.df['tissue_id'].nunique()),
            "num_databases": int(self.df['database'].nunique()),
        }

        # Database contribution (for treemap)
        db_agg = self.df.groupby('database').agg(
            count=('gene', 'size'),
            unique_genes=('gene', 'nunique'),
            unique_cells=('cell_id', 'nunique'),
            unique_tissues=('tissue_id', 'nunique')
        ).sort_values('count', ascending=False)

        total = len(self.df)
        self.db_contribution = [
            {
                "name": db,
                "count": int(row['count']),
                "genes": int(row['unique_genes']),
                "cells": int(row['unique_cells']),
                "tissues": int(row['unique_tissues']),
                "pct": round(row['count'] / total * 100, 1)
            }
            for db, row in db_agg.iterrows()
        ]

        # Source type distribution (for donut)
        src = self.df.groupby('source_type').size().sort_values(ascending=False)
        self.source_type_dist = [
            {"name": st, "count": int(c), "pct": round(c / total * 100, 1)}
            for st, c in src.items()
        ]

        # Cross-database agreement
        agreement = self.df.groupby(['gene', 'cell_id'])['database'].nunique()
        hist = agreement.value_counts().sort_index()
        multi_db = int((agreement > 1).sum())
        self.cross_db_agreement = {
            "labels": [str(n) for n in hist.index],
            "values": [int(v) for v in hist.values],
            "total_pairs": int(len(agreement)),
            "multi_db_pairs": multi_db,
            "multi_db_pct": round(multi_db / len(agreement) * 100, 1)
                            if len(agreement) > 0 else 0
        }

        self.logger.info(f"Dashboard: {self.dashboard_stats['num_databases']} databases, "
                         f"{self.cross_db_agreement['multi_db_pct']}% multi-DB agreement.")

    # =================================================================
    # 3. CELL ONTOLOGY LINEAGES
    # =================================================================

    def _build_cell_lineages(self):
        """Build lineage maps, pruning ancestors above 'cell' root."""
        self.logger.info("Building ancestor lineages...")
        lineage_map = {}
        unique_cell_ids = self.df['cell_id'].unique()

        for cell_id in unique_cell_ids:
            if cell_id not in self.ontology.cl_id_to_name:
                continue
            if self.ontology.ontology_parser.is_term_deprecated(cell_id):
                continue

            try:
                ancestors = self.ontology.ontology_parser.get_term_ancestors_with_distances(
                    cell_id)
                ancestors[cell_id] = 0

                # Prune above root 'cell' node
                root_cell_id = 'CL:0000000'
                if root_cell_id in ancestors:
                    root_dist = ancestors[root_cell_id]
                    ancestors = {a: d for a, d in ancestors.items()
                                 if d <= root_dist}

                lineage_for_cell = {}
                for ancestor_id, dist in ancestors.items():
                    if self.ontology.ontology_parser.is_term_deprecated(ancestor_id):
                        continue
                    ancestor_name = self.ontology.cl_id_to_name.get(
                        ancestor_id, ancestor_id)
                    lineage_for_cell[ancestor_id] = {
                        'dist': dist, 'name': ancestor_name
                    }

                lineage_map[cell_id] = lineage_for_cell
            except Exception:
                pass

        self.cell_lineage_map = lineage_map
        self.logger.info(f"Built lineages for {len(lineage_map)} cell types.")

    # =================================================================
    # 4. ONTOLOGY SUNBURST (DAG → TREE)
    # =================================================================

    def _build_ontology_sunburst(self):
        """Build CL ontology tree for interactive sunburst chart.

        Converts the DAG to a tree using majority-vote parent assignment:
        for each node, the most frequently observed parent across all cell
        lineages wins.
        """
        self.logger.info("Building ontology sunburst...")

        if not self.cell_lineage_map:
            self.logger.warning("No lineage data — sunburst will be empty.")
            self.sunburst_data = {'ids': [], 'labels': [], 'parents': [],
                                  'values': [], 'in_db': []}
            return

        cell_marker_counts = self.df.groupby('cell_id')['gene'].nunique().to_dict()
        cells_in_db = set(self.df['cell_id'].unique())

        # Collect all nodes from lineage data
        all_nodes = {}  # node_id -> name
        for cell_id, lineage in self.cell_lineage_map.items():
            for node_id, data in lineage.items():
                if node_id not in all_nodes:
                    all_nodes[node_id] = data['name']

        # Vote on parent-child relationships using CL-only lineage chains
        parent_votes = defaultdict(Counter)
        root_id = 'CL:0000000'

        for cell_id, lineage in self.cell_lineage_map.items():
            # Filter to CL-only terms to avoid non-CL ancestors causing cycles
            cl_items = [(nid, data) for nid, data in lineage.items()
                        if nid.startswith('CL:')]
            # Sort by distance descending: root first, cell last
            chain = sorted(cl_items, key=lambda x: x[1]['dist'],
                           reverse=True)
            for i in range(1, len(chain)):
                child = chain[i][0]
                parent = chain[i - 1][0]
                # Root is never a child; guard against self-loops
                if child != parent and child != root_id:
                    parent_votes[child][parent] += 1

        # Pick most-voted parent per node
        node_parents = {}
        for node_id, votes in parent_votes.items():
            best_parent = votes.most_common(1)[0][0]
            node_parents[node_id] = best_parent

        # Force root parent to '' (MUST be after vote loop)
        node_parents[root_id] = ''

        # Break any remaining cycles via parent-chain walk
        for node_id in list(node_parents.keys()):
            visited = set()
            cur = node_id
            while cur and node_parents.get(cur, '') != '':
                if cur in visited:
                    # Cycle — attach this node to root to break it
                    node_parents[cur] = root_id
                    break
                visited.add(cur)
                cur = node_parents.get(cur)

        # Cells without lineage → attach to root
        for cell_id in cells_in_db:
            if cell_id not in node_parents and cell_id in all_nodes:
                node_parents[cell_id] = root_id

        # Build Plotly sunburst arrays
        ids, labels, parents, values, in_db = [], [], [], [], []

        for node_id in all_nodes:
            if node_id not in node_parents:
                continue
            ids.append(node_id)
            labels.append(all_nodes[node_id])
            parents.append(node_parents[node_id])
            values.append(cell_marker_counts.get(node_id, 0))
            in_db.append(node_id in cells_in_db)

        self.sunburst_data = {
            'ids': ids, 'labels': labels, 'parents': parents,
            'values': values, 'in_db': in_db
        }
        self.logger.info(f"Sunburst: {len(ids)} nodes "
                         f"({sum(in_db)} in database, "
                         f"{len(ids) - sum(in_db)} ancestors).")

    # =================================================================
    # 5. TISSUE × CELL COVERAGE MATRIX
    # =================================================================

    def _compute_tissue_cell_matrix(self):
        """Compute top tissue × cell type marker count matrix for heatmap."""
        self.logger.info("Computing tissue-cell coverage matrix...")

        top_n_tissues = 20
        top_n_cells = 25

        tissue_counts = (self.df.groupby('tissue_name')['gene']
                         .nunique().nlargest(top_n_tissues))
        cell_counts = (self.df.groupby('cell_name')['gene']
                       .nunique().nlargest(top_n_cells))

        subset = self.df[
            self.df['tissue_name'].isin(tissue_counts.index) &
            self.df['cell_name'].isin(cell_counts.index)
        ]

        if subset.empty:
            self.tissue_cell_matrix = {
                'tissues': [], 'cell_types': [], 'values': []
            }
            return

        matrix = (subset.groupby(['tissue_name', 'cell_name'])['gene']
                  .nunique().unstack(fill_value=0))
        matrix = matrix.reindex(
            index=tissue_counts.index, columns=cell_counts.index, fill_value=0)

        self.tissue_cell_matrix = {
            'tissues': matrix.index.tolist(),
            'cell_types': matrix.columns.tolist(),
            'values': matrix.values.tolist()
        }
        self.logger.info(f"Coverage matrix: {len(matrix)} tissues × "
                         f"{len(matrix.columns)} cell types.")

    # =================================================================
    # 6. CELL MARKER EVIDENCE
    # =================================================================

    def _build_cell_marker_evidence(self):
        """Build top markers per cell type with evidence quality metrics."""
        self.logger.info("Building cell marker evidence data...")

        df = self.df[self.df['gene'] != 'Unknown Gene']

        # Aggregate: cell × gene → db_count, total_mentions
        agg = df.groupby(['cell_id', 'gene']).agg(
            db_count=('database', 'nunique'),
            total_mentions=('gene', 'size')
        ).reset_index()

        self.cell_markers = {}
        for cell_id, grp in agg.groupby('cell_id'):
            top = grp.nlargest(20, 'db_count')
            self.cell_markers[cell_id] = [
                {
                    "gene": row['gene'],
                    "db_count": int(row['db_count']),
                    "mentions": int(row['total_mentions'])
                }
                for _, row in top.iterrows()
            ]

        self.logger.info(f"Built marker evidence for "
                         f"{len(self.cell_markers)} cell types.")

    # =================================================================
    # 7. CELL → TISSUE SUMMARY
    # =================================================================

    def _build_cell_tissue_summary(self):
        """Build simplified cell → tissue marker count for bar charts."""
        self.logger.info("Building cell-tissue summary...")

        df = self.df[self.df['gene'] != 'Unknown Gene']

        ct_summary = (df.groupby(['cell_id', 'tissue_name'])['gene']
                      .nunique().reset_index(name='marker_count'))

        self.cell_tissue_summary = {}
        for cell_id, grp in ct_summary.groupby('cell_id'):
            tissues = (grp.sort_values('marker_count', ascending=False)
                       .head(20)
                       .apply(lambda r: {
                           "name": r['tissue_name'],
                           "count": int(r['marker_count'])
                       }, axis=1)
                       .tolist())
            self.cell_tissue_summary[cell_id] = tissues

        self.logger.info(f"Built tissue summaries for "
                         f"{len(self.cell_tissue_summary)} cell types.")

    # =================================================================
    # 8. GENE-CENTRIC DATA (optimized for JSON size)
    # =================================================================

    def _build_gene_centric_data(self):
        """Build gene-centric summary with size-optimized two-tier approach.

        All genes get bar chart data (top tissues/cells by source count).
        Genes supported by GENE_DETAIL_THRESHOLD+ databases also get
        source detail rows for the filterable table.
        """
        self.logger.info("Building gene-centric summary...")

        df = self.df[self.df['gene'] != 'Unknown Gene']

        # Pre-aggregate to unique combos
        agg = (df.groupby(['gene', 'tissue_name', 'cell_name',
                           'database', 'source_type'])
               .size().reset_index(name='count'))

        # Identify genes eligible for source detail table
        gene_db_counts = df.groupby('gene')['database'].nunique()
        detail_genes = set(
            gene_db_counts[gene_db_counts >= GENE_DETAIL_THRESHOLD].index)

        # All genes get bar chart data (single-DB genes still have tissue/cell info)
        bar_chart_genes = set(gene_db_counts.index)

        self.gene_to_summary = {}
        for gene, grp in agg.groupby('gene'):
            n_dbs = int(grp['database'].nunique())

            entry = {
                "total_tissues": int(grp['tissue_name'].nunique()),
                "total_cells": int(grp['cell_name'].nunique()),
                "total_databases": n_dbs,
            }

            # Bar chart data for genes with 2+ databases
            if gene in bar_chart_genes:
                t_agg = (grp.groupby('tissue_name')
                         .agg(mentions=('count', 'sum'),
                              source_count=('database', 'nunique'))
                         .reset_index()
                         .sort_values('source_count', ascending=False)
                         .head(GENE_TOP_N_BARS))
                c_agg = (grp.groupby('cell_name')
                         .agg(mentions=('count', 'sum'),
                              source_count=('database', 'nunique'))
                         .reset_index()
                         .sort_values('source_count', ascending=False)
                         .head(GENE_TOP_N_BARS))
                entry["tissue_plot"] = t_agg[['tissue_name', 'mentions',
                                              'source_count']].to_dict('records')
                entry["cell_plot"] = c_agg[['cell_name', 'mentions',
                                            'source_count']].to_dict('records')
                # Narrative summary fields: top tissue and cell by source count
                top_tissue = t_agg.iloc[0]['tissue_name'] if len(t_agg) > 0 else ''
                top_cell = c_agg.iloc[0]['cell_name'] if len(c_agg) > 0 else ''
                entry["top_tissue"] = top_tissue
                entry["top_cell"] = top_cell

            # Source detail rows for genes with THRESHOLD+ databases
            if gene in detail_genes:
                sources = (grp.nlargest(GENE_TOP_N_SOURCES, 'count')
                           [['tissue_name', 'cell_name', 'database',
                             'source_type']]
                           .to_dict('records'))
                entry["sources"] = sources

            self.gene_to_summary[gene] = entry

        self.logger.info(f"Built summaries for {len(self.gene_to_summary)} genes "
                         f"({len(detail_genes)} with source details).")

    # =================================================================
    # 9a. HIERARCHY TREE  (System → Tissue → Cell Type)
    # =================================================================

    def _build_hierarchy_tree(self):
        """Build 3-level hierarchy: organ system → tissue → cell type with marker counts.

        Used by Tab 1 Atlas zoomable sunburst. Each leaf carries:
          - marker_count  : unique genes
          - cell_count    : unique cell types
          - nlp_flag      : majority NLP tier for this tissue
          - top_cells     : top 5 cell types by marker count
          - db_count      : number of databases contributing
        """
        self.logger.info("Building hierarchy tree...")

        # Load CellxGene system list and tissue_descendants for system assignment
        _cxg_data = os.path.join(os.path.dirname(__file__), '..', '..',
            'SapiensOntoCellMap_env', 'lib', 'python3.10',
            'site-packages', 'cellxgene_ontology_guide', 'data')
        try:
            with open(os.path.join(_cxg_data, 'system_list.json')) as f:
                system_ids = set(json.load(f))
        except FileNotFoundError:
            system_ids = set()

        try:
            from cellxgene_ontology_guide.ontology_parser import OntologyParser
            onto_parser = OntologyParser()
        except Exception:
            onto_parser = None

        def _get_system_label(tissue_id):
            """Walk UBERON ancestors to find best system label."""
            if onto_parser is None:
                return 'Other'
            try:
                ancestors = set(onto_parser.get_term_ancestors(tissue_id))
            except Exception:
                ancestors = set()
            matching = ancestors & system_ids
            for sys_id in self._SYSTEM_PRIORITY:
                if sys_id in matching:
                    return onto_parser.get_term_label(sys_id).title()
            fb = self._SYSTEM_FALLBACK.get(tissue_id)
            if fb:
                return fb
            try:
                return onto_parser.get_term_label(tissue_id).title()
            except Exception:
                return 'Other'

        df = self.df[self.df['gene'] != 'Unknown Gene'].copy()

        # NLP quality map: tissue_name → dominant nlp_tissue_flag
        nlp_map = (df.groupby('tissue_name')['nlp_tissue_flag']
                   .agg(lambda x: x.mode().iloc[0] if len(x) > 0 else 'exact')
                   .to_dict())
        self.nlp_quality_map = nlp_map

        # Per-tissue aggregates
        tis_agg = df.groupby(['tissue_id', 'tissue_name']).agg(
            marker_count=('gene', 'nunique'),
            cell_count=('cell_name', 'nunique'),
            db_count=('database', 'nunique'),
        ).reset_index()

        # Top cell types per tissue
        top_cells_by_tissue = {}
        for (tid, tname), grp in df.groupby(['tissue_id', 'tissue_name']):
            top = (grp.groupby('cell_name')['gene']
                   .nunique()
                   .nlargest(5)
                   .reset_index())
            top_cells_by_tissue[(tid, tname)] = [
                {"name": r['cell_name'], "markers": int(r['gene'])}
                for _, r in top.iterrows()
            ]

        # Assign system label per tissue
        tree = {}  # system_label → { tissue_id → tissue_node }
        for _, row in tis_agg.iterrows():
            sys_label = _get_system_label(row['tissue_id'])
            if sys_label not in tree:
                tree[sys_label] = {}
            key = (row['tissue_id'], row['tissue_name'])
            tree[sys_label][row['tissue_name']] = {
                'tissue_id': row['tissue_id'],
                'tissue_name': row['tissue_name'],
                'marker_count': int(row['marker_count']),
                'cell_count': int(row['cell_count']),
                'db_count': int(row['db_count']),
                'nlp_flag': nlp_map.get(row['tissue_name'], 'exact'),
                'top_cells': top_cells_by_tissue.get(key, []),
            }

        # Flatten to sunburst-ready lists (ids, labels, parents, values, meta)
        ids, labels, parents, values, meta = [], [], [], [], []

        # Root node
        ids.append('root')
        labels.append('Human Body')
        parents.append('')
        values.append(0)
        meta.append({'type': 'root'})

        for sys_label, tissues in tree.items():
            sys_id = f'sys_{sys_label}'
            sys_marker_total = sum(t['marker_count'] for t in tissues.values())
            ids.append(sys_id)
            labels.append(sys_label)
            parents.append('root')
            values.append(sys_marker_total)
            meta.append({'type': 'system', 'tissue_count': len(tissues)})

            for tname, tnode in tissues.items():
                tis_id = f'tis_{tnode["tissue_id"]}'
                ids.append(tis_id)
                labels.append(tname)
                parents.append(sys_id)
                values.append(tnode['marker_count'])
                meta.append({
                    'type': 'tissue',
                    'tissue_id': tnode['tissue_id'],
                    'cell_count': tnode['cell_count'],
                    'db_count': tnode['db_count'],
                    'nlp_flag': tnode['nlp_flag'],
                    'top_cells': tnode['top_cells'],
                })

        self.hierarchy_tree = {
            'ids': ids, 'labels': labels, 'parents': parents,
            'values': values, 'meta': meta,
        }
        self.logger.info(
            f"Hierarchy tree: {len(tree)} systems, "
            f"{sum(len(v) for v in tree.values())} tissues.")

    # =================================================================
    # 9b. DOT MATRIX  (Tissue × Cell Type, with db_count)
    # =================================================================

    def _build_dot_matrix(self):
        """Build 30×40 tissue × cell type matrix for Evidence Matrix tab.

        Each cell carries:
          marker_count : unique genes
          db_count     : number of databases confirming
        Stored as parallel arrays for compact JSON.
        """
        self.logger.info("Building dot matrix (30 tissues × 40 cell types)...")

        TOP_TISSUES = 30
        TOP_CELLS = 40

        df = self.df[self.df['gene'] != 'Unknown Gene']

        # Exclude "All Tissues" — it's a catch-all not a real tissue
        df = df[df['tissue_name'] != 'All Tissues']

        tissue_ranks = (df.groupby('tissue_name')['gene']
                        .nunique().nlargest(TOP_TISSUES).index.tolist())
        cell_ranks = (df.groupby('cell_name')['gene']
                      .nunique().nlargest(TOP_CELLS).index.tolist())

        subset = df[
            df['tissue_name'].isin(tissue_ranks) &
            df['cell_name'].isin(cell_ranks)
        ]

        # marker_count matrix
        marker_mat = (subset.groupby(['tissue_name', 'cell_name'])['gene']
                      .nunique().unstack(fill_value=0)
                      .reindex(index=tissue_ranks, columns=cell_ranks, fill_value=0))

        # db_count matrix
        db_mat = (subset.groupby(['tissue_name', 'cell_name'])['database']
                  .nunique().unstack(fill_value=0)
                  .reindex(index=tissue_ranks, columns=cell_ranks, fill_value=0))

        self.dot_matrix = {
            'tissues': tissue_ranks,
            'cell_types': cell_ranks,
            'marker_values': marker_mat.values.tolist(),
            'db_values': db_mat.values.tolist(),
        }
        self.logger.info(
            f"Dot matrix: {len(tissue_ranks)} tissues × {len(cell_ranks)} cell types.")

    # =================================================================
    # 9c. SOURCE TYPE PER TISSUE
    # =================================================================

    def _build_source_type_per_tissue(self):
        """Per-tissue source_type breakdown for detail cards."""
        self.logger.info("Building source_type per tissue...")

        df = self.df[self.df['gene'] != 'Unknown Gene']
        agg = (df.groupby(['tissue_name', 'source_type'])
               .size().reset_index(name='count'))

        result = {}
        for tname, grp in agg.groupby('tissue_name'):
            total = grp['count'].sum()
            result[tname] = [
                {
                    'source_type': row['source_type'],
                    'count': int(row['count']),
                    'pct': round(row['count'] / total * 100, 1)
                }
                for _, row in grp.sort_values('count', ascending=False).iterrows()
            ]

        self.source_type_per_tissue = result
        self.logger.info(f"Source type per tissue: {len(result)} tissues.")

    # =================================================================
    # 9. ANATOMY DATA
    # =================================================================

    # Region → (primary_uberon_for_system_lookup, sex, uberon_parent_ids_for_DB_matching)
    # system_label is derived at runtime from CellxGene OntologyParser — not hardcoded.
    # sex: None=both, 'female'=female-only, 'male'=male-only
    # Emojis are visual metadata only — no ontology field exists for these.
    _REGION_DEF = {
        'region-brain':     ('UBERON:0001016', None,     ['UBERON:0001016', 'UBERON:0001017', 'UBERON:0000955', 'UBERON:0002240'], '🧠'),
        'region-eye':       ('UBERON:0000970', None,     ['UBERON:0000970'],                                                       '👁️'),
        'region-lung':      ('UBERON:0002048', None,     ['UBERON:0001004', 'UBERON:0002048'],                                     '🫁'),
        'region-heart':     ('UBERON:0000948', None,     ['UBERON:0004535', 'UBERON:0000948'],                                     '❤️'),
        'region-liver':     ('UBERON:0002107', None,     ['UBERON:0002107'],                                                       '🫀'),
        'region-stomach':   ('UBERON:0000945', None,     ['UBERON:0000945'],                                                       '🔵'),
        'region-intestine': ('UBERON:0000160', None,     ['UBERON:0001007', 'UBERON:0000160', 'UBERON:0001155', 'UBERON:0002108', 'UBERON:0000059'], '🌀'),
        'region-pancreas':  ('UBERON:0001264', None,     ['UBERON:0001264'],                                                       '🔶'),
        'region-kidney':    ('UBERON:0002113', None,     ['UBERON:0001008', 'UBERON:0002113', 'UBERON:0018707'],                   '🫘'),
        'region-skin':      ('UBERON:0002097', None,     ['UBERON:0002097'],                                                       '🧬'),
        'region-blood':     ('UBERON:0002390', None,     ['UBERON:0002390', 'UBERON:0000178', 'UBERON:0002371'],                   '🩸'),
        'region-lymph':     ('UBERON:0002405', None,     ['UBERON:0002405', 'UBERON:0000029', 'UBERON:0002106', 'UBERON:0002370'], '🛡️'),
        'region-breast':    ('UBERON:0000310', 'female', ['UBERON:0000310'],                                                       '🎀'),
        'region-adipose':   ('UBERON:0001013', None,     ['UBERON:0001013'],                                                       '💛'),
        'region-thyroid':   ('UBERON:0002046', None,     ['UBERON:0000949', 'UBERON:0002046', 'UBERON:0002369'],                   '⚗️'),
        'region-muscle':    ('UBERON:0000383', None,     ['UBERON:0000383', 'UBERON:0001434'],                                     '💪'),
        'region-uterus':    ('UBERON:0000995', 'female', ['UBERON:0000990', 'UBERON:0000995', 'UBERON:0000992', 'UBERON:0003889'], '♀️'),
        'region-placenta':  ('UBERON:0001987', 'female', ['UBERON:0001987'],                                                       '🤰'),
        'region-prostate':  ('UBERON:0002367', 'male',   ['UBERON:0002367'],                                                       '♂️'),
    }
    # Keep backward-compat alias
    _REGION_UBERON = {k: v[2] for k, v in _REGION_DEF.items()}
    # System-lookup priority: when an organ belongs to multiple systems (e.g. liver → digestive+endocrine+exocrine),
    # prefer the most biologically primary one. Order = highest priority first.
    # These IDs come directly from CellxGene system_list.json — resolved via OntologyParser at runtime.
    _SYSTEM_PRIORITY = [
        'UBERON:0001016',  # nervous system
        'UBERON:0001017',  # central nervous system
        'UBERON:0004535',  # cardiovascular system
        'UBERON:0001004',  # respiratory system
        'UBERON:0001007',  # digestive system
        'UBERON:0000949',  # endocrine system
        'UBERON:0002330',  # exocrine system
        'UBERON:0001008',  # renal system
        'UBERON:0000990',  # reproductive system
        'UBERON:0002405',  # immune system
        'UBERON:0002390',  # hematopoietic system
        'UBERON:0000383',  # musculature of body
        'UBERON:0001434',  # skeletal system
        'UBERON:0001032',  # sensory system
        'UBERON:0001009',  # circulatory system
        'UBERON:0000010',  # peripheral nervous system
    ]
    # For organs whose ancestors don't intersect system_list, use tissue_general fallback label
    _SYSTEM_FALLBACK = {
        'UBERON:0002097': 'Integumentary System',   # skin — not in UBERON system hierarchy
        'UBERON:0001013': 'Connective / Adipose',   # adipose tissue
        'UBERON:0000383': 'Musculoskeletal System', # muscle (is a system itself)
        'UBERON:0000310': 'Reproductive System',    # breast — no system ancestor in UBERON
        'UBERON:0001264': 'Digestive / Endocrine',  # pancreas — dual-system organ
        'UBERON:0001987': 'Reproductive System',    # placenta
        'UBERON:0002367': 'Reproductive System',    # prostate
    }

    # Legacy alias for _REGION_META (only the tuple layout changed)
    @property
    def _REGION_META(self):
        return {k: (None, v[3], v[1], v[2]) for k, v in self._REGION_DEF.items()}

    def _compute_anatomy_data(self):
        """Build per-region stats for anatomy SVG using CellxGene OntologyParser + tissue_descendants."""
        self.logger.info("Computing anatomy region data...")
        df = self.df

        # Load CellxGene resources
        _cxg_data = os.path.join(os.path.dirname(__file__), '..', '..',
            'SapiensOntoCellMap_env', 'lib', 'python3.10',
            'site-packages', 'cellxgene_ontology_guide', 'data')
        try:
            with open(os.path.join(_cxg_data, 'tissue_descendants.json')) as f:
                tissue_desc = json.load(f)
            with open(os.path.join(_cxg_data, 'system_list.json')) as f:
                system_ids = set(json.load(f))
            self.logger.info(f"Loaded tissue_descendants.json ({len(tissue_desc)} entries)")
        except FileNotFoundError:
            tissue_desc = {}
            system_ids = set()
            self.logger.warning("CellxGene data files not found — falling back to empty")

        # Load OntologyParser for label + ancestor resolution
        try:
            from cellxgene_ontology_guide.ontology_parser import OntologyParser
            onto_parser = OntologyParser()
        except Exception:
            onto_parser = None
            self.logger.warning("OntologyParser unavailable — system labels will be blank")

        def _norm(uid):
            return str(uid).split(' ')[0].strip()

        def _system_label(primary_uid):
            """Derive the most relevant system label for a given organ UBERON ID."""
            if onto_parser is None:
                return self._SYSTEM_FALLBACK.get(primary_uid, '')
            # If the primary_uid IS a system, return it directly
            if primary_uid in system_ids:
                return onto_parser.get_term_label(primary_uid).title()
            # Walk ancestors, intersect with system_list, pick by priority order
            try:
                ancestors = set(onto_parser.get_term_ancestors(primary_uid))
            except Exception:
                ancestors = set()
            matching = ancestors & system_ids
            for sys_id in self._SYSTEM_PRIORITY:
                if sys_id in matching:
                    return onto_parser.get_term_label(sys_id).title()
            # Fallback for organs not connected to a system in UBERON hierarchy
            return self._SYSTEM_FALLBACK.get(primary_uid, onto_parser.get_term_label(primary_uid).title())

        df = df.copy()
        df['tissue_id_norm'] = df['tissue_id'].apply(_norm)

        # Build tissue ID sets per region (parents + all CellxGene descendants)
        region_tissue_ids = {}
        for region, (primary_uid, sex, uberon_parents, emoji) in self._REGION_DEF.items():
            ids = set(uberon_parents)
            for parent in uberon_parents:
                for desc in tissue_desc.get(parent, []):
                    ids.add(_norm(desc))
            region_tissue_ids[region] = ids

        region_stats = {}
        for region, tissue_id_set in region_tissue_ids.items():
            primary_uid, sex, uberon_parents, emoji = self._REGION_DEF[region]
            sys_label = _system_label(primary_uid)
            sub = df[df['tissue_id_norm'].isin(tissue_id_set)]
            if sub.empty:
                region_stats[region] = {
                    "system": sys_label, "emoji": emoji, "sex": sex,
                    "cell_type_count": 0, "marker_count": 0,
                    "top_cell_types": [], "tissue_names": [], "tissue_breakdown": [],
                }
                continue
            top_cells = (sub.groupby("cell_name")["gene"]
                           .nunique()
                           .sort_values(ascending=False)
                           .head(12)
                           .index.tolist())
            tis_grp = (sub.groupby("tissue_name")["gene"]
                          .nunique()
                          .sort_values(ascending=False))
            tissue_breakdown = [
                {"name": name, "markers": int(cnt)}
                for name, cnt in tis_grp.head(20).items()
            ]
            region_stats[region] = {
                "system":           sys_label,
                "emoji":            emoji,
                "sex":              sex,
                "cell_type_count":  int(sub["cell_name"].nunique()),
                "marker_count":     int(sub["gene"].nunique()),
                "top_cell_types":   top_cells,
                "tissue_names":     [t["name"] for t in tissue_breakdown],
                "tissue_breakdown": tissue_breakdown,
            }

        self.anatomy_data = region_stats
        self.logger.info(f"Anatomy data built for {len(region_stats)} regions "
                         f"({sum(1 for v in region_stats.values() if v['cell_type_count'] > 0)} non-empty).")

    # =================================================================
    # 10. ASSEMBLE & EMBED
    # =================================================================

    def _assemble_and_embed(self):
        """Assemble final JSON and inject into HTML template."""
        self.logger.info("Assembling final data structure...")

        self.final_data = {
            "stats": self.dashboard_stats,
            "db_contribution": self.db_contribution,
            "source_types": self.source_type_dist,
            "cross_db_agreement": self.cross_db_agreement,
            "sunburst": self.sunburst_data,
            "coverage_matrix": self.tissue_cell_matrix,
            "cell_markers": self.cell_markers,
            "cell_tissue_summary": self.cell_tissue_summary,
            "cell_lineage_map": self.cell_lineage_map,
            "gene_to_summary": self.gene_to_summary,
            "all_genes_list": self.all_genes,
            "all_tissues_list": self.all_tissues_list,
            "all_cells_list": self.all_cells_list,
            "anatomy_data": self.anatomy_data,
            # New data structures for 5-tab architecture
            "hierarchy_tree": self.hierarchy_tree,
            "dot_matrix": self.dot_matrix,
            "nlp_quality_map": self.nlp_quality_map,
            "source_type_per_tissue": self.source_type_per_tissue,
        }

        json_str = json.dumps(self.final_data)
        json_mb = len(json_str) / (1024 * 1024)
        self.logger.info(f"JSON payload: {json_mb:.1f} MB")

        # Embed into template
        self.logger.info(f"Embedding data into {self.template_html_path}...")
        try:
            with open(self.template_html_path, 'r', encoding='utf-8') as f:
                template = f.read()

            placeholder = 'const embeddedData = "%%__DATA_PLACEHOLDER__%%";'
            if placeholder not in template:
                self.logger.error("Placeholder not found in HTML template.")
                return

            html = template.replace(
                placeholder, f'const embeddedData = {json_str};')
            with open(self.output_html_path, 'w', encoding='utf-8') as f:
                f.write(html)

            html_mb = len(html) / (1024 * 1024)
            self.logger.info(f"Generated: {self.output_html_path} ({html_mb:.1f} MB)")
        except Exception as e:
            self.logger.error(f"Embedding failed: {e}", exc_info=True)

    # =================================================================
    # 10. STATIC PUBLICATION PLOT
    # =================================================================

    def generate_summary_plot(self):
        """Generate 4-panel publication-ready static summary plot."""
        self.logger.info("Generating static summary plot...")

        # Aggregate
        agg_df = self.df.groupby('database').agg(
            total_markers=('gene', 'size'),
            unique_cell_ids=('cell_id', 'nunique'),
            unique_tissue_ids=('tissue_id', 'nunique')
        ).reset_index()

        source_comp = (self.df.groupby(['database', 'source_type'])
                       .size().reset_index(name='count'))

        order = agg_df.sort_values(
            by='total_markers', ascending=False)['database']

        # Cross-database agreement for panel C
        agreement = self.df.groupby(
            ['gene', 'cell_id'])['database'].nunique()
        agreement_hist = agreement.value_counts().sort_index()

        # Style
        sns.set_theme(style="ticks", context="paper", font_scale=1.4)
        plt.rcParams.update({
            "figure.dpi": 300,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "legend.title_fontsize": 11,
            "font.family": "sans-serif",
        })

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        palette = sns.color_palette("deep")

        def format_k(value, _):
            if value >= 1_000_000:
                return f'{value / 1_000_000:.1f}M'
            if value >= 1_000:
                return f'{int(value / 1000)}K'
            return f'{int(value)}'

        # Panel A: Markers per database
        sns.barplot(data=agg_df, x="database", y="total_markers",
                    order=order, ax=axes[0], color=palette[0])
        axes[0].set_title("Marker Entries per Database")
        axes[0].set_ylabel("Total Entries")
        axes[0].yaxis.set_major_formatter(FuncFormatter(format_k))

        # Panel B: Unique cell types per database
        sns.barplot(data=agg_df, x="database", y="unique_cell_ids",
                    order=order, ax=axes[1], color=palette[2])
        axes[1].set_title("Unique Cell Types per Database")
        axes[1].set_ylabel("Cell Type Count")

        # Panel C: Cross-database agreement (replaces tissue count)
        axes[2].bar(agreement_hist.index, agreement_hist.values,
                    color=palette[3], edgecolor='white', linewidth=0.5)
        axes[2].set_title("Cross-Database Marker Agreement")
        axes[2].set_xlabel("Number of Supporting Databases")
        axes[2].set_ylabel("Gene–Cell Type Pairs")
        axes[2].yaxis.set_major_formatter(FuncFormatter(format_k))
        multi_pct = round((agreement > 1).mean() * 100, 1)
        axes[2].annotate(f'{multi_pct}% supported by 2+ databases',
                         xy=(0.95, 0.92), xycoords='axes fraction',
                         ha='right', fontsize=10, color='#555555',
                         fontstyle='italic')

        # Panel D: Source type composition (stacked)
        pivot_src = source_comp.pivot(
            index='database', columns='source_type', values='count').fillna(0)
        pivot_src = pivot_src.reindex(order)
        unique_sources = source_comp['source_type'].unique()
        src_palette = sns.color_palette("tab10", len(unique_sources))
        color_map = dict(zip(unique_sources, src_palette))
        pivot_src.plot(
            kind='bar', stacked=True, ax=axes[3], width=0.8,
            color=[color_map.get(col) for col in pivot_src.columns])
        axes[3].set_title("Source Type Composition")
        axes[3].set_ylabel("Total Entries")
        axes[3].set_xlabel("Database")
        axes[3].yaxis.set_major_formatter(FuncFormatter(format_k))
        axes[3].legend(title="Source Type", bbox_to_anchor=(1.05, 1),
                       loc="upper left", frameon=False, fontsize=9)

        # Polish
        panel_labels = ['A', 'B', 'C', 'D']
        for i, ax in enumerate(axes):
            sns.despine(ax=ax)
            ax.grid(False)
            ax.text(-0.1, 1.05, panel_labels[i], transform=ax.transAxes,
                    fontsize=20, fontweight='bold', va='top', ha='right')

        for ax in [axes[0], axes[1], axes[3]]:
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        plt.tight_layout(rect=[0, 0, 0.9, 1])
        plt.savefig(PLOT_OUTPUT_FILENAME, dpi=300, bbox_inches='tight',
                    facecolor='white', transparent=False)
        self.logger.info(f"Saved static plot: '{PLOT_OUTPUT_FILENAME}'")
        plt.close(fig)

    # =================================================================
    # 11. RUN PIPELINE
    # =================================================================

    def run(self):
        """Execute the full generation pipeline."""
        self.logger.info("=== SapiensOntoCellMap Generator START ===")

        if not self._load_and_clean_data():
            return

        # Static publication plot
        self.generate_summary_plot()

        # Interactive visualizer data
        self._compute_dashboard_data()
        self._build_cell_lineages()
        self._build_ontology_sunburst()
        self._compute_tissue_cell_matrix()
        self._build_cell_marker_evidence()
        self._build_cell_tissue_summary()
        self._build_gene_centric_data()
        self._build_hierarchy_tree()
        self._build_dot_matrix()
        self._build_source_type_per_tissue()
        self._compute_anatomy_data()
        self._assemble_and_embed()

        self.logger.info("=== SapiensOntoCellMap Generator END ===")


if __name__ == "__main__":
    script_directory = os.path.dirname(os.path.abspath(__file__))
    template_file = os.path.join(script_directory, 'visualizer_template.html')

    try:
        generator = SapiensMapGenerator(
            input_csv=PROCESSED_COMBINED_DATABASE_FILE,
            template_html=template_file,
            output_html=PROCESSED_COMBINED_DATABASE_FILE_HTML
        )
        generator.run()
    except NameError:
        print("Configuration Error: Config constants not found.")
    except Exception as e:
        print(f"Initialization Error: {e}")
        import traceback
        traceback.print_exc()
