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

        # Global gene list
        all_genes_set = set(df['gene'].unique()) - {'Unknown Gene'}
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

        # Genes with only 1 DB get summary stats only (saves ~50% of JSON)
        bar_chart_genes = set(
            gene_db_counts[gene_db_counts >= 2].index)

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
    # 9. ASSEMBLE & EMBED
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
