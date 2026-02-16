#!/usr/bin/env python3

"""
SapiensOntoCellMap Data Processor & Embedder (v8.5)

Key Features:
1. **HTML VISUALIZER**: Generates the interactive SapiensOntoCellMap.
2. **STATIC PLOTS**: Generates the "database_summary_2x2.png" publication figure.
3. **DATA FIXES**: Includes Hierarchy Pruning, Source Counts, and Cell ID Aggregation.
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
from typing import Dict, Any, Set, List
from collections import defaultdict

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

# --- Configuration for Plots ---
PLOT_OUTPUT_FILENAME = "database_summary_2x2_publishable_v4.png"

class SapiensMapGenerator:
    """
    Encapsulates the pipeline for generating the SapiensOntoCellMap visualizer
    AND the static publication summary plots.
    """
    
    LOG_LEVEL = logging.INFO

    def __init__(self, input_csv: str, template_html: str, output_html: str):
        self.input_csv_path = input_csv
        self.template_html_path = template_html
        self.output_html_path = output_html
        
        # Initialize Logging
        logging.basicConfig(level=self.LOG_LEVEL, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

        # Data Containers
        self.df: pd.DataFrame = pd.DataFrame()
        self.all_genes: List[str] = []
        self.all_tissues_list: List[Dict] = []
        self.all_cells_list: List[Dict] = []
        self.cell_lineage_map: Dict[str, Any] = {}
        self.gene_to_summary: Dict[str, Any] = {} 
        self.final_data: Dict[str, Any] = {}

        # Initialize Ontology Parser
        try:
            self.ontology = CellxGeneOntologyParser()
            self.logger.info("CellxGeneOntologyParser initialized successfully.")
        except Exception as e:
            self.logger.error(f"Failed to initialize CellxGeneOntologyParser: {e}", exc_info=True)
            raise

    def _load_and_clean_data(self) -> bool:
        """
        Loads raw CSV, cleans IDs, and populates helper lists.
        """
        self.logger.info(f"Loading raw data from {self.input_csv_path}...")
        try:
            # Low_memory=False suppresses DtypeWarning on large files
            df = pd.read_csv(self.input_csv_path, low_memory=False)
        except FileNotFoundError:
            self.logger.error(f"Input file not found: {self.input_csv_path}")
            return False

        self.logger.info(f"Loaded {len(df)} total rows.")
        
        # Initial cleanup
        df.dropna(subset=['tissue_id', 'cell_id'], inplace=True)
        
        # 1. Standardize separators (CL_123 -> CL:123)
        df['cell_id'] = df['cell_id'].astype(str).str.replace('_', ':')
        df['tissue_id'] = df['tissue_id'].astype(str).str.replace('_', ':')
        
        # 2. Extract valid ontology patterns only
        df['cell_id'] = df['cell_id'].str.extract(r'(CL:\d+)', expand=False)
        df['tissue_id'] = df['tissue_id'].str.extract(r'(UBERON:\d+)', expand=False)
        
        # 3. Drop rows that became invalid
        df.dropna(subset=['tissue_id', 'cell_id'], inplace=True)
        
        # Fill Metadata NaNs
        df['tissue_name'] = df['tissue_name'].fillna(df['db_tissue_name']).fillna('Unknown Tissue')
        df['cell_name'] = df['cell_name'].fillna(df['db_cell_name']).fillna('Unknown Cell')
        df['gene'] = df['gene'].fillna('Unknown Gene').astype(str)
        df['database'] = df['database'].fillna('N/A').astype(str)
        df['source_type'] = df['source_type'].fillna('N/A').astype(str)
        df['source_info'] = df['source_info'].fillna('N/A').astype(str)
        
        # Populate Global Lists
        all_genes_set = set(df['gene'].unique())
        if 'Unknown Gene' in all_genes_set:
            all_genes_set.remove('Unknown Gene')
        self.all_genes = sorted(list(all_genes_set))
            
        tissue_df = df[['tissue_id', 'tissue_name']].drop_duplicates().sort_values(by='tissue_name')
        self.all_tissues_list = tissue_df.to_dict('records')
        
        cell_name_map = df.groupby('cell_id')['cell_name'].apply(lambda x: x.mode().iloc[0]).to_dict()
        self.all_cells_list = sorted(
            [{'cell_id': cid, 'cell_name': cname} for cid, cname in cell_name_map.items()],
            key=lambda x: x['cell_name']
        )
        
        self.logger.info(f"Cleaned data: {len(df)} rows retained.")
        self.df = df
        return True

    def generate_summary_plot(self):
        """
        Generates the 4-panel publication-ready static plot.
        Uses the internal self.df so no reload is necessary.
        """
        self.logger.info("📊 Generating static summary plots...")
        
        # --- 1. Aggregate data for plots ---
        # Note: We use self.df directly
        agg_df = self.df.groupby('database').agg(
            total_markers=('gene', 'size'),
            unique_cell_ids=('cell_id', 'nunique'),
            unique_tissue_ids=('tissue_id', 'nunique')
        ).reset_index()

        source_comp = (
            self.df.groupby(['database', 'source_type'])
            .size()
            .reset_index(name='count')
        )

        order = agg_df.sort_values(by='total_markers', ascending=False)['database']

        # --- 2. Set style and figure layout ---
        sns.set_theme(style="ticks", context="paper", font_scale=1.4)
        plt.rcParams.update({
            "figure.dpi": 300,
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11,
            "legend.title_fontsize": 12
        })

        fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True)
        axes = axes.flatten()

        # --- Helper Functions (Inner Scope) ---
        def add_labels(ax):
            for p in ax.patches:
                height = p.get_height()
                if height > 0:
                    ax.annotate(f"{int(height):,}",
                        xy=(p.get_x() + p.get_width() / 2, height),
                        xytext=(0, 5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=9, fontweight='bold')

        def format_k(value, tick_number):
            if value >= 1_000_000: return f'{value / 1_000_000:.1f}M'
            if value >= 1_000: return f'{int(value / 1000)}K'
            return f'{int(value)}'

        # --- 3. Plot 1: Number of Markers ---
        sns.barplot(data=agg_df, x="database", y="total_markers", order=order, ax=axes[0], color="#4C72B0")
        axes[0].set_title("Number of Markers")
        axes[0].set_ylabel("Total Markers")
        axes[0].yaxis.set_major_formatter(FuncFormatter(format_k))
        add_labels(axes[0])

        # --- 4. Plot 2: Unique Cell IDs ---
        sns.barplot(data=agg_df, x="database", y="unique_cell_ids", order=order, ax=axes[1], color="#55A868")
        axes[1].set_title("Unique Cell Types")
        axes[1].set_ylabel("Unique Cell ID Count")
        add_labels(axes[1])

        # --- 5. Plot 3: Unique Tissue IDs ---
        sns.barplot(data=agg_df, x="database", y="unique_tissue_ids", order=order, ax=axes[2], color="#C44E52")
        axes[2].set_title("Unique Tissues")
        axes[2].set_ylabel("Unique Tissue ID Count")
        axes[2].set_xlabel("Database")
        add_labels(axes[2])

        # --- 6. Plot 4: Source type per database (stacked bar) ---
        pivot_src = source_comp.pivot(index='database', columns='source_type', values='count').fillna(0)
        pivot_src = pivot_src.reindex(order)

        unique_sources = source_comp['source_type'].unique()
        palette = sns.color_palette("tab10", len(unique_sources))
        color_map = {source: color for source, color in zip(unique_sources, palette)}

        pivot_src.plot(
            kind='bar', stacked=True, ax=axes[3], width=0.8,
            color=[color_map.get(col) for col in pivot_src.columns]
        )
        axes[3].set_title("Source Type Composition")
        axes[3].set_ylabel("Total Markers")
        axes[3].set_xlabel("Database")
        axes[3].yaxis.set_major_formatter(FuncFormatter(format_k))
        axes[3].legend(title="Source Type", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

        # --- 7. Final Polish ---
        panel_labels = ['A', 'B', 'C', 'D']
        for i, ax in enumerate(axes):
            sns.despine(ax=ax)
            ax.grid(False)
            ax.text(-0.1, 1.05, panel_labels[i], transform=ax.transAxes,
                    fontsize=20, fontweight='bold', va='top', ha='right')

        for ax in [axes[2], axes[3]]:
            labels = ax.get_xticklabels()
            ax.set_xticklabels(labels, rotation=45, ha='right')

        plt.setp(axes[0].get_xticklabels(), visible=False)
        plt.setp(axes[1].get_xticklabels(), visible=False)

        plt.tight_layout(rect=[0, 0, 0.9, 1])

        plt.savefig(PLOT_OUTPUT_FILENAME, dpi=300, bbox_inches='tight', facecolor='white', transparent=False)
        self.logger.info(f"✅ Saved static plot: '{PLOT_OUTPUT_FILENAME}'")
        # plt.show() # Commented out for pipeline usage

    def _simulate_umap_coordinates(self):
        """Simulates UMAP coordinates for tissues and cells."""
        self.logger.info("Simulating UMAP coordinates...")
        unique_tissues = self.df['tissue_id'].unique()
        tissue_coords = pd.DataFrame({
            'tissue_id': unique_tissues,
            'tissue_umap_x': np.random.randn(len(unique_tissues)),
            'tissue_umap_y': np.random.randn(len(unique_tissues))
        })
        unique_cells_in_tissue = self.df[['tissue_id', 'cell_id']].drop_duplicates()
        cell_coords = pd.DataFrame({
            'tissue_id': unique_cells_in_tissue['tissue_id'],
            'cell_id': unique_cells_in_tissue['cell_id'],
            'cell_umap_x': np.random.randn(len(unique_cells_in_tissue)),
            'cell_umap_y': np.random.randn(len(unique_cells_in_tissue))
        })
        self.df = pd.merge(self.df, tissue_coords, on='tissue_id', how='left')
        self.df = pd.merge(self.df, cell_coords, on=['tissue_id', 'cell_id'], how='left')

    def _build_cell_lineages(self):
        """Builds lineage maps, pruning ancestors above 'cell' and deprecated terms."""
        self.logger.info("Building ancestor lineages...")
        lineage_map = {}
        unique_cell_ids = self.df['cell_id'].unique()

        for cell_id in unique_cell_ids:
            if cell_id not in self.ontology.cl_id_to_name: continue
            if self.ontology.ontology_parser.is_term_deprecated(cell_id): continue
            
            try:
                ancestors_with_dist = self.ontology.ontology_parser.get_term_ancestors_with_distances(cell_id)
                ancestors_with_dist[cell_id] = 0
                
                # Pruning logic
                root_cell_id = 'CL:0000000'
                if root_cell_id in ancestors_with_dist:
                    root_dist = ancestors_with_dist[root_cell_id]
                    ancestors_with_dist = {a: d for a, d in ancestors_with_dist.items() if d <= root_dist}

                lineage_for_cell = {}
                for ancestor_id, dist in ancestors_with_dist.items():
                    if self.ontology.ontology_parser.is_term_deprecated(ancestor_id): continue
                    ancestor_name = self.ontology.cl_id_to_name.get(ancestor_id, ancestor_id)
                    lineage_for_cell[ancestor_id] = {'dist': dist, 'name': ancestor_name}
                
                lineage_map[cell_id] = lineage_for_cell
            except Exception: pass

        self.cell_lineage_map = lineage_map
        
    def _build_gene_centric_data(self):
        """Builds gene_to_summary for Tab 3 with source counts."""
        self.logger.info("Building gene-centric summary...")
        self.gene_to_summary = {}
        
        for gene, gene_group in self.df.groupby('gene'):
            if gene == 'Unknown Gene': continue
            
            gene_data = {
                "tissues": defaultdict(lambda: {'mentions': 0, 'databases': set()}),
                "cells": defaultdict(lambda: {'mentions': 0, 'databases': set()}),
                "sources": []
            }
            
            for _, row in gene_group.iterrows():
                gene_data['tissues'][row['tissue_name']]['mentions'] += 1
                gene_data['tissues'][row['tissue_name']]['databases'].add(row['database'])
                gene_data['cells'][row['cell_name']]['mentions'] += 1
                gene_data['cells'][row['cell_name']]['databases'].add(row['database'])
                
                gene_data['sources'].append({
                    "tissue_name": row['tissue_name'], "cell_name": row['cell_name'],
                    "database": row['database'], "source_type": row['source_type'],
                    "source_info": row['source_info']
                })
            
            # Sort for plots
            for key in ['tissues', 'cells']:
                plot_key = 'tissue_plot' if key == 'tissues' else 'cell_plot'
                gene_data[plot_key] = sorted(
                    [{"name": k, "mentions": v['mentions'], "source_count": len(v['databases'])} 
                     for k, v in gene_data[key].items()],
                    key=lambda x: x['source_count'], reverse=True
                )
            del gene_data['tissues']
            del gene_data['cells']
            self.gene_to_summary[gene] = gene_data

    def _aggregate_data_to_json(self):
        """Aggregates data grouping by cell_id."""
        self.logger.info("Aggregating data for visualizer...")
        final_data = {
            "tissues_umap": [], "cell_to_tissue_summary": {}, 
            "gene_to_summary": self.gene_to_summary,
            "cell_lineage_map": self.cell_lineage_map, 
            "all_genes_list": self.all_genes,
            "all_tissues_list": self.all_tissues_list,
            "all_cells_list": self.all_cells_list
        }
        
        cell_centric_agg = {} 
        for tissue_id, tissue_group in self.df.groupby('tissue_id'):
            tissue_name = tissue_group['tissue_name'].mode().iloc[0]
            tissue_obj = {
                "tissue_id": tissue_id, "tissue_name": tissue_name,
                "umap_x": tissue_group['tissue_umap_x'].iloc[0],
                "umap_y": tissue_group['tissue_umap_y'].iloc[0],
                "cells": [], "features_in_tissue": list(tissue_group['gene'].unique()),
            }
            
            for cell_id, cell_group in tissue_group.groupby('cell_id'):
                cell_name = cell_group['cell_name'].mode().iloc[0]
                cell_obj = {
                    "cell_id": cell_id, "cell_name": cell_name,
                    "umap_x": cell_group['cell_umap_x'].iloc[0],
                    "umap_y": cell_group['cell_umap_y'].iloc[0], "features": []
                }
                
                if cell_id not in cell_centric_agg: cell_centric_agg[cell_id] = {}
                cell_centric_agg[cell_id][tissue_id] = {"tissue_name": tissue_name, "markers": set()}

                for gene, feature_group in cell_group.groupby('gene'):
                    if gene == 'Unknown Gene': continue
                    cell_centric_agg[cell_id][tissue_id]["markers"].add(gene)
                    
                    feature_obj = { "name": gene, "sources": [] }
                    for _, row in feature_group.iterrows():
                        feature_obj["sources"].append({
                            "database": row['database'], "source_type": row['source_type'],
                            "source_info": row['source_info']
                        })
                    cell_obj["features"].append(feature_obj)
                tissue_obj["cells"].append(cell_obj)
            
            tissue_obj["cell_count"] = len(tissue_obj["cells"])
            tissue_obj["feature_count"] = len(tissue_obj["features_in_tissue"])
            final_data["tissues_umap"].append(tissue_obj)
        
        for cell_id, tissues in cell_centric_agg.items():
            for tissue_id, data in tissues.items():
                data['markers'] = sorted(list(data['markers']))
        final_data['cell_to_tissue_summary'] = cell_centric_agg
        self.final_data = final_data

    def _embed_data_in_template(self):
        """Injects JSON into HTML."""
        self.logger.info(f"Embedding data into {self.template_html_path}...")
        try:
            with open(self.template_html_path, 'r', encoding='utf-8') as f: template_content = f.read()
            json_data_string = json.dumps(self.final_data)
            
            placeholder_line = 'const embeddedData = "%%__DATA_PLACEHOLDER__%%";'
            if placeholder_line not in template_content:
                self.logger.error("Placeholder not found in HTML template.")
                return
            
            final_html = template_content.replace(placeholder_line, f'const embeddedData = {json_data_string};')
            with open(self.output_html_path, 'w', encoding='utf-8') as f: f.write(final_html)
            self.logger.info(f"✅ Generated Visualizer: {self.output_html_path}")
        except Exception as e:
            self.logger.error(f"Embedding failed: {e}")

    def run(self):
        """Executes full pipeline."""
        self.logger.info("--- SapiensOntoCellMap Generator START ---")
        if not self._load_and_clean_data(): return

        # 1. Generate Static Plots (Using loaded data)
        self.generate_summary_plot()

        # 2. Process Data for Visualizer
        self._simulate_umap_coordinates()
        self._build_cell_lineages()
        self._build_gene_centric_data()
        self._aggregate_data_to_json()
        self._embed_data_in_template()
        
        self.logger.info("--- SapiensOntoCellMap Generator END ---")

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