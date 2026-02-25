#!/usr/bin/env python3
"""
Author : Sonal Rashmi
Date   : 2026-02-17
Description :
Hierarchical Annotation Engine for SapiensOntoCellMap.

After flat enrichment produces significant hits for each cluster, this engine
traverses the CL (Cell Ontology) to produce multi-resolution annotations with
confidence scoring at each hierarchy depth.

Key differentiator: No existing cell annotation tool performs ontology-aware
hierarchical marker enrichment analysis.
"""

import logging
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


class HierarchicalAnnotator:
    """
    Produces multi-resolution annotations by traversing the CL ontology graph
    and aggregating enrichment evidence at each hierarchy level.

    For each cluster:
    1. Maps significant cell type hits to CL IDs
    2. Gets ancestor chains for each CL ID
    3. Aggregates evidence at each ancestor node
    4. Computes confidence at each depth
    5. Recommends optimal resolution
    """

    def __init__(self, ontology_parser, master_db_df, confidence_threshold=0.5):
        """
        Args:
            ontology_parser: CellxGeneOntologyParser instance with
                             get_cell_ontology_graph() and
                             get_term_ancestors_with_distances() methods
            master_db_df: The master marker DB DataFrame (used to build
                          cell_name → cell_id mapping)
            confidence_threshold: Minimum confidence to recommend a resolution
                                  level (default: 0.5)
        """
        self.ontology_parser = ontology_parser
        self.confidence_threshold = confidence_threshold

        # Build cell_name → cell_id mapping from master DB
        self.cell_name_to_id = self._build_cell_name_to_id_map(master_db_df)
        logger.info(f"Built cell_name→CL_ID map with {len(self.cell_name_to_id)} entries")

        # Build CL ID → name mapping from ontology
        self.cl_id_to_name = dict(ontology_parser.cl_id_to_name)

        # Build children map (parent → set of children) for counting descendants
        self._children_map = defaultdict(set)
        self._ancestor_cache = {}
        self._build_children_map()

    def _build_cell_name_to_id_map(self, master_db_df):
        """Build {cell_name_upper: cell_id} from the master DB."""
        name_to_id = {}
        if 'cell_name' in master_db_df.columns and 'cell_id' in master_db_df.columns:
            pairs = master_db_df[['cell_name', 'cell_id']].dropna().drop_duplicates()
            for _, row in pairs.iterrows():
                name = str(row['cell_name']).strip()
                cid = str(row['cell_id']).strip()
                if name and cid and cid.startswith('CL:'):
                    name_to_id[name.upper()] = cid
        return name_to_id

    def _build_children_map(self):
        """Build parent→children relationships from ontology ancestor queries.

        Obsolete CL terms are excluded from both the child and parent roles.
        CL retires deprecated cell type concepts (e.g. CL:0000215 'obsolete
        barrier cell') whose names begin with 'obsolete' in the cl_id_to_name
        map. Including them would propagate deprecated biology into the ancestor
        tree and produce misleading broad/intermediate resolution labels.
        """
        for cl_id in self.cl_id_to_name:
            # Skip obsolete CL terms as children
            if self.cl_id_to_name.get(cl_id, '').lower().startswith('obsolete'):
                continue
            ancestors = self._get_ancestors(cl_id)
            if ancestors:
                # Direct parents (distance=1)
                for anc_id, dist in ancestors.items():
                    if dist == 1:
                        # Skip obsolete CL terms as parents
                        if self.cl_id_to_name.get(anc_id, '').lower().startswith('obsolete'):
                            continue
                        self._children_map[anc_id].add(cl_id)

    def _get_ancestors(self, cl_id):
        """Get ancestors with distances, with caching."""
        if cl_id in self._ancestor_cache:
            return self._ancestor_cache[cl_id]
        try:
            ancestors = self.ontology_parser.ontology_parser.get_term_ancestors_with_distances(cl_id)
            self._ancestor_cache[cl_id] = ancestors
            return ancestors
        except Exception:
            self._ancestor_cache[cl_id] = {}
            return {}

    def _get_all_descendants(self, cl_id, max_depth=10):
        """Get all descendants of a CL ID by traversing children map."""
        descendants = set()
        frontier = {cl_id}
        for _ in range(max_depth):
            next_frontier = set()
            for node in frontier:
                children = self._children_map.get(node, set())
                new_children = children - descendants - {cl_id}
                next_frontier.update(new_children)
                descendants.update(new_children)
            if not next_frontier:
                break
            frontier = next_frontier
        return descendants

    def _get_db_descendants(self, cl_id):
        """
        Get descendants of cl_id that exist in our marker database
        (i.e., have entries in cell_name_to_id).
        """
        all_descendants = self._get_all_descendants(cl_id)
        db_cl_ids = set(self.cell_name_to_id.values())
        return all_descendants.intersection(db_cl_ids)

    def annotate_cluster(self, cluster_id, cluster_sig_results_df):
        """
        Produce hierarchical annotations for a single cluster.

        Args:
            cluster_id: e.g., "Cluster 0"
            cluster_sig_results_df: DataFrame of significant results for this
                                    cluster (must have Cell_type, adj_p_value columns)

        Returns:
            List of dicts, one per hierarchy level, sorted from broad to fine:
            [
                {'Cluster': ..., 'Depth': 0, 'CL_ID': ..., 'Cell_Type': ...,
                 'N_Supporting': ..., 'Combined_Score': ..., 'Confidence': ...,
                 'Resolution': ..., 'Supporting_Types': ...},
                ...
            ]
        """
        if cluster_sig_results_df.empty:
            return []

        # Step 1: Map cell type names to CL IDs
        sig_cl_ids = {}  # cl_id → (cell_type_name, adj_p_value)
        for _, row in cluster_sig_results_df.iterrows():
            cell_type = str(row['Cell_type']).strip()
            p_val = float(row['adj_p_value'])
            cl_id = self.cell_name_to_id.get(cell_type.upper())
            if cl_id:
                # Keep best p-value if multiple entries for same CL ID
                if cl_id not in sig_cl_ids or p_val < sig_cl_ids[cl_id][1]:
                    sig_cl_ids[cl_id] = (cell_type, p_val)

        if not sig_cl_ids:
            return []

        # Step 2: Get ancestor chains and aggregate evidence
        ancestor_evidence = defaultdict(lambda: {
            'supporting_ids': set(),
            'supporting_names': set(),
            'combined_score': 0.0,
            'min_depth': float('inf'),
        })

        for cl_id, (cell_name, p_val) in sig_cl_ids.items():
            # The significant cell type itself counts
            score = -np.log10(max(p_val, 1e-300))
            ancestor_evidence[cl_id]['supporting_ids'].add(cl_id)
            ancestor_evidence[cl_id]['supporting_names'].add(cell_name)
            ancestor_evidence[cl_id]['combined_score'] += score
            ancestor_evidence[cl_id]['min_depth'] = 0  # relative depth

            # Walk up ancestors
            ancestors = self._get_ancestors(cl_id)
            for anc_id, dist in ancestors.items():
                if anc_id == cl_id:
                    continue
                # Skip very generic root terms
                if anc_id in ('CL:0000000', 'CL:0000001'):
                    continue
                # Skip obsolete CL terms — retired concepts should not appear
                # as resolution levels in annotations (e.g. CL:0000493
                # 'obsolete regulatory T cell', CL:0000215 'obsolete barrier cell').
                if self.cl_id_to_name.get(anc_id, '').lower().startswith('obsolete'):
                    continue
                ancestor_evidence[anc_id]['supporting_ids'].add(cl_id)
                ancestor_evidence[anc_id]['supporting_names'].add(cell_name)
                ancestor_evidence[anc_id]['combined_score'] += score

        # Step 3: Compute confidence and depth for each ancestor node
        results = []
        for node_id, evidence in ancestor_evidence.items():
            node_name = self.cl_id_to_name.get(node_id, node_id)
            n_supporting = len(evidence['supporting_ids'])

            # N_possible = descendants of this node that exist in our DB
            db_descendants = self._get_db_descendants(node_id)
            # Include the node itself if it's in the DB
            if node_id in set(self.cell_name_to_id.values()):
                db_descendants.add(node_id)

            n_possible = max(len(db_descendants), n_supporting)  # at least as many as supporting
            confidence = n_supporting / n_possible if n_possible > 0 else 0.0

            # Compute depth: distance from root concept to this node
            ancestors = self._get_ancestors(node_id)
            # Depth is max ancestor distance (how far from root)
            if ancestors:
                depth = max(ancestors.values())
            else:
                depth = 0

            results.append({
                'Cluster': cluster_id,
                'Depth': depth,
                'CL_ID': node_id,
                'Cell_Type': node_name,
                'N_Supporting': n_supporting,
                'Combined_Score': round(evidence['combined_score'], 2),
                'Confidence': round(confidence, 3),
                'Supporting_Types': ', '.join(sorted(evidence['supporting_names'])),
            })

        if not results:
            return []

        # Sort by depth (broad to fine), then by combined score (desc)
        results.sort(key=lambda x: (x['Depth'], -x['Combined_Score']))

        # Step 4: Assign resolution labels
        # Walk from broadest to finest, label based on confidence
        for r in results:
            if r['Confidence'] >= self.confidence_threshold:
                if r['Depth'] <= 3:
                    r['Resolution'] = 'broad'
                elif r['Depth'] <= 6:
                    r['Resolution'] = 'intermediate'
                else:
                    r['Resolution'] = 'fine'
            else:
                r['Resolution'] = 'uncertain'

        # Filter: keep only nodes with >=2 supporting types OR nodes that ARE
        # a significant hit themselves (leaf nodes)
        filtered = []
        sig_cl_id_set = set(sig_cl_ids.keys())
        for r in results:
            if r['N_Supporting'] >= 2 or r['CL_ID'] in sig_cl_id_set:
                filtered.append(r)

        return filtered

    def annotate_all_clusters(self, sig_results_df):
        """
        Run hierarchical annotation for all clusters in the significant results.

        Args:
            sig_results_df: DataFrame with columns [Cluster, Cell_type, adj_p_value, ...]

        Returns:
            DataFrame with hierarchical annotations for all clusters
        """
        if sig_results_df.empty:
            return pd.DataFrame()

        all_annotations = []
        clusters = sig_results_df['Cluster'].unique()
        logger.info(f"Running hierarchical annotation for {len(clusters)} clusters...")

        for cluster_id in sorted(clusters, key=lambda s: [int(t) if t.isdigit() else t.lower() for t in __import__('re').split(r'(\d+)', str(s))]):
            cluster_df = sig_results_df[sig_results_df['Cluster'] == cluster_id]
            annotations = self.annotate_cluster(cluster_id, cluster_df)
            all_annotations.extend(annotations)

        if not all_annotations:
            return pd.DataFrame()

        result_df = pd.DataFrame(all_annotations)
        col_order = ['Cluster', 'Depth', 'CL_ID', 'Cell_Type', 'N_Supporting',
                      'Combined_Score', 'Confidence', 'Resolution', 'Supporting_Types']
        result_df = result_df[[c for c in col_order if c in result_df.columns]]

        logger.info(f"Hierarchical annotation complete: {len(result_df)} entries across {len(clusters)} clusters")
        return result_df

    def get_best_resolution(self, cluster_annotations_df, cluster_id):
        """
        Get the recommended annotation at the best confidence-supported resolution.

        Returns the finest-grained annotation that still meets the confidence threshold.

        Args:
            cluster_annotations_df: Output of annotate_all_clusters
            cluster_id: e.g., "Cluster 0"

        Returns:
            Dict with best annotation, or None if no confident annotation exists
        """
        if cluster_annotations_df.empty:
            return None

        cluster_df = cluster_annotations_df[cluster_annotations_df['Cluster'] == cluster_id]
        confident = cluster_df[cluster_df['Confidence'] >= self.confidence_threshold]

        if confident.empty:
            return None

        # Return the deepest (most specific) confident annotation
        best = confident.loc[confident['Depth'].idxmax()]
        return best.to_dict()

    def get_broad_type(self, cluster_annotations_df, cluster_id, top_cell_type=None):
        """
        Get the broadest confident ancestor type for a cluster.

        When top_cell_type is provided, the search is constrained to ancestors
        of that cell type's CL ID. This prevents false-positive broad labels
        caused by unrelated significant hits sharing housekeeping genes
        (e.g., keratinocyte subtypes contaminating a T cell cluster's Broad_Type).

        The deepest (most specific) qualifying ancestor in the top annotation's
        lineage is returned. Falls back to the shallowest confident ancestor
        across all results if no lineage ancestor is found.

        Args:
            cluster_annotations_df: Output of annotate_all_clusters
            cluster_id: e.g., "Cluster 0"
            top_cell_type: Cell_type name string of the top-ranked annotation
                           (used to constrain the ancestor search)

        Returns:
            Cell type name string, or None
        """
        if cluster_annotations_df.empty:
            return None

        cluster_df = cluster_annotations_df[cluster_annotations_df['Cluster'] == cluster_id]
        # Candidate ancestors: confident nodes with true aggregation support
        broad = cluster_df[
            (cluster_df['Confidence'] >= self.confidence_threshold) &
            (cluster_df['N_Supporting'] >= 2)
        ]

        if broad.empty:
            return None

        # --- Constrained search: only ancestors of top annotation's lineage ---
        if top_cell_type:
            top_cl_id = self.cell_name_to_id.get(str(top_cell_type).strip().upper())
            if top_cl_id:
                # Get all ancestors of the top annotation (dict cl_id → distance)
                top_ancestors = self._get_ancestors(top_cl_id)
                ancestor_cl_ids = set(top_ancestors.keys())
                # Exclude the top annotation itself, root cells, and terms too
                # shallow to be informative (depth < 2 from root)
                ancestor_cl_ids.discard(top_cl_id)
                ancestor_cl_ids.discard('CL:0000000')
                ancestor_cl_ids.discard('CL:0000001')

                lineage_broad = broad[broad['CL_ID'].isin(ancestor_cl_ids)]
                # Non-obsolete only (belt-and-suspenders: annotate_cluster already
                # filters, but this guards against stale cached DataFrames)
                lineage_broad = lineage_broad[
                    ~lineage_broad['Cell_Type'].str.lower().str.startswith('obsolete')
                ]

                if not lineage_broad.empty:
                    # Return the deepest (most specific) lineage ancestor, tiebreak
                    # on Combined_Score descending
                    lineage_broad = lineage_broad.sort_values(
                        by=['Depth', 'Combined_Score'], ascending=[False, False]
                    )
                    return lineage_broad.iloc[0]['Cell_Type']

        # --- Fallback: shallowest confident ancestor across all results ---
        broadest = broad.loc[broad['Depth'].idxmin()]
        return broadest['Cell_Type']
