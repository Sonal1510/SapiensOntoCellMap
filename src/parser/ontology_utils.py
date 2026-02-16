#!/usr/bin/python3
"""
Author      : Sonal Rashmi (reviewed and improved)
Date        : 01/10/2025 (updated)
Description : Robust CellxGene ontology helper that builds label<->id maps
              and exposes resilient name-to-id lookups for CL and UBERON.
              This version includes improved normalization, matching logic,
              and an optional 3-stage hybrid search (Exact, Fuzzy, Semantic).
             
              FIXED: _normalize_text is now a method of the class,
              so it can be accessed by BaseParser.
"""

from __future__ import annotations
import re
import unicodedata
import os
import pickle
from typing import Optional, Dict, Iterable, Tuple, Any, List

# primary dependency
from cellxgene_ontology_guide.ontology_parser import OntologyParser

# --- NEW: Dependencies for Semantic Search (Step 3) ---
try:
    from sentence_transformers import SentenceTransformer, util
    import torch
    _USE_SEMANTIC = True
except ImportError:
    _USE_SEMANTIC = False
    print("Warning: 'sentence-transformers' or 'torch' not installed. Semantic search will be disabled.")
    print("To enable it, run: pip install sentence-transformers torch")
    class SentenceTransformer: pass # type: ignore
    class torch: # type: ignore
        class Tensor: pass

# --- Dependencies for Fuzzy Search (Step 2) ---
try:
    from rapidfuzz import process as rf_process, fuzz as rf_fuzz
    _USE_RAPIDFUZZ = True
except Exception:
    from thefuzz import process as fuzz_process, fuzz as a_fuzz # type: ignore
    _USE_RAPIDFUZZ = False
    print("Warning: 'rapidfuzz' not found, falling back to 'thefuzz'.")
    print("For faster fuzzy matching, run: pip install rapidfuzz")

# --- Stop Words and Generic IDs ---
STOP_WORDS = {"cell", "cells", "tissue", "of", "the", "and", "or", "in", "a"}
GENERIC_TERM_IDS = {
    "CL:0000000",      # cell
    "CL:0000001",      # primary cultured cell
    "UBERON:0000479",  # tissue
    "UBERON:0000062",  # organ
    "UBERON:0001062",  # anatomical entity
}


class CellxGeneOntologyParser:
    """
    Robust wrapper around cellxgene_ontology_guide.OntologyParser that:
     - builds label->id and id->label maps for CL and UBERON
     - includes synonyms in the search space
     - supports fuzzy matching with score and match metadata
     - supports an optional, 3-stage hybrid search (Exact, Fuzzy, Semantic)
    """

    def __init__(
        self,
        fuzzy_threshold: int = 85,
        semantic_model_name: Optional[str] = None,  # e.g., 'all-MiniLM-L6-v2'
        semantic_index_path: Optional[str] = None # e.g., './data/ontology_semantic_index.pkl'
    ):
        self.ontology_parser = OntologyParser()
        self.fuzzy_threshold = int(fuzzy_threshold)

        # maps
        self.cl_name_to_id: Dict[str, str] = {}
        self.cl_id_to_name: Dict[str, str] = {}
        self.uberon_name_to_id: Dict[str, str] = {}
        self.uberon_id_to_name: Dict[str, str] = {}

        # build maps
        self._build_maps()

        # --- Load Semantic Model and Index ---
        self.semantic_model: Optional[SentenceTransformer] = None
        self.cl_index: Optional[Dict[str, Any]] = None
        self.uberon_index: Optional[Dict[str, Any]] = None

        if _USE_SEMANTIC and semantic_model_name and semantic_index_path:
            try:
                print(f"Loading semantic model '{semantic_model_name}'...")
                self.semantic_model = SentenceTransformer(semantic_model_name)
                print(f"Loading semantic index from '{semantic_index_path}'...")
                if os.path.exists(semantic_index_path):
                    with open(semantic_index_path, 'rb') as f:
                        index_data = pickle.load(f)
                        self.cl_index = index_data.get('cl')
                        self.uberon_index = index_data.get('uberon')
                    print("Semantic model and index loaded successfully.")
                else:
                    print(f"Warning: Semantic index file not found at '{semantic_index_path}'.")
                    self.semantic_model = None # Disable semantic search
            except Exception as e:
                print(f"Error loading semantic model/index: {e}")
                self.semantic_model = None
        elif (semantic_model_name or semantic_index_path) and not _USE_SEMANTIC:
                print("Semantic search is disabled because 'sentence-transformers' or 'torch' are not installed.")

    def _normalize_text(self, s: str) -> str:
        """
        Improved helper to normalize input strings for consistent matching.
        """
        s = "" if s is None else str(s)
        s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("utf-8")
        s = re.sub(r"[^\w\s+]+", " ", s)
        s = re.sub(r"[-_]+", " ", s)
        s = s.lower()
        s_parts = [word for word in s.split() if word not in STOP_WORDS]
        s = " ".join(s_parts)
        s = re.sub(r"\s+", " ", s).strip()
        return s

    def _get_term_label_to_id_map(self, ontology_prefix: str) -> Dict[str, str]:
        # ... (this method is unchanged) ...
        candidates = [
            "get_term_label_to_id_map", "get_label_to_id_map",
            "get_term_label_to_id", "get_term_map",
        ]
        for method in candidates:
            if hasattr(self.ontology_parser, method):
                try:
                    fn = getattr(self.ontology_parser, method)
                    result = fn(ontology_prefix)
                    if isinstance(result, dict):
                        return result
                except Exception:
                    continue
        # ... (fallback logic) ...
        ontology_obj = getattr(self.ontology_parser, "ontology", None)
        if isinstance(ontology_obj, dict):
            label_to_id: Dict[str, str] = {}
            for tid, term in ontology_obj.items():
                label = getattr(term, "label", None) or getattr(term, "name", None) or str(term)
                if isinstance(label, str) and tid.startswith(ontology_prefix + ":"):
                    label_to_id[label] = tid
            if label_to_id:
                return label_to_id
        if hasattr(self.ontology_parser, "get_terms_by_ontology"):
            try:
                result = {}
                for t in self.ontology_parser.get_terms_by_ontology(ontology_prefix):
                    tid = getattr(t, "id", None)
                    lab = getattr(t, "label", None) or getattr(t, "name", None) or str(t)
                    if tid and lab:
                        result[lab] = tid
                if result:
                    return result
            except Exception:
                pass
        raise RuntimeError(f"Unable to build label->id map for {ontology_prefix}")

    def _get_term_synonyms(self, term_id: str) -> Iterable[str]:
        # ... (this method is unchanged) ...
        candidates = ["get_term_synonyms", "get_synonyms", "term_synonyms"]
        for method in candidates:
            if hasattr(self.ontology_parser, method):
                try:
                    fn = getattr(self.ontology_parser, method)
                    syn = fn(term_id)
                    if syn: return syn
                except Exception:
                    continue
        ont = getattr(self.ontology_parser, "ontology", None)
        if isinstance(ont, dict) and term_id in ont:
            term = ont[term_id]
            if hasattr(term, "synonyms"):
                s = getattr(term, "synonyms")
                if s: return s
        return []

    def _build_maps(self) -> None:
        cl_map_raw = self._get_term_label_to_id_map("CL")
        uberon_map_raw = self._get_term_label_to_id_map("UBERON")
        self.cl_id_to_name = {v: k for k, v in cl_map_raw.items()}
        self.uberon_id_to_name = {v: k for k, v in uberon_map_raw.items()}

        self.cl_name_to_id = {self._normalize_text(k): v for k, v in cl_map_raw.items()}
        self.uberon_name_to_id = {self._normalize_text(k): v for k, v in uberon_map_raw.items()}

        for cl_id in list(self.cl_id_to_name.keys()):
            for s in self._get_term_synonyms(cl_id) or []:
                if s and isinstance(s, str):
                    self.cl_name_to_id[self._normalize_text(s)] = cl_id
        for uber_id in list(self.uberon_id_to_name.keys()):
            for s in self._get_term_synonyms(uber_id) or []:
                if s and isinstance(s, str):
                    self.uberon_name_to_id[self._normalize_text(s)] = uber_id

    def _find_best_match(
        self,
        query_name: str,
        name_to_id_map: Dict[str, str],
        id_to_name_map: Dict[str, str],
        semantic_index: Optional[Dict[str, Any]],
        approx_match_threshold: int,
    ) -> Optional[Dict[str, Any]]:
        # ... (this method is unchanged) ...
        q_norm = self._normalize_text(query_name)
        if not q_norm:
            return None
        if q_norm in name_to_id_map:
            match_id = name_to_id_map[q_norm]
            if match_id not in GENERIC_TERM_IDS:
                return {"query": query_name, "match_id": match_id, "match_label": id_to_name_map.get(match_id), "score": 100, "type": "exact"}
        choices = list(name_to_id_map.keys())
        if choices: 
            limit = 5
            if _USE_RAPIDFUZZ:
                scorer = rf_fuzz.token_set_ratio
                best_hits = rf_process.extract(q_norm, choices, scorer=scorer, limit=limit)
            else:
                scorer = a_fuzz.token_set_ratio # type: ignore
                best_hits = fuzz_process.extract(q_norm, choices, scorer=scorer, limit=limit) # type: ignore
            for match_text, score, _ in best_hits:
                score = int(round(score))
                if score < approx_match_threshold:
                    continue 
                match_id = name_to_id_map.get(match_text)
                if not match_id or match_id in GENERIC_TERM_IDS:
                    continue 
                len_ratio = len(match_text) / len(q_norm) if q_norm else 0
                if len_ratio < 0.3:
                    score = int(score * (len_ratio + 0.5))
                if score >= approx_match_threshold:
                    return {"query": query_name, "match_id": match_id, "match_label": id_to_name_map.get(match_id), "score": score, "type": "fuzzy"}
        if self.semantic_model and semantic_index:
            try:
                query_embedding = self.semantic_model.encode(q_norm, convert_to_tensor=True)
                hits = util.semantic_search(query_embedding, semantic_index['embeddings'], top_k=1)
                if hits and hits[0]:
                    best_hit = hits[0][0]
                    score = int(round(best_hit['score'] * 100))
                    semantic_threshold = 60
                    if score >= semantic_threshold:
                        match_id = semantic_index['ids'][best_hit['corpus_id']]
                        if match_id not in GENERIC_TERM_IDS:
                            return {
                                "query": query_name, "match_id": match_id,
                                "match_label": id_to_name_map.get(match_id),
                                "score": score, "type": "semantic"
                            }
            except Exception as e:
                print(f"Error during semantic search for '{query_name}': {e}")
        return None

    def get_cell_name_given_id(self, cell_id: str) -> Optional[str]:
        return self.cl_id_to_name.get(cell_id)

    def find_best_cell_name_match(self, cell_name: str, approx_match_threshold: Optional[int] = None) -> Optional[Dict[str, Any]]:
        threshold = approx_match_threshold if approx_match_threshold is not None else self.fuzzy_threshold
        return self._find_best_match(cell_name, self.cl_name_to_id, self.cl_id_to_name, self.cl_index, threshold)

    def get_cell_id_given_name(self, cell_name: str, approx_match_threshold: Optional[int] = None) -> Optional[str]:
        info = self.find_best_cell_name_match(cell_name, approx_match_threshold)
        return info.get("match_id") if info else None

    def get_pairwise_cell_distances(self, cell_ids: List[str]) -> Dict[str, Dict[str, int]]:
        # ... (this method is unchanged) ...
        if not hasattr(self.ontology_parser, "get_term_ancestors_with_distances"):
            print("❌ Error: The 'get_term_ancestors_with_distances' method is not available...")
            return {}
        unique_ids = sorted(list(set(cell_id.replace("_", ":") for cell_id in cell_ids)))
        if len(unique_ids) < 2:
            return {}
        ancestor_cache: Dict[str, Optional[Dict[str, int]]] = {}
        for cell_id in unique_ids:
            if cell_id not in self.cl_id_to_name:
                print(f"⚠️ Warning: Cell ID '{cell_id}' not found in Cell Ontology. Skipping.")
                ancestor_cache[cell_id] = None
                continue
            try:
                ancestors = self.ontology_parser.get_term_ancestors_with_distances(cell_id)
                ancestors[cell_id] = 0 
                ancestor_cache[cell_id] = ancestors
            except Exception as e:
                print(f"⚠️ Warning: Could not get ancestors for {cell_id}: {e}")
                ancestor_cache[cell_id] = None
        results: Dict[str, Dict[str, int]] = {uid: {} for uid in unique_ids if ancestor_cache[uid]}
        for i in range(len(unique_ids)):
            cell_id_1 = unique_ids[i]
            ancestors_1 = ancestor_cache.get(cell_id_1)
            if ancestors_1 is None:
                continue 
            set_1 = set(ancestors_1.keys())
            for j in range(i + 1, len(unique_ids)):
                cell_id_2 = unique_ids[j]
                ancestors_2 = ancestor_cache.get(cell_id_2)
                if ancestors_2 is None:
                    continue 
                set_2 = set(ancestors_2.keys())
                common_ancestors = set_1 & set_2
                if not common_ancestors:
                    results[cell_id_1][cell_id_2] = 999
                    results[cell_id_2][cell_id_1] = 999
                    continue
                min_distance = 999
                for ancestor_id in common_ancestors:
                    dist = ancestors_1[ancestor_id] + ancestors_2[ancestor_id]
                    if dist < min_distance:
                        min_distance = dist
                results[cell_id_1][cell_id_2] = min_distance
                results[cell_id_2][cell_id_1] = min_distance
        return results
    
    # === **** START OF BUGFIX **** ===
    # This method was missing from your file, causing the error.
    # It provides the data needed for the hierarchy.
    def get_cell_ontology_graph(self) -> Dict[str, Dict[str, Any]]:
        """
        Exposes the CL (Cell Ontology) graph from the wrapped parser
        using robust, public methods.
        
        This is required by the SapiensMapGenerator to build the hierarchy.
        
        Returns:
            Dict[str, Dict[str, Any]]: 
                A dictionary where keys are CL IDs and values are dicts
                containing 'name' (str) and 'parents' (List[str]).
                e.g., {'CL:0000084': {'name': 'T cell', 'parents': ['CL:0000086']}}
        """
        print("Building CL graph using public methods...")
        
        if not hasattr(self.ontology_parser, "get_term_ancestors_with_distances"):
            print("❌ Error: The 'get_term_ancestors_with_distances' method is not available...")
            return {}
            
        cl_graph: Dict[str, Dict[str, Any]] = {}

        if not self.cl_id_to_name:
             print("❌ Error: `cl_id_to_name` map is empty. `_build_maps` may have failed.")
             return {}
             
        for cl_id, cl_name in self.cl_id_to_name.items():
            cl_graph[cl_id] = {
                'name': cl_name,
                'parents': []
            }

        for cl_id in cl_graph.keys():
            try:
                ancestors = self.ontology_parser.get_term_ancestors_with_distances(cl_id)
                direct_parents = [
                    ancestor_id for ancestor_id, distance in ancestors.items()
                    if distance == 1
                ]
                cl_graph[cl_id]['parents'] = direct_parents
            except Exception:
                cl_graph[cl_id]['parents'] = []

        print(f"Built CL graph with {len(cl_graph)} nodes.")
        return cl_graph
        
    # ** NEW METHOD TO FIX HIERARCHY BUG **
    def is_term_obsolete(self, term_id: str) -> bool:
        """
        Checks if a term ID is obsolete/deprecated.
        """
        if not hasattr(self.ontology_parser, "ontology"):
             print("❌ Error: `ontology_parser.ontology` attribute not found.")
             return False
             
        term = self.ontology_parser.ontology.get(term_id)
        if term:
            return getattr(term, 'deprecated', False)
        return False # Assume not obsolete if not found
    # === **** END OF BUGFIX **** ===

    # ---------- UBERON (tissue) methods ----------
    def get_tissue_name_given_id(self, tissue_id: str) -> Optional[str]:
        return self.uberon_id_to_name.get(tissue_id)

    def find_best_tissue_name_match(self, tissue_name: str, approx_match_threshold: Optional[int] = None) -> Optional[Dict[str, Any]]:
        threshold = approx_match_threshold if approx_match_threshold is not None else self.fuzzy_threshold
        return self._find_best_match(tissue_name, self.uberon_name_to_id, self.uberon_id_to_name, self.uberon_index, threshold)

    def get_tissue_id_given_name(self, tissue_name: str, approx_match_threshold: Optional[int] = None) -> Optional[str]:
        info = self.find_best_tissue_name_match(tissue_name, approx_match_threshold)
        return info.get("match_id") if info else None