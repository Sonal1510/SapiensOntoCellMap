#!/usr/bin/python3
"""
Author      : Sonal Rashmi (reviewed and improved)
Date        : 01/10/2025 (updated)
Description : Robust CellxGene ontology helper that builds label<->id maps
              and exposes resilient name-to-id lookups for CL and UBERON.
              This version includes improved normalization and matching logic
              to avoid overly generic terms and increase specificity.
"""

from __future__ import annotations
import re
import unicodedata
from typing import Optional, Dict, Iterable, Tuple, Any, List

# primary dependency
from cellxgene_ontology_guide.ontology_parser import OntologyParser

# fuzzy matching: prefer rapidfuzz if available (faster), otherwise fall back to thefuzz
try:
    from rapidfuzz import process as rf_process, fuzz as rf_fuzz
    _USE_RAPIDFUZZ = True
except Exception:
    from thefuzz import process as fuzz_process, fuzz as a_fuzz # type: ignore
    _USE_RAPIDFUZZ = False

# --- IMPROVEMENT 1: Enhanced Normalization and Stop Words ---
# Define common, low-information words to remove during normalization
STOP_WORDS = {"cell", "cells", "tissue", "of", "the", "and", "or", "in", "a"}
# Define overly generic ontology IDs that we want to avoid matching to
GENERIC_TERM_IDS = {
    "CL:0000000",       # cell
    "CL:0000001",       # primary cultured cell
    "UBERON:0000479",   # tissue
    "UBERON:0000062",   # organ
    "UBERON:0001062",   # anatomical entity
}

def _normalize_text(s: str) -> str:
    """
    Improved helper to normalize input strings for consistent matching.
    - Strips accents and normalizes unicode.
    - Removes complex punctuation, parentheses, and brackets.
    - Replaces hyphens and underscores with spaces.
    - Removes common biological stop words.
    - Converts to lowercase and collapses whitespace.
    """
    s = "" if s is None else str(s)
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("utf-8")
    # Keep alphanumeric, plus sign (for markers like CD4+), and spaces
    s = re.sub(r"[^\w\s+]+", " ", s)
    s = re.sub(r"[-_]+", " ", s)
    s = s.lower()
    # Remove stop words
    s_parts = [word for word in s.split() if word not in STOP_WORDS]
    s = " ".join(s_parts)
    s = re.sub(r"\s+", " ", s).strip()
    return s

class CellxGeneOntologyParser:
    """
    Robust wrapper around cellxgene_ontology_guide.OntologyParser that:
     - builds label->id and id->label maps for CL and UBERON
     - includes synonyms in the search space
     - supports fuzzy matching with score and match metadata
     - tolerates multiple versions of the underlying API
     - **IMPROVED**: Avoids matching to generic terms and increases specificity.

    Example:
        parser = CellxGeneOntologyParser()
        parser.get_cell_id_given_name("b cell")
        parser.find_best_tissue_name_match("skin of body")
    """

    def __init__(self, fuzzy_threshold: int = 85):
        self.ontology_parser = OntologyParser()
        self.fuzzy_threshold = int(fuzzy_threshold)

        # maps (kept lowercase on keys where needed)
        self.cl_name_to_id: Dict[str, str] = {}
        self.cl_id_to_name: Dict[str, str] = {}
        self.uberon_name_to_id: Dict[str, str] = {}
        self.uberon_id_to_name: Dict[str, str] = {}

        # build maps (robust to API differences)
        self._build_maps()

    # ---------- Internal helpers (API fallbacks) ----------
    def _get_term_label_to_id_map(self, ontology_prefix: str) -> Dict[str, str]:
        # ... (no changes to this method) ...
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
        # ... (no changes to this method) ...
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

    # ---------- Build maps ----------
    def _build_maps(self) -> None:
        cl_map_raw = self._get_term_label_to_id_map("CL")
        uberon_map_raw = self._get_term_label_to_id_map("UBERON")
        self.cl_id_to_name = {v: k for k, v in cl_map_raw.items()}
        self.uberon_id_to_name = {v: k for k, v in uberon_map_raw.items()}
        # --- IMPROVEMENT: Use normalized keys for lookup ---
        self.cl_name_to_id = {_normalize_text(k): v for k, v in cl_map_raw.items()}
        self.uberon_name_to_id = {_normalize_text(k): v for k, v in uberon_map_raw.items()}

        for cl_id in list(self.cl_id_to_name.keys()):
            for s in self._get_term_synonyms(cl_id) or []:
                if s and isinstance(s, str):
                    self.cl_name_to_id[_normalize_text(s)] = cl_id
        for uber_id in list(self.uberon_id_to_name.keys()):
            for s in self._get_term_synonyms(uber_id) or []:
                if s and isinstance(s, str):
                    self.uberon_name_to_id[_normalize_text(s)] = uber_id

    # ---------- Generic search method (New and Improved) ----------
    def _find_best_match(
        self,
        query_name: str,
        name_to_id_map: Dict[str, str],
        id_to_name_map: Dict[str, str],
        approx_match_threshold: int,
    ) -> Optional[Dict[str, Any]]:
        """
        A generic, improved search function for finding the best ontology match.
        """
        q_norm = _normalize_text(query_name)
        if not q_norm:
            return None

        # 1. Exact match on normalized text (fastest and best)
        if q_norm in name_to_id_map:
            match_id = name_to_id_map[q_norm]
            if match_id not in GENERIC_TERM_IDS:
                return {"query": query_name, "match_id": match_id, "match_label": id_to_name_map.get(match_id), "score": 100, "type": "exact"}

        # 2. Fuzzy match if no exact match is found
        choices = list(name_to_id_map.keys())
        if not choices:
            return None

        # --- IMPROVEMENT 2: Use a more robust scorer and get multiple candidates ---
        limit = 5  # Get top 5 candidates to filter through
        if _USE_RAPIDFUZZ:
            # token_set_ratio is good at handling different word orders and extra words
            scorer = rf_fuzz.token_set_ratio
            best_hits = rf_process.extract(q_norm, choices, scorer=scorer, limit=limit)
        else:
            scorer = a_fuzz.token_set_ratio
            best_hits = fuzz_process.extract(q_norm, choices, scorer=scorer, limit=limit)

        # --- IMPROVEMENT 3: Filter candidates to find the best non-generic match ---
        for match_text, score, _ in best_hits:
            score = int(round(score))
            if score < approx_match_threshold:
                continue

            match_id = name_to_id_map.get(match_text)
            if not match_id or match_id in GENERIC_TERM_IDS:
                continue # Skip this generic or invalid match

            # --- IMPROVEMENT 4: Penalize matches that are much shorter than the query ---
            # This helps prevent "progenitor cell" from strongly matching "cell"
            len_ratio = len(match_text) / len(q_norm) if q_norm else 0
            if len_ratio < 0.5: # If match is less than half the query length
                score = int(score * (len_ratio + 0.5)) # Apply a penalty

            if score >= approx_match_threshold:
                return {"query": query_name, "match_id": match_id, "match_label": id_to_name_map.get(match_id), "score": score, "type": "fuzzy"}

        return None # No suitable match found

    # ---------- CL (cell type) methods ----------
    def get_cell_name_given_id(self, cell_id: str) -> Optional[str]:
        return self.cl_id_to_name.get(cell_id)

    def find_best_cell_name_match(self, cell_name: str, approx_match_threshold: Optional[int] = None) -> Optional[Dict[str, Any]]:
        threshold = approx_match_threshold if approx_match_threshold is not None else self.fuzzy_threshold
        return self._find_best_match(cell_name, self.cl_name_to_id, self.cl_id_to_name, threshold)

    def get_cell_id_given_name(self, cell_name: str, approx_match_threshold: Optional[int] = None) -> Optional[str]:
        info = self.find_best_cell_name_match(cell_name, approx_match_threshold)
        return info.get("match_id") if info else None

    # ---------- UBERON (tissue) methods ----------
    def get_tissue_name_given_id(self, tissue_id: str) -> Optional[str]:
        return self.uberon_id_to_name.get(tissue_id)

    def find_best_tissue_name_match(self, tissue_name: str, approx_match_threshold: Optional[int] = None) -> Optional[Dict[str, Any]]:
        threshold = approx_match_threshold if approx_match_threshold is not None else self.fuzzy_threshold
        return self._find_best_match(tissue_name, self.uberon_name_to_id, self.uberon_id_to_name, threshold)

    def get_tissue_id_given_name(self, tissue_name: str, approx_match_threshold: Optional[int] = None) -> Optional[str]:
        info = self.find_best_tissue_name_match(tissue_name, approx_match_threshold)
        return info.get("match_id") if info else None