# SapiensOntoCellMap — Feature & Fix Changelog

## 2026-04-28: Annotation Quality, CL ID Normalisation, Output Enrichment

### Fix 1: CL_ → CL: ID Normalisation in cell_name→CL_ID Map (hierarchical_annotation.py)
**Issue:** The master DB stores CL IDs as `CL_0000235` (underscore). The ontology
parser expects `CL:0000235` (colon). The guard `cid.startswith('CL:')` silently
dropped all 386 underscore-format entries, leaving only 1152/1547 unique cell
types mapped to the CL ontology. This caused get_broad_type() to fail for ~25%
of cell types and fall through to the generic fallback.
**Root Cause:** Format inconsistency between the DB build pipeline (uses underscores)
and the CellxGene ontology parser (uses colons).
**Fix:** In `_build_cell_name_to_id_map()`, normalise on read:
`if cid.startswith('CL_'): cid = 'CL:' + cid[3:]`. Map now covers all 1206 entries.

### Fix 2: get_broad_type() Fallback Returns None When top_cell_type Given (hierarchical_annotation.py)
**Issue:** When the constrained ancestor search (lineage ancestors of Top_Cell_Type)
found no results — either due to unmapped CL ID (Fix 1) or no confident ancestors —
the code fell through to a generic fallback: shallowest confident ancestor across
ALL hierarchical results for that cluster. This picked the ancestor of an unrelated
significant hit (e.g., fibroblast cluster → "hematopoietic oligopotent progenitor
cell" because an immune cell type also scored in that cluster).
**Root Cause:** Missing early return inside the `if top_cell_type:` block.
**Fix:** After the constrained search, return `None` instead of falling through when
`top_cell_type` is provided. The generic fallback is now only reached when
`top_cell_type=None` (unconstrained call). Lineage_Conflict dropped from 8–12/20
clusters to 0/20 after both Fix 1 and Fix 2.
**Note:** Broad_Type=None for some clusters (macrophage, fibroblast, endothelial) is
correct — their CL ancestors (leukocyte, myeloid cell) have many CL descendants in
the DB, so Confidence = n_supporting/n_possible stays below 0.5 unless many related
cell types are significant in the cluster. This is honest uncertainty.

### Fix 3: min_cluster_degs Parameter — WE Inflation for Small-n Clusters (get_marker_enrichment_test.py)
**Issue:** Weighted Enrichment formula: `WE = (overlap_w_sum / n) / (ref_w_sum / N)`.
When cluster DEG count `n` is small (e.g., 6–14), the `N/n` ratio becomes very
large, artificially inflating WE and Combined_Score for any cell type capturing
even 2–3 DEGs. This caused off-target cell types to win on small clusters.
**Root Cause:** The WE formula has no guard against small n; min_db_markers guards
the reference size but not the cluster size.
**Fix:** New parameter `min_cluster_degs` (default 0). When `n < min_cluster_degs`
in weighted mode, WE is not computed for that cluster; Combined_Score = NaN; ranking
falls back to adj_p_value. Added CLI arg `--min_cluster_degs` (recommended: 10
for spatial data). Logged as WARNING when triggered.

### Fix 4: Score-Gated Tissue Priority — tissue_priority_ratio (get_cluster_annotation.py)
**Issue:** Hard unconditional priority: selected_tissue Level 2 always wins over
all_tissue regardless of score magnitude. When the tissue-specific DB is sparse for
a cell type (e.g., immune cells in skin), a weak selected_tissue hit suppresses a
strong all_tissue annotation.
**Fix:** New parameter `--tissue_priority_ratio` (default 0.0 = original hard
priority). When `selected_tissue_score < ratio × all_tissue_score`, all_tissue is
used and the override is logged at INFO level.

### Enhancement: Runner_Up_Cell_Type and Score_Gap Output Columns (get_cluster_annotation.py)
**Purpose:** Expose annotation uncertainty directly in the top_annotation_summary.csv.
**Implementation:** When extracting the top-ranked hit per cluster, also extract
the second-ranked hit. Compute `Score_Gap = top_score − runner_up_score`. Add
`Runner_Up_Cell_Type` and `Score_Gap` to the output row.
**Usage:** Score_Gap=0 indicates a tied winner — treat annotation as ambiguous.
Large Score_Gap (>100) indicates a confident, unambiguous annotation.

### Enhancement: Broad_Type_Consensus Output Column (hierarchical_annotation.py, get_cluster_annotation.py)
**Purpose:** Check whether the runner-up annotation agrees on the same broad lineage,
providing a consensus signal for Broad_Type reliability.
**Implementation:** New method `get_broad_type_with_consensus(top_cell_type,
runner_up_cell_type)`. Calls `get_broad_type()` for both. Returns
`{'Broad_Type': str, 'Broad_Type_Consensus': bool}`.

### Enhancement: Proliferating_Lineage Output Column (get_cluster_annotation.py)
**Purpose:** Make Proliferative_Flag interpretable across all cell types. A
proliferating macrophage and a proliferating keratinocyte both flag identically —
this disambiguates them.
**Implementation:** When `Proliferative_Flag=True` and `Broad_Type` is not None,
output `Proliferating_Lineage = f"{Broad_Type}_proliferating"` (e.g.,
`epidermal cell_proliferating`, `myeloid cell_proliferating`).

### Enhancement: Lineage_Conflict Output Column (get_cluster_annotation.py)
**Purpose:** Flag rows where Top_Cell_Type is not a CL ontology descendant of
Broad_Type — indicating an inconsistency between the two fields.
**Implementation:** After computing Broad_Type, check if Top_Cell_Type's CL ID
appears in the ancestor set of Broad_Type's CL ID. If not, set
`Lineage_Conflict=True`. Fires 0 times in current test suite — serves as a
diagnostic guard for future edge cases and DB changes.

### Fix: min_db_markers CLI Default 0→5 (get_cluster_annotation.py)
**Issue:** Default was 0 (no filter), allowing cell types with 1–4 reference
markers to compete via WE. Raising to 10 caused regressions — legitimate cell
types (mature B cell, malignant cell, contractile cell) had fewer than 10 markers
in the tissue-specific DB background intersection and were silently filtered,
allowing tissue-inappropriate winners.
**Fix:** Default set to 5. Filters trivially small reference sets while preserving
legitimate cell types with moderate marker coverage.
**Note:** Raise explicitly via `--min_db_markers` when handling datasets where
WE inflation from small reference sets is a known issue.

### CLI Changes (get_cluster_annotation.py)
- Added `--min_cluster_degs` (int, default=0): Skip WE for clusters below this DEG count
- Added `--tissue_priority_ratio` (float, default=0.0): Score-gated tissue priority fallback
- Changed `--min_db_markers` default: 0→5

### Validation
Tested against 9-sample suite (2 scRNA-seq, 3 Visium HD, 4 Xenium). Zero regressions
in scRNA and Visium samples. Broad_Type artefacts (wrong lineage ancestors) eliminated
across all samples. Lineage_Conflict=0 everywhere. 7 annotation changes in Xenium:
5 improvements (biogenic amine secreting cell → more specific immune/stromal types),
2 pre-existing ambiguous clusters correctly flagged by Score_Gap=0.


## 2026-02-16: Database Schema Validation Layer

### Fix: CL IDs in tissue_id (cellmarkerdb_parser.py)
**Issue:** 90 rows had Cell Ontology IDs (CL_0002248, CL_0002322) in the
`db_tissue_id` column, which should only contain UBERON IDs.
**Root Cause:** CellMarkerDB source `uberonongology_id` column contains
erroneous CL IDs. Parser checked for 'CL:' prefix but source uses 'CL_'.
**Fix:** Positive filter — only keep values starting with 'UBERON:' or 'UBERON_'.
**Justification:** The column's semantic definition requires UBERON ontology
identifiers only. Positive filtering is safer than negative filtering.

### Fix: LMHA IDs in cell_id (cellmarkerdb_parser.py)
**Issue:** 13 rows had LMHA:00142 in `db_cell_id`.
**Root Cause:** CellMarkerDB source `cellontology_id` contains non-CL IDs.
**Fix:** Positive filter — only keep values starting with 'CL:' or 'CL_'.

### Fix: Obsolete Cell Ontology Terms (base_parser.py)
**Issue:** 10,520 rows have cell_name prefixed with "obsolete".
**Root Cause:** cellxgene_ontology_guide maps deprecated CL IDs to their
obsolete labels without filtering.
**Fix:** After name resolution, detect obsolete prefix and attempt replacement
via OBO get_term_replacement(). Fallback: strip prefix.
**Justification:** OBO Foundry FP-009 principle — deprecated terms should
be replaced by their successors.

### Fix: NaN tissue_name Fallback (base_parser.py)
**Issue:** 65,934 rows have NaN tissue_name despite having db_tissue_name.
**Root Cause:** CellxGene "All Tissues" entries not resolved by NLP pipeline.
**Fix:** Fallback to db_tissue_name when tissue_name remains NaN.
**Justification:** NaN makes rows invisible in tissue-filtered queries.

### Enhancement: DatabaseValidator (database_validator.py)
**Purpose:** Strict schema enforcement during DB creation.
**Rules enforced:**
- db_tissue_id, tissue_id: Must start with UBERON: or UBERON_ (or NaN)
- db_cell_id, cell_id: Must start with CL: or CL_ (or NaN)
- cell_name: Must not start with "obsolete"
- gene: Must be non-null
- source_type: Must be in controlled vocabulary (title-cased)
- database: Must be non-empty

### Enhancement: Deduplication
**Issue:** 32,720 duplicate rows on (tissue_id, cell_id, gene, database).
**Fix:** drop_duplicates(keep='first') after concat.
**Justification:** Duplicates inflate weighted enrichment scores.

## 2026-02-16: P0 Enrichment Engine Statistical Correctness Fixes

### Fix 1: Global FDR Correction (get_marker_enrichment_test.py)
**Issue:** BH-FDR was applied per-cluster inside the loop, correcting over ~800
tests per cluster instead of the true burden of ~9,600 total tests (12 clusters x
800 cell types). This inflated significance rates: cSCC showed 31.8% significant,
Visium HD showed 45.2%.
**Fix:** Removed per-cluster FDR. After the enrichment loop, all raw p-values
across all clusters are collected and a single global BH-FDR correction is applied.
**Expected impact:** Significance rates for spatial data should decrease to ~10%,
matching scRNAseq/Xenium levels.

### Fix 2: Robust N / Background Gene Count (get_marker_enrichment_test.py)
**Issue:** Background gene count N was always derived from unique genes in the DEG
file. This is correct for Flex panels (N=2,000) and spatial (N=18,000) but will
fail for standard whole-transcriptome scRNA-seq where FindAllMarkers only exports
significant DEGs (N=3,000-8,000 instead of true ~20,000).
**Fix:** Added `background_gene_count` parameter (CLI: `--background_gene_count`)
to override N. When not provided and deg_type is scrna, a heuristic checks if
N < 15,000 AND >80% of adjusted p-values are < 0.05 — if so, logs a warning
that N may be underestimated.
**Usage:** `--background_gene_count 20000` for standard scRNA-seq with Seurat
FindAllMarkers output.

### Fix 3: Gene Name Normalization (get_marker_enrichment_test.py)
**Issue:** Gene name matching between DEG file and marker DB was exact string
comparison. Case mismatches (e.g., "Cd3e" vs "CD3E") caused silent false negatives.
**Fix:** Both marker DB gene names and DEG Feature Names are normalized to
uppercase + stripped whitespace before intersection. A reverse display-name map
preserves original case for output reporting.
**Scope:** Addresses case and whitespace mismatches. Alias/synonym resolution
(e.g., TNFRSF1B vs TNFR2) is deferred to a future HGNC integration.

### Fix 4: Minimum Overlap Filter (get_marker_enrichment_test.py)
**Issue:** Single-gene overlaps (k=1) dominated spatial results — 50% of cSCC
results and 64% of Xenium results had k=1. These are biologically noisy and
often based on housekeeping or broadly expressed genes.
**Fix:** Added `min_overlap` parameter (default=2, CLI: `--min_overlap`). Results
with fewer than `min_overlap` overlapping genes are excluded.
**Usage:** `--min_overlap 1` to restore old behavior if needed.

### CLI Changes (get_cluster_annotation.py)
- Added `--min_overlap` (int, default=2): Minimum gene overlap to report
- Added `--background_gene_count` (int, default=None): Override N for hypergeometric test
