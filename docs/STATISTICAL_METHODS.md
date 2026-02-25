# SapiensOntoCellMap: Statistical Methods Reference

**Version:** 0.1.0 | **Date:** 2026-02-23 | **Author:** Sonal Rashmi

This document describes the complete statistical strategy used in
SapiensOntoCellMap for cell type annotation (hypergeometric enrichment),
evidence-weighted scoring, and annotation-derived cell type composition.

---

## Table of Contents

1. [Notation](#1-notation)
2. [Hypergeometric Enrichment Test](#2-hypergeometric-enrichment-test)
3. [Multiple Testing Correction](#3-multiple-testing-correction)
4. [Evidence-Weighted Scoring System](#4-evidence-weighted-scoring-system)
5. [Combined Score (Ranking Metric)](#5-combined-score-ranking-metric)
6. [Hierarchical Annotation Confidence](#6-hierarchical-annotation-confidence)
7. [Cell Type Composition Scoring](#7-cell-type-composition-scoring)
8. [Input Filtering and Gene Matching](#8-input-filtering-and-gene-matching)
9. [Known Limitations and Caveats](#9-known-limitations-and-caveats)

---

## 1. Notation

| Symbol | Definition |
|--------|-----------|
| `N` | Total number of background genes (all genes tested; see §9) |
| `K` | Number of reference genes for cell type `c` (from marker DB, intersected with background) |
| `n` | Number of DEGs passing filters for cluster `i` |
| `k` | Number of genes in the overlap: DEGs ∩ reference markers |
| `W_g` | Evidence weight of gene `g` for cell type `c` (see §4) |
| `W_overlap` | Sum of evidence weights over overlapping genes: Σ W_g for g ∈ DEGs ∩ markers |
| `W_ref` | Sum of evidence weights over all reference genes: Σ W_g for g ∈ markers |
| `adj_p` | BH-corrected p-value |

---

## 2. Hypergeometric Enrichment Test

### Null hypothesis

Under the null, the `n` DEGs for cluster `i` are a uniform random draw
(without replacement) from the `N` background genes. The probability of
observing an overlap of ≥ `k` genes with the `K` marker genes of cell type
`c` follows the hypergeometric distribution:

```
P(X ≥ k) = Σ_{x=k}^{min(n,K)}  C(K,x) * C(N-K, n-x) / C(N, n)
```

In code (`scipy.stats.hypergeom.sf`):

```python
p_value = hypergeom.sf(k - 1, N, K, n)
```

This is the one-tailed survival function (upper tail): probability of observing
at least `k` overlapping genes by chance.

### Background gene set N

`N` is determined as follows (in priority order):

1. **User-specified** `--background_gene_count` — exact user override.
2. **Inferred from DEG file** — number of unique genes present in the input
   (all clusters combined), after HGNC alias resolution.

> **Critical for scRNA-seq:** Seurat `FindAllMarkers` returns only genes that
> passed Wilcoxon filtering, typically 3,000–8,000 genes. For whole-transcriptome
> 10x 3′/5′ data, `N` should be ~20,000. Use `--background_gene_count 20000`.
> A warning is emitted when `N < 15,000` and >80% of adjusted p-values are
> significant (pattern indicating DEG-only input).

### Minimum overlap filter

Results with `k < 2` (default, configurable via `--min_overlap`) are discarded
before p-value calculation to reduce false discoveries from single-gene coincidences.

### Unweighted enrichment ratio

```
Enrichment_ratio = (k/n) / (K/N)
```

Interpretation: observed overlap fraction divided by expected overlap fraction
under the null. Values > 1 indicate enrichment; equals the ratio of observed
to expected overlap count.

---

## 3. Multiple Testing Correction

Global Benjamini–Hochberg (BH) FDR correction is applied across **all
(cluster, cell_type) pairs simultaneously** — not per-cluster.

```python
adj_p_value = statsmodels.stats.multitest.multipletests(
    all_p_values, method="fdr_bh"
)[1]
```

**Why global?** Per-cluster correction is anticonservative when clusters have
very few significant hits (inflated family size within one cluster masks the
real FDR). Global correction maintains the intended 5% FDR across the entire
experiment.

**Reporting threshold:** Results with `adj_p_value ≥ --pval` (default 0.05)
are retained in the full results table but flagged. Significant results used
in heatmaps and summaries require `adj_p_value < 0.05`.

---

## 4. Evidence-Weighted Scoring System

### Motivation

Not all marker–cell type associations are equally credible. An `Experiment`-
supported marker (validated by flow cytometry or immunostaining) is more
reliable than a `Computational` prediction. Markers reported by multiple
independent databases are more reliable than those from a single source.

The weighted mode (`Level 2`) applies evidence weights as follows:

### Source-type weight

| `source_type` | Weight |
|---------------|--------|
| Experiment | 4.0 |
| Single-Cell Sequencing | 3.0 |
| Literature | 2.0 |
| Review | 2.0 |
| Company | 1.0 |
| Computational | 0.5 |

Weights reflect the hierarchy: experimental evidence > transcriptomic profiling
> curated literature > computational inference.

### Cross-database agreement multiplier

For each (gene, cell_type) pair, if `D` independent databases report the
association:

```
cross_db_multiplier = 1 + log2(D)
```

This rewards consistency across sources: a gene reported by 4 databases
gets a multiplier of 3.0× vs. a single-database gene at 1.0×.

### Final evidence weight W_g

Within the annotation pipeline (`get_cluster_annotation.py`), weights are
aggregated per (cell_type, gene) across unique databases:

```
W_g = Σ_{d reporting (g, c)} source_type_weight(d)
```

This is the sum of source-type weights across all unique databases reporting
gene `g` for cell type `c`. It naturally encodes both per-source credibility
and cross-database agreement.

The `sum` aggregation rewards genes supported by multiple databases and
higher-quality source types simultaneously.

### Weighted Recall

```
Weighted_Recall = W_overlap / W_ref
                = Σ W_g [g ∈ DEGs ∩ markers] / Σ W_g [g ∈ markers ∩ background]
```

**Interpretation:** The fraction of total evidence weight for this cell type
that is "explained" by the cluster's DEGs. Values in [0, 1]. A value of 0.7
means the cluster's DEGs account for 70% of the total marker evidence weight.

### Weighted Enrichment Ratio

```
Weighted_Enrichment = (W_overlap / n) / (W_ref / N)
```

This is the weighted analog of `Enrichment_ratio`. Decomposed:

- `W_overlap / n` = mean evidence weight per DEG (observed evidence density)
- `W_ref / N` = mean evidence weight per background gene (expected density under null)

**Interpretation:** How much denser (in evidence weight) are the cluster's
DEGs compared to what is expected if DEGs were a random draw from the
background? Values > 1 indicate enrichment.

**Algebraic identity:**

```
Weighted_Enrichment = Weighted_Recall / (n/N)
```

This shows the relationship to `Weighted_Recall` directly: it is the
observed weighted recall normalized by the expected recall under the null.

---

## 5. Combined Score (Ranking Metric)

After FDR correction:

```
Combined_Score = Weighted_Enrichment * (-log10(adj_p_value))
```

**Purpose:** Single-column ranking that integrates:
- **Statistical significance** via `-log10(adj_p_value)` — penalizes high p-values
- **Effect size** via `Weighted_Enrichment` — penalizes low evidence density

**Properties:**
- Combined_Score = 0 when adj_p_value = 1 (no statistical signal)
- Combined_Score scales linearly with effect size and logarithmically with significance
- Suitable for ranking but **not a formal test statistic** — do not interpret as a p-value

**Use:** The top annotation summary and the heatmap sort by `adj_p_value`
ascending and `Combined_Score` descending within each cluster.

---

## 6. Hierarchical Annotation Confidence

After flat enrichment (Level 2), significant hits are mapped onto the
Cell Ontology (CL) graph via ancestor traversal.

### Confidence score

For each CL node `v` with evidence from `k_v` supporting descendant cell types
out of `K_v` total descendants represented in the marker database:

```
Confidence(v) = k_v / K_v
```

where:
- `k_v` = number of significant cell type hits that are `v` or its descendants
- `K_v` = number of cell types in the marker DB that are `v` or its descendants
  (intersected with the set of background genes to ensure testability)

**Interpretation:** Confidence ≈ 0.8 means 80% of the known subtypes of
cell type `v` (that exist in the database) were found significant. This is
conceptually analogous to a sensitivity measure at this resolution level.

### Resolution labels

| Confidence | Ontology depth | Resolution label |
|------------|---------------|-----------------|
| ≥ 0.5 | 1–3 | `broad` |
| ≥ 0.5 | 4–6 | `intermediate` |
| ≥ 0.5 | ≥ 7 | `fine` |
| < 0.5 | any | `uncertain` |

**Reporting threshold:** Default `confidence_threshold = 0.5`. Nodes with
`N_Supporting < 2` are suppressed unless they are direct significant hits
(leaf nodes), to avoid over-interpreting single-subtype support at broad levels.

### Combined hierarchical score

For each ancestor node, an aggregated `-log10(p)` score is summed across all
supporting child hits:

```
Combined_Score(v) = Σ_{c ∈ supporting descendants} -log10(p_c)
```

This rewards nodes supported by many highly-significant subtypes.

---

## 7. Cell Type Composition Scoring

### Rationale

Cell type composition and annotation address complementary questions:

- **Annotation** asks: *What is the most likely cell type identity of this cluster?*
- **Composition** asks: *What is the relative cell type evidence balance across
  all enriched types in this cluster?*

Prior versions used NNLS deconvolution on the marker database treated as a
signature expression matrix. This approach is fundamentally inappropriate:
the marker database is a sparse, evidence-weighted binary presence matrix
(most cell types have 4–42 genes), not a dense expression profile. L1
normalization of sparse signature columns amplifies cell types with few
markers — a cell type with 1 marker gene receives unit weight per gene
versus a cell type with 42 markers at 0.024 weight per gene — producing
biologically implausible dominance of rare, narrowly-supported cell types.

The annotation-derived composition approach uses the database as intended
(via enrichment testing) and derives composition from enrichment scores.

### CL Ontology Ancestor Pruning

Raw Level 2 enrichment results include both specific cell types (e.g.
"keratinocyte", "helper T cell") and their CL ontological ancestors
(e.g. "epithelial cell", "lymphocyte", "T cell"). Because ancestors share
all the marker gene sets of their descendants, they are enriched in the same
clusters and receive comparable Combined_Scores. Without pruning, this dilutes
the specific cell type signal: "keratinocyte" appears at only ~12% composition
when swamped by "epithelial cell", "keratin accumulating cell", "barrier cell",
and other ancestors.

Pruning uses the hierarchical annotation output (`N_Supporting` column):

| `N_Supporting` | Interpretation | Action |
|---------------|----------------|--------|
| 1 | Leaf node: no more specific significant descendant exists | **Keep** |
| > 1 | Ancestor node: more specific descendants are also significant | **Remove** |

Two additional filters prevent garbage passing the `N_Supporting == 1` check:

- **OBSOLETE terms**: CL retires old classes with `obsolete` prefix in the
  Cell_Type name (e.g. "obsolete barrier cell"). Their database-facing names
  would pollute the leaf set.
- **Depth < 2**: root-adjacent terms ("eukaryotic cell", "nucleate cell") are
  not integrated into the functional CL sub-graph.

**Top-1 safety net**: if the highest-scoring significant cell type for a
cluster was pruned as an ancestor (e.g. "macrophage" is pruned because
"inflammatory macrophage" and "alternatively activated macrophage" are also
significant), it is re-added to the pool. This preserves the dominant
biological label while still removing lower-ranked generic ancestors.

### Composition formula

After ancestor pruning, each remaining cell type `c` contributes:

```
raw_score(c) = Weighted_Enrichment(c) × −log₁₀(adj_p_value(c))
```

Scores are clipped to ≥ 0 (Combined_Score is non-negative by construction
when adj_p_value ≤ 1), then normalized to sum to 1.0:

```
composition(c) = raw_score(c) / Σ_c raw_score(c)
```

Cell types not significantly enriched in a cluster receive composition = 0
for that cluster. **Each cluster's composition scores sum exactly to 1.0.**

### Properties

- **No expression matrix required**: works identically for scRNA-seq and
  spatial data, using only the enrichment test output.
- **Scale-invariant**: the normalization removes absolute magnitude differences
  between clusters with few vs. many significant cell types.
- **Consistent across clusters**: all clusters are tested against the same
  cell type universe, so composition = 0 unambiguously means "not
  significantly enriched" rather than "not tested" (as was ambiguous in the
  per-cluster NNLS approach).
- **Ancestor-pruned**: CL ontological ancestor types are excluded; only the
  most specific (leaf) cell types contribute, making the composition
  biologically interpretable.
- **Priority**: selected_tissue Level 2 results are used when a `--tissue`
  filter was specified; otherwise all_tissue Level 2 results are used.

### Output

The composition DataFrame (clusters × cell types) is saved as
`{job_name}_composition_scores.csv`. The HTML **Composition** tab displays:
1. A stacked bar chart (top 15 cell types by mean composition + "Other" segment)
2. A full per-cluster table of all non-zero composition scores

---

## 8. Input Filtering and Gene Matching

### DEG filtering thresholds

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--pval` | 0.05 | Adj. p-value cutoff for DEG inclusion |
| `--log2fc` | 1.0 | Minimum log2 fold-change (upregulated genes only) |
| `--mean` | 0.0 (auto for spatial) | Minimum mean counts per gene |
| `--topgenes` | None | Restrict to top-N genes by log2FC per cluster |

For spatial data, `mean_counts_thresh` is auto-calibrated to the 75th
percentile of positive mean-counts values (opt-out with `--no_auto_spatial_filter`).

### HGNC gene alias resolution

Gene symbols in the DEG file and marker database are normalized to
UPPERCASE, then mapped to canonical HGNC approved symbols using the HGNC
complete gene set (`hgnc_complete_set.txt`). Both `alias_symbol` and
`prev_symbol` columns are parsed, with first-mapping-wins priority.

This resolves common synonyms: e.g., `TNFRSF1B` ↔ `TNFR2`,
`CD3` → `CD3D/CD3E/CD3G`, historical aliases for rebranded gene symbols.

### Two-level annotation

| Level | Key | Enrichment mode | Purpose |
|-------|-----|-----------------|---------|
| Level 1 | `database + "_" + db_cell_name` | Conventional (unweighted) | Fine-grained per-database hits |
| Level 2 | `cell_name` (CL-normalized) | Weighted | Cross-database consensus |

Level 2 is used for the Top Annotation Summary, hierarchical annotation,
and Combined_Score ranking. Level 1 provides per-database resolution for
audit and interpretability.

---

## 9. Known Limitations and Caveats

### Hypergeometric test assumption

The hypergeometric test assumes:
- DEGs are drawn uniformly at random from the background (null)
- Genes are independent (no co-regulation structure)

In practice, scRNA-seq DEGs violate both: high-variance genes are
over-represented, and co-expressed genes are correlated. The result is
that p-values are anti-conservative (too small). Global FDR correction
partially mitigates this but does not eliminate it.

### N underestimation for standard scRNA-seq

Seurat `FindAllMarkers` typically returns 3,000–8,000 filtered DEGs, not all
~20,000 tested genes. If `N` is set to the DEG count (default behavior when
no `--background_gene_count` is supplied), the hypergeometric test has a
smaller denominator than the true experiment, making p-values too small.

**Mitigation:** Always supply `--background_gene_count 20000` for standard
10x 3′/5′ Chromium data. The pipeline emits a warning when `N < 15,000`
and the DEG pattern suggests a filtered-only file.

### Weighted enrichment distributional theory

`Weighted_Enrichment` and `Combined_Score` are **descriptive ranking metrics**,
not formal test statistics. Their distributions under the null are not
analytically characterized. Do not assign p-values to these metrics.
The only inferential quantity is `adj_p_value` (hypergeometric, BH-corrected).

A formal weighted test (e.g., Wallenius' noncentral hypergeometric or
permutation-based) is planned for a future release.

### Spatial mean_counts_thresh auto-calibration

For spatial data, the 75th percentile mean-counts threshold may suppress
signal in low-expression regions (e.g., necrotic tissue, sparse Xenium panels).
Use `--no_auto_spatial_filter` to disable and set `--mean` manually if
auto-calibration suppresses biologically relevant hits.

### Composition scores reflect enrichment, not absolute proportions

`composition(c)` measures the relative share of enrichment evidence attributed
to each cell type within a cluster. It is **not** a quantitative cell count
fraction. Clusters with a single dominant hit (e.g., one cell type at adj_p < 1e-10
and all others at adj_p > 0.01) will show near-100% composition for that type
even if the underlying biology is mixed. Interpretation should account for the
number of significantly enriched types and their relative Combined_Score spread.

### Composition is sensitive to gene panel size

For targeted assays (Xenium panels, 50–500 genes), fewer marker genes overlap
the database, reducing the power of the hypergeometric test. This may result
in fewer significant cell types and composition scores dominated by the most
broadly expressed markers in the panel.

---

## References

1. Boyle EI et al. GO::TermFinder. *Bioinformatics.* 2004;20(18):3710.
   — Hypergeometric test for gene set enrichment.

2. Benjamini Y, Hochberg Y. Controlling the FDR. *J Royal Stat Soc B.*
   1995;57(1):289–300.

3. Hu C et al. CellMarker 2.0. *Nucleic Acids Res.* 2023;51(D1):D1091.
4. Franzen O et al. PanglaoDB. *Database.* 2019;2019:baz046.
