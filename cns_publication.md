# SapiensOntoCellMap: Publication Strategy & Competitive Landscape

## Target Venue
**Primary:** Nature Methods or Genome Biology (tool/resource paper)
**Alternative:** Bioinformatics (application note) or Nucleic Acids Research (database/web server)

---

## 1. The Problem Statement

Cell type annotation remains a critical bottleneck in single-cell and spatial
transcriptomics. Despite 30+ published tools, three fundamental gaps persist:

1. **No tool unifies markers across databases with formal ontology normalization.**
   Existing marker databases (CellMarker 2.0, PanglaoDB, CellTalkDB, etc.) each
   use different nomenclature for the same cell types and tissues. A user must
   manually reconcile "CD8+ T cell", "cytotoxic T lymphocyte", and
   "CL:0000625 CD8-positive, alpha-beta T cell" as the same entity.

2. **No marker-based tool provides ontology-aware hierarchical annotation.**
   Current tools return flat ranked lists. When a cluster matches both "T cell"
   and "CD8+ alpha-beta T cell", existing tools treat these as independent
   annotations. The user gets no guidance on annotation resolution — whether
   to call a cluster "T cell" (safe but vague) or "CD8+ effector memory T cell"
   (specific but risky).

3. **No tool handles both scRNA-seq and spatial transcriptomics natively.**
   Reference-based tools (Azimuth, CellTypist) require matched reference atlases.
   Marker-based tools (scType, SCSA) are designed only for scRNA-seq. Spatial
   platforms (Visium HD, Xenium, MERFISH) have distinct DEG formats, smaller gene
   panels, and different statistical properties that existing tools ignore.

---

## 2. Competitive Landscape

### 2.1 Marker-Based Enrichment Tools (Direct Competitors)

| Tool | Databases | Ontology | Spatial | Hierarchical | Multi-DB Evidence | Year | Journal |
|------|-----------|----------|---------|-------------|-------------------|------|---------|
| **SapiensOntoCellMap** | **14+** | **CL + UBERON** | **Yes** | **Planned** | **Yes (weighted)** | 2026 | — |
| CellMarker 2.0 | 1 (own) | Partial (MeSH) | No | No | No | 2023 | NAR |
| scType | 1 (curated) | No | No | No | No | 2022 | Nat Comms |
| SCSA | 3 (CellMarker, CancerSEA, PanglaoDB) | No | No | No | No | 2020 | Front Genet |
| CellAssign | User-supplied | No | No | No | No | 2019 | Nat Methods |
| Garnett | User-supplied | No | No | No | No | 2019 | Nat Methods |
| ClusterMole | 3 (CellMarker, PanglaoDB, PANTHER) | Partial | No | No | No | 2024 | Bioinformatics |
| EasyCellType | 4 (CellMarker, PanglaoDB, SaVanT, SCINA) | No | No | No | Majority vote | 2024 | Bioinformatics |

**Key differentiators of SapiensOntoCellMap:**
- 14+ databases vs. 1-4 in competitors (4-14x broader coverage)
- Formal CL/UBERON normalization (no other tool does both)
- Native spatial support (Visium HD, Xenium, Space Ranger DEG format)
- Cross-database evidence weighting (source-type + agreement scoring)

### 2.2 Reference-Based / Transfer Learning Tools (Indirect Competitors)

| Tool | Method | Requires Reference | Spatial | Scalability | Year | Journal |
|------|--------|-------------------|---------|-------------|------|---------|
| Azimuth/Seurat v5 | Label transfer (RPCA) | Yes (atlas) | Limited | ~100K cells | 2024 | Nat Biotech |
| CellTypist | Logistic regression | Yes (atlas) | Yes | ~1M cells | 2022 | Science |
| SingleR | Correlation-based | Yes (reference) | No | ~50K cells | 2019 | Nat Immunol |
| scArches | Transfer learning (VAE) | Yes (model) | Yes | ~500K cells | 2022 | Nat Biotech |
| scGPT | Foundation model | Pre-trained | Yes | ~100K cells | 2024 | Nat Methods |

**Why these are NOT direct competitors (but complementary):**
- Reference-based tools require a matched atlas for each tissue/disease context.
  For rare tissues, developmental stages, or disease-perturbed states, suitable
  references often don't exist.
- They are "black box" — the user cannot inspect which markers drove the annotation.
- Marker-based tools are interpretable, auditable, and work without reference data.
- **Complementary use case:** SapiensOntoCellMap validates/explains reference-based
  annotations by showing which known markers support each call.

### 2.3 Ontology-Aware Tools (Closest Conceptual Competitors)

| Tool | Ontology Use | Method | Hierarchical Output | Status |
|------|-------------|--------|-------------------|--------|
| **SapiensOntoCellMap** | **CL + UBERON for normalization + planned hierarchy** | **Enrichment** | **Planned (multi-resolution)** | **Active** |
| OnClass | CL graph embedding | ML classifier + ontology propagation | Yes (ancestors scored) | 2021, Nat Comms |
| CellO | CL hierarchy | Hierarchical classifier (one per CL node) | Yes (multi-level) | 2021, iScience |
| Cell Ontology (OBO) | — | Standard vocabulary | — | Ongoing |

**OnClass** (Sheng et al., 2021) is the closest conceptual competitor. It embeds CL
terms in a graph and propagates predictions to ancestor nodes. However:
- It requires reference data (expression profiles mapped to CL terms)
- It cannot use multiple marker databases
- It does not work with spatial DEG data
- The ontology propagation is a post-hoc embedding, not a formal statistical test

**CellO** (Kimmel & Kelley, 2021) trains one classifier per CL node and enforces
hierarchy consistency. Limitations:
- Requires pre-trained models per tissue
- No marker interpretability
- No spatial support

**SapiensOntoCellMap's planned hierarchical engine** would be unique in combining:
1. Marker-based enrichment (interpretable, no reference needed)
2. CL ontology graph traversal (formal ancestor/descendant propagation)
3. Multi-resolution output with confidence scoring
4. Evidence from 14+ databases aggregated at each ontology level

---

## 3. Knowledge Gap & Unique Contribution

### The Gap
```
           Marker databases exist (CellMarker, PanglaoDB, ...)
                        |
                        v
        BUT they are fragmented, use inconsistent nomenclature,
        and no tool unifies them with ontology normalization
                        |
                        v
     Enrichment tools exist (scType, SCSA, ClusterMole, ...)
                        |
                        v
        BUT they use 1-4 databases, no ontology hierarchy,
        no spatial support, no cross-database evidence weighting
                        |
                        v
     Ontology tools exist (OnClass, CellO, ...)
                        |
                        v
        BUT they require reference data, are not marker-based,
        and cannot leverage multi-database evidence
                        |
                        v
    ┌─────────────────────────────────────────────────────────┐
    │  SapiensOntoCellMap fills the intersection:              │
    │  Multi-DB + Ontology-normalized + Enrichment-based       │
    │  + Spatial-native + Hierarchical (planned)               │
    └─────────────────────────────────────────────────────────┘
```

### Unique Contributions (Publishable Claims)

1. **Largest unified cell marker database** with formal CL/UBERON normalization
   (14+ sources, >200K marker-cell type associations, all mapped to CL terms)

2. **First marker-based tool with ontology-aware hierarchical annotation**
   (multi-resolution output: "definitely T cell, probably CD8+, possibly effector memory")

3. **First enrichment tool natively supporting spatial transcriptomics**
   (Visium HD, Xenium, MERFISH differential expression formats)

4. **Source-aware evidence weighting** across databases
   (experimental > scRNA-seq > literature > computational)

5. **Interactive clinical-grade reports** with lazy-loaded heatmaps, violin
   distributions, and exportable result tables in a single self-contained HTML

---

## 4. Proposed Features for Publication Readiness

### P0: Must-Have for Submission

| Feature | Status | Impact | Effort |
|---------|--------|--------|--------|
| Global FDR correction | DONE | Correctness | Complete |
| Robust N detection (scRNA-seq) | DONE (warning) | Correctness | Complete |
| Gene name normalization (case + alias) | DONE (case) | Recall | Alias mapping: 1 week |
| Minimum overlap filter (k >= 2) | DONE | Precision | Complete |
| Hierarchical annotation engine | NOT STARTED | **Key differentiator** | 2-3 weeks |
| Benchmarking suite (vs scType, SCSA, CellTypist) | NOT STARTED | Mandatory | 2 weeks |
| Unit test suite | NOT STARTED | Reviewer requirement | 1 week |

### P1: Strongly Recommended

| Feature | Why | Effort |
|---------|-----|--------|
| Source-aware weighted enrichment | Methodology claim | 1 week |
| Scanpy/MAST/edgeR input support | Broadens user base 10x | 3 days |
| Confidence scoring per annotation | "How sure are you?" | 1 week |
| Ambiguity detection ("T cell vs NK cell") | Clinical relevance | 1 week |

### P2: Nice-to-Have (Impact Multipliers)

| Feature | Why | Effort |
|---------|-----|--------|
| pip installable package | Citations | 3 days |
| Scanpy integration (`sc.tl.annotate_sapiensonto()`) | Discovery | 1 week |
| Web interface (Streamlit/Shiny) | Accessibility for wet-lab users | 1 week |
| Batch mode with multi-sample summary | Core facility use case | 3 days |
| Export to CellxGene / HCA format | Standards compliance | 3 days |

---

## 5. Benchmarking Strategy

### 5.1 Datasets

| Dataset | Source | Why |
|---------|--------|-----|
| PBMC 10x (Zheng et al.) | Well-characterized ground truth | Standard benchmark in every paper |
| Human Lung Cell Atlas | Tabula Sapiens | Complex tissue, many cell types |
| Mouse brain (Allen Brain Atlas) | Cross-species generalization | Tests HGNC normalization |
| Skin scRNA-seq (in-house) | Domain expertise + spatial complement | Unique to this paper |
| Visium HD skin (in-house) | Spatial benchmarking | No other tool benchmarks spatial |
| Xenium skin panel (in-house) | Panel-based spatial | Tests small-N robustness |

### 5.2 Metrics

| Metric | What it measures |
|--------|-----------------|
| Accuracy (exact match) | Fraction of clusters correctly annotated at finest level |
| Hierarchical accuracy | Correct at any ancestor level (CD8 T cell matches if ground truth is effector CD8) |
| Resolution score | How specific are annotations? (T cell = 1, CD8+ effector memory = 4) |
| Cross-tool agreement | How often does SapiensOntoCellMap agree with Azimuth/CellTypist? |
| Interpretability | Can a biologist understand WHY each annotation was made? |
| Runtime | Wall-clock time per dataset |

### 5.3 Comparators

| Tool | Why include |
|------|------------|
| scType | Most popular marker-based tool (1,500+ citations) |
| CellTypist | State-of-the-art reference-based (800+ citations) |
| Azimuth | Gold standard for PBMC/common tissues |
| SingleR | Classic reference-based baseline |
| Manual expert annotation | Ground truth for in-house datasets |

---

## 6. Figures Plan (8-10 figures)

### Main Figures
1. **Pipeline overview schematic** — Database aggregation, ontology normalization,
   enrichment testing, hierarchical annotation, report generation
2. **Database coverage** — UpSet plot showing overlap of cell types across 14+ databases;
   bar chart of unique markers per source
3. **Ontology normalization impact** — Sankey diagram: raw terms -> CL-normalized terms;
   quantify how many synonyms collapse to canonical CL IDs
4. **Benchmarking: PBMC** — Confusion matrix vs. manual annotation; comparison
   with scType, CellTypist, Azimuth
5. **Hierarchical annotation showcase** — Tree visualization showing multi-resolution
   annotation for a complex tissue (lung or brain); confidence score at each depth
6. **Spatial annotation** — Visium HD and Xenium results overlaid on tissue images;
   comparison with manual pathologist annotation

### Supplementary Figures
7. **Source-aware weighting** — Effect of evidence weighting on annotation quality
8. **Sensitivity analysis** — Effect of k-threshold, FDR method, N estimation
9. **Runtime comparison** — Benchmarking speed across dataset sizes
10. **Interactive report screenshots** — Heatmap, violin plots, results table

---

## 7. Manuscript Outline

### Title Options
- "SapiensOntoCellMap: Ontology-aware cell type annotation from unified marker
  databases for single-cell and spatial transcriptomics"
- "Multi-database ontology-normalized marker enrichment for hierarchical cell
  type annotation in scRNA-seq and spatial transcriptomics"

### Abstract (150 words)
Cell type annotation is essential for interpreting single-cell and spatial
transcriptomics. Existing marker-based tools rely on 1-4 databases with
inconsistent nomenclature and lack ontological structure. Reference-based tools
require matched atlases and offer no marker-level interpretability. We present
SapiensOntoCellMap, which unifies 14+ marker databases under Cell Ontology (CL)
and Uberon (UBERON) normalization, enabling the largest harmonized marker
knowledge base for human cell types. Our hierarchical annotation engine traverses
the CL ontology graph to provide multi-resolution annotations with confidence
scores, answering "Is this a CD8+ T cell, or can I only say T cell?" We natively
support Visium HD, Xenium, and scRNA-seq inputs with statistically rigorous
enrichment testing (global FDR, gene normalization, source-aware weighting).
Benchmarking on PBMC, lung, and skin datasets demonstrates superior accuracy and
interpretability compared to scType, CellTypist, and Azimuth.

### Sections
1. Introduction (problem + gap)
2. Results
   - 2.1 Unified marker database with ontology normalization
   - 2.2 Enrichment-based annotation with global FDR
   - 2.3 Hierarchical ontology-aware annotation
   - 2.4 Source-aware evidence weighting
   - 2.5 Benchmarking against existing tools
   - 2.6 Application to spatial transcriptomics
   - 2.7 Interactive reporting for clinical workflows
3. Discussion
4. Methods

---

## 8. Critical Path to Submission

### Phase 1: Core Engine (Weeks 1-3)
- [ ] Implement hierarchical annotation engine (CL graph traversal)
- [ ] Add confidence scoring per annotation level
- [ ] Implement HGNC alias resolution for gene normalization
- [ ] Support Scanpy rank_genes_groups input format

### Phase 2: Benchmarking (Weeks 4-5)
- [ ] Download and process PBMC 10x benchmark dataset
- [ ] Run SapiensOntoCellMap + scType + CellTypist + Azimuth on same data
- [ ] Compute accuracy, hierarchical accuracy, resolution score
- [ ] Generate comparison figures

### Phase 3: Writing & Figures (Weeks 6-7)
- [ ] Draft manuscript (Methods first, then Results, then Intro/Discussion)
- [ ] Generate all main + supplementary figures
- [ ] Package as pip-installable
- [ ] Create documentation site

### Phase 4: Polish & Submit (Week 8)
- [ ] Internal review
- [ ] Prepare supplementary materials
- [ ] Submit to Nature Methods / Genome Biology

---

## 9. Literature References

### Marker Databases
- Hu C, et al. CellMarker 2.0. *Nucleic Acids Res.* 2023;51(D1):D1091-D1097.
- Franzen O, et al. PanglaoDB. *Database.* 2019;2019:baz046.
- Zhang X, et al. CellMarker. *Nucleic Acids Res.* 2019;47(D1):D721-D728.

### Cell Annotation Tools
- Ianevski A, et al. scType. *Nat Commun.* 2022;13:4027.
- Cao Y, et al. SCSA. *Front Genet.* 2020;11:490.
- Zhang AW, et al. CellAssign. *Nat Methods.* 2019;16:1007-1015.
- Pliner HA, et al. Garnett. *Nat Methods.* 2019;16:983-986.
- Bhuva DD, et al. ClusterMole. *Bioinformatics.* 2024;40(1):btad726.

### Reference-Based Tools
- Hao Y, et al. Seurat v5 / Azimuth. *Nat Biotechnol.* 2024;42:293-304.
- Dominguez Conde C, et al. CellTypist. *Science.* 2022;376:eabl5197.
- Aran D, et al. SingleR. *Nat Immunol.* 2019;20:163-172.
- Lotfollahi M, et al. scArches. *Nat Biotechnol.* 2022;40:235-244.

### Ontology-Aware Tools
- Sheng J, et al. OnClass. *Nat Commun.* 2021;12:5556.
- Kimmel JC, Kelley DR. CellO. *iScience.* 2021;24(1):101913.

### Foundation Models
- Cui H, et al. scGPT. *Nat Methods.* 2024;21:1470-1480.

### Statistical Methods
- Boyle EI, et al. GO::TermFinder. *Bioinformatics.* 2004;20(18):3710-3715.
- Benjamini Y, Hochberg Y. FDR. *J Royal Stat Soc B.* 1995;57(1):289-300.

### Reviews
- Pasquini G, et al. Automated methods for cell type annotation. *Nat Methods.* 2021;18:1124-1126.
- Abdelaal T, et al. Comparison of automatic cell identification methods. *Nat Biotechnol.* 2019;37:411-419.
- Clarke ZA, et al. Tutorial: cell type annotation guidelines. *Nat Protoc.* 2021;16:3212-3232.
