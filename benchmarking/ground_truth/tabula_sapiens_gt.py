"""
Tabula Sapiens Ground Truth Definitions
========================================
Benchmarking subsets: Skin and Blood compartments.
Ground truth labels are from expert-curated atlas annotations in:
  The Tabula Sapiens Consortium, Science 2022 (doi:10.1126/science.abl4896).

Data access: CELLxGENE (https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)
"""

# Ground truth: atlas cell type label → (canonical_name, CL_ID, broad_lineage)
# Blood compartment
TABULA_SAPIENS_BLOOD_GT: dict[str, tuple[str, str, str]] = {
    "cd4-positive, alpha-beta t cell":  ("CD4-positive, alpha-beta T cell",  "CL:0000624", "T cell"),
    "cd8-positive, alpha-beta t cell":  ("CD8-positive, alpha-beta T cell",  "CL:0000625", "T cell"),
    "b cell":                           ("B cell",                           "CL:0000236", "B cell"),
    "classical monocyte":               ("classical monocyte",               "CL:0000860", "Monocyte"),
    "non-classical monocyte":           ("non-classical monocyte",           "CL:0000875", "Monocyte"),
    "natural killer cell":              ("natural killer cell",              "CL:0000623", "NK cell"),
    "plasmablast":                      ("plasmablast",                      "CL:0000946", "B cell"),
    "platelet":                         ("platelet",                         "CL:0000233", "Megakaryocyte"),
    "erythrocyte":                      ("erythrocyte",                      "CL:0000232", "Erythrocyte"),
    "neutrophil":                       ("neutrophil",                       "CL:0000775", "Granulocyte"),
}

# Skin compartment
TABULA_SAPIENS_SKIN_GT: dict[str, tuple[str, str, str]] = {
    "keratinocyte":                     ("keratinocyte",                     "CL:0000312", "Keratinocyte"),
    "fibroblast of dermis":             ("fibroblast",                       "CL:0000057", "Fibroblast"),
    "endothelial cell":                 ("endothelial cell",                 "CL:0000115", "Endothelial"),
    "cd4-positive, alpha-beta t cell":  ("CD4-positive, alpha-beta T cell",  "CL:0000624", "T cell"),
    "cd8-positive, alpha-beta t cell":  ("CD8-positive, alpha-beta T cell",  "CL:0000625", "T cell"),
    "macrophage":                       ("macrophage",                       "CL:0000235", "Macrophage"),
    "mast cell":                        ("mast cell",                        "CL:0000097", "Mast cell"),
    "melanocyte":                       ("melanocyte",                       "CL:0000148", "Melanocyte"),
    "pericyte cell":                    ("pericyte",                         "CL:0000669", "Pericyte"),
    "schwann cell":                     ("Schwann cell",                     "CL:0002573", "Neural"),
}

# Combined ground truth (blood + skin) used by benchmark_tabula_sapiens.py
TABULA_SAPIENS_GROUND_TRUTH: dict[str, tuple[str, str, str]] = {
    **TABULA_SAPIENS_BLOOD_GT,
    **TABULA_SAPIENS_SKIN_GT,
}
