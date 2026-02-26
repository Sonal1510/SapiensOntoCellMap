"""
PBMC3k Ground Truth Definitions
================================
Standard benchmark: Zheng et al. 2017, 2,638 PBMCs, 9 Louvain clusters, 8 cell types.
Ground truth labels are taken from the Seurat PBMC3k tutorial (Hao et al. 2021),
which is the de-facto standard for scRNA-seq annotation benchmarks.

Reference: Hao et al., Nat Biotech 2021 (doi:10.1038/s41587-021-01033-z)
"""

# Ground truth: cluster label → (canonical_name, CL_ID, broad_lineage)
# Keys are exact Louvain cluster label strings from scanpy.datasets.pbmc3k_processed()
PBMC3K_GROUND_TRUTH: dict[str, tuple[str, str, str]] = {
    "CD4 T cells":       ("CD4-positive, alpha-beta T cell",  "CL:0000624", "T cell"),
    "CD14+ Monocytes":   ("classical monocyte",               "CL:0000860", "Monocyte"),
    "B cells":           ("B cell",                           "CL:0000236", "B cell"),
    "CD8 T cells":       ("CD8-positive, alpha-beta T cell",  "CL:0000625", "T cell"),
    "NK cells":          ("natural killer cell",              "CL:0000623", "NK cell"),
    "FCGR3A+ Monocytes": ("non-classical monocyte",           "CL:0000875", "Monocyte"),
    "Dendritic cells":   ("myeloid dendritic cell",           "CL:0000451", "Dendritic cell"),
    "Megakaryocytes":    ("megakaryocyte",                    "CL:0000556", "Megakaryocyte"),
}

# Acceptable fuzzy aliases for Top-1 string matching.
# Keys: ground truth labels (lowercased); values: list of acceptable substrings.
GT_ALIASES: dict[str, list[str]] = {
    "cd4 t cells":       ["cd4", "helper t", "helper cd4", "alpha-beta"],
    "cd14+ monocytes":   ["monocyte", "cd14", "classical monocyte"],
    "b cells":           ["b cell", "b-cell"],
    "cd8 t cells":       ["cd8", "cytotoxic", "killer"],
    "fcgr3a+ monocytes": ["monocyte", "non-classical", "fcgr3a", "cd16"],
    "nk cells":          ["natural killer", "nk cell"],
    "dendritic cells":   ["dendritic", " dc"],
    "megakaryocytes":    ["megakaryocyte", "platelet"],
}

# CL IDs that constitute each broad lineage (used for broad-type accuracy)
CL_LINEAGE_MAP: dict[str, set[str]] = {
    "T cell":          {"CL:0000084", "CL:0000624", "CL:0000625", "CL:0000798"},
    "B cell":          {"CL:0000236", "CL:0000785", "CL:0000946"},
    "NK cell":         {"CL:0000623", "CL:0000825"},
    "Monocyte":        {"CL:0000576", "CL:0000860", "CL:0000875"},
    "Dendritic cell":  {"CL:0000451", "CL:0000990", "CL:0000782"},
    "Megakaryocyte":   {"CL:0000556"},
}
