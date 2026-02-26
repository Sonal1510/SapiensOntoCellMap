#!/usr/bin/env Rscript
# =============================================================================
# scType Wrapper Script
# =============================================================================
# Runs scType annotation on a Seurat FindAllMarkers DEG CSV.
# scType has no CRAN/Bioconductor package; it is sourced from GitHub at runtime.
#
# Reference: Ianevski et al., Nat Commun 2022 (doi:10.1038/s41467-022-28803-w)
#
# Usage: Rscript run_sctype.R <deg_csv> <output_csv> <tissue_type>
#
# Arguments:
#   deg_csv      — Seurat FindAllMarkers CSV (must have: cluster, gene, avg_log2FC, p_val_adj)
#   output_csv   — output path for predictions (cluster, predicted_cell_type)
#   tissue_type  — scType tissue string, e.g. "Immune system", "Skin", "Brain"
#
# Outputs:
#   output_csv   — two-column CSV: cluster, predicted_cell_type
#
# Requirements (install in R before running):
#   install.packages(c("HGNChelper", "openxlsx"))
# =============================================================================

suppressMessages({
  library(HGNChelper)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_sctype.R <deg_csv> <output_csv> <tissue_type>")
}

deg_csv     <- args[1]
output_csv  <- args[2]
tissue_type <- args[3]

cat("[scType] Loading DEG file:", deg_csv, "\n")
markers <- read.csv(deg_csv, stringsAsFactors = FALSE)

# Validate required columns
required_cols <- c("cluster", "gene")
missing <- setdiff(required_cols, colnames(markers))
if (length(missing) > 0) {
  stop(paste("[scType] Missing columns:", paste(missing, collapse = ", ")))
}

# Source scType functions from GitHub
cat("[scType] Sourcing scType from GitHub...\n")
source_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/"
suppressMessages({
  source(paste0(source_url, "gene_sets_prepare.r"))
  source(paste0(source_url, "sctype_score_.r"))
})

# Prepare gene sets (uses built-in scType marker database)
cat("[scType] Preparing gene sets for tissue:", tissue_type, "\n")
gs_list <- gene_sets_prepare(
  path_to_db_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
  cell_type = tissue_type
)

# Build a cluster × gene matrix from DEG (top genes by log2FC per cluster)
clusters <- unique(markers$cluster)
es_max_rows <- list()
for (cl in clusters) {
  cl_genes <- markers[markers$cluster == cl, ]
  cl_genes <- cl_genes[order(-cl_genes$avg_log2FC), ]
  top_genes <- head(cl_genes$gene, 200)
  score_vec <- sapply(names(gs_list$gs_positive), function(ct) {
    pos <- length(intersect(gs_list$gs_positive[[ct]], top_genes))
    neg <- length(intersect(gs_list$gs_negative[[ct]], top_genes))
    pos - neg
  })
  es_max_rows[[as.character(cl)]] <- score_vec
}
es_max <- do.call(rbind, es_max_rows)

# Pick best cell type per cluster
best_ct <- apply(es_max, 1, function(row) {
  if (max(row) <= 0) "Unknown" else names(which.max(row))
})

# Write output CSV
out_df <- data.frame(
  cluster            = names(best_ct),
  predicted_cell_type = as.character(best_ct),
  stringsAsFactors   = FALSE
)
write.csv(out_df, output_csv, row.names = FALSE)
cat("[scType] Wrote predictions to:", output_csv, "\n")
