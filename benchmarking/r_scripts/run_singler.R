#!/usr/bin/env Rscript
# =============================================================================
# SingleR Wrapper Script
# =============================================================================
# Runs SingleR annotation on a .h5ad file via zellkonverter/BPCells.
# Produces per-cluster majority-vote predictions.
#
# Reference: Aran et al., Nat Immunol 2019 (doi:10.1038/s41590-018-0276-y)
#
# Usage: Rscript run_singler.R <h5ad_path> <output_csv> <reference> <cluster_key>
#
# Arguments:
#   h5ad_path    — path to AnnData .h5ad file
#   output_csv   — output path for predictions (cluster, predicted_cell_type)
#   reference    — celldex reference name:
#                  HumanPrimaryCellAtlasData | MonacoImmuneData | BlueprintEncodeData
#   cluster_key  — colData (obs) column for cluster assignment (default: louvain)
#
# Outputs:
#   output_csv   — two-column CSV: cluster, predicted_cell_type
#
# Requirements:
#   BiocManager::install(c("SingleR", "celldex", "zellkonverter"))
# =============================================================================

suppressMessages({
  library(SingleR)
  library(celldex)
  library(zellkonverter)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_singler.R <h5ad_path> <output_csv> <reference> [cluster_key]")
}

h5ad_path   <- args[1]
output_csv  <- args[2]
ref_name    <- args[3]
cluster_key <- if (length(args) >= 4) args[4] else "louvain"

cat("[SingleR] Loading h5ad:", h5ad_path, "\n")
sce <- readH5AD(h5ad_path, use_hdf5 = TRUE)

# Load reference dataset
cat("[SingleR] Loading reference:", ref_name, "\n")
ref <- switch(ref_name,
  HumanPrimaryCellAtlasData = HumanPrimaryCellAtlasData(),
  MonacoImmuneData           = MonacoImmuneData(),
  BlueprintEncodeData        = BlueprintEncodeData(),
  stop(paste("[SingleR] Unknown reference:", ref_name))
)

# Run SingleR (per-cell)
cat("[SingleR] Running per-cell annotation...\n")
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)

# Aggregate to cluster level (majority vote)
if (!cluster_key %in% colnames(colData(sce))) {
  stop(paste("[SingleR] Cluster key not found:", cluster_key,
             "\nAvailable:", paste(colnames(colData(sce)), collapse = ", ")))
}
cluster_labels <- as.character(colData(sce)[[cluster_key]])
df_cell <- data.frame(
  cluster   = cluster_labels,
  predicted = pred$labels,
  stringsAsFactors = FALSE
)
majority <- tapply(df_cell$predicted, df_cell$cluster, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})

# Write output CSV
out_df <- data.frame(
  cluster             = names(majority),
  predicted_cell_type = as.character(majority),
  stringsAsFactors    = FALSE
)
write.csv(out_df, output_csv, row.names = FALSE)
cat("[SingleR] Wrote predictions to:", output_csv, "\n")
