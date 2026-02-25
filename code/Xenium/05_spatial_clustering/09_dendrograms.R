#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(dreamlet)
})

# Parameters
lambda <- 0.8
k_geom <- 50
resolutions <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

input_dir <- here("processed-data", "Xenium", "04_clustering", "Banksy")
output_dir <- here("plots", "Xenium", "06_dendrograms")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over resolutions
for (res in resolutions) {
  cluster_col <- paste0("clust_HARMONY_M0_lam", lambda, "_k", k_geom, "_res", res)
  file_rds <- file.path(input_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res, "_gex.rds"))

  if (!file.exists(file_rds)) {
    warning("Missing file: ", file_rds)
    next
  }

  spe <- readRDS(file_rds)
  if (!cluster_col %in% colnames(colData(spe))) {
    warning("Missing cluster column: ", cluster_col)
    next
  }

12

  # aggregate by cluster + sample
  pb <- aggregateToPseudoBulk(
    spe,
    assay = "counts",
    cluster_id = cluster_col,
    sample_id = "sample_id",
    verbose = FALSE
  )

  # build dendrogram
  hcl <- buildClusterTreeFromPB(pb, method = "ward.D")

  # save plot
  pdf_path <- file.path(output_dir, paste0("Dendrogram_Banksy_lambda_", lambda, "_res", res, ".pdf"))
  pdf(pdf_path, width = 8, height = 4)
  plot(hcl, hang = -1, main = paste0("res = ", res))
  dev.off()

  message("Saved dendrogram: ", pdf_path)
}
