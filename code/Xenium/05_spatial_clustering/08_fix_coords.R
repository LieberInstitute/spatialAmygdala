#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
})

# ======= Parameters =======
lambda <- 0.8
resolutions <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
input_dir <- here("processed-data", "Xenium", "04_clustering", "Banksy")
ref_path <- here("processed-data", "Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds")

# ======= Load original coordinates =======
spe_ref <- readRDS(ref_path)
coords_orig <- spatialCoords(spe_ref)

# ======= Loop over files =======
for (res in resolutions) {
  res_label <- gsub("\\.", "_", as.character(res))
  rds_path <- file.path(input_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res_label, "_gex.rds"))

  if (!file.exists(rds_path)) {
    warning("File not found: ", rds_path)
    next
  }

  cat("Fixing:", rds_path, "\n")
  spe <- readRDS(rds_path)

  # reset spatial coordinates to original values
  spatialCoords(spe) <- coords_orig[colnames(spe), ]

  # Overwrite with fixed version
  saveRDS(spe, file = rds_path)
}
