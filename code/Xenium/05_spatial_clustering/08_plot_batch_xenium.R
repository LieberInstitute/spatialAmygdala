#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(ggspavis)
  library(scCustomize)
  library(patchwork)
})

# ======== Parameters ========
lambda <- 0.8
k_geom <- 50
resolutions <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

input_dir <- here("processed-data", "Xenium", "04_clustering", "Banksy")
coord_reset_path <- here("processed-data", "Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds")
output_dir <- here("plots", "Xenium", "05_spatial_clustering")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load coordinates from original unmodified SPE
spe_orig <- readRDS(coord_reset_path)
coords_orig <- spatialCoords(spe_orig)

file_rds <- file.path(input_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", 0.6, "_gex.rds"))

# Loop over resolutions
for (res in resolutions) {
  res_label <- gsub("\\.", "_", as.character(res))  # for filenames
  clust_col <- paste0("clust_HARMONY_M0_lam", lambda, "_k", k_geom, "_res", res_label)
  file_rds <- file.path(input_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res_label, "_gex.rds"))

  if (!file.exists(file_rds)) {
    warning("File does not exist: ", file_rds)
    next
  }

  message("Processing resolution: ", res)

  # Load and reset spatialCoords
  spe <- readRDS(file_rds)
  spatialCoords(spe) <- coords_orig[colnames(spe), ]

  # Skip if cluster column missing
  if (!clust_col %in% colnames(colData(spe))) {
    warning("Missing cluster column: ", clust_col)
    next
  }

  # Get unique cluster count and colors
  num_groups <- length(unique(spe[[clust_col]]))
  colors <- scCustomize::scCustomize_Palette(
    num_colors = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )

  # Plot with ggspavis
  p <- ggspavis::plotCoords(spe,
    annotate = clust_col,
    in_tissue = NULL,
    sample_id = "brnum"
  ) +
    scale_color_manual(values = colors) +
    ggtitle(paste0("BANKSY λ=", lambda, ", res=", res)) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

  # Save
  pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res, ".pdf"))
  ggsave(pdf_path, plot = p, width = 40, height = 10, useDingbats = FALSE)

  message("Saved to: ", pdf_path)
}
