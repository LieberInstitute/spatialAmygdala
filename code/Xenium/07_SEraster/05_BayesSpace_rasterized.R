suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("BayesSpace")
    library("RColorBrewer")
    library("ggplot2")
    library("gridExtra")
    library("patchwork")
})

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

processed_dir <- here("processed-data", "Xenium", "07_SEraster")
plots_dir <- here("plots", "Xenium", "07_SEraster", "BayesSpace")

# load 9017, 9192, 9206, 9280
spe <- readRDS(here(processed_dir, "spe_rasterized_harmony_259.rds"))

# ====== Introducing offset ======

# Get unique sample IDs
sample_ids <- sort(unique(spe$sample_id))
n_samples <- length(sample_ids)

# Define how many samples per row in the grid
samples_per_row <- ceiling(sqrt(n_samples))  # e.g., 3 samples => 2x2 grid

# Grid spacing (padding between samples)
row_spacing <- 25000
col_spacing <- 15000

spe$array_col <- spatialCoords(spe)[, 1]
spe$array_row <- spatialCoords(spe)[, 2]

# Offset coordinates into a 2D grid
for (i in seq_along(sample_ids)) {
  id <- sample_ids[i]
  
  row_idx <- (i - 1) %/% samples_per_row
  col_idx <- (i - 1) %% samples_per_row
  
  sel <- spe$sample_id == id
  colData(spe)$row[sel] <- spe$array_row[sel] + row_idx * row_spacing
  colData(spe)$col[sel] <- spe$array_col[sel] + col_idx * col_spacing
}

spatial_coords <- data.frame(
    row = colData(spe)$row,
    col = colData(spe)$col,
    sample_id = colData(spe)$sample_id
)

pdf(here(plots_dir, "offset_check.pdf"), width = 10, height = 8)
ggplot(spatial_coords, aes(x = col, y = row, color = sample_id)) +
    geom_point() +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Spatial Coordinates with Sample Offsets") +
    theme_minimal() +
    coord_fixed()
dev.off()

metadata(spe)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

message("Running spatialCluster()")
Sys.time()
set.seed(2)
spe <- spatialCluster(spe, use.dimred = "HARMONY", q = k, nrep=10000, burn.in=100, gamma=3)
Sys.time()

saveRDS(spe, here(processed_dir, paste0("spe_BayesSpace_clustered_k", k, ".rds")))

bayesSpace_name <- paste0("BayesSpace_", k)
colnames(colData(spe))[ncol(colData(spe))] <- bayesSpace_name

cluster_export(
    spe,
    bayesSpace_name,
    cluster_dir = here::here("processed-data","Xenium", "07_SEraster", "BayesSpace")
)

