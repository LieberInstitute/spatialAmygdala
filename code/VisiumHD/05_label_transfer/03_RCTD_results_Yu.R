library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggspavis)
library(spacexr)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")
plot_dir <- here("plots", "VisiumHD", "05_label_transfer")

# load xenium data
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_cells.rds"))
spe

spe <- logNormCounts(spe)

# load
myRCTD <- readRDS(here(processed_dir, "rctd_results_HDcells.rds"))


## doublet mode ##
print('Examining doublet mode results')

# Extract only what we need
results_df <- myRCTD@results$results_df
coords <- myRCTD@spatialRNA@coords
celltypes <- levels(myRCTD@reference@cell_types)

# Build the index in one step without copying spe
idx <- match(rownames(results_df), colnames(spe))
idx <- idx[!is.na(idx)]  # drop any non-matches
spe <- spe[, idx]

# Align results_df to the subsetted spe
results_df <- results_df[colnames(spe), ]

# Add doublet results to colData
colData(spe)$first_type <- results_df$first_type
colData(spe)$second_type <- results_df$second_type
colData(spe)$spot_class <- results_df$spot_class
colData(spe)$first_class_weight <- results_df$first_class  # weight for first type

# If you want to also reconstruct a full weight matrix (useful for per-celltype spatial plots):
# Create a weight matrix with a column per cell type, filling in doublet weights
celltypes <- levels(myRCTD@reference@cell_types)
weight_mat <- matrix(0, nrow = nrow(results_df), ncol = length(celltypes),
                     dimnames = list(rownames(results_df), celltypes))

for (i in seq_len(nrow(results_df))) {
  ct1 <- as.character(results_df$first_type[i])
  ct2 <- as.character(results_df$second_type[i])
  w1 <- results_df$first_class[i]  # weight of first type
  weight_mat[i, ct1] <- w1
  weight_mat[i, ct2] <- 1 - w1
}

weight_df <- as.data.frame(weight_mat)
colnames(weight_df) <- gsub(" ", ".", colnames(weight_df))
colData(spe) <- cbind(colData(spe), weight_df)

celltypes_clean <- gsub(" ", ".", celltypes)

plot_dir <- here("plots", "VisiumHD", "05_label_transfer", "RCTD_human_Yu")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

samples <- unique(spe$sample_id)


# Build a plotting data frame once
plot_df <- as.data.frame(spatialCoords(spe))
plot_df <- cbind(plot_df, as.data.frame(colData(spe)[, c("sample_id", celltypes_clean)]))

for (s in samples) {
  print(s)
  sub_df <- plot_df[plot_df$sample_id == s, ]
  
  plot_list <- lapply(celltypes_clean, function(ct) {
    print(ct)
    ggplot(sub_df, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = .data[[ct]])) +
      geom_point(size = 0.3) +
      scale_color_gradientn(colors = viridisLite::rocket(10, direction = -1)) +
      labs(title = ct, color = "Weight") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  pdf(file.path(plot_dir, paste0(s, "_doublet_weights_Yu.pdf")), width = 7, height = 7)
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]])
  }
  dev.off()
}