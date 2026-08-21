library(here)
library(ggplot2)
library(patchwork)
library(spacexr)
library(SpatialExperiment)
library(scales)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
myRCTD <- readRDS(here(processed_dir, "rctd_results_human.rds"))

## multi mode ##
print('Examining multi mode all weights results')

results <- myRCTD@results
weights <- lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
coords <- myRCTD@spatialRNA@coords

# restore barcode rownames (rbind drops them) and normalize
rownames(weights_df) <- rownames(coords)
norm_weights <- normalize_weights(weights_df)
norm_weights <- as.data.frame(as.matrix(norm_weights))
rownames(norm_weights) <- rownames(coords)

# align spe to coords order
spe <- spe[ , colnames(spe) %in% rownames(coords)]
spe <- spe[ , match(rownames(coords), colnames(spe))]

# drop "Human " prefix from cell type names to match reference image labels
colnames(norm_weights) <- gsub("^Human_", "", colnames(norm_weights))

stopifnot(identical(colnames(spe), rownames(norm_weights)))

# central cell types -> color channels (matches reference image)
central_types <- c("DRD2.ISL1", "DRD2.PAX6", "PRKCD")
channel_cols <- c(
  DRD2.ISL1 = "#E69F00",  # orange
  DRD2.PAX6 = "#0072B2",  # blue
  PRKCD     = "#458224"   # reddish purple / magenta
)
stopifnot(all(central_types %in% colnames(norm_weights)))

# assemble plotting data frame: colData + central weights + coords
cd <- as.data.frame(colData(spe))
cd <- cbind(cd, norm_weights[, central_types, drop = FALSE])  # positional, aligned

xy <- spatialCoords(spe)
cd$.x <- xy[, 1]
cd$.y <- xy[, 2]

# build composite RGB color per spot
#  - each type drives one channel proportional to its weight
#  - alpha encodes total central content so background spots fade to grey
# build color per spot
#  - grey if max central weight < threshold
#  - otherwise: color of the highest-weighted central type, scaled by that weight
threshold <- 0.10   # from the distribution; try 0.05 / 0.15 to taste

w <- as.matrix(cd[, central_types])
w[is.na(w)] <- 0




library(spatialLIBD)

# winning central type as a factor; sub-threshold -> "none"
win_idx  <- max.col(w, ties.method = "first")
win_val  <- w[cbind(seq_len(nrow(w)), win_idx)]
win_type <- central_types[win_idx]

pass <- win_val >= threshold
central_call <- ifelse(pass, win_type, "none")
central_call <- factor(central_call, levels = c(central_types, "none"))

# attach to spe (vis_clus reads from colData)
spe$central_call <- central_call

# palette: real types get their channel color, "none" is grey
clus_colors <- c(channel_cols[central_types], none = "#D1D1D1")
names(clus_colors)[seq_along(central_types)] <- central_types

plot_dir <- here("plots", "Visium", "09_deconvolution", "RCTD_human")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

samples <- unique(spe$sample_id)

for (s in samples) {
  print(s)
  p <- vis_clus(
    spe,
    sampleid    = s,
    clustervar  = "central_call",
    colors      = clus_colors,
    is_stitched = TRUE,
    ... = " RCTD central types",      # title suffix; or drop this arg
    spatial     = TRUE                  # draw histology background
  )
  ggsave(file.path(plot_dir, paste0(s, "_central_call_hist.pdf")),
         p, width = 8, height = 7)
}



for (s in samples) {
  plot_list <- lapply(central_types, function(ct) {
    grad <- colorRampPalette(c("#D1D1D1", channel_cols[[ct]]))(10)
    vis_gene(spe, sampleid = s, geneid = ct, is_stitched = TRUE,
             cont_colors = grad, spatial = TRUE)
  })
  pdf(file.path(plot_dir, paste0(s, "_central_perType_hist.pdf")), width = 7, height = 7)
  for (p in plot_list) print(p)
  dev.off()
}