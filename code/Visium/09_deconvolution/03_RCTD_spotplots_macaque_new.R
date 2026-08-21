library(here)
library(ggplot2)
library(patchwork)
library(spacexr)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))


# load
myRCTD <- readRDS(here(processed_dir, "rctd_results_macaque.rds"))


## multi mode ## 
print('Examining multi mode all weights results')

results = myRCTD@results
weights = lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
norm_weights <- normalize_weights(weights_df)
coords = myRCTD@spatialRNA@coords

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights_df)


celltypes <- levels(myRCTD@reference@cell_types)
plot_dir <- here("plots", "Visium", "09_deconvolution", "RCTD_macaque")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

samples <- unique(spe$sample_id)

for (s in samples) {
  print(s)
  plot_list <- lapply(celltypes, function(i) {
    print(i)
    spatialLIBD::vis_gene(spe, sampleid = s, geneid = i, is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
  })
  
  pdf(file.path(plot_dir, paste0(s, "_multi_allWeight.pdf")), width = 7, height = 7)
  for (i in seq_along(plot_list)) {
    print(plot_list[[i]])
  }
  dev.off()
}
