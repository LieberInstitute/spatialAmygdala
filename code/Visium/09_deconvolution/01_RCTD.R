library(here)
library(SpatialFeatureExperiment)
library(spacexr)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

# load
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce


# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="fine_celltype")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="multi", max_multi_types=5)

# save
saveRDS(res, here(processed_dir, "rctd_results.rds"))