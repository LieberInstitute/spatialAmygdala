library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)
library(spacexr)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_relabel_ITC_smoothed.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

# drop
load(here("processed-data", "snRNAseq", "yu_sce_gtf.rda"))
sce <- sce.amy

# remove Human_ prefix
sce$ident <- gsub("Human_","", sce$ident)

# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="ident")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="multi", max_multi_types=5)

# save
saveRDS(res, here(processed_dir, "rctd_results_Yu_etal.rds"))