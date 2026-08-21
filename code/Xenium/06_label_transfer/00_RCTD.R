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

processed_dir <- here("processed-data", "Xenium", "06_label_transfer")

# load xenium data
spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds"))
spe

# make colnames unique
colnames(spe) <- make.unique(colnames(spe))

# load the snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce

# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="fine_celltype")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="doublet")

# save
saveRDS(res, here(processed_dir, "rctd_results.rds"))