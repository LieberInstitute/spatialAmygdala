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
plots_dir <- here("plots", "VisiumHD", "04_clustering")
processed_dir <- here("processed-data", "VisiumHD", "04_clustering")

# load
sfe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "sfe_qc.rds"))
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))


# ==== conver to spe ====
# get assays and such

 spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = as(counts(sfe), "dgCMatrix")),
    rowData = rowData(sfe),
    colData = cbind(colData(sfe), spatialCoords(sfe)),
    spatialCoordsNames = spatialCoordsNames(sfe),
    spatialCoords = spatialCoords(sfe)
  )

# ========== RCTD ==========
rctd_data <- createRctd(spe, sce, cell_type_col="fine_celltype")
res <- runRctd(rctd_data, max_cores=30, rctd_mode="doublet")

# save
saveRDS(res, here(processed_dir, "rctd_results.rds"))