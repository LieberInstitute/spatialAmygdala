library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)
library(SEraster)
library(here)

# set directories
plots_dir <- here("plots", "Xenium", "06_label_transfer")
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")

# load xenium data
spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe

# load
res <- readRDS(here(processed_dir, "rctd_results.rds"))

ws <- assay(res, "weights")
# table(colSums(ws) == 0)
#  FALSE   TRUE 
# 421643  60876 

table(res$first_type)



# ========= Plottig first predicted cell type ==========
imgData(spe) <- NULL # subsetting is throwing an error due to imgData

# subset to common cells (colData) between spe and res
common_cells <- intersect(colnames(spe), rownames(colData(res)))
spe <- spe[, common_cells]
res <- res[, common_cells]

spe$first_type <- colData(res)$first_type

colnames(colData(spe))


# =======================
# SEraster 
# =======================


rastCt <- SEraster::rasterizeCellType(spe,
                                      col_name = "first_type",
                                      resolution = 100,
                                      fun = "sum",
                                      square = FALSE)
