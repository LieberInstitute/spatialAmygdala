#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("SpatialExperiment")
  library("scuttle")
  library("scran")
  library("scater")
  library("here")
  library("dplyr")
  library("patchwork")
  library("Banksy")
  library("harmony")
  library("data.table")
})

# ========================
# Parse resolution argument
# ========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 02_banksy_many_harmony_many_lambda_resSweep.R <resolution>")
res <- as.numeric(args[[1]])
res_label <- gsub("\\.", "_", as.character(res))  # for filenames

# ============== Setup ==============

processed_dir <- here("processed-data", "Xenium", "04_clustering", "Banksy")
plot_dir <- here("plots", "Xenium", "05_spatial_clustering")

spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds"))

# ======= Feature Selection =======
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe.gex <- spe[gene_expression_idx,]

# ======= Spatial Coord adjustment =======
locs <- spatialCoords(spe.gex)
locs <- cbind(locs, sample_id = factor(spe.gex$brnum))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe.gex) <- locs

# ======= BANKSY Parameters =======
lambda <- 0.8
k_geom <- 36
npcs <- 50
aname <- "nucleus_normcounts"

# ======= Run BANKSY and Harmony =======
set.seed(1000)
spe.gex <- Banksy::computeBanksy(spe.gex, assay_name = aname, k_geom = k_geom)
spe.gex <- Banksy::runBanksyPCA(spe.gex, lambda = lambda, npcs = npcs)

# rename PCA
reducedDim(spe.gex, "PCA") <- reducedDim(spe.gex, "PCA_M0_lam0.8")
reducedDim(spe.gex, "PCA_M0_lam0.8") <- NULL

# Harmony
spe.gex <- RunHarmony(spe.gex, "brnum", reduction.save = "HARMONY_M0_lam0.8")

# ========== Clustering ==========
spe.gex <- Banksy::clusterBanksy(
  spe.gex,
  lambda = lambda,
  npcs = npcs,
  resolution = res,
  dimred = "HARMONY_M0_lam0.8"
)

# Transfer cluster labels to original spe
clust_col <- paste0("clust_Banksy_k", k_geom, "_res", res_label, "_lam", lambda)
spe[[clust_col]] <- spe.gex[[clust_col]]

# ========== Save ==========
saveRDS(spe, file = here(processed_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res_label, ".rds")))
saveRDS(spe.gex, file = here(processed_dir, paste0("Banksy_integrated_lambda_", lambda, "_res", res_label, "_gex.rds")))
