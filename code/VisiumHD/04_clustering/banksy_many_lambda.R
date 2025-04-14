#!/usr/bin/env Rscript

# Libraries
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(patchwork)
library(scran)
library(scater)

# install banksy is not installed
if (!requireNamespace("Banksy", quietly = TRUE)) {
  BiocManager::install("Banksy")
}

# Get lambda value from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
lambda <- as.numeric(args[1])
cat("Running with lambda =", lambda, "\n")

# Define directories and load data
plot_dir <- here("plots", "VisiumHD", "04_clustering", "banksy_16um")
spe <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "AmyHD_016_pilot.rds"))

# QC Metrics
rownames(spe) <- rowData(spe)$symbol
is.mito <- grepl("^MT-", rownames(spe))
df <- scuttle::perCellQCMetrics(spe, subsets = list(Mito = is.mito))
spe$sum <- df$sum
spe$detected <- df$detected
spe$subsets_mito_percent <- df$subsets_Mito_percent
spe$sum_discard <- isOutlier(spe$sum, nmads = 3, type = "lower", log = TRUE)

# Banksy clustering
spe <- spe[!grepl("^MT-", rownames(spe)), ]
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
dec <- modelGeneVar(spe)
chosen <- getTopHVGs(dec, n = 500)
spe.hvg <- spe[chosen, ]
spe.hvg <- Banksy::computeBanksy(spe.hvg, assay_name = "logcounts", k_geom = 30)
spe.hvg <- Banksy::runBanksyPCA(spe.hvg, lambda = lambda, npcs = 30)
spe.hvg <- Banksy::clusterBanksy(spe.hvg, lambda = lambda, npcs = 30, resolution = 0.8)

# Save results
saveRDS(spe.hvg, here("processed-data", "VisiumHD", "04_clustering", "banksy_16um",
                      paste0("AmyHD_016_pilot_banksy_hvg_lambda_", lambda, ".rds")))
spe$cluster <- spe.hvg[[paste0("clust_M0_lam", lambda, "_k50_res0.8")]]
saveRDS(spe, here("processed-data", "VisiumHD", "04_clustering", "banksy_16um",
                  paste0("AmyHD_016_pilot_banksy_lambda_", lambda, ".rds")))
