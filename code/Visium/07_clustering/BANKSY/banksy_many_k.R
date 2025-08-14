#!/usr/bin/env Rscript

library("here")
library("SpatialExperiment")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("spatialLIBD")
library("ggplot2")
library("patchwork")
library("Banksy")
library("Seurat")
library("scater")
library("scran")
library("harmony")
library("escheR")
library("data.table")

# Parse k from command line arguments
args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])
if (is.na(k)) stop("You must provide a value for k (number of clusters).")

# Save directories
plot_dir <- here("plots", "07_clustering", "BANKSY")
processed_dir <- here("processed-data","07_clustering")

# Load data
load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))

# Drop unwanted samples
spe <- spe[, !colData(spe)$sample_id %in% c("Br9469", "Br9017", "Br9206")]

# Re-normalize
spe <- logNormCounts(spe)

# ====== Stagger spatial coordinates ========
locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = factor(spe$sample_id))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs

# Load SVGs and subset to top 2000
SVGs.df <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
genes <- SVGs.df$gene_name[1:2000]
spe <- spe[genes, ]

# Parameters
lambda <- 0.4
k_geom <- 18
npcs <- 20
aname <- "logcounts"

# Run Banksy steps
spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)
set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs, group = "capture_area")
spe <- RunHarmony(
  spe,
  group.by.vars = c("sample_id"),
  reduction = "PCA_M0_lam0.4",
  reduction.save = "HARMONY_M0_lam0.4"
)

# Cluster with kmeans
set.seed(1000)
spe <- Banksy::clusterBanksy(
  spe, lambda = lambda, npcs = npcs,
  dimred = "HARMONY_M0_lam0.4",
  algo = "kmeans",
  kmeans.centers = k
)

# Drop duplicated colData columns (safety)
colData(spe) <- colData(spe)[, !duplicated(colnames(colData(spe)))]

# Save
outfile <- here("processed-data", "Visium", "07_clustering", "BANKSY", paste0("spe_banksy_kmeans_slideOnly_", k, ".rds"))
saveRDS(spe, file = outfile)


