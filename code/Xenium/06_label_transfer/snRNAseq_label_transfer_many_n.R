#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
de_n <- as.numeric(args[1])

library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")
library("BiocParallel")

# Load data
spe <- readRDS(here("processed-data", "Xenium", "04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
sce <- readRDS(here("processed-data", "snRNAseq", "sce.human_all_genes.rds"))

# Normalize
logcounts(sce) <- logNormCounts(sce)
logcounts(spe) <- assay(spe, "cell_normcounts")

# Subset to common genes
common_genes <- intersect(rownames(spe), rownames(sce))
spe <- spe[common_genes, ]
sce <- sce[common_genes, ]

# Run SingleR
pred.broad <- SingleR(test = logcounts(spe), ref = logcounts(sce), labels = sce$broad_celltype,
                      de.method = "wilcox", BPPARAM = MulticoreParam(20), de.n = de_n)

pred.fine <- SingleR(test = logcounts(spe), ref = logcounts(sce), labels = sce$fine_celltype,
                     de.method = "wilcox", BPPARAM = MulticoreParam(20), de.n = de_n)

# Save output
out_dir <- here("processed-data", "Xenium", "06_label_transfer")
write.csv(pred.broad, file = file.path(out_dir, paste0("xenium_5um_pred_broad_celltype_n", de_n, ".csv")), row.names = FALSE)
write.csv(pred.fine, file = file.path(out_dir, paste0("xenium_5um_pred_fine_celltype_n", de_n, ".csv")), row.names = FALSE)
