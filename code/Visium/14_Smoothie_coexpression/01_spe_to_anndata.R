#!/usr/bin/env Rscript
## =============================================================================
## 01_spe_to_anndata.R
##
## Convert per-donor SpatialExperiment subsets directly to .h5ad files
## using zellkonverter. Each donor gets its own AnnData with raw counts
## and spatial coordinates ready for Smoothie.
## =============================================================================

library(here)
library(SpatialExperiment)
library(zellkonverter)
library(SingleCellExperiment)

## ---- Load SPE ----
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

## Reorder donors: anterior → posterior
spe$sample_id <- factor(spe$sample_id,
                         levels = c("Br9192", "Br9280", "Br2743", "Br6471",
                                    "Br8325", "Br6423", "Br6660"))

## ---- Create output directory ----
out_dir <- here("processed-data", "Visium", "14_smoothie_coexpression", "h5ad")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Inspect ----
cat("SPE dimensions:", dim(spe), "\n")
cat("Assays:", assayNames(spe), "\n")
cat("spatialCoords columns:", colnames(spatialCoords(spe)), "\n")
cat("Donors:", levels(spe$sample_id), "\n\n")

## ---- Export per-donor h5ad ----
donors <- levels(spe$sample_id)

for (donor in donors) {
    cat("Processing", donor, "...\n")

    ## Subset
    spe_sub <- spe[, spe$sample_id == donor]

    ## Keep only the raw counts assay (Smoothie will normalize itself)
    ## Convert to a basic SCE so zellkonverter writes it cleanly
    sce <- SingleCellExperiment(
        assays = list(X = assay(spe_sub, "counts")),
        colData = colData(spe_sub)
    )

    ## Store spatial coordinates in obsm$spatial (AnnData convention)
    ## spatialCoords returns pxl_col_in_fullres, pxl_row_in_fullres
    coords <- spatialCoords(spe_sub)
    reducedDim(sce, "spatial") <- coords

    ## Write h5ad
    h5ad_path <- file.path(out_dir, paste0(donor, ".h5ad"))
    writeH5AD(sce, file = h5ad_path)

    cat("  Saved:", h5ad_path,
        " (", ncol(sce), "spots x", nrow(sce), "genes)\n")
}

cat("\nDone! All .h5ad files in:", out_dir, "\n")
cat("Next: run 02_smoothie_multi_dataset.py\n")