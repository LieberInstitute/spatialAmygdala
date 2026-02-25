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
spe <- readRDS(here("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
spe

# make colnames unique
colnames(spe) <- make.unique(colnames(spe))

# load the snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce

# --------------------------------------------------------------------------
# Convert SPE to SpatialRNA (spacexr format)
# --------------------------------------------------------------------------
## Remove any remaining duplicate gene names
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

## Extract raw counts
spatial_counts <- counts(spe)

## Final safety check on the extracted matrix
stopifnot("Duplicate rownames still present in spatial_counts" =
              !any(duplicated(rownames(spatial_counts))))

## Extract spatial coordinates as a data.frame with columns x, y
coords <- as.data.frame(spatialCoords(spe))
if (ncol(coords) == 2) {
    colnames(coords) <- c("x", "y")
} else {
    stop("Expected 2 coordinate columns in spatialCoords(spe)")
}

## Ensure rownames on coords match colnames of counts
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

## nUMI per spot/cell
nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes

## Sanity check
stopifnot(
    "Barcode mismatch" =
        identical(rownames(coords), colnames(spatial_counts)) &&
        identical(names(nUMI_spatial), colnames(spatial_counts))
)

## Create SpatialRNA object
puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --------------------------------------------------------------------------
# Convert SCE to Reference (spacexr format)
# --------------------------------------------------------------------------
## Remove duplicate genes from reference
dup_ref <- duplicated(rownames(sce))
if (any(dup_ref)) {
    message("Removing ", sum(dup_ref), " duplicated gene names from sce")
    sce <- sce[!dup_ref, ]
}

ref_counts <- counts(sce)

## Cell type labels
cell_types <- factor(sce$fine_celltype)
names(cell_types) <- colnames(sce)

## nUMI per cell
nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- colnames(sce)

## Create Reference object
reference <- Reference(ref_counts, cell_types, nUMI_ref)

# ========== RCTD ==========
myRCTD <- create.RCTD(puck, reference, UMI_min = 10, MAX_MULTI_TYPES = 5, max_cores = 30)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# save
saveRDS(myRCTD, here(processed_dir, "rctd_results_xenium.rds"))