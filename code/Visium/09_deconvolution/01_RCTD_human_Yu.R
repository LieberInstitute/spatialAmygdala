library(here)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(spacexr)
library(patchwork)
library(ggplot2)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

# load
load(here("processed-data", "snRNAseq", "yu_sce_gtf.rda"))
sce.amy

counts(sce.amy) <- round(expm1(counts(sce.amy)))

sce <- sce.amy

# --------------------------------------------------------------------------
# Convert SPE to SpatialRNA (spacexr format)
# --------------------------------------------------------------------------
## Remove any remaining duplicate gene names (belt-and-suspenders)
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
## Ensure column names are exactly "x" and "y"
if (ncol(coords) == 2) {
    colnames(coords) <- c("x", "y")
} else {
    stop("Expected 2 coordinate columns in spatialCoords(spe)")
}

## Ensure rownames on coords match colnames of counts
## (spatialCoords may or may not carry rownames depending on the SPE version)
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

## nUMI per spot — names must also match
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
## Drop cell types with < 25 cells to avoid issues with RCTD
cell_type_counts <- table(sce$ident)
keep_types <- names(cell_type_counts)[cell_type_counts >= 25]
sce <- sce[, sce$ident %in% keep_types]

## Remove duplicate genes from reference as well
dup_ref <- duplicated(rownames(sce))
if (any(dup_ref)) {
    message("Removing ", sum(dup_ref), " duplicated gene names from sce")
    sce <- sce[!dup_ref, ]
}

ref_counts <- counts(sce)

## Cell type labels — adjust "ident" to match your colData column
cell_types <- factor(sce$ident)
names(cell_types) <- colnames(sce)

## nUMI per cell
nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- colnames(sce)

## Create Reference object
reference <- Reference(ref_counts, cell_types, nUMI_ref)


# ========== RCTD ==========
myRCTD <- create.RCTD(puck, reference, MAX_MULTI_TYPES = 5, max_cores = 30)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "multi")


# save
saveRDS(myRCTD, here(processed_dir, "rctd_results_human.rds"))