library(here)
library(SpatialExperiment)
library(spacexr)
library(SingleCellExperiment)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")

# --------------------------------------------------------------------------
# Load VisiumHD data
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_cells.rds"))
colnames(spe) <- make.unique(colnames(spe))

# --------------------------------------------------------------------------
# Load snRNAseq reference
# --------------------------------------------------------------------------
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
counts(sce.amy) <- round(expm1(counts(sce.amy)))
sce <- sce.amy
rm(sce.amy)

sce$ident <- sub("^Human_", "", sce$ident)

# --------------------------------------------------------------------------
# Replace TSHZ1+ idents with MetaNeighbor cluster labels
# --------------------------------------------------------------------------
sce_itc <- readRDS(here("processed-data", "snRNAseq", "BICCN", "sce_ITC_harmony_transferred.rds"))

# Subset to Yu cells only
sce_itc_yu <- sce_itc[, sce_itc$dataset == "Yu"]
message("Yu ITC cells in ITC SCE: ", ncol(sce_itc_yu))

# Find matching cells in sce and sce_itc $original_barcode
sce_barcodes <- colnames(sce)
sce_itc_barcodes <- sce_itc_yu$original_barcode
yu_match <- match(sce_itc_barcodes, sce_barcodes)
matched <- which(!is.na(yu_match))
message("Matched ", length(matched), " Yu ITC cells to SCE barcodes")
# Matched 3803 Yu ITC cells to SCE barcodes

# Sanity check: what idents are these cells currently?
message("Original idents being replaced:")
print(table(sce$ident[yu_match[matched]]))
# TSHZ1 CALCRL TSHZ1 SEMA3C 
#         2751         1052 

# Overwrite ident with mn_cluster for matched cells
sce$ident <- as.character(sce$ident)
sce$ident[yu_match[matched]] <- as.character(sce_itc_yu$mn_cluster_transfer[matched])
sce$ident <- factor(sce$ident)

message("Updated idents:")
print(table(sce$ident))

rm(sce_itc, sce_itc_yu)

# --------------------------------------------------------------------------
# Convert SPE to SpatialRNA (spacexr format)
# --------------------------------------------------------------------------
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

spatial_counts <- counts(spe)

stopifnot("Duplicate rownames still present in spatial_counts" =
              !any(duplicated(rownames(spatial_counts))))

coords <- as.data.frame(spatialCoords(spe))
if (ncol(coords) == 2) {
    colnames(coords) <- c("x", "y")
} else {
    stop("Expected 2 coordinate columns in spatialCoords(spe)")
}

barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes

stopifnot(
    "Barcode mismatch" =
        identical(rownames(coords), colnames(spatial_counts)) &&
        identical(names(nUMI_spatial), colnames(spatial_counts))
)

puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --------------------------------------------------------------------------
# Convert SCE to Reference (spacexr format)
# --------------------------------------------------------------------------
cell_type_counts <- table(sce$ident)
keep_types <- names(cell_type_counts)[cell_type_counts >= 25]
dropped <- names(cell_type_counts)[cell_type_counts < 25]
if (length(dropped) > 0) {
    message("Dropping cell types with < 25 cells: ", paste(dropped, collapse = ", "))
}
sce <- sce[, sce$ident %in% keep_types]
sce$ident <- droplevels(sce$ident)

dup_ref <- duplicated(rownames(sce))
if (any(dup_ref)) {
    message("Removing ", sum(dup_ref), " duplicated gene names from sce")
    sce <- sce[!dup_ref, ]
}
ref_counts <- counts(sce)

cell_types <- factor(sce$ident)
names(cell_types) <- colnames(sce)

nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- colnames(sce)

reference <- Reference(ref_counts, cell_types, nUMI_ref)

# --------------------------------------------------------------------------
# Run RCTD
# --------------------------------------------------------------------------
myRCTD <- create.RCTD(puck, reference, UMI_min = 25, MAX_MULTI_TYPES = 2, max_cores = 30)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# Save
saveRDS(myRCTD, here(processed_dir, "rctd_results_HDcells_ITCmn.rds"))
message("Done! RCTD results saved.")