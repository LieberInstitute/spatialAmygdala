library(here)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(SingleCellExperiment)
library(spacexr)

plots_dir     <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Load Visium SPE
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                    "MarkerGenes",
                    "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe <- spe[!duplicated(rownames(spe)), ]

# --------------------------------------------------------------------------
# Load snRNAseq reference (Yu)
# --------------------------------------------------------------------------
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
counts(sce.amy) <- round(expm1(counts(sce.amy)))
sce <- sce.amy
rm(sce.amy)

sce$ident <- sub("^Human_", "", sce$ident)

# --------------------------------------------------------------------------
# Replace TSHZ1+ idents with MetaNeighbor ITC subtype labels
# --------------------------------------------------------------------------
sce_itc <- readRDS(here("processed-data", "snRNAseq", "BICCN",
                        "sce_ITC_harmony_transferred.rds"))
sce_itc_yu <- sce_itc[, sce_itc$dataset == "Yu"]
message("Yu ITC cells in ITC SCE: ", ncol(sce_itc_yu))

yu_match <- match(sce_itc_yu$original_barcode, colnames(sce))
matched  <- which(!is.na(yu_match))
message("Matched ", length(matched), " Yu ITC cells to SCE barcodes")

message("Original idents being replaced:")
print(table(sce$ident[yu_match[matched]]))

sce$ident <- as.character(sce$ident)
sce$ident[yu_match[matched]] <- as.character(
    sce_itc_yu$mn_cluster_transfer[matched]
)
sce$ident <- factor(sce$ident)

message("Updated idents:")
print(table(sce$ident))

rm(sce_itc, sce_itc_yu)

# --------------------------------------------------------------------------
# Convert SPE -> SpatialRNA
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
# Convert SCE -> Reference
# --------------------------------------------------------------------------
cell_type_counts <- table(sce$ident)
keep_types <- names(cell_type_counts)[cell_type_counts >= 25]
dropped <- setdiff(names(cell_type_counts), keep_types)
if (length(dropped) > 0) {
    message("Dropping cell types with < 25 cells: ",
            paste(dropped, collapse = ", "))
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
# Run RCTD (multi mode, as before)
# --------------------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 30))
message("Running RCTD with ", n_cores, " cores on ", ncol(spe), " spots")

myRCTD <- create.RCTD(puck, reference,
                      MAX_MULTI_TYPES = 5,
                      max_cores = n_cores)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "multi")

saveRDS(myRCTD, file.path(processed_dir, "rctd_results_human_ITCmn.rds"))
message("Done.")
sessioninfo::session_info()