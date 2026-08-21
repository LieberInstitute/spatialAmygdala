library(here)
library(dplyr)
library(SpatialExperiment)
library(spacexr)

processed_dir <- here("processed-data", "Xenium", "06_label_transfer")

# --- Get sample index from command line ---
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

# --- Load xenium data ---
spe <- readRDS(here("processed-data", "Xenium", "04_clustering", "Banksy",
                     "Banksy_integrated_res2.0_collapsed_v6.rds"))
colnames(spe) <- make.unique(colnames(spe))

# --- Get unique sample IDs and subset to this task's sample ---
sample_ids <- sort(unique(spe$brnum))
stopifnot("Task ID out of range" = task_id >= 1 && task_id <= length(sample_ids))

current_sample <- sample_ids[task_id]
message("Processing sample: ", current_sample, " (task ", task_id, " of ", length(sample_ids), ")")

spe <- spe[, spe$brnum == current_sample]
message("Cells in this sample: ", ncol(spe))

# --- Load and prepare reference ---
load(here("processed-data", "snRNAseq", "sce_FINAL_human.rda"))
sce <- rda.human
sce$ident <- sce$fine_celltype

## Drop cell types with < 25 cells
cell_type_counts <- table(sce$ident)
keep_types <- names(cell_type_counts)[cell_type_counts >= 25]
sce <- sce[, sce$ident %in% keep_types]

## Remove duplicate genes from reference
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

# --- Prepare SpatialRNA for this sample ---
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

spatial_counts <- counts(spe)
stopifnot(!any(duplicated(rownames(spatial_counts))))

coords <- as.data.frame(spatialCoords(spe))
colnames(coords) <- c("x", "y")
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes

puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --- Run RCTD ---
myRCTD <- create.RCTD(puck, reference, UMI_min = 10, MAX_MULTI_TYPES = 2, max_cores = 30)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# --- Save per-sample result ---
out_file <- here(processed_dir, paste0("rctd_results_", current_sample, "_totty.rds"))
saveRDS(myRCTD, out_file)
message("Saved: ", out_file)