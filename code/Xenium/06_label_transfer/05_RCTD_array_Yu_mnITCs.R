library(here)
library(dplyr)
library(SpatialExperiment)
library(SingleCellExperiment)
library(spacexr)

processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# --- Get sample index from SLURM array task id (fall back to argv) ---
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                 if (length(args)) args[1] else NA))
stopifnot(!is.na(task_id))

# --- Load xenium data ---
spe <- readRDS(here("processed-data", "Xenium", "04_clustering", "Banksy",
                    "Banksy_integrated_res2.0_collapsed_v6.rds"))
colnames(spe) <- make.unique(colnames(spe))

sample_ids <- sort(unique(spe$brnum))
stopifnot("Task ID out of range" = task_id >= 1 && task_id <= length(sample_ids))

current_sample <- sample_ids[task_id]
message("Processing sample: ", current_sample,
        " (task ", task_id, " of ", length(sample_ids), ")")

spe <- spe[, spe$brnum == current_sample]
message("Cells in this sample: ", ncol(spe))

out_file <- file.path(processed_dir,
                      paste0("rctd_results_ITCmn_", current_sample, ".rds"))
if (file.exists(out_file)) {
    message("Output already exists, skipping: ", out_file)
    quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# Load / build Reference (cached across array tasks)
# --------------------------------------------------------------------------
ref_cache <- file.path(processed_dir, "rctd_reference_yu_itcmn.rds")

if (file.exists(ref_cache)) {
    message("Loading cached Reference: ", ref_cache)
    reference <- readRDS(ref_cache)
} else {
    message("Building Reference from scratch")

    load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
    counts(sce.amy) <- round(expm1(counts(sce.amy)))
    sce <- sce.amy
    rm(sce.amy)

    sce$ident <- sub("^Human_", "", sce$ident)

    # ---- Replace TSHZ1+ idents with MetaNeighbor ITC subtype labels ----
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

    # ---- Drop small clusters + dedup genes ----
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
    saveRDS(reference, ref_cache)
    message("Reference cached to ", ref_cache)
    rm(sce)
}

# --------------------------------------------------------------------------
# Prepare SpatialRNA for this sample
# --------------------------------------------------------------------------
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) {
    message("Removing ", sum(dup_genes), " duplicated gene names from spe")
    spe <- spe[!dup_genes, ]
}

spatial_counts <- counts(spe)
stopifnot(!any(duplicated(rownames(spatial_counts))))

coords <- as.data.frame(spatialCoords(spe))
stopifnot(ncol(coords) == 2)
colnames(coords) <- c("x", "y")
barcodes <- colnames(spatial_counts)
rownames(coords) <- barcodes

nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- barcodes

puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --------------------------------------------------------------------------
# Run RCTD
# --------------------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 30))
message("Running RCTD with ", n_cores, " cores on ", ncol(spe), " spots")

myRCTD <- create.RCTD(puck, reference,
                      UMI_min = 10,
                      MAX_MULTI_TYPES = 2,
                      max_cores = n_cores)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

saveRDS(myRCTD, out_file)
message("Saved: ", out_file)

sessioninfo::session_info()