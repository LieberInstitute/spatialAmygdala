library(here)
library(SpatialExperiment)
library(spacexr)
library(SingleCellExperiment)

# --------------------------------------------------------------------------
# Parse sample index from SLURM array task id
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID",
                                 if (length(args)) args[1] else NA))
stopifnot(!is.na(task_id))

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Load VisiumHD data and pick this sample
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control",
                    "spe_qc_cells.rds"))
colnames(spe) <- make.unique(colnames(spe))

sample_ids <- sort(unique(as.character(unlist(spe$sample_id))))
stopifnot(task_id >= 1, task_id <= length(sample_ids))
this_sample <- sample_ids[task_id]

message("Task ", task_id, " of ", length(sample_ids), ": ", this_sample)

spe <- spe[, as.character(unlist(spe$sample_id)) == this_sample]
message("SPE subset to ", ncol(spe), " spots")

out_file <- file.path(processed_dir,
                      paste0("rctd_results_HDcells_ITCmn_", this_sample, ".rds"))
if (file.exists(out_file)) {
    message("Output already exists, skipping: ", out_file)
    quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# Load snRNAseq reference (cached after first task builds it)
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

    # ---- Replace TSHZ1+ idents with MetaNeighbor cluster labels ----
    sce_itc <- readRDS(here("processed-data", "snRNAseq", "BICCN",
                            "sce_ITC_harmony_transferred.rds"))
    sce_itc_yu <- sce_itc[, sce_itc$dataset == "Yu"]
    yu_match <- match(sce_itc_yu$original_barcode, colnames(sce))
    matched <- which(!is.na(yu_match))
    message("Matched ", length(matched), " Yu ITC cells")

    sce$ident <- as.character(sce$ident)
    sce$ident[yu_match[matched]] <- as.character(
        sce_itc_yu$mn_cluster_transfer[matched]
    )
    sce$ident <- factor(sce$ident)
    rm(sce_itc, sce_itc_yu)

    # ---- Filter small clusters and dedup ----
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
    if (any(dup_ref)) sce <- sce[!dup_ref, ]

    cell_types <- factor(sce$ident)
    names(cell_types) <- colnames(sce)
    nUMI_ref <- Matrix::colSums(counts(sce))
    names(nUMI_ref) <- colnames(sce)

    reference <- Reference(counts(sce), cell_types, nUMI_ref)
    saveRDS(reference, ref_cache)
    message("Reference cached to ", ref_cache)
    rm(sce)
}

# --------------------------------------------------------------------------
# Convert SPE -> SpatialRNA
# --------------------------------------------------------------------------
dup_genes <- duplicated(rownames(spe))
if (any(dup_genes)) spe <- spe[!dup_genes, ]

spatial_counts <- counts(spe)
coords <- as.data.frame(spatialCoords(spe))
stopifnot(ncol(coords) == 2)
colnames(coords) <- c("x", "y")
rownames(coords) <- colnames(spatial_counts)

nUMI_spatial <- Matrix::colSums(spatial_counts)
names(nUMI_spatial) <- colnames(spatial_counts)

puck <- SpatialRNA(coords, spatial_counts, nUMI_spatial)

# --------------------------------------------------------------------------
# Run RCTD
# --------------------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 8))
message("Running RCTD with ", n_cores, " cores on ", ncol(spe), " spots")

myRCTD <- create.RCTD(puck, reference,
                      UMI_min = 25,
                      MAX_MULTI_TYPES = 2,
                      max_cores = n_cores)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

saveRDS(myRCTD, out_file)
message("Done: ", out_file)

sessioninfo::session_info()