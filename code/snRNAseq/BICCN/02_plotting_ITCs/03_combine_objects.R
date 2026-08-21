## 01_join_sce_objects.R
## Join three SCE objects (BICCN, LIBD, Yu) for ITC integration
## - Subset LIBD and Yu to ITC (TSHZ1) populations
## - Harmonize colData across datasets
## - Intersect genes, combine, and save

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    library(scater)
})

# ---- Load data ----
sce.biccn <- readRDS(here("processed-data/snRNAseq/BICCN/WHB_amygdala_Final_ITC_subsets_v2.rds"))
sce.libd  <- readRDS(here("processed-data", "snRNAseq", "sce_amy_human.rds"))
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce.yu <- sce.amy
rm(sce.amy)

message("BICCN: ", ncol(sce.biccn), " cells, ", nrow(sce.biccn), " genes")
message("LIBD:  ", ncol(sce.libd),  " cells, ", nrow(sce.libd),  " genes")
message("Yu:    ", ncol(sce.yu),    " cells, ", nrow(sce.yu),    " genes")

# ---- 1. Assign broad ITC populations in BICCN ----
# Fix: include EMSN_224 (left branch) and EMSN_232 (right branch)
# based on dendrogram hierarchy
EMSN_1 <- c("EMSN_222", "EMSN_223", "EMSN_224", "EMSN_225", "EMSN_226", "EMSN_227")
EMSN_2 <- c("EMSN_228", "EMSN_229", "EMSN_230", "EMSN_231", "EMSN_232",
             "EMSN_233", "EMSN_234", "EMSN_426")

sce.biccn$celltype_broad <- ifelse(sce.biccn$cluster %in% EMSN_1, "ITC_1",
                             ifelse(sce.biccn$cluster %in% EMSN_2, "ITC_2", NA))

stopifnot(sum(is.na(sce.biccn$celltype_broad)) == 0)
message("BICCN ITC grouping: ", paste(names(table(sce.biccn$celltype_broad)),
        table(sce.biccn$celltype_broad), sep = "=", collapse = ", "))

# ---- 2. Subset LIBD to ITC populations ----
itc_libd <- c("TSHZ1_PRKG1", "TSHZ1_CPNE4")
sce.libd <- sce.libd[, sce.libd$fine_celltype %in% itc_libd]

# Map TSHZ1 subtypes to ITC_1 / ITC_2
# TSHZ1_PRKG1 -> ITC_1, TSHZ1_CPNE4 -> ITC_2
sce.libd$celltype_broad <- ifelse(sce.libd$fine_celltype == "TSHZ1_PRKG1", "ITC_1", "ITC_2")

message("LIBD ITC: ", ncol(sce.libd), " cells")
message("  ", paste(names(table(sce.libd$celltype_broad)),
        table(sce.libd$celltype_broad), sep = "=", collapse = ", "))

# ---- 3. Subset Yu to ITC populations ----
itc_yu <- c("Human_TSHZ1 CALCRL", "Human_TSHZ1 SEMA3C")
sce.yu <- sce.yu[, sce.yu$ident %in% itc_yu]

# Map TSHZ1 subtypes to ITC_1 / ITC_2
# TSHZ1 CALCRL -> ITC_1, TSHZ1 SEMA3C -> ITC_2
sce.yu$celltype_broad <- ifelse(sce.yu$ident == "Human_TSHZ1 CALCRL", "ITC_1", "ITC_2")

message("Yu ITC: ", ncol(sce.yu), " cells")
message("  ", paste(names(table(sce.yu$celltype_broad)),
        table(sce.yu$celltype_broad), sep = "=", collapse = ", "))

# ---- 4. Harmonize colData ----
# Keep a minimal, consistent set of columns across datasets

# BICCN: use library_label as sample
sce.biccn$dataset     <- "BICCN"
sce.biccn$sample_id   <- sce.biccn$library_label
sce.biccn$celltype_fine <- sce.biccn$cluster
sce.biccn$donor <- sce.biccn$donor_label
# celltype_broad already set above

# LIBD: use Sample as sample
sce.libd$dataset      <- "LIBD"
sce.libd$sample_id    <- sce.libd$Sample
sce.libd$celltype_fine <- as.character(sce.libd$fine_celltype)
sce.libd$donor <- sce.libd$Sample
# celltype_broad already set above

# Yu: use library_id as sample
sce.yu$dataset        <- "Yu"
sce.yu$sample_id      <- sce.yu$library_id
sce.yu$celltype_fine  <- as.character(sce.yu$ident)
sce.yu$donor <- sce.yu$sample
# celltype_broad already set above

# Save original cell barcodes before renaming
sce.biccn$original_barcode <- colnames(sce.biccn)
sce.libd$original_barcode  <- colnames(sce.libd)
sce.yu$original_barcode    <- colnames(sce.yu)

# Define the columns to keep
cols_keep <- c("dataset", "sample_id", "celltype_broad", "celltype_fine", "original_barcode", "donor")

# Rebuild colData with only harmonized columns
colData(sce.biccn) <- colData(sce.biccn)[, cols_keep]
colData(sce.libd)  <- colData(sce.libd)[, cols_keep]
colData(sce.yu)    <- colData(sce.yu)[, cols_keep]

# ---- 5. Intersect genes ----
# Use gene symbols as the common identifier
# Make sure rownames are gene symbols
message("Intersecting genes across datasets...")

genes_common <- Reduce(intersect, list(
    rownames(sce.biccn),
    rownames(sce.libd),
    rownames(sce.yu)
))
message("Common genes: ", length(genes_common))

sce.biccn <- sce.biccn[genes_common, ]
sce.libd  <- sce.libd[genes_common, ]
sce.yu    <- sce.yu[genes_common, ]

# ---- 6. Rebuild as clean SCEs ----
# Reconstruct from scratch to eliminate GRanges in rowRanges,
# differing metadata, altExps, reducedDims, etc.
# This guarantees cbind will work cleanly.

rd <- DataFrame(gene_symbol = genes_common, row.names = genes_common)

rebuild_sce <- function(sce_in, prefix) {
    new_colnames <- paste0(prefix, "_", seq_len(ncol(sce_in)))
    SingleCellExperiment(
        assays  = list(counts = counts(sce_in)),
        colData = colData(sce_in)[, cols_keep],
        rowData = rd
    ) |> (\(x) { colnames(x) <- new_colnames; x })()
}

sce.biccn <- rebuild_sce(sce.biccn, "BICCN")
sce.libd  <- rebuild_sce(sce.libd,  "LIBD")
sce.yu    <- rebuild_sce(sce.yu,    "Yu")

# ---- 8. Combine ----
message("Combining SCE objects...")
sce.combined <- cbind(sce.biccn, sce.libd, sce.yu)

message("Combined: ", ncol(sce.combined), " cells x ", nrow(sce.combined), " genes")
message("Datasets: ", paste(names(table(sce.combined$dataset)),
        table(sce.combined$dataset), sep = "=", collapse = ", "))
message("ITC types: ", paste(names(table(sce.combined$celltype_broad)),
        table(sce.combined$celltype_broad), sep = "=", collapse = ", "))
message("Samples: ", length(unique(sce.combined$sample_id)))

# ---- 9. Save ----
out_path <- here("processed-data", "snRNAseq", "sce_ITC_combined.rds")
saveRDS(sce.combined, out_path)
message("Saved to: ", out_path)

sessioninfo::session_info()