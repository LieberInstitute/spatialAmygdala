suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(data.table)
})

# -----------------------------
# Inputs
# -----------------------------
sce_path <- here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50.rds")

# MapMyCells output files (edit these paths)
tsv_path  <- here("processed-data","snRNAseq","Zhou_cocaine", "MapMyCells_outs",
                  "MapMyCells_HVG3000_countscsvgz_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1769098702954.csv")

# -----------------------------
# Load SCE
# -----------------------------
sce <- readRDS(sce_path)

# -----------------------------
# Read MapMyCells annotations
# -----------------------------
ann <- fread(tsv_path)  # handles .gz automatically
stopifnot("cell_id" %in% names(ann))

# Make sure IDs are character
ann[, cell_id := as.character(cell_id)]


# -----------------------------
# Align rows to sce colnames and add to colData
# -----------------------------
# Match MapMyCells rows to SCE columns
m <- match(colnames(sce), ann$cell_id)

# Quick diagnostics
message("SCE cells: ", ncol(sce))
message("MapMyCells rows: ", nrow(ann))
message("Matched cells: ", sum(!is.na(m)))
if (anyNA(m)) {
  warning("Some SCE colnames were not found in MapMyCells cell_id. ",
          "They will get NA annotations in colData.")
}

# Create an aligned annotation table (same order/length as SCE columns)
ann_aligned <- ann[m]

# Drop the redundant cell_id column before binding into colData
ann_aligned[, cell_id := NULL]

# Add columns into colData (keeps existing colData)
# data.table -> data.frame to play nicely with S4Vectors
new_cols <- as.data.frame(ann_aligned, stringsAsFactors = FALSE, check.names = FALSE)

# Add with a prefix to avoid name collisions (optional but recommended)
prefix <- "mmc_"
names(new_cols) <- paste0(prefix, names(new_cols))

# Bind into colData
for (nm in names(new_cols)) {
  colData(sce)[[nm]] <- new_cols[[nm]]
}

# -----------------------------
# Save updated SCE
# -----------------------------
out_sce <- here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50_with_MapMyCells.rds")
saveRDS(sce, out_sce)
message("Saved updated SCE: ", out_sce)
