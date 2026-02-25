library(zellkonverter)
library(SingleCellExperiment)
library(here)
library(data.table)
library(Matrix)

# --- Load SCEs ---
sce <- readH5AD(here("processed-data/snRNAseq/01_build_sce/WHB_amygdala_ITC_subset.h5ad"))
sce
# class: SingleCellExperiment 
# dim: 59357 66072 
# metadata(0):
# assays(1): X
# rownames(59357): ENSG00000000003 ENSG00000000005 ... ENSG00000288643
#   ENSG00000288645
# rowData names(1): gene_symbol
# colnames(66072): 10X386_2:CCCTCAAAGTCCCAAT 10X383_5:ATGCCTCCATCCTATT
#   ... 10X218_5:CGACAGCAGATTAGTG 10X204_3:TTTCACATCGACGACC
# colData names(3): cell_barcode library_label anatomical_division_label
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

# load metadata
metadata_path <- here("processed-data/snRNAseq/abc_atlas/abc_cell_metadata_extended.csv")
metadata <- fread(metadata_path)
metadata <- as.data.frame(metadata)

# Sanity check: is the matching column present?
stopifnot("cell_barcode" %in% colnames(metadata))
stopifnot("cell_barcode" %in% colnames(colData(sce)))

# Match metadata to the order of columns in SCE
merged_meta <- metadata[match(sce$cell_barcode, metadata$cell_barcode), ]

# Sanity check
stopifnot(all(sce$cell_barcode == merged_meta$cell_barcode))

# Combine existing colData with extended metadata
colData(sce) <- cbind(colData(sce), merged_meta)

# relabel X assay to counts
assayNames(sce)[which(assayNames(sce) == "X")] <- "counts"

# save
saveRDS(sce, here("processed-data/snRNAseq/01_build_sce/WHB_amygdala_ITC_subset_SCE.rds"))
