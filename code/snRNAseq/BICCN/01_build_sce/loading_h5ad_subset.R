## -------------------------------
## SCRIPT TO LOAD H5AD SUBSET INTO SEURAT
## -------------------------------

# Load required packages
library(rhdf5)
library(data.table)
library(Matrix)
library(Seurat)
library(here)

# -------------------------------
# 1. Define file paths
# -------------------------------

# Adjust this to point to the H5AD expression matrix
h5_path <- here("processed-data/snRNAseq/abc_atlas/expression_matrices/WHB-10Xv3/20240330/WHB-10Xv3-Neurons-raw.h5ad")

# Adjust this to point to the extended metadata file
metadata_path <- here("processed-data/snRNAseq/abc_atlas/metadata/abc_cell_metadata_extended.csv")

# -------------------------------
# 2. Read metadata
# -------------------------------
metadata <- fread(metadata_path)
metadata <- as.data.frame(metadata)

# -------------------------------
# 3. Define subset: e.g., Amygdala excitatory neurons
# -------------------------------
subset_label <- "Amygdala excitatory"
subset_cells <- metadata$cell_label[metadata$supercluster == subset_label]
length(subset_cells)

# -------------------------------
# 4. Read sample (cell) and gene names from HDF5
# -------------------------------

# Cell barcodes / sample names
sample_names <- h5read(h5_path, "/obs/cell_label")

# Gene symbols (e.g., SLC17A7, GAD1, etc.)
gene_names <- h5read(h5_path, "/var/gene_symbol/categories")

# Optional: If you want Ensembl IDs instead
# gene_ids <- h5read(h5_path, "/var/gene_identifier")

# -------------------------------
# 5. Match subsetted cells to sample_names
# -------------------------------
read_samples <- match(subset_cells, sample_names)
read_samples <- read_samples[!is.na(read_samples)]
length(read_samples)  # Final number of cells to read

# -------------------------------
# 6. Read expression matrix subset
# -------------------------------
system.time({
  counts_subset <- h5read(
    h5_path,
    name = "/X/data",
    index = list(NULL, read_samples)  # rows = genes, cols = selected cells
  )
  counts_subset <- t(counts_subset)  # now cells × genes
  subcounts <- as(counts_subset, "dgCMatrix")  # convert to sparse matrix
})

# -------------------------------
# 7. Add gene and sample names
# -------------------------------
rownames(subcounts) <- subset_cells[match(rownames(subcounts), sample_names[read_samples])]
colnames(subcounts) <- gene_names

# -------------------------------
# 8. Match metadata
# -------------------------------
meta <- metadata[match(rownames(subcounts), metadata$cell_label), ]
stopifnot(nrow(meta) == nrow(subcounts))  # Sanity check

# -------------------------------
# 9. Create Seurat object
# -------------------------------
seu <- CreateSeuratObject(counts = t(subcounts))  # Seurat expects genes × cells
seu <- AddMetaData(seu, meta)

# -------------------------------
# 10. Save output or inspect
# -------------------------------
saveRDS(seu, file = here("processed-data/seurat_objects/WHB_amygdala_excitatory_subset.rds"))
print(seu)
