library(NMFscape)
library(here)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

# set random seed for reproducibility
set.seed(12345)

# load data
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy

sce <- as.SingleCellExperiment(rat.amy)


# ====== run NMF ========

# Run NMF with 50 factors
sce <- runNMFscape(sce, k = 50, verbose = FALSE)

# Save NMF results
saveRDS(sce, here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50.rds"))