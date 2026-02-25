library(here)
library(SingleCellExperiment)
library(scater)
library(scran)
library(RcppICA)
library(Seurat)

# set random seed for reproducibility
set.seed(12345)

# load data
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy

sce <- as.SingleCellExperiment(rat.amy)



# run ICA with 50 components
system.time(
ica_result <- RcppICA::fastICA(t(logcounts(sce)), n.comp = 50)
)

# save ica_result
saveRDS(ica_result, here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_ICA50_result.rds"))