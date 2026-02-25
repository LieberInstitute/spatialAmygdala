library("here")
library("SpatialExperiment")
library("SummarizedExperiment")
library("SingleCellExperiment")
library("spatialLIBD")
library("ggplot2")
library("patchwork")
library("Banksy")
library("Seurat")
library("scater")
library("scran")
library("harmony")
library("escheR")
library("data.table")
library("PCAtools")


# Save directories
plot_dir = here("plots","Visium", "07_clustering","BANKSY")
processed_dir = here("processed-data","07_clustering")

spe <- readRDS(here("processed-data", "Visium","07_clustering", "BANKSY", "spe_banksy_lambda_0.8_harmony.rds"))
colnames(colData(spe))

mat <- getExplanatoryPCs(spe, dimred="PCA", variables = c("sample_id","slide_id","sum_umi","capture_area"), n_dimred=50)
mat

pdf(file.path(plot_dir, "PCA_variance_explained.pdf"), width = 8, height = 6)
plotExplanatoryPCs(mat)
dev.off()