suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("harmony"))

processed_dir <- here("processed-data", "Visium", "07_batch_correction")

load(here("processed-data","Visium", "06_dim_reduction", "spe_stitched_pca.Rdata"))
spe
# Loading required package: SpatialExperiment
# class: SpatialExperiment 
# dim: 36601 227302 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
#   ENSG00000277196
# rowData names(1): symbol
# colnames(227302): AAACAAGTATCTCCCA-1_V13Y24-346_A1
#   AAACAATCTACTAGCA-1_V13Y24-346_A1 ... TTGTTTCATTAGTCTA-1_V13B23-407_C1
#   TTGTTTCCATACAACT-1_V13B23-407_C1
# colData names(37): in_tissue array_row ... subsets_mito_percent_z
#   sizeFactor
# reducedDimNames(1): PCA
# mainExpName: NULL
# altExpNames(0):

reducedDimNames(spe)
#  [1] "PCA"

#run harmony
spe <- RunHarmony(spe, "brnum")

#remove PCA for future use
reducedDim(spe, "PCA") <- NULL

#run UMAP
spe <- runUMAP(spe, dimred = "PCA", name = "UMAP-PCA", min_dist=0.3)
spe <- runUMAP(spe, dimred = "HARMONY", name = "UMAP-HARMONY", , min_dist=0.3)

#explore UMAP results
pdf(file = here::here("plots","Visium", "07_batch_correction", "UMAP_uncorrected.pdf"))
plotReducedDim(spe, dimred="UMAP-PCA", colour_by="brnum")
dev.off()

pdf(file = here::here("plots","Visium", "07_batch_correction", "UMAP_harmony.pdf"))
plotReducedDim(spe, dimred="UMAP-HARMONY", colour_by="brnum")
dev.off()

save(spe, file = here::here("processed-data","Visium", "07_batch_correction", "spe_harmony.Rdata"))