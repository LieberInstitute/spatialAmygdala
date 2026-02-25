suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("patchwork"))

load(here("processed-data", "Visium", "03_qc_metrics", "spe_stitched_local_outliers.Rdata"))

# drop local outliers
spe <- spe[, !colData(spe)$local_outliers]
dim(spe)

markers <- readRDS(here("processed-data","Visium","08_marker_genes", "markers_bs_manual_ITC_smoothed.rds"))
markers
# List of length 12
# names(12): BADL BL BLVM BM Ce HPC ITC LA Me PCo Vascular WM

# get top 100 markers per cluster
top_n <- 100
top_markers <- unlist(lapply(markers, function(x) {
    head(rownames(x), top_n)
}))

# get only nuique genes
top_markers <- unique(top_markers)
length(top_markers) 

# normalize data
spe <- computeLibraryFactors(spe)
spe <- spe[, sizeFactors(spe) > 0]
spe <- logNormCounts(spe)

rownames(spe) <- rowData(spe)$gene_name

set.seed(195)
message("running PCA - ", Sys.time())
spe <- scater::runPCA(spe, 
                      subset_row=top_markers,
                      ncomponents = 50,
                      exprs_values='logcounts',
                      scale = TRUE, name = "PCA")

# save
save(spe, file = here::here("processed-data", "Visium", "05_dim_reduction", "spe_stitched_pca_markers.Rdata"))
