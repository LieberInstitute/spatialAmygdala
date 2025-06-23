suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("BiocSingular"))
suppressPackageStartupMessages(library("PCAtools"))
suppressPackageStartupMessages(library("patchwork"))

load(here("processed-data", "Visium", "03_qc_metrics", "spe_stitched_local_outliers.Rdata"))

SVGs.df <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
top_svg <- SVGs.df$gene_name[1:2000]

# normalize data
spe <- computeLibraryFactors(spe)
spe <- spe[, sizeFactors(spe) > 0]
spe <- logNormCounts(spe)

rownames(spe) <- rowData(spe)$gene_name

set.seed(195)
message("running PCA - ", Sys.time())
spe <- scater::runPCA(spe, 
                      subset_row=top_svg,
                      ncomponents = 50,
                      exprs_values='logcounts',
                      scale = TRUE, name = "PCA")

# save
save(spe, file = here::here("processed-data", "Visium", "05_dim_reduction", "spe_stitched_pca.Rdata"))

