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


# Save directories
plot_dir = here("plots", "07_clustering", "BANKSY")
processed_dir = here("processed-data","07_clustering")

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)
 
# class: SpatialExperiment 
# dim: 25801 308044 
# metadata(0):
# assays(2): counts logcounts
# rownames(25801): MIR1302-2HG AL627309.1 ... AC007325.4 AC007325.2
# rowData names(6): source type ... gene_name Symbol.uniq
# colnames(308044): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
#   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
# colData names(36): sample_id in_tissue ... local_outliers sizeFactor
# reducedDimNames(4): 10x_pca 10x_tsne 10x_umap PCA
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor



colnames(colData(spe))
#unique(spe$exclude_overlapping)
#[1] FALSE  TRUE    NA

# renormalization
spe <- logNormCounts(spe)



# ====== Stagger spatial coordinates ========
# per: https://prabhakarlab.github.io/Banksy/articles/batch-correction.html

locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = factor(spe$sample_id))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs

# ===== Load SVGs ======
#load nnSVG results
SVGs.df <- read.csv(here("processed-data", "Visium", "04_feature_selection", "nnSVG_summary.csv"))
genes <- SVGs.df$gene_name[1:2000]


# subset to hvgs
spe <- spe[genes, ]


# == Banksy ==sq
lambda <- 0.4
k_geom <- 18
npcs <- 20
aname <- "logcounts"

spe <- Banksy::computeBanksy(spe, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs, group = "capture_area")


spe <- RunHarmony(spe,  
                group.by.vars =c("sample_id","slide_id"),  
                reduction = "PCA_M0_lam0.4",
                reduction.save="HARMONY_M0_lam0.4"
                )

set.seed(1000)
spe$clust_HARMONY_k50_res0.8 <- NULL
spe <- Banksy::clusterBanksy(spe, lambda = lambda, npcs = npcs, resolution = 0.4, dimred = "HARMONY_M0_lam0.4")

# drop duplicate columns
colData(spe) <- colData(spe)[, !duplicated(colnames(colData(spe)))]
colnames(colData(spe))

saveRDS(spe, here("processed-data", "Visium","07_clustering", "BANKSY", "spe_banksy_lambda_0.8_harmony.rds"))












