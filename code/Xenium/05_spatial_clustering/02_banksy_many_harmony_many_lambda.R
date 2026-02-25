library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("Banksy")
library("harmony")
library("data.table")

# save directories
processed_dir <- here("processed-data", "Xenium", "05_spatial_clustering")
plot_dir <- here("plots", "Xenium", "05_spatial_clustering")

spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "sce_combined_harmonized_singlecell.rds"))
spe
# class: SpatialExperiment 
# dim: 541 1018069 
# metadata(0):
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018069): aaaadgkh-1 aaaadlfe-1 ... oihobmbk-1 oihoeehh-1
# colData names(19): cell_id transcript_counts ... cell_area.sf
#   nucleus_area.sf
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# ======= Feature Selection =======

# subset to only target genes probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe.gex <- spe[gene_expression_idx,]
dim(spe.gex)
# [1]     366 1018069

# ======== Banksy clustering ========
# per: https://prabhakarlab.github.io/Banksy/articles/batch-correction.html


# aajifchf-1   130.469330  9324.461914
# aajilhii-1   128.405975  9308.610352
# aajjaoeb-1   145.384537  9299.897461
# aajjghhf-1    97.161903  9324.034180
# aajjnefj-1    92.433006  9351.603516
# aajjoggg-1   165.425476  9412.909180
# aajjpabb-1   157.702469  9439.750977
# aajkccda-1   160.587875  9420.091797

locs <- spatialCoords(spe.gex)
locs <- cbind(locs, sample_id = factor(spe.gex$brnum))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe.gex) <- locs

# print min/max of coords per sample
locs_dt[, .(min_sdimx = min(sdimx), max_sdimx = max(sdimx)), by = group]
#    group min_sdimx max_sdimx
#    <num>     <num>     <num>
# 1:     3  51885.92  63411.03
# 2:     1  17295.31  28823.90
# 3:     2  34590.61  46120.82
# 4:     4  69181.23  80701.34


# == Banksy ==sq
lambda <- 0.8
k_geom <- 36
npcs <- 50
aname <- "nucleus_normcounts"

set.seed(1000)
spe.gex <- Banksy::computeBanksy(spe.gex, assay_name = aname, k_geom = k_geom)


spe.gex <- Banksy::runBanksyPCA(spe.gex, lambda = lambda, npcs = npcs)
spe.gex


# drop PCA, rename PCA_M0_lam0.8 to PCA
reducedDim(spe.gex, "PCA") <- reducedDim(spe.gex, "PCA_M0_lam0.8")
reducedDim(spe.gex, "PCA_M0_lam0.8") <- NULL
spe.gex <- RunHarmony(spe.gex, "brnum", reduction.save="HARMONY_M0_lam0.8")

# Cluster
spe.gex <- Banksy::clusterBanksy(spe.gex, lambda = lambda, npcs = npcs, resolution = 0.8, , dimred = "HARMONY_M0_lam0.8")

# clusters to original spe
spe$clust_Banksy_k36_res0.6_lam0.8 <- spe.gex$clust_Banksy_k36_res0.6_lam0.8


saveRDS(spe, here("processed-data", "Xenium","04_clustering", "Banksy", "Banksy_integrated_lambda_0.8_res0.6.rds"))
saveRDS(spe.gex, here("processed-data", "Xenium","04_clustering", "Banksy", "Banksy_integrated_lambda_0.8_res0.6_gex.rds"))