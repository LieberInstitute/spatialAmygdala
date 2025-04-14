library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library(Banksy)
library(harmony)

# save directories
processed_dir <- here("processed-data", "Xenium", "04_clusterig")
plot_dir <- here("plots", "Xenium", "04_clustering")

load(here("processed-data","Xenium", "03_quality_control", "spe_normcounts.Rdata"))
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


# == Banksy ==sq
lambda <- 0.8
k_geom <- 36
npcs <- 50
aname <- "nucleus_normcounts"

set.seed(1000)
spe.gex <- Banksy::computeBanksy(spe.gex, assay_name = aname, k_geom = k_geom)


spe.gex <- Banksy::runBanksyPCA(spe.gex, lambda = lambda, npcs = npcs)

# harmony
spe.gex <- RunHarmony(spe.gex, "brnum")

# drop PCA, rename PCA_M0_lam0.4 to PCA
reducedDim(spe.gex, "PCA") <- reducedDim(spe.gex, "PCA_M0_lam0.8")
reducedDim(spe.gex, "PCA_M0_lam0.4") <- NULL
spe.gex <- RunHarmony(spe.gex, "capture_area")

# Cluster
spe.gex <- Banksy::clusterBanksy(spe.gex, lambda = lambda, npcs = npcs, resolution = 0.6)

saveRDS(spe.gex, here("processed-data", "Xenium","04_clustering", "Banksy", "Banksy_integrated_lambda_0.8_res0.6.rds"))