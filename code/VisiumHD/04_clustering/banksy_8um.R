library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","04_clustering", "banksy_8um")


spe <- readRDS(here("processed-data", "VisiumHD","02_build_spe", "AmyHD_008_pilot.rds"))
spe
# class: SpatialExperiment 
# dim: 18085 605189 
# metadata(0):
# assays(2): '' counts
# rownames(18085): ENSG00000187634 ENSG00000188976 ... ENSG00000198695
#   ENSG00000198727
# rowData names(1): symbol
# colnames(605189): s_008um_00301_00321-1 s_008um_00602_00290-1 ...
#   s_008um_00353_00477-1 s_008um_00595_00611-1
# colData names(5): barcode in_tissue array_row array_col sample_id


# ===== Add QC Metrics =====

# drop out of tissue
#spe <- spe[, spe$in_tissue]

# get mito genes
rownames(spe) <- rowData(spe)$symbol
is.mito <- grepl("^MT-", rownames(spe))

# subset to SPE to is mito
spe.mito <- spe[rownames(spe) %in% is.mito,]

# get qc metrics
df <- scuttle::perCellQCMetrics(spe, subsets=list(Mito=is.mito))

# add to colData
spe$sum <- df$sum
spe$detected <- df$detected
spe$subsets_mito_percent <- df$subsets_Mito_percent
spe$subsets_mito_sum <- df$subsets_Mito_sum

spe$sum_discard <- isOutlier(spe$sum, nmads=3, type="lower", log=TRUE)
spe$detected_discard <- isOutlier(spe$detected, nmads=3, type="lower", log=TRUE)
spe$subsets_mito_percent_discard <- isOutlier(spe$subsets_mito_percent, nmads=3, type="higher")

# ======= SpotSweeper =========
# spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
# spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
# spe <- localOutliers(spe, metric = "subsets_mito_percent", direction = "higher", log = TRUE)

# ========== Banksy clustering ===========
library(Banksy)

# == Normalization ==
# drop MT genes
spe <- spe[!grepl("^MT-", rownames(spe)),]

# calculate library size factors
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))

spe <- spe[, sizeFactors(spe) > 0]
dim(spe)

spe <- logNormCounts(spe)

# == Feature Selection ==
# get HGs
dec <- modelGeneVar(spe)
chosen <- getTopHVGs(dec, n=500)

spe.hvg <- spe[chosen,]


# == Banksy ==sq
lambda <- 0.8
k_geom <- 30
npcs <- 30
aname <- "logcounts"
spe.hvg <- Banksy::computeBanksy(spe.hvg, assay_name = aname, k_geom = k_geom)

set.seed(1000)
spe.hvg <- Banksy::runBanksyPCA(spe.hvg, lambda = lambda, npcs = npcs)

set.seed(1000)
spe.hvg <- Banksy::clusterBanksy(spe.hvg, lambda = lambda, npcs = npcs, resolution = 0.8)

saveRDS(spe.hvg, here("processed-data", "VisiumHD","04_clustering", "banksy_8um","AmyHD_008_pilot_banksy_hvg.rds"))

# add cluster names back to original SPE
spe$cluster <- spe.hvg$clust_M0_lam0.8_k50_res0.8

saveRDS(spe, here("processed-data", "VisiumHD","04_clustering", "banksy_8um","AmyHD_008_pilot_banksy.rds"))
