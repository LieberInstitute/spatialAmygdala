library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library(Banksy)

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

lambda <- 0.8
k_geom <- 36
npcs <- 50
aname <- "nucleus_normcounts"
resolution <- 0.8

donors <- sort(unique(spe.gex$brnum))
donors

for (d in donors) {
  message("=== Donor: ", d, " ===")
  spe.gex.subset <- spe.gex[, spe.gex$brnum == d]

  # compute / PCA / cluster
  spe.gex.subset <- Banksy::computeBanksy(spe.gex.subset, assay_name = aname, k_geom = k_geom)

  set.seed(1000)
  spe.gex.subset <- Banksy::runBanksyPCA(spe.gex.subset, lambda = lambda, npcs = npcs)

  set.seed(1000)
  spe.gex.subset <- Banksy::clusterBanksy(spe.gex.subset, lambda = lambda, npcs = npcs, resolution = resolution)

  # save RDS (using your existing 04_clustering path for outputs)
  rds_fn <- here("processed-data", "Xenium","04_clustering", "Banksy",
                 sprintf("%s_Banksy_lambda_%.1f_subset.rds", d, lambda))
  saveRDS(spe.gex.subset, rds_fn)

}