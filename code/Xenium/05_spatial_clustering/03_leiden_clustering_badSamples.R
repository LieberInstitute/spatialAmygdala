library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("bluster")


# load xenium data
spe <- readRDS(here("processed-data","Xenium", "04_dim_reduction", "spe_proseg_5um_badSamples_harmonized_singlecell.rds"))
spe


# ======= Clustering ========

clust.louvain <- clusterCells(spe, use.dimred="HARMONY", 
    BLUSPARAM=NNGraphParam(k = 25, cluster.fun="leiden"))

spe$leiden_k25 <- clust.louvain

# save the clustering results
output_file <- here("processed-data", "Xenium", "05_spatial_clustering", "spe_proseg_5um_badSamples_harmonized_singlecell_leiden_k25.rds")
saveRDS(spe, file=output_file)
message("Clustering results saved to: ", output_file)