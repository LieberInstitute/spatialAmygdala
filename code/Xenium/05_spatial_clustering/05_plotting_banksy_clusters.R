library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("bluster")



spe.9280 <- readRDS(here("processed-data", "Xenium", "04_clustering","Banksy", "Br9280_Banksy_lambda_0.8_subset.rds"))
spe.9206 <- readRDS(here("processed-data", "Xenium", "04_clustering","Banksy",  "Br9206_Banksy_lambda_0.8_subset.rds"))
spe.9192 <- readRDS(here("processed-data", "Xenium", "04_clustering", "Banksy", "Br9192_Banksy_lambda_0.8_subset.rds"))
spe.9017 <- readRDS(here("processed-data", "Xenium", "04_clustering","Banksy", "Br9017_Banksy_lambda_0.8_subset.rds"))

colnames(colData(spe.9280))

# ======= Plotting =======
library(patchwork)

pdf(file = here("plots", "Xenium", "05_spatial_clustering", "Banksy_clustering_all_samples.pdf"),
    width = 40, height = 10)

p1 <- ggspavis::plotSpots(spe.9280, annotate="clust_M0_lam0.8_k50_res0.8", in_tissue=NULL, sample_id="brnum")

p2 <- ggspavis::plotSpots(spe.9206, annotate="clust_M0_lam0.8_k50_res0.8", in_tissue=NULL, sample_id="brnum")

p3 <- ggspavis::plotSpots(spe.9192, annotate="clust_M0_lam0.8_k50_res0.8", in_tissue=NULL, sample_id="brnum")

p4 <- ggspavis::plotSpots(spe.9017, annotate="clust_M0_lam0.8_k50_res0.8", in_tissue=NULL, sample_id="brnum")

cowplot::plot_grid(p1,p2,p3,p4, ncol = 4)
dev.off()



