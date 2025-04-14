library(SpatialExperiment)
library(here)
library(spatialLIBD)
library(scater)
library(patchwork)
library(ggpubr)
library(SpotSweeper)
library(escheR)


# Save directories
plot_dir = here("plots", "Visium","03_qc_metrics")
processed_dir = here("processed-data","Visium","03_qc_metrics")


# ===== Load brains Br2743, 6423, 6471, 6660, 8325, 9017, 9192, 9206, 9280, 9469 =====
spe.9280 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br9280.rds"))
spe.2743 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br2743.rds"))
spe.6423 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br6423.rds"))
spe.6471 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br6471.rds"))
spe.6660 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br6660.rds"))
spe.8325 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br8325.rds"))
spe.9017 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br9017.rds"))
spe.9192 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br9192.rds"))
spe.9206 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br9206.rds"))
spe.9469 <- readRDS(here("processed-data", "visium_stitching", "NacUtils", "spe_inputs", "spe_Br9469.rds"))

# ===== Second round of samples (#2-3) =====


# =========== Merging SCE objects ===========
common_genes <- Reduce(intersect, list(rownames(spe.2743), rownames(spe.6423),
                                         rownames(spe.6471), rownames(spe.6660), 
                                         rownames(spe.8325), rownames(spe.9017), 
                                         rownames(spe.9192), rownames(spe.9206), 
                                         rownames(spe.9280), rownames(spe.9469))
                                         )

spe.2743 <- spe.2743[common_genes,]
spe.6423 <- spe.6423[common_genes,]
spe.6471 <- spe.6471[common_genes,]
spe.6660 <- spe.6660[common_genes,]
spe.8325 <- spe.8325[common_genes,]
spe.9017 <- spe.9017[common_genes,]
spe.9192 <- spe.9192[common_genes,]
spe.9206 <- spe.9206[common_genes,]
spe.9280 <- spe.9280[common_genes,]
spe.9469 <- spe.9469[common_genes,]



# combine
spe <- Reduce(cbind, list(spe.2743, spe.6423, spe.6471, spe.6660, 
                           spe.8325, spe.9017, spe.9192, spe.9206, 
                           spe.9280, spe.9469)
                           )


# save combined object as rds
saveRDS(spe, here(processed_dir, "spe_stitched_combined_noQC.Rds"))
