
library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","04_clustering", "banksy_16um")


spe_0.2 <- readRDS(here("processed-data", "VisiumHD","04_clustering", "banksy_16um","AmyHD_016_pilot_banksy_lambda_0.2.rds"))
spe_0.2 

spe_0.4 <- readRDS(here("processed-data", "VisiumHD","04_clustering", "banksy_16um","AmyHD_016_pilot_banksy_lambda_0.4.rds"))
spe_0.4 

spe_0.6 <- readRDS(here("processed-data", "VisiumHD","04_clustering", "banksy_16um","AmyHD_016_pilot_banksy_lambda_0.6.rds"))
spe_0.6

spe_0.8 <- readRDS(here("processed-data", "VisiumHD","04_clustering", "banksy_16um","AmyHD_016_pilot_banksy_lambda_0.8.rds"))
spe_0.8

spe_0.2$banksy_0.2 <- spe_0.2$cluster
spe_0.2$banksy_0.4 <- spe_0.4$cluster
spe_0.2$banksy_0.6 <- spe_0.6$cluster
spe_0.2$banksy_0.8 <- spe_0.8$cluster

spe <- spe_0.2
colnames(colData(spe))
#  [1] "barcode"              "in_tissue"            "array_row"           
#  [4] "array_col"            "sample_id"            "sum"                 
#  [7] "detected"             "subsets_mito_percent" "sum_discard"         
# [10] "sizeFactor"           "cluster"              "banksy_0.2"          
# [13] "banksy_0.4"           "banksy_0.6"           "banksy_0.8" 


# ======== spot plots ========

# lambda 0.2
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$banksy_0.2)))
png(here(plot_dir, "AmyHD_banksy_spotplots_016_lambda_0.2.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="banksy_0.2", point_size=.3) +
    scale_color_manual(values=pal) 
dev.off()

# lambda 0.4
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$banksy_0.4)))
png(here(plot_dir, "AmyHD_banksy_spotplots_016_lambda_0.4.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="banksy_0.4", point_size=.3) +
    scale_color_manual(values=pal)
dev.off()

# lambda 0.6
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$banksy_0.6)))
png(here(plot_dir, "AmyHD_banksy_spotplots_016_lambda_0.6.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="banksy_0.6", point_size=.3) +
    scale_color_manual(values=pal)
dev.off()

# lambda 0.8
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$banksy_0.8)))
png(here(plot_dir, "AmyHD_banksy_spotplots_016_lambda_0.8.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="banksy_0.8", point_size=.3) +
    scale_color_manual(values=pal)
dev.off()


# save final banksy spe
saveRDS(spe, here("processed-data", "VisiumHD","04_clustering", "banksy_16um","AmyHD_016_pilot_banksy_lambda_many.rds"))