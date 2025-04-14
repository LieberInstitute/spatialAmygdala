
library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)


plot_dir <- here("plots","VisiumHD","04_clustering")


spe <- readRDS(here("processed-data", "VisiumHD","04_clustering", "banksy_8um","AmyHD_008_pilot_banksy.rds"))
spe 


colnames(colData(spe))


# ======== spot plots ========
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(spe$cluster)))
png(here(plot_dir, "AmyHD_banksy_spotplots.png"), width=10, height=10, units="in", res=300)
make_escheR(spe) |>
  add_ground(var="cluster", point_size=.3) +
    scale_color_manual(values=pal) 
dev.off()

