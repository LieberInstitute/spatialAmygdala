library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("bluster")



spe.good <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering", "spe_proseg_5um_badSamples_harmonized_singlecell_leiden_k25.rds"))
spe.bad <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering", "spe_proseg_5um_goodSamples_harmonized_singlecell_leiden_k25.rds"))

# join spe
spe <- cbind(spe.good, spe.bad)

# ====== load transferred labels ======

# predictions on xenium segmentations
pred.x.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "pred_broad_celltype.csv"))
pred.x.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "pred_fine_celltype.csv"))

# predictions on proseg segmentations
pred.p.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "proseg_pred_broad_celltype.csv"))
pred.p.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "proseg_pred_fine_celltype.csv"))

# ====== add labels to spe ======
# split SPE again
spe.good <- spe[, spe$brnum %in% c("Br9280")]
spe.bad <- spe[, spe$brnum %in% c("Br9017", "Br9206")]

spe.good$labels_proseg_broad <- pred.p.broad$pruned.labels
spe.good$labels_proseg_fine <- pred.p.fine$pruned.labels

# spe$labels_xenium_broad <- pred.x.broad$pruned.labels
# spe$labels_xenium_fine <- pred.x.fine$pruned.labels




# ======= Plotting =======
library(patchwork)

pdf(file = here("plots", "Xenium", "05_spatial_clustering", "UMAP_goodSamples_labels.pdf"),
    width = 40, height = 10)

p1 <- plotReducedDim(spe.good, dimred = "UMAP_HARMONY", colour_by = "brnum") +
    ggtitle("Proseg UMAP: Good Samples")

p2 <- plotReducedDim(spe.good, dimred ="UMAP_HARMONY",  colour_by = "leiden_k25") +
    ggtitle("Proseg UMAP: Leiden Clusters (k=25)")

p3 <- plotReducedDim(spe.good, dimred = "UMAP_HARMONY", colour_by = "labels_proseg_broad") +
    ggtitle("Proseg Broad Cell Types")

p4 <- plotReducedDim(spe.good, dimred = "UMAP_HARMONY", colour_by = "labels_proseg_fine") +
    ggtitle("Proseg Fine Cell Types")

cowplot::plot_grid(p1,p2,p3,p4, ncol = 4)
dev.off()



