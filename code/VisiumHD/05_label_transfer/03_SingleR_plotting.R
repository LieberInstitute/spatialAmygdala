library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")

# save directories
processed_dir <- here("processed-data", "VisiumHD", "06_label_transfer", "SingleR")
plot_dir <- here("plots", "VisiumHD", "05_label_transfer", "SingleR")


# load VisiumHD data
spe <- readRDS(here("processed-data/VisiumHD/02_build_spe/Br9280_CeA_Spatial.Polygons_spe.rds"))
spe

# colnames(spe) <- make.unique(colnames(spe), sep = "-")
# rownames(spatialCoords(spe)) <- colnames(spe)

# load label predictions
pred.broad <- read.csv(here("processed-data", "VisiumHD", "05_label_transfer", "SingleR", "HDsegmentations_pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "VisiumHD", "05_label_transfer", "SingleR", "HDsegmentations_pred_fine_celltype.csv"))


# ======== Adding labels back in ========

# add labels
spe$pred_broad_celltype <- pred.broad$pruned.labels
spe$pred_fine_celltype <- pred.fine$pruned.labels


# ==== Plotting the label transfer results ====

mycolors <- scCustomize::DiscretePalette_scCustomize(length(unique(spe$pred_fine_celltype)), palette = "polychrome")
celltype_colors <- setNames(mycolors, unique(spe$pred_fine_celltype))

# Spot Plots
pdf(file=here(plot_dir, "SpotPlots_VisiumHD_fine_labels.pdf"), width=10, height=8)
ggspavis::plotSpots(spe, annotate="pred_fine_celltype", in_tissue=NULL, point_size=.5) +
    scale_color_manual(values = celltype_colors)
dev.off()

pdf(file=here(plot_dir, "SpotPlots_VisiumHD_broad_labels.pdf"), width=10, height=8)
ggspavis::plotSpots(spe, annotate="pred_broad_celltype", in_tissue=NULL, point_size=.5)
dev.off()




# ====== Subset to only broad excitatory =====

spe.broad <- spe[, spe$pred_broad_celltype == "Excitatory"]

pdf(file=here(plot_dir, "SpotPlots_VisiumHDfine_labels_excitatory.pdf"), width=10, height=8)
ggspavis::plotSpots(spe.broad, annotate="pred_fine_celltype", in_tissue=NULL, point_size=.5) +
    scale_color_manual(values = celltype_colors)
dev.off()


# ====== subst to only fine names with "TSHZ1_" ======

spe.itc <- spe[, grepl("TSHZ1_", spe$pred_fine_celltype)]

pdf(file=here(plot_dir, "SpotPlots_VisiumHD_fine_labels_ITCs.pdf"), width=10, height=8)
ggspavis::plotSpots(spe.itc, annotate="pred_fine_celltype", in_tissue=NULL, point_size=.5) +
    scale_color_manual(values = celltype_colors)
dev.off()







# Drop NA labels once
spe <- spe[, !is.na(spe$pred_fine_celltype)]

# reset levels after dropping NAs
spe$pred_fine_celltype <- factor(spe$pred_fine_celltype)

fine_levels <- sort(unique(spe$pred_fine_celltype))
for (ct in fine_levels) {
  spe.subset <- spe[, spe$pred_fine_celltype == ct]

  p <- ggspavis::plotSpots(
    spe.subset,
    annotate = "pred_fine_celltype",
    in_tissue = NULL,
    point_size = 0.5
  ) +
    scale_color_manual(values = celltype_colors[ct, drop = FALSE]) +
    ggtitle(ct) +
    guides(color = "none")

  fn <- file.path(plot_dir, paste0("SpotPlot_", ct, ".png"))
  ggsave(fn, p, width = 10, height = 8, dpi = 300)
}