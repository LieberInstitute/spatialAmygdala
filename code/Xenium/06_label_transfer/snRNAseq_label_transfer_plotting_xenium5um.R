library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")

# save directories
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
plot_dir <- here("plots", "Xenium", "06_label_transfer")


# load xenium data
spe <- readRDS(here("processed-data","Xenium","04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe


# subset to only Br9280
#spe <- spe[, spe$brnum == "Br9280"]

colnames(spe) <- make.unique(colnames(spe), sep = "-")
rownames(spatialCoords(spe)) <- colnames(spe)

# load label predictions
pred.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype.csv"))


# ======== Adding labels back in ========

# add labels
spe$pred_broad_celltype <- pred.broad$pruned.labels
spe$pred_fine_celltype <- pred.fine$pruned.labels


# ==== Plotting the label transfer results ====

mycolors <- scCustomize::DiscretePalette_scCustomize(length(unique(spe$pred_fine_celltype)), palette = "polychrome")
celltype_colors <- setNames(mycolors, unique(spe$pred_fine_celltype))

# UMAP
png(file=here(plot_dir, "UMAP_xenium_5um_fine_labels.png"), width=8, height=8, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", colour_by = "pred_fine_celltype", point_size=0.3) +
    scale_color_manual(values = celltype_colors)
dev.off()

png(file=here(plot_dir, "UMAP_xenium_5um_broad_labels.png"), width=8, height=8, units="in", res=300)
plotReducedDim(spe, dimred="UMAP", colour_by = "pred_broad_celltype", point_size=0.3)
dev.off()

# Spot Plots
png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_new.png"), width=15, height=15, units="in", res=300)
ggspavis::plotSpots(spe, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()



# ====== Subset to only broad excitatory =====

spe.broad <- spe[, spe$pred_broad_celltype == "Excitatory"]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_excitatory.png"), width=15, height=15, units="in", res=300)
ggspavis::plotSpots(spe.broad, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()


# ====== subst to only fine names with "TSHZ1_" ======

spe.itc <- spe[, grepl("TSHZ1_", spe$pred_fine_celltype)]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_ITCs.png"), width=15, height=15, units="in", res=300)
ggspavis::plotSpots(spe.itc, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()

# ====== subset to broad inhibitory ======

spe.inh <- spe[, spe$pred_broad_celltype == "Inhibitory"]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_inhibitory.png"), width=15, height=15, units="in", res=300)
ggspavis::plotSpots(spe.inh, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()