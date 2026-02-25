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
ggspavis::plotCoords(spe, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()



# ====== Subset to only broad excitatory =====

spe.broad <- spe[, spe$pred_broad_celltype == "Excitatory"]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_excitatory.png"), width=15, height=15, units="in", res=300)
ggspavis::plotCoords(spe.broad, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()


# ====== subst to only fine names with "TSHZ1_" ======

spe.itc <- spe[, grepl("TSHZ1_", spe$pred_fine_celltype)]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_ITCs.png"), width=15, height=15, units="in", res=300)
ggspavis::plotCoords(spe.itc, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()


# ====== subset to broad inhibitory ======

spe.inh <- spe[, spe$pred_broad_celltype == "Inhibitory"]

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_inhibitory.png"), width=15, height=15, units="in", res=300)
ggspavis::plotSpots(spe.inh, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.01, sample_id="brnum") +
    scale_color_manual(values = celltype_colors)
dev.off()






# ======= ITC plotting for Warren Alpert grant =======

# subset to only excitatory and inhibitory
spe.neurons <- spe[, spe$pred_broad_celltype %in% c("Excitatory", "Inhibitory")]

library(scales)

tshz1_groups <- c("TSHZ1_CPNE4", "TSHZ1_PRKG1")
all_types <- names(celltype_colors)
other_types <- setdiff(all_types, tshz1_groups)

# subtle but clean base palette (greys, blue-grey, olive, lavender)
nice_base <- c(
  "#a26d6dff", # medium grey
  "#C7CAD0", # silver-blue
  "#9FAEC8", # slate lavender
  "#B0C4A3", # dusty sage
  "#CFA9A9", # rose sand
  "#8DA7A1", # sea-foam grey
  "#B7B6D9"  # soft violet
)

# recycle as many as needed
nice_base <- rep(nice_base, length.out = length(other_types))
names(nice_base) <- other_types

# Tshz1 highlight greens
highlight <- c(
  "TSHZ1_CPNE4"  = "#00B300", # rich emerald
  "TSHZ1_PRKG1"  = "#1a902eff"  # mint-neon (tasteful, not tacky)
)

celltype_colors2 <- c(nice_base, highlight)

png(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_neuronal_only.png"), width=15, height=15, units="in", res=300)
ggspavis::plotCoords(spe.neurons, annotate="pred_fine_celltype", in_tissue=NULL, point_size=.5, sample_id="brnum") +
    scale_color_manual(values = celltype_colors2)
dev.off()

# plot only last sample
spe.neurons.last <- spe.neurons[, spe.neurons$brnum == "Br9280"]

pdf(file=here(plot_dir, "SpotPlots_xenium_5um_fine_labels_neuronal_9280.pdf"),
    width=7.5, height=7.5) 
ggspavis::plotCoords(spe.neurons.last, annotate="pred_fine_celltype", in_tissue=NULL, point_size=.5, sample_id="brnum") +
    scale_color_manual(values = celltype_colors2)
dev.off()

