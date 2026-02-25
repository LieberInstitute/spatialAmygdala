suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("dendextend")
    library("pheatmap")
    library("dreamlet")
    library("SingleCellExperiment")  # needed for dreamlet
    library("escheR")
})

# Load your data
spe <-readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_relabel.rds"))
spe.bs <- spe


load(here("processed-data","Visium", "98_NMF", "spe_NMF_Yu.rda"))
spe.nmf <- spe
spe.nmf

# drop low quality samples
spe.nmf <- spe.nmf[, !colData(spe.nmf)$sample_id %in% c("Br9469", "Br9017", "Br9206")]

# subset to same spots as in spe.bs
spe.nmf <- spe.nmf[, colnames(spe.bs)]
all(colnames(spe.nmf) == colnames(spe.bs))
# [1] TRUE

# in spe.bs, set spots > 0 in NMF 44 + 66 as ITC in spe.bs clusters
spe.bs$ITC <- "Non-ITC"
spe.bs$ITC[ (reducedDims(spe.nmf)$NMF_proj[,44] > 1)] <- "ITC"
table(spe.bs$ITC)

# plot ITC vs Non-ITC, make ITC red and Non-ITC grey
pal_itc <- c("Non-ITC" = "lightgrey", "ITC" = "red")

pdf(here("plots", "Visium", "07_clustering","BayesSpace", "MarkerGenes","ITCs", "Spatial_ITC_vs_NonITC.pdf"), width = 5 * length(unique(spe.bs$sample_id)), height = 5)
plots <- lapply(unique(spe.bs$sample_id), function(sample_id) {
    spe.subset <- spe.bs[, spe.bs$sample_id == sample_id]
    make_escheR(spe.subset) |>
        add_fill(var = "ITC", point_size = 1) +
        scale_fill_manual(values = pal_itc) +
        ggtitle(sample_id) 
})

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)

dev.off()

# okay > 1 isn't great. Let's try to GMM of k2, where the high value cluster is ITC
library(mclust)
gmm_44 <- Mclust(reducedDims(spe.nmf)$NMF_proj[,44], G=2)
gmm_65 <- Mclust(reducedDims(spe.nmf)$NMF_proj[,65], G=2)


# okay now just do 44 + 65
spe.bs$ITC_gmm2 <- "Non-ITC"
spe.bs$ITC_gmm2[ (gmm_44$classification == 2) | (gmm_65$classification == 2) ] <- "ITC"
table(spe.bs$ITC_gmm2)

# make ITC in BS_k16_relabel
spe.bs$BS_k16_relabel_ITC <- as.factor(spe.bs$BS_k16_relabel)
spe.bs$BS_k16_relabel_ITC <- factor(
  spe.bs$BS_k16_relabel,
  levels = c(levels(spe.bs$BS_k16_relabel), "ITC")
)
spe.bs$BS_k16_relabel_ITC[ as.factor(spe.bs$ITC_gmm2) == "ITC" ] <- "ITC"
table(spe.bs$BS_k16_relabel_ITC)


WM1_pal <- c("#D3D3D3")
WM2_pal <- c("#808080ff")

green_pal <- colorRampPalette(c("#0fdb71", "#05a150"))(2) # Shades of green
Ce_pal <- green_pal[1]
Me_pal <- green_pal[2]
HPC_pal <- c("#FFA500")
Other_pal <- c("#ff00eaff")
ITC_pal <- c("red")

blue_pal <- colorRampPalette(c("#03dffc", "#038cfc"))(6) # Shades of blue
CoA_pal <- blue_pal[1]
LA_pal <- blue_pal[2]
aBA_pal <- blue_pal[3]
BA_pal <- blue_pal[4]
BLVM1_pal <- blue_pal[5]
BLVM2_pal <- blue_pal[6]
Meninges_pal <- c("#000000ff")

pal <- c(WM1_pal, WM2_pal, Ce_pal, Me_pal, HPC_pal, Other_pal, CoA_pal, LA_pal, aBA_pal, BA_pal, BLVM1_pal, BLVM2_pal, Meninges_pal, ITC_pal)
names(pal) <- c("WM.1", "WM.2", "CeA", "MeA", "HPC", "Other", "CoA", "LA", "aBA", "BA", "BLVM.1", "BLVM.2", "Meninges", "ITC")


# replot with custom color sclae +Red for itc
pdf(here("plots", "Visium",  "07_clustering","BayesSpace", "MarkerGenes","ITCs", "Spatial_Clusters_Grid_BS_manual_ITC_custom_palette.pdf"), width = 5 * length(unique(spe.bs$sample_id)), height = 5)
plots <- lapply(unique(spe.bs$sample_id), function(sample_id) {
    spe.subset <- spe.bs[, spe.bs$sample_id == sample_id]
    make_escheR(spe.subset) |>
        add_fill(var = "BS_k16_relabel_ITC", point_size = .8) +
        scale_fill_manual(values = pal) +
        ggtitle(sample_id) 
})

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)

dev.off()





library(FNN)

coords <- spatialCoords(spe.bs)
sample_ids <- spe.bs$sample_id
initial_labels <- spe.bs$ITC_gmm2  # "ITC" or "Non-ITC"

# Initialize vector for smoothed labels
smoothed_labels <- rep(NA, nrow(spe.bs))

# Loop through each sample
for (s in unique(sample_ids)) {
  idx <- which(sample_ids == s)
  coords_d <- coords[idx, , drop = FALSE]
  labels_d <- initial_labels[idx]
  binary_labels <- as.integer(labels_d == "ITC")
  
  # Check for enough spots for kNN
  if (length(idx) <= 6) {
    warning("Sample ", s, " has too few spots (<=6); skipping smoothing for this sample.")
    smoothed_labels[idx] <- labels_d  # fallback to original
    next
  }

  # Get 6-nearest neighbors
  nn <- get.knn(coords_d, k = 6)$nn.index
  
  # Count how many of each spot's neighbors are ITC
  itc_counts <- rowSums(matrix(binary_labels[nn], nrow = nrow(nn)))

  # If ≥ 4 neighbors are ITC, relabel as ITC
  smoothed_labels[idx] <- ifelse(itc_counts >= 3, "ITC", "Non-ITC")
}

# Save smoothed labels to object
spe.bs$ITC_gmm2_smooth <- smoothed_labels

# Update the BS_k16_relabel_ITC label accordingly
spe.bs$BS_k16_relabel_ITC_smooth <- as.character(spe.bs$BS_k16_relabel)
spe.bs$BS_k16_relabel_ITC_smooth[smoothed_labels == "ITC"] <- "ITC"
spe.bs$BS_k16_relabel_ITC_smooth <- factor(spe.bs$BS_k16_relabel_ITC_smooth, levels = c(levels(spe.bs$BS_k16_relabel), "ITC"))

# Sanity check
table(spe.bs$ITC_gmm2, spe.bs$ITC_gmm2_smooth)
table(spe.bs$BS_k16_relabel_ITC_smooth)



# plot
pal_itc <- c("Non-ITC" = "lightgrey", "ITC" = "red")
pdf(here("plots", "Visium",  "07_clustering","BayesSpace", "MarkerGenes","ITCs","Spatial_ITC_vs_NonITC_GMM_smoothed.pdf"), width = 5 * length(unique(spe.bs$sample_id)), height = 5)
plots <- lapply(unique(spe.bs$sample_id), function(sample_id) {
    spe.subset <- spe.bs[, spe.bs$sample_id == sample_id]
    make_escheR(spe.subset) |>
        add_fill(var = "BS_k16_relabel_ITC_smooth", point_size = 1) +
        scale_fill_manual(values = pal_itc) +
        ggtitle(sample_id) 
})

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)
dev.off()

# make ITC in BS_manual
spe.bs$BS_k16_relabel_ITC <- spe.bs$BS_k16_relabel
spe.bs$BS_k16_relabel_ITC <- factor(
  spe.bs$BS_k16_relabel,
  levels = c(levels(spe.bs$BS_k16_relabel), "ITC")
)
spe.bs$BS_k16_relabel_ITC[ spe.bs$ITC_gmm2_smooth == "ITC" ] <- "ITC"
table(spe.bs$BS_k16_relabel_ITC)

# replot with custom color sclae +Red for itc
pdf(here("plots", "Visium", "07_clustering","BayesSpace", "MarkerGenes","ITCs","Spatial_Clusters_Grid_BS_manual_ITC_smoothed_custom_palette.pdf"), width = 5 * length(unique(spe.bs$sample_id)), height = 5)
plots <- lapply(unique(spe.bs$sample_id), function(sample_id) {
    spe.subset <- spe.bs[, spe.bs$sample_id == sample_id]
    make_escheR(spe.subset) |>
        add_fill(var = "BS_k16_relabel_ITC", point_size = .8) +
        scale_fill_manual(values = pal) +
        ggtitle(sample_id) 
})

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)
dev.off()




# save
saveRDS(spe.bs, file=here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_relabel_ITC_smoothed.rds"))
