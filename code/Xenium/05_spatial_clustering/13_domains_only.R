suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("RColorBrewer")
    library("ggplot2")
    library("scran")
    library("scater")
    library(patchwork)
})


output_dir <- here("plots", "Xenium", "05_spatial_clustering")

# load
spe <- readRDS(file.path("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
spe



# === Subset to Cluster 16 using grep - Mea + ITCs ===
clusters_of_interest <- c("16")
spe_subset <- spe[, grepl(paste(clusters_of_interest, collapse = "|"), spe$Banksy_res2.0_collapsed_v6)]
spe_subset
# class: SpatialExperiment 
# dim: 366 53027 
# metadata(9): polygons technology ... technology BANKSY_params
# assays(4): counts nucleus_normcounts cell_normcounts H0
# rownames(366): ABCC9 ADAMTS12 ... ZIC2 ZNF536
# rowData names(3): ID Symbol Type
# colnames(53027): aadahfem-1 aadanoag-1 ... ohhlhklp-1 ohhlibfb-1
# colData names(32): cell_id transcript_counts ...
#   clust_HARMONY_M0_lam0.8_k50_res2 Banksy_res2.0_collapsed_v6
# reducedDimNames(4): PCA HARMONY UMAP HARMONY_M0_lam0.8
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# Computer new PCA and UMAP on subsetted data
spe_subset <- runPCA(spe_subset, exprs_values = "nucleus_normcounts", ncomponents = 10)
spe_subset <- runUMAP(spe_subset, dimred = "PCA")

# plot PCA and UMAP colored by brnum
pdf(file.path(output_dir, "PCA_UMAP_Cluster16.pdf"), width = 12, height = 6)
p1 <- plotReducedDim(spe_subset, "HARMONY_M0_lam0.8", colour_by = "brnum") + ggtitle("PCA - Cluster 16")
p2 <- plotReducedDim(spe_subset, "HARMONY_M0_lam0.8", colour_by = "brnum") + ggtitle("UMAP - Cluster 16")
p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()

pdf(file.path(output_dir, "UMAP_Cluster16_TSHZ1_SLC17A6.pdf"), width = 12, height = 6)
p1 <- plotReducedDim(spe_subset, "UMAP", colour_by = "TSHZ1", by_exprs_values = "nucleus_normcounts") + ggtitle("TSHZ1")
p2 <- plotReducedDim(spe_subset, "UMAP", colour_by = "SLC17A6", by_exprs_values = "nucleus_normcounts") + ggtitle("SLC17A6")
p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()