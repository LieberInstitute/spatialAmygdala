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
})

# Load your data
spe <-readRDS(here("processed-data","Visium", "06_batch_correction", "spe_harmony_markers.rds"))

# Load cluster assignments for k=16
bs_folders <- list.files(here::here("processed-data","Visium", "07_clustering", "BayesSpace","MarkerGenes","cluster_csv_markers"), full.names = TRUE)

# Find the folder for k=16
k16_folder <- bs_folders[grepl("BayesSpace_16$", bs_folders)]

if (length(k16_folder) == 0) {
    stop("No BayesSpace_16 folder found")
}

# Load the CSV file from k=16 folder
bs_csv <- list.files(k16_folder, pattern = "csv$", full.names = TRUE)[1]
bs_df <- read.csv(bs_csv)

# Add clusters to spe colData
colData(spe)$BS_k16 <- factor(bs_df$cluster)

# Verify clusters loaded correctly
cat("Number of clusters:", length(unique(colData(spe)$BS_k16)), "\n")
cat("Spots per cluster:\n")
print(table(colData(spe)$BS_k16))

# Create pseudobulk using dreamlet::aggregateToPseudoBulk
pb <- aggregateToPseudoBulk(spe,
                            assay = "counts",  
                            cluster_id = "BS_k16",
                            sample_id = "sample_id") 

# Now pb contains pseudobulk profiles: one column per cluster per sample
# The columns are named like "cluster.sample"


hcl <- buildClusterTreeFromPB(pb)


pdf(here("plots", "Visium", "07_clustering", "BayesSpace", "MarkerGenes", "BS_k16_dendrogram.pdf"), 
    width = 12, height = 8)

plot(hcl, main = "Hierarchical Clustering of k=16 Clusters\n(Based on Pseudobulk Gene Expression)",
     xlab = "Clusters", ylab = "Height")

dev.off()




library(ggplot2)

# Subset to sample Br8316
spe_sample <- spe[, spe$sample_id == "Br8325"]

# Get spatial coordinates and cluster info
plot_df <- data.frame(
    x = spatialCoords(spe_sample)[, "pxl_col_in_fullres"],
    y = spatialCoords(spe_sample)[, "pxl_row_in_fullres"],
    cluster = spe_sample$BS_k16
)

# Simple faceted plot
pdf(here("plots", "Visium", "07_clustering","BayesSpace", "MarkerGenes", "BS_k16_faceted_Br8325.pdf"), 
    width = 20, height = 20)

ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(aes(color = cluster), size = 1.5) +
    facet_wrap(~cluster, ncol = 5) +
    scale_y_reverse() +  # Flip y-axis for proper spatial orientation
    coord_fixed() +  # Keep aspect ratio
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Clusters faceted - Sample Br8316")

dev.off()

# Or if you want to highlight each cluster individually in its facet
plot_df$is_current_cluster <- plot_df$cluster

pdf(here("plots", "Visium", "07_clustering", "BayesSpace", "MarkerGenes",  "BS_k16_faceted_highlighted_Br8325.pdf"), 
    width = 20, height = 20)

ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(data = transform(plot_df, cluster = NULL), color = "grey90", size = 1) +
    geom_point(aes(color = is_current_cluster), size = 1.5) +
    facet_wrap(~cluster, ncol = 5) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Clusters highlighted - Sample Br8316")

dev.off()


# ========= Re labeling ========

# make a new cluster ID with these from BS k16. For multiples use WM.1, WM.2, etc
colData(spe)$BS_k16_relabel <- factor(dplyr::case_when(
    colData(spe)$BS_k16 == 8 ~ "Meninges",
    colData(spe)$BS_k16 == 10 ~ "WM.2",
    colData(spe)$BS_k16 == 15 ~ "BLVM.2",
    colData(spe)$BS_k16 == 14 ~ "WM.1",
    colData(spe)$BS_k16 == 7 ~ "WM.1",
    colData(spe)$BS_k16 == 13 ~ "WM.1",
    colData(spe)$BS_k16 == 5 ~ "CeA",
    colData(spe)$BS_k16 == 2 ~ "MeA",
    colData(spe)$BS_k16 == 11 ~ "BLVM.1",
    colData(spe)$BS_k16 == 9 ~ "BLVM.2",
    colData(spe)$BS_k16 == 3 ~ "LA",
    colData(spe)$BS_k16 == 12 ~ "BA",
    colData(spe)$BS_k16 == 16 ~ "HPC",
    colData(spe)$BS_k16 == 1 ~ "CoA",
    colData(spe)$BS_k16 == 6 ~ "aBA",
    colData(spe)$BS_k16 == 4 ~ "Other",
))


# re-plot all same samples with new labels

library("escheR")
pdf(file = here::here("plots", "Visium", "07_clustering", "BayesSpace", "MarkerGenes", "markers", "BS_k16_relabeled.pdf"), width = 5 * length(unique(spe$sample_id)), height = 5)

# Create one plot per sample, then combine in a row
plots <- lapply(unique(spe$sample_id), function(sample_id) {
    spe.subset <- spe[, spe$sample_id == sample_id]
    make_escheR(spe.subset) |>
        add_fill(var = "BS_k16_relabel", point_size = 1) +
        scale_fill_manual(values = pal) +
        ggtitle(sample_id) 
})

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)

dev.off()

# save
saveRDS(spe, here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_relabel.rds"))

