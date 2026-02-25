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
load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))

# Drop low quality samples
spe <- spe[, !colData(spe)$sample_id %in% c("Br9469", "Br9017", "Br9206")]

# Load cluster assignments for k=25
bs_folders <- list.files(here::here("processed-data","Visium", "07_clustering", "BayesSpace","SVGs","cluster_csv_new"), full.names = TRUE)

# Find the folder for k=25
k25_folder <- bs_folders[grepl("BayesSpace_25$", bs_folders)]

if (length(k25_folder) == 0) {
    stop("No BayesSpace_25 folder found")
}

# Load the CSV file from k=25 folder
bs_csv <- list.files(k25_folder, pattern = "csv$", full.names = TRUE)[1]
bs_df <- read.csv(bs_csv)

# Add clusters to spe colData
colData(spe)$BS_k25 <- factor(bs_df$cluster)

# Verify clusters loaded correctly
cat("Number of clusters:", length(unique(colData(spe)$BS_k25)), "\n")
cat("Spots per cluster:\n")
print(table(colData(spe)$BS_k25))

# Create pseudobulk using dreamlet::aggregateToPseudoBulk
pb <- aggregateToPseudoBulk(spe,
                            assay = "counts",  
                            cluster_id = "BS_k25",
                            sample_id = "sample_id") 

# Now pb contains pseudobulk profiles: one column per cluster per sample
# The columns are named like "cluster.sample"


hcl <- buildClusterTreeFromPB(pb)


pdf(here("plots", "Visium", "07_clustering", "BayesSpace", "SVGs", "BS_k25_dendrogram.pdf"), 
    width = 12, height = 8)

plot(hcl, main = "Hierarchical Clustering of k=25 Clusters\n(Based on Pseudobulk Gene Expression)",
     xlab = "Clusters", ylab = "Height")

dev.off()




library(ggplot2)

# Subset to sample Br8325
spe_sample <- spe[, spe$sample_id == "Br8325"]

# Get spatial coordinates and cluster info
plot_df <- data.frame(
    x = spatialCoords(spe_sample)[, "pxl_col_in_fullres"],
    y = spatialCoords(spe_sample)[, "pxl_row_in_fullres"],
    cluster = spe_sample$BS_k25
)

# Simple faceted plot
pdf(here("plots", "Visium", "07_clustering", "BayesSpace", "SVGs", "BS_k25_faceted_Br8325.pdf"), 
    width = 20, height = 20)

ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(aes(color = cluster), size = 1.5) +
    facet_wrap(~cluster, ncol = 5) +
    scale_y_reverse() +  # Flip y-axis for proper spatial orientation
    coord_fixed() +  # Keep aspect ratio
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Clusters faceted - Sample Br8325")

dev.off()

# Or if you want to highlight each cluster individually in its facet
plot_df$is_current_cluster <- plot_df$cluster

pdf(here("plots", "Visium", "07_clustering", "BayesSpace", "SVGs", "BS_k25_faceted_highlighted_Br8325.pdf"), 
    width = 20, height = 20)

ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(data = transform(plot_df, cluster = NULL), color = "grey90", size = 1) +
    geom_point(aes(color = is_current_cluster), size = 1.5) +
    facet_wrap(~cluster, ncol = 5) +
    scale_y_reverse() +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Clusters highlighted - Sample Br8325")

dev.off()