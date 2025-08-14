library("here")
library("SpatialExperiment")
library("SingleCellExperiment")
library("ggplot2")
library("Banksy")
library("escheR")
library("data.table")
library("RColorBrewer")
library("patchwork")

# Directories
plot_dir <- here("plots", "Visium", "07_clustering", "BANKSY")
rds_dir <- here("processed-data", "Visium", "07_clustering", "BANKSY")

# Make sure output directory exists
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Loop through k=2 to 20
for (k in 2:20) {
    message("Processing k = ", k)

    # File path
    rds_file <- file.path(rds_dir, paste0("spe_banksy_kmeans_", k, ".rds"))
    if (!file.exists(rds_file)) {
        warning("Missing file: ", rds_file)
        next
    }

    # Read RDS
    spe <- readRDS(rds_file)

    # Get column name of clustering result (should be last added or predictable)
    clust_col <- grep("^clust_.*kmeans.*", colnames(colData(spe)), value = TRUE)
    if (length(clust_col) != 1) {
        warning("Expected 1 clustering column, found ", length(clust_col), " in k = ", k)
        next
    }

    # Make color palette
    n_clusters <- length(unique(colData(spe)[[clust_col]]))
    pal <- colorRampPalette(brewer.pal(9, "Set1"))(n_clusters)

    # Open PDF
    pdf_file <- file.path(plot_dir, paste0("banksy_k", k, ".pdf"))
    pdf(pdf_file, width = 10, height = 10)

    # Plot per sample
    for (id in unique(colData(spe)$sample_id)) {
        spe.sub <- spe[, spe$sample_id == id]
        p <- make_escheR(spe.sub) |>
            add_fill(var = clust_col, point_size = 1.5) +
            scale_fill_manual(values = pal) +
            ggtitle(paste0("Sample: ", id, " | k = ", k))
        print(p)
    }

    dev.off()
    message("Saved: ", pdf_file)
}
