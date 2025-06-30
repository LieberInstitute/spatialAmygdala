setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')

library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("BayesSpace")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("patchwork")

# Load the harmonized dataset
load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)

# Get sample list
sample_list <- unique(colData(spe)$sample_id)

# Get cluster number from SLURM array task ID
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Loop over each sample_id
for (sample_id in sample_list) {
    message(paste0("Processing sample: ", sample_id))

    # Subset to current sample
    spe_sub <- spe[, colData(spe)$sample_id == sample_id]

    # Preprocessing
    spe_sub <- spatialPreprocess(spe_sub, platform = "Visium", skip.PCA = TRUE)
    colData(spe_sub)$row <- spe_sub$array_row
    colData(spe_sub)$col <- spe_sub$array_col
    metadata(spe_sub)$BayesSpace.data <- list(platform = "Visium", is.enhanced = FALSE)

    # Clustering
    message("Running spatialCluster()")
    Sys.time()
    set.seed(2)
    spe_sub <- spatialCluster(
        spe_sub,
        use.dimred = "PCA-HARMONY_sample",
        q = k,
        nrep = 10000,
        burn.in = 100,
        gamma = 3
    )
    Sys.time()

    # Rename cluster column
    bayesSpace_name <- paste0("BayesSpace_", k)
    colnames(colData(spe_sub))[ncol(colData(spe_sub))] <- bayesSpace_name
    clustV <- bayesSpace_name

    # Define output directories
    out_dir <- here("processed-data", "Visium", "07_clustering", "BayesSpace", "SVGs", sample_id)
    plot_dir <- here("plots", "Visium", "07_clustering", "BayesSpace", "SVGs", sample_id)

    # Create directories if they don't exist
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    # Export cluster object
    cluster_export(
        spe_sub,
        bayesSpace_name,
        cluster_dir = out_dir
    )

    message(paste0("Finished processing sample: ", sample_id))
}

