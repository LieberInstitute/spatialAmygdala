library("SpatialFeatureExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("ggspavis")
library("scCustomize")

# ======= Load base spe =======
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_016.rds"))

# drop out of tissue for all except Br8325_MeA
spe <- spe[, spe$in_tissue | spe$sample_id == "Br8325_MeA"]

# ======= Find all Banksy output files =======
banksy_dir <- here("processed-data", "VisiumHD", "04_clustering", "016_wAI_markers")
banksy_files <- list.files(banksy_dir, pattern = "Banksy_integrated_lambda_0\\.8_res.*\\.rds", full.names = TRUE)
cat("Found", length(banksy_files), "Banksy output files:\n")
cat(paste(" ", basename(banksy_files), collapse = "\n"), "\n\n")

# ======= Extract resolution from filename and sort =======
resolutions <- gsub(".*res([0-9.]+).*\\.rds", "\\1", basename(banksy_files))
file_order <- order(as.numeric(resolutions))
banksy_files <- banksy_files[file_order]
resolutions <- resolutions[file_order]

# ======= Loop and plot =======
plots_dir <- here("plots", "VisiumHD", "04_clustering", "016_wAI_markers")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(banksy_files)) {
    res <- resolutions[i]
    cat("Processing resolution:", res, "\n")

    # Load clustered subset
    sfe.sub <- readRDS(banksy_files[i])

    # Find the cluster column (matches pattern clust_HARMONY_M0_lam0.8_k50_resX.X)
    clust_cols <- grep("clust_HARMONY_M1_lam0\\.8", colnames(colData(sfe.sub)), value = TRUE)
    if (length(clust_cols) == 0) {
        cat("  WARNING: no cluster column found, skipping\n")
        next
    }
    clust_col <- clust_cols[1]
    cat("  Cluster column:", clust_col, "\n")

    # Transfer labels to full spe
    spe[[clust_col]] <- factor(sfe.sub[[clust_col]])

    num_groups <- length(unique(spe[[clust_col]]))
    cat("  Number of clusters:", num_groups, "\n")

    colors <- scCustomize_Palette(
        num_groups,
        ggplot_default_colors = FALSE,
        color_seed = 123
    )

    # Plot
    pdf_path <- file.path(plots_dir, paste0("Banksy_integrated_lambda_0.8_res", res, ".pdf"))
    pdf(file = pdf_path, width = 40, height = 10)

    p <- plotCoords(spe, annotate = clust_col, in_tissue = NULL, sample_id = "sample_id") +
        scale_color_manual(values = colors) +
        ggtitle(paste0("Banksy Integrated Clustering (lambda=0.8, res=", res, ", k=", num_groups, " clusters)")) +
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 20),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE))

    print(p)
    dev.off()
    cat("  Saved:", pdf_path, "\n")

    # Clean up cluster column from spe to avoid accumulation
    spe[[clust_col]] <- NULL
}

cat("\nDone! All plots saved to:", plots_dir, "\n")