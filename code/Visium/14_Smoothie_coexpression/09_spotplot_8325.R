library(here)
library(SpatialExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(matrixStats)

# ==============================================================================
# Configuration
# ==============================================================================
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

modules_df <- read.csv(here("processed-data", "Visium", "14_smoothie_coexpression",
                             "results", "modules_df_0.8_9.csv"))

plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Subset to Br8325
spe <- spe[, spe$sample_id == "Br8325"]
cat("Spots in Br8325:", ncol(spe), "\n")

# drop overlapping
spe <- spe[, !spe$exclude_overlapping]
# ==============================================================================
# Modules to plot (module_label -> domain label)
# ==============================================================================
plot_modules <- c(
    "6"  = "LA",
    "27" = "BLD",
    "5"  = "PL",
    "10" = "CoA",
    "13" = "CeA",
    "19" = "MeA",
    "20" = "AI",
    "16" = "HPC"
)

# Domain colors
pal <- c(
    AI  = "#D62728",
    BLD = "#9B59B6",
    PL  = "#f1e438ff",
    LA  = "#F4B400",
    CoA = "#5DA5DA",
    CeA = "#197d43ff",
    MeA = "#baf739ff",
    HPC = "#d6a8f8ff"
)

# ==============================================================================
# Compute module scores (Smoothie method)
# ==============================================================================
lc <- logcounts(spe)
coords <- data.frame(
    x = spatialCoords(spe)[, "pxl_col_in_fullres"],
    y = -spatialCoords(spe)[, "pxl_row_in_fullres"]
)

score_list <- list()
for (mod in names(plot_modules)) {
    mod_genes <- modules_df$name[modules_df$module_label == as.integer(mod)]
    present <- mod_genes[mod_genes %in% rownames(spe)]

    if (length(present) == 0) {
        coords[[mod]] <- 0
        next
    }

    mat <- as.matrix(lc[present, , drop = FALSE])
    row_maxs <- matrixStats::rowMaxs(mat)
    row_maxs[row_maxs == 0] <- 1
    coords[[mod]] <- colSums(mat / row_maxs)
}

# ==============================================================================
# Plot each module
# ==============================================================================
plot_list <- list()

for (mod in names(plot_modules)) {
    domain <- plot_modules[mod]

    p <- ggplot(coords, aes(x = x, y = y, color = .data[[mod]])) +
        geom_point(size = 0.3, stroke = 0) +
        scale_color_viridis_c(option = "magma", direction = 1, name = NULL) +
        labs(title = paste0("M", mod, " - ", domain)) +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 12, face = "bold",
                                      color = pal[domain]),
            legend.key.height = unit(0.3, "cm"),
            legend.key.width = unit(0.15, "cm"),
            legend.text = element_text(size = 6),
            plot.margin = margin(2, 2, 2, 2)
        )

    plot_list[[mod]] <- p
}
# ==============================================================================
# Assemble patchwork (2 rows x 4 columns)
# ==============================================================================
combined <- wrap_plots(plot_list, nrow = 2)

ggsave(file.path(plots_dir, "module_scores_Br8325_patchwork.pdf"),
       combined, width = 8, height = 4)
ggsave(file.path(plots_dir, "module_scores_Br8325_patchwork.png"),
       combined, width = 8, height = 4, dpi = 300)

cat("Saved: module_scores_Br8325_patchwork.pdf / .png\n")