library(here)
library(SpatialExperiment)
library(ggplot2)
library(patchwork)
library(matrixStats)
library(scales)
library(viridisLite)

# ==============================================================================
# Load
# ==============================================================================
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                    "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

modules_df <- read.csv(here("processed-data", "Visium", "14_smoothie_coexpression",
                            "results", "modules_df_0.8_9.csv"))

plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

spe <- spe[, spe$sample_id == "Br8325"]
spe <- spe[, !spe$exclude_overlapping]
cat("Spots in Br8325:", ncol(spe), "\n")

# ==============================================================================
# Colormap: pale cream (low) -> orange -> red -> near-black (high)
# ==============================================================================
rctd_pal <- viridisLite::magma(20, direction = -1)

# ==============================================================================
# Modules
# ==============================================================================
plot_modules <- c(
    "6"  = "LA",
    "27" = "BLD",
    "5"  = "PL",
    "10" = "CoA"
)

pal <- c(
     BLD = "#9B59B6", PL  = "#f1e438ff", LA  = "#F4B400", CoA = "#5DA5DA"
)

# ==============================================================================
# Module scores (Smoothie method), then rescale each to 0-1
# ==============================================================================
lc <- logcounts(spe)
coords <- data.frame(
    x = spatialCoords(spe)[, "pxl_col_in_fullres"],
    y = -spatialCoords(spe)[, "pxl_row_in_fullres"]
)

for (mod in names(plot_modules)) {
    mod_genes <- modules_df$name[modules_df$module_label == as.integer(mod)]
    present   <- mod_genes[mod_genes %in% rownames(spe)]

    if (length(present) == 0) {
        coords[[paste0("M", mod)]] <- 0
        next
    }

    mat <- as.matrix(lc[present, , drop = FALSE])
    row_maxs <- matrixStats::rowMaxs(mat)
    row_maxs[row_maxs == 0] <- 1
    s <- colSums(mat / row_maxs)

    # scale to 0-1 against the 99th pct so all panels share one colorbar
    hi <- stats::quantile(s, 0.99, na.rm = TRUE)
    if (hi <= 0) hi <- 1
    coords[[paste0("M", mod)]] <- s / hi
}

# ==============================================================================
# RCTD-style panel
# ==============================================================================
plot_module_rctd <- function(df, col, title, title_col = "black",
                             ylimit = c(0, 1), size = 0.0025, sort_points = TRUE,
                             legend_name = "Module\nscore") {

    if (sort_points) df <- df[order(df[[col]]), ]   # high spots drawn last

    ggplot(df, aes(x = x, y = y, color = .data[[col]])) +
        geom_point(size = size, shape = 19) +
        scale_color_gradientn(colours = rctd_pal,
                              limits  = ylimit,
                              oob     = scales::squish,
                              name    = legend_name) +
        coord_fixed() +
        ggtitle(title) +
        scale_x_continuous(breaks = NULL) +
        scale_y_continuous(breaks = NULL) +
        theme_bw() +
        theme(
            panel.grid       = element_blank(),
            panel.background = element_blank(),
            panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
            axis.title       = element_blank(),
            axis.text        = element_blank(),
            axis.ticks       = element_blank(),
            plot.title   = element_text(hjust = 0.5, size = 11, face = "bold",
                                        colour = title_col),
            legend.position   = "right",
            legend.title      = element_text(size = 7, hjust = 0.5, lineheight = 0.9),
            legend.key.height = unit(0.55, "cm"),
            legend.key.width  = unit(0.18, "cm"),
            legend.text       = element_text(size = 6),
            plot.margin  = margin(2, 2, 2, 2)
        )
}

# ==============================================================================
# Build + assemble
# ==============================================================================
plot_list <- list()

for (mod in names(plot_modules)) {
    domain <- plot_modules[mod]
    plot_list[[mod]] <- plot_module_rctd(
        coords,
        col       = paste0("M", mod),
        title     = paste0("M", mod, " - ", domain),
        title_col = pal[domain]
    )
}

# single shared colorbar for all panels
combined <- wrap_plots(plot_list, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

ggsave(file.path(plots_dir, "module_scores_Br8325_RCTDstyle.pdf"),
       combined, width = 10, height = 2.5)
ggsave(file.path(plots_dir, "module_scores_Br8325_RCTDstyle.png"),
       combined, width = 10, height = 2.5, dpi = 300)

cat("Saved: module_scores_Br8325_RCTDstyle.pdf / .png\n")