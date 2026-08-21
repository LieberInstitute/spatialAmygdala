library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(scattermore)

# set directories
plots_dir <- here("plots", "VisiumHD", "03_quality_control")
processed_dir <- here("processed-data", "VisiumHD", "03_quality_control")

# load
spe <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "spe_016_combined.rds"))
spe


# ======= Standardize coordinates for each samples =======
# Center coordinates per sample so each tissue is centered at (0,0)
cd <- as.data.frame(colData(spe))
coords <- spatialCoords(spe)
cd$x <- coords[, 1]
cd$y <- coords[, 2]

# Per-sample centering
cd <- cd %>%
  group_by(sample_id) %>%
  mutate(
    x_centered = x - median(x),
    y_centered = y - median(y)
  ) %>%
  ungroup()

# Write back into the spe object
spatialCoords(spe)[, 1] <- cd$x_centered
spatialCoords(spe)[, 2] <- cd$y_centered

# ======= Add perCellQC metrics =======

gradient_theme <- theme(plot.title = element_text(hjust = 0.5),
                        legend.title = element_blank(),
                        legend.key.width = grid::unit(0.5, "lines"), 
                        legend.key.height = grid::unit(1, "lines"))

# add perCellQC metrics - match "MT-" with grep
spe <- addPerCellQC(spe, subsets=list(Mito=grep("^MT-", rowData(spe)$Symbol)))
spe$log_sum <- log1p(spe$sum)
spe$log_detected <- log1p(spe$detected)


# spot all plots across donors
pdf(here(plots_dir, "spotplots_016_sum_by_donor_all.pdf"), width = 9, height = 6)
plotCoords(spe, sample_id = "sample_id", annotate="log_sum", point_size=0.02, point_shape=15, in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()

# spot all plots across donors, but only in tissue
pdf(here(plots_dir, "spotplots_016_sum_by_donor_in_tissue.pdf"), width = 9, height = 6)
plotCoords(spe, sample_id = "sample_id", annotate="log_sum", point_size=0.02, point_shape=15) +
    scale_color_viridis_c()
dev.off()

# drop out of tissue for all except Br8325_MeA
spe_in_tissue <- spe[, spe$in_tissue == 1 | spe$sample_id == "Br8325_MeA"]

# plot all plots across donors, but only in tissue (with Br8325_MeA)
pdf(here(plots_dir, "spotplots_016_sum_by_donor_in_tissue_with_8325_MeA.pdf"), width = 9, height = 6)
plotCoords(spe_in_tissue, sample_id = "sample_id", annotate="log_sum", point_size=0.005, point_shape=15, , in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()


# plot detected genes across donors, but only in tissue (with Br8325_MeA)
pdf(here(plots_dir, "spotplots_016_detected_by_donor_in_tissue_with_8325_MeA.pdf"), width = 9, height = 6) 
plotCoords(spe_in_tissue, sample_id = "sample_id", annotate="log_detected", point_size=0.005, point_shape=15, , in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()

# plot mito % across donors, but only in tissue (with Br8325_MeA)
pdf(here(plots_dir, "spotplots_016_mito_by_donor_in_tissue_with_8325_MeA.pdf"), width = 9, height = 6)
plotCoords(spe_in_tissue, sample_id = "sample_id", annotate="subsets_Mito_percent", point_size=0.005, point_shape=15, , in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()


# visualize QC metrics with violins
cd <- as.data.frame(colData(spe))

make_plot <- function(df, y_col, thresh) {
  ggplot(df, aes(x = sample_id, y = .data[[y_col]], colour = sample_id, fill = sample_id)) +
    geom_violin(alpha = 0.3, scale = "width") +
    geom_scattermore(
      position = position_jitter(width = 0.4, height = 0),
      pointsize = 0.5,
      pixels = c(1024, 1024)
    ) +
    labs(x = "sample_id", y = y_col, colour = "sample_id", fill = "sample_id") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # horizontal line for QC thresholds with "tresh" variable
    if (!is.null(thresh))
        geom_hline(data = data.frame(yintercept = thresh), aes(yintercept = yintercept), linetype = "dashed", color = "red")
}

p1 <- make_plot(cd, "sum", thresh = 10) + scale_y_log10()
p2 <- make_plot(cd, "detected", thresh = 10) + scale_y_log10()
p3 <- make_plot(cd, "subsets_Mito_percent", thresh = NULL)

pdf(here(plots_dir, "violins_qc_metrics_016.pdf"), width = 9, height = 3.5)
p1 + p2 + p3
dev.off()



# ======= Dropping low quality spots =======

# drop bins with < 10 UMIs, < 10 detected genes
qc_filter <- (spe$sum > 10) & (spe$detected > 10) 
table(qc_filter)
#  FALSE   TRUE 
#  24169 853636 



# drop
spe$qc_filter <- qc_filter
spe_qc <- spe[, spe$qc_filter]
spe_qc

# save
saveRDS(spe_qc, here(processed_dir, "spe_qc_016.rds"))