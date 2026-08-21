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
spe <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "spe_cells_combined.rds"))
spe
# class: SpatialExperiment 
# dim: 18085 284986 
# metadata(20): resouces spatialList ... cellseg boundary
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(2): ID Symbol
# colnames(284986): cellid_000000004-1 cellid_000000007-1 ...
#   cellid_000063230-1 cellid_000063232-1
# colData names(2): sample_id map
# reducedDimNames(0):
# mainExpName: Gene Expression
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


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


# add perCellQC metrics - match "MT-" with grep
spe <- addPerCellQC(spe, subsets=list(Mito=grep("^MT-", rowData(spe)$Symbol)))
spe$log_sum <- log1p(spe$sum)
spe$log_detected <- log1p(spe$detected)


# spot plots of log_sum across donors
pdf(here(plots_dir, "spotplots_cells_sum_by_donor.pdf"), width = 9, height = 6)
plotCoords(spe, sample_id = "sample_id", annotate="log_sum", point_size=0.1, in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()

# spot plots of log_detected across donors
pdf(here(plots_dir, "spotplots_cells_detected_by_donor.pdf"), width = 9, height = 6)
plotCoords(spe, sample_id = "sample_id", annotate="log_detected", point_size=0.1, in_tissue=NULL) +
    scale_color_viridis_c()
dev.off()

# spot plots of subsets_Mito_percent across donors
pdf(here(plots_dir, "spotplots_cells_mito_percent_by_donor.pdf"), width = 9, height = 6)
plotCoords(spe, sample_id = "sample_id", annotate="subsets_Mito_percent", point_size=0.1, in_tissue=NULL) +
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

p1 <- make_plot(cd, "sum", thresh = 25) + scale_y_log10()
p2 <- make_plot(cd, "detected", thresh = 25) + scale_y_log10()
p3 <- make_plot(cd, "subsets_Mito_percent", thresh = 25)

pdf(here(plots_dir, "violins_cells_qc_metrics.pdf"), width = 9, height = 3.5)
p1 + p2 + p3
dev.off()



# ======= Dropping low quality spots =======

# drop bins with < 25 UMIs, < 25 detected genes
qc_filter <- (spe$sum > 25) & (spe$detected > 25)  & (spe$subsets_Mito_percent < 25)
table(qc_filter)
#  FALSE   TRUE 
#   2825 282161 

# drop
spe$qc_filter <- qc_filter
spe_qc <- spe[, spe$qc_filter]
spe_qc
# class: SpatialExperiment 
# dim: 18085 282161 
# metadata(20): resouces spatialList ... cellseg boundary
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(2): ID Symbol
# colnames(282161): cellid_000000032-1 cellid_000000033-1 ...
#   cellid_000063230-1 cellid_000063232-1
# colData names(11): sample_id map ... log_detected qc_filter
# reducedDimNames(0):
# mainExpName: Gene Expression
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# save
saveRDS(spe_qc, here(processed_dir, "spe_qc_cells.rds"))