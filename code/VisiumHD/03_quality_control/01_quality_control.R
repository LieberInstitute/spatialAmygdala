library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)

# set directories
plots_dir <- here("plots", "VisiumHD", "03_quality_control")
processed_dir <- here("processed-data", "VisiumHD", "03_quality_control")

# load
sfe <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_combined.rds"))
sfe
# class: SpatialFeatureExperiment 
# dim: 18085 284986 
# metadata(0):
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(3): ID Symbol Type
# colnames(284986): 4 7 ... 63230-2 63232-2
# colData names(5): Barcode sample_id pxl_col_in_hires pxl_row_in_hires
#   Sample
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# ======= Add perCellQC metrics =======

gradient_theme <- theme(plot.title = element_text(hjust = 0.5),
                        legend.title = element_blank(),
                        legend.key.width = grid::unit(0.5, "lines"), 
                        legend.key.height = grid::unit(1, "lines"))

sfe <- vhd_all
# add perCellQC metrics - match "MT-" with grep
sfe <- addPerCellQC(sfe, subsets=list(Mito=grep("^MT-", rowData(sfe)$Symbol)))
sfe$log_sum <- log1p(sfe$sum)
sfe$log_detected <- log1p(sfe$detected)

# visualize QC metrics
pdf(here(plots_dir, "violins_qc_metrics.pdf"), width=9, height=3.5)
p1 <- plotColData(sfe, x="Sample", y="sum", colour="Sample") + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- plotColData(sfe, x="Sample", y="detected", colour="Sample") + scale_y_log10() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3 <- plotColData(sfe, x="Sample", y="subsets_Mito_percent", colour="Sample") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1 + p2 + p3
dev.off()


pdf(here(plots_dir, "spotplots_sum_new.pdf"), width=9, height=6)
plotSpatialFeature(sfe, features = "log_sum", colGeometryName = "cellSeg") + 
    ggtitle("log library size") +
    gradient_theme +
    scale_fill_gradientn(colors = pals::jet())
dev.off()

pdf(here(plots_dir, "spotplots_qc_metrics_mito.pdf"), width=9, height=6)
plotSpatialFeature(sfe, features = "subsets_Mito_percent", colGeometryName = "cellSeg") + 
    ggtitle("Mito %") +
    gradient_theme +
    scale_fill_gradientn(colors = pals::jet())
dev.off()


pdf(here(plots_dir, "spotplots_qc_metrics_detected.pdf"), width=9, height=6)
plotSpatialFeature(sfe, features = "log_detected", colGeometryName = "cellSeg") + 
    ggtitle("log detected genes") +
    gradient_theme +
    scale_fill_gradientn(colors = pals::jet())
dev.off()

# ======= Filter out low quality spots =======

# filter out low quality spots based on QC metrics
# library size and number of detected genes: log(sum) or log(detected) <
# mito percent: subsets_Mito_percent > 20%
qc_filter <- (sfe$sum > 50) & (sfe$detected > 50) & (sfe$subsets_Mito_percent < 20)  
table(qc_filter)

sfe$qc_filter <- qc_filter

# visualize filtered spots
pdf(here(plots_dir, "spotplots_qc_filter.pdf"), width=9, height=6)
plotSpatialFeature(sfe, features = "qc_filter", colGeometryName = "cellSeg") + 
    ggtitle("QC filter") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.width = grid::unit(0.5, "lines"), 
          legend.key.height = grid::unit(1, "lines"))
dev.off()


# save
saveRDS(sfe, here(processed_dir, "sfe_qc.rds"))