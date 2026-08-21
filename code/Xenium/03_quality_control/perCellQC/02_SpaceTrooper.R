library(here)
library(SpaceTrooper)
library(scater)
library(scran)
library(ggspavis)
library(patchwork)

processed_dir <- here("processed-data", "Xenium", "03_quality_control")
plot_dir <- here("plots", "Xenium", "03_quality_control")

load(here("processed-data","Xenium", "02_build_spe", "spe_xenium_5um_spacetrooper.Rdata"))
spe
# class: SpatialExperiment 
# dim: 541 1018068 
# metadata(8): polygons technology ... polygons technology
# assays(1): counts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0344
#   DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018068): aaaaeomf-1 aaaajkhp-1 ... oimbboka-1 oimbcgpk-1
# colData names(13): cell_id transcript_counts ... polygons brnum
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id



# ======= QC with SpaceTrooper =======

spe <- spatialPerCellQC(spe, rmZeros=TRUE,
        negProbList=c("NegPrb", "Negative", "SystemControl"))
# Removing 389 cells with 0 counts!

colnames(colData(spe))
#  [1] "cell_id"                    "transcript_counts"         
#  [3] "control_probe_counts"       "control_codeword_counts"   
#  [5] "unassigned_codeword_counts" "deprecated_codeword_counts"
#  [7] "total_counts"               "nucleus_area"              
#  [9] "sample_id"                  "AspectRatio"               
# [11] "Area_um"                    "polygons"                  
# [13] "brnum"                      "sum"                       
# [15] "detected"                   "total"                     
# [17] "control_sum"                "control_detected"          
# [19] "target_sum"                 "target_detected"           
# [21] "x_centroid"                 "y_centroid"                
# [23] "ctrl_total_ratio"           "log2Ctrl_total_ratio"      
# [25] "log2AspectRatio"            "SignalDensity"             
# [27] "log2SignalDensity"   

# QC by area
spe <- computeSpatialOutlier(spe, computeBy="Area_um", method="both")
spe <- computeQCScore(spe)
spe <- computeQCScoreFlags(spe, qsThreshold=0.5)

table(spe$low_qcscore, spe$brnum)
#         Br9017 Br9192 Br9206 Br9280
#   FALSE 242411 266993 250289 186559
#   TRUE   10354  12469  36944  11660

# plot QC score as below
pdf(file = here(plot_dir, "QC_score_all_FOVs.pdf"), width = 12, height = 12)
plot_list_qc <- list()
for (sample in unique(spe$sample_id)) {
  spe_subset <- spe[, spe$sample_id == sample]
  p <- plotCentroids(spe_subset, colourBy="QC_score", size=0.1) + ggtitle(sample) + coord_flip() +  scale_x_reverse() 
  plot_list_qc[[sample]] <- p
}
wrap_plots(plot_list_qc, ncol=2)
dev.off()


# QC score violins by sample
pdf(file = here(plot_dir, "QC_score_violin_by_sample.pdf"), width = 10, height = 7)
p1 <- plotColData(spe, y="QC_score", x="sample_id", colour_by="sample_id") +
      theme(axis.text.x = element_blank()) +
      ggtitle("QC Score by Sample ID")
p1
dev.off()

# Sum violins by sample
pdf(file = here(plot_dir, "Sum_violin_by_sample.pdf"), width = 10, height = 7)
p2 <- plotColData(spe, y="sum", x="sample_id", colour_by="sample_id") +
      theme(axis.text.x = element_blank()) +
      ggtitle("Sum by Sample ID")
p2
dev.off()

# detected violins by sample
pdf(file = here(plot_dir, "Detected_violin_by_sample.pdf"), width = 10, height = 7)
p3 <- plotColData(spe, y="detected", x="sample_id", colour_by="sample_id") +
      theme(axis.text.x = element_blank()) +
      ggtitle("Detected Genes by Sample ID")
p3
dev.off()









library(scattermore)
library(patchwork)
library(here)

cd <- as.data.frame(colData(spe))

make_violin <- function(df, y_col, thresh = NULL, log_y = FALSE, title = NULL) {
  p <- ggplot(df, aes(x = brnum, y = .data[[y_col]], colour = brnum, fill = brnum)) +
    geom_violin(alpha = 0.3, scale = "width") +
    geom_scattermore(
      position = position_jitter(width = 0.4, height = 0),
      pointsize = 0.5,
      pixels = c(1024, 1024)
    ) +
    labs(x = "brnum", y = y_col, colour = "brnum", fill = "brnum", title = title) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # red dashed cutoff line
  if (!is.null(thresh)) {
    p <- p + geom_hline(yintercept = thresh, linetype = "dashed", colour = "red")
  }
  if (log_y) p <- p + scale_y_log10()
  p
}

p1 <- make_violin(cd, "sum",       thresh = 10,  log_y = TRUE,  title = "Counts per cell")
p2 <- make_violin(cd, "detected",  thresh = 10,  log_y = TRUE,  title = "Detected genes per cell")
p3 <- make_violin(cd, "QC_score",  thresh = 0.5, log_y = FALSE, title = "QC score")
p4 <- make_violin(cd, "Area_um",   thresh = NULL, log_y = TRUE, title = "Cell area (um^2)")

png(here(plot_dir, "violins_xenium_qc_metrics.png"), width = 12, height = 3.5, units = "in", res = 300)
(p1 + p2 + p3 + p4) + plot_layout(ncol = 4, guides = "collect")
dev.off()



















# ====== Normalization =======

# make cell and nucleus area scaling factors
spe$cell_area.sf <- spe$Area_um / median(spe$Area_um)
spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area)

# histograms of scaling factor
png(file = here("plots","Xenium", "03_quality_control", "cell_area_scaling_factors.png"), width = 10, height = 5, units = "in", res = 300)
hist(spe$cell_area.sf, breaks = 50, main = "Cell area", xlab = "Scaling factor")
dev.off()


png(file = here("plots","Xenium", "03_quality_control", "nucleus_area_scaling_factors.png"), width = 10, height = 5, units = "in", res = 300)
hist(spe$nucleus_area.sf, breaks = 50, main = "Nucleus area", xlab = "Scaling factor")
dev.off()

# normalize the counts by the nucleus and cell area scaling factors
assay(spe, "nucleus_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$nucleus_area.sf, transform="log", assay.type="counts")
assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$cell_area.sf, transform="log", assay.type="counts")


# ===== Drop low QC cells and save =====

spe_filtered <- spe[, !spe$low_qcscore]
spe_filtered
# class: SpatialExperiment 
# dim: 541 946252 
# metadata(8): polygons technology ... polygons technology
# assays(1): counts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0344
#   DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(946252): aaaaeomf-1 aaaajkhp-1 ... oimbboka-1 oimbcgpk-1
# colData names(31): cell_id transcript_counts ... QC_score low_qcscore
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# save
saveRDS(spe_filtered, file = here(processed_dir, "spe_spacetrooper_QCed.rds"))