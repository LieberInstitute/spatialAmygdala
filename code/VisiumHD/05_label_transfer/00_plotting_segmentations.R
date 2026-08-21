library(SpotSweeper)
library(RANN)
library(here)
library(SpatialExperiment)
library(spatialLIBD)
library(escheR)
library(patchwork)
library(scran)
library(scater)
library(ggspavis)

plot_dir <- here("plots","VisiumHD","05_label_transfer", "segmentations")


spe <- readRDS(here("processed-data/VisiumHD/02_build_spe/Br9280_CeA_Spatial.Polygons_spe.rds"))
spe
# dim: 18085 50908 
# metadata(0):
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(1): gene_id
# colnames(50908): cellid_000000004-1 cellid_000000007-1 ...
#   cellid_000060097-1 cellid_000060106-1
# colData names(9): orig.ident nCount_Spatial.008um ... sample_id
#   in_tissue
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x y
# imgData names(4): sample_id image_id data scaleFactor


# ======= Plot segmentations ========


spe <- logNormCounts(spe)

# CeA marekrs
pdf(file.path(plot_dir, "spotplot_segmentations_cea.pdf"), width = 12, height = 18)
p1 <- plotSpots(spe, annotate="PENK")
p2 <- plotSpots(spe, annotate="SST")
p3 <- plotSpots(spe, annotate="PPP1R1B")
p4 <- plotSpots(spe, annotate="PRKCD")
p5 <- plotSpots(spe, annotate="CARTPT")
p6 <- plotSpots(spe, annotate="CRH")
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol=2))
dev.off()


# itc markers
pdf(file.path(plot_dir, "spotplot_segmentations_ITCs.pdf"), width = 12, height = 18)
p1 <- plotSpots(spe, annotate="TSHZ1")
p2 <- plotSpots(spe, annotate="FOXP2")

p3 <- plotSpots(spe, annotate="CPNE4")
p4 <- plotSpots(spe, annotate="DSCAM")

p5 <- plotSpots(spe, annotate="PRKG1")
p6 <- plotSpots(spe, annotate="DRD3")
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol=2))
dev.off()

# other markers
pdf(file.path(plot_dir, "spotplot_segmentations_gaba.pdf"), width = 12, height = 18)
p1 <- plotSpots(spe, annotate="PVALB")
p2 <- plotSpots(spe, annotate="LAMP5")

p3 <- plotSpots(spe, annotate="VIP")
p4 <- plotSpots(spe, annotate="CCK")

p5 <- plotSpots(spe, annotate="GAD1")
p6 <- plotSpots(spe, annotate="GAD2")
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol=2))
dev.off()

