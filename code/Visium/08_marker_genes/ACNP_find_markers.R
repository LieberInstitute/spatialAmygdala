suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("DeconvoBuddies")
})

plot_dir <- here("plots", "Visium", "08_marker_genes")

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)

#subset to brnum 9280
spe <- spe[, colData(spe)$sample_id == "Br8325"]
spe

# get folders in cluster_Csv
bs_folders <- list.files(here::here("processed-data","Visium", "08_clustering", "BayesSpace","HVGs","Br8325"), full.names = TRUE)
bs_folders
# [1] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/HVGs/cluster_csv/BayesSpace_10"
# [2] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/HVGs/cluster_csv/BayesSpace_12"

# get the number of clusters at the end of the folder name
bs_k <- gsub(".*BayesSpace_", "", bs_folders)

# loop through each folder, open the csv inside, and add clusters to spe colData
for (i in seq_along(bs_folders)) {
    bs_folder <- bs_folders[i]
    bs_csv <- list.files(bs_folder, full.names = TRUE)
    bs_csv <- bs_csv[grepl("csv", bs_csv)]
    bs_csv <- bs_csv[1]
    bs_df <- read.csv(bs_csv)

    colData(spe)[[paste0("BS_k", bs_k[i])]] <- factor(bs_df$cluster)
}

colnames(colData(spe))



# ========== Plotting 



# ========== Manually collapse clusters to get desired AMY regions  ===========

unique(spe$BS_k20) 
# [1] 2  7  8  6  9  12 20 17 3  14 4  13 16 11 15 10 5  18 19 1 
# Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

# make a grid plot of the different spatial clusters with BS_k10
# Create a dataframe of spatial coordinates and clusters
df <- data.frame(
    x = spatialCoords(spe)[, 1],
    y = spatialCoords(spe)[, 2],
    cluster = spe$BS_k20
)

# rotate the coordinates 180
df$y <- -df$y

# Plot each cluster in a grid
pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(20)
p <- ggplot(df, aes(x = x, y = y, color = factor(cluster))) +
    geom_point(size = 0.5) +
    facet_wrap(~ cluster, ncol = 5) +
    scale_color_manual(values = pal) +
    #theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Spatial Clusters", x = "X Coordinate", y = "Y Coordinate") 


# Save the plot
ggsave(here("plots", "Visium", "09_marker_genes", "Spatial_Clusters_Grid.png"), plot = p, width = 15, height = 10, units = "in", dpi = 300)

# New annotations fomr BS_k20
# 1 = Ce
# 2 = part of BLVM
# 3 = WM
# 4 = Vascular
# 5 = BADL
# 6 = WM
# 7 = WM
# 8 = WM
# 9 = BLVM
# 10 = HPC
# 11 = PCo
# 12 = LA
# 13 = WM
# 14 = Me
# 15 = PCo
# 16 = WM
# 17 = Vascular
# 18 = BM
# 19 = BL
# 20 = WM

# new $BS_manual clusters with teh above reassignments
spe$BS_manual <- 1
spe$BS_manual[spe$BS_k20 == 1] <- "Ce"
spe$BS_manual[spe$BS_k20 == 2] <- "BLVM"
spe$BS_manual[spe$BS_k20 == 3] <- "WM"
spe$BS_manual[spe$BS_k20 == 4] <- "Vascular"
spe$BS_manual[spe$BS_k20 == 5] <- "BADL"
spe$BS_manual[spe$BS_k20 == 6] <- "WM"
spe$BS_manual[spe$BS_k20 == 7] <- "WM"
spe$BS_manual[spe$BS_k20 == 8] <- "WM"
spe$BS_manual[spe$BS_k20 == 9] <- "BLVM"
spe$BS_manual[spe$BS_k20 == 10] <- "HPC"
spe$BS_manual[spe$BS_k20 == 11] <- "PCo"
spe$BS_manual[spe$BS_k20 == 12] <- "LA"
spe$BS_manual[spe$BS_k20 == 13] <- "WM"
spe$BS_manual[spe$BS_k20 == 14] <- "Me"
spe$BS_manual[spe$BS_k20 == 15] <- "PCo"
spe$BS_manual[spe$BS_k20 == 16] <- "WM"
spe$BS_manual[spe$BS_k20 == 17] <- "Vascular"
spe$BS_manual[spe$BS_k20 == 18] <- "BM"
spe$BS_manual[spe$BS_k20 == 19] <- "BL"
spe$BS_manual[spe$BS_k20 == 20] <- "WM"

# replotting
library(escheR)
pdf(here("plots", "Visium", "09_marker_genes", "Spatial_Clusters_Grid_BS_manual.pdf"), width = 10, height = 10)
p <- make_escheR(spe) |>
    add_fill(var="BS_manual", point_size=1.75) +
    scale_fill_manual(values = pal)
print(p)
dev.off()









# Designing a custom color palette using colorRampPalette for shades of green and blue
# WM = grey
# Ce and Me = different shades of green
# HPC = orange
# PCo, LA, BLVM, BL, BM, BADL = different shades of blue
# vascular = dark grey

WM_pal <- c("#D3D3D3")
green_pal <- colorRampPalette(c("#0fdb71", "#05a150"))(2) # Shades of green
Ce_pal <- green_pal[1]
Me_pal <- green_pal[2]
HPC_pal <- c("#FFA500")
blue_pal <- colorRampPalette(c("#03dffc", "#038cfc"))(6) # Shades of blue
PCo_pal <- blue_pal[1]
LA_pal <- blue_pal[2]
BLVM_pal <- blue_pal[3]
BL_pal <- blue_pal[4]
BM_pal <- blue_pal[5]
BADL_pal <- blue_pal[6]
Vascular_pal <- c("#808080")

pal <- c(WM_pal, Ce_pal, Me_pal, HPC_pal, PCo_pal, LA_pal, BLVM_pal, BL_pal, BM_pal, BADL_pal, Vascular_pal)
names(pal) <- c("WM", "Ce", "Me", "HPC", "PCo", "LA", "BLVM", "BL", "BM", "BADL", "Vascular")

# replotting
pdf(here("plots", "Visium", "09_marker_genes", "Spatial_Clusters_Grid_BS_manual_custom_palette.pdf"), width = 10, height = 10)
p <- make_escheR(spe) |>
    add_fill(var="BS_manual", point_size=1.75) +
    scale_fill_manual(values = pal)
print(p)

dev.off()




# ======== MARKER GENES =========

# renormalize
spe <- scuttle::logNormCounts(spe)

rownames(spe) <- rowData(spe)$symbol

markers <- scran::findMarkers(
  spe,
  groups=spe$BS_manual,
  test.type = c("t"),
  pval.type = c("all"),
  full.stats = TRUE,
  sorted = TRUE,
  direction="up"
)

# get top 10 gene names for each cluster
top_markers <- lapply(markers, function(x) {
  top10 <- x[1:10,]
  top10_genes <- rownames(top10)
  top10_genes
})
# $`1`
#  [1] "PENK"    "SYNPR"   "TMEM272" "GAD2"    "SLC32A1" "SLC35F1" "GPR88"  
#  [8] "NLRP1"   "TSHZ1"   "GAD1"   


# plot violins of the top 10 markers for each cluster in top_markers
for (i in seq_along(top_markers)) {
    cluster_name <- names(top_markers[i])
    png(here("plots", "Visium", "09_marker_genes", paste0("Markers_cluster_", cluster_name, "_violin.png")), width = 5, height = 10, units = "in", res = 300)
        p <- scater::plotExpression(spe, features=top_markers[[i]], x="BS_manual", colour_by="BS_manual", ncol=2)
        print(p)
    dev.off()
}



# BADL = PEX5L/EDIL3/STXBP6/COL25A1
# BL = MTPN/COL25A1/NRXN2
# BLVM = LAMP5/NPTX1/STMN4/ATP2B4
# BM = NCAM2/CNR1/CCK/GABRD
# PCo = PDYN/CDH13/ESR1
# Ce = PENK/SYNPR/GAD2/TSHZ1/SST/PRKCD
# LA = CAMK2N1/TTC9B/CYP26B1/SLC30A3/ARPP19
# Me = CALB2/CALB1/ CARTPT/SLC17A6/GABRE/OTP

features <- c("PEX5L", "EDIL3", "STXBP6", "COL25A1",
                "MTPN",  "NRXN2",
                "LAMP5", "NPTX1", "STMN4", "ATP2B4",
                "NCAM2", "CNR1", "CCK", "GABRD",
                "PDYN", "CDH13", "ESR1",
                "PENK", "SYNPR", "GAD2", "TSHZ1", "SST", "PRKCD",
                "CAMK2N1", "TTC9B", "CYP26B1", "SLC30A3", "ARPP19",
                "CALB2", "CALB1", "CARTPT", "SLC17A6", "GABRE", "OTP")


spe$BS_manual <- factor(spe$BS_manual)

# grouped heatmap
pdf(here("plots", "Visium", "09_marker_genes", "grouped_heatmap.pdf"), width = 10, height = 10)
p <- scater::plotDots(spe, features, group="BS_manual", center=TRUE, scale=TRUE)
print(p)
dev.off()



# subset to only amygdala regions
amygdala_regions <- c("Ce", "BLVM", "LA", "Me", "PCo", "BADL", "BL", "BM")
spe.amy<- spe[, spe$BS_manual %in% amygdala_regions]
spe.amy
# class: SpatialExperiment 
# dim: 36601 19868 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(1): symbol
# colnames(19868): AAACAAGTATCTCCCA-1_V13M06-387_A1
#   AAACAGAGCGACTCCT-1_V13M06-387_A1 ... TTGTTTCATTAGTCTA-1_V13M06-388_D1
#   TTGTTTCCATACAACT-1_V13M06-388_D1
# colData names(57): in_tissue array_row ... BS_k9 BS_manual
# reducedDimNames(4): PCA HARMONY UMAP-PCA UMAP-HARMONY
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor


# heatmap of known marker genes
known <- c("GULP1", "SATB1","COL25A1","PEX5L", "ESR1")

# blue to red coloscale
png(here("plots", "Visium", "09_marker_genes", "amygdala_heatmap_known_markers.png"), width = 4, height = 3, units = "in", res = 300)
p <- scater::plotGroupedHeatmap(spe.amy, known, group="BS_manual", center=TRUE, scale=TRUE, treeheight_row = 0, treeheight_col = 0)
print(p)
dev.off()


novel <- c("STXBP6", "EDIL3", "NRXN2", "MTPN", "STMN4", "LAMP5", "CNR1", "NCAM2", "PENK", "SYNPR", "CYP26B1", "TTC9B", "CALB2", "CARTPT", "PDYN", "CDH13")

# blue to red coloscale
png(here("plots", "Visium", "09_marker_genes", "amygdala_heatmap_novel_markers.png"), width = 4, height = 5, units = "in", res = 300)
p <- scater::plotGroupedHeatmap(spe.amy, novel, group="BS_manual", center=TRUE, scale=TRUE, cluster_cols=FALSE, cluster_rows=FALSE, legend=TRUE)
print(p)
dev.off()