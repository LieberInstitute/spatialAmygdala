library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("harmony")

# save directories
processed_dir <- here("processed-data", "Xenium", "04_dim_reduction")
plot_dir <- here("plots", "Xenium", "04_dim_reduction")

spe <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_proseg_5um_harmonized_singlecell.rds"))
spe
# dim: 541 1018069 
# metadata(0):
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018069): aaaadgkh-1 aaaadlfe-1 ... oihobmbk-1 oihoeehh-1
# colData names(19): cell_id transcript_counts ... cell_area.sf
#   nucleus_area.sf
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# ======= QC metrics ========

spe <- scuttle::addPerCellQC(spe)

#scatter sum by volume
png(file.path(plot_dir, "proseg_counts_by_volume.png"), width = 5, height = 5, units = "in", res = 300)
scater::plotColData(spe, x = "volume", y = "sum", colour_by = "brnum", scattermore=TRUE) +
    ggtitle("Sum counts") +
    xlab("Volume)") +
    ylab("Sum of counts") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# violin of sum
png(file.path(plot_dir, "proseg_counts_by_donors.png"), width = 5, height = 5, units = "in", res = 300)
scater::plotColData(spe, x = "brnum", y = "sum", colour_by = "brnum", scattermore=TRUE) +
    ggtitle("Donor") +
    ylab("Sum of counts") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# detected
png(file.path(plot_dir, "proseg_detected_by_violin.png"), width = 5, height = 5, units = "in", res = 300)
scater::plotColData(spe, x = "brnum", y = "detected", colour_by = "brnum", scattermore=TRUE) +
    ggtitle("Donor") +
    ylab("Number of detected genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ======== Dim reduction and Harmony batch correction ========

# split good and bad samples into different objects
# good: Br9192 and 9280
# bad: Br9017 and 9206
spe.good <- spe[, spe$brnum %in% c("Br9192", "Br9280")]
spe.bad <- spe[, spe$brnum %in% c("Br9017", "Br9206")]

unique(spe.good$brnum)
# [1] "Br9192" "Br9280"

unique(spe.bad$brnum)
# [1] "Br9017" "Br9206"


#  ==== renormalize =====

# Compute scaling factors
spe.bad$volume.sf <- spe.bad$volume / median(spe.bad$volume, na.rm = TRUE)
spe.good$volume.sf <- spe.good$volume / median(spe.good$volume, na.rm = TRUE)

# normalize
assay(spe.good, "volume_normcounts") <- scuttle::normalizeCounts(
        spe.good, size.factors = spe.good$volume.sf, transform = "log", assay.type = "counts"
    )

assay(spe.bad, "volume_normcounts") <- scuttle::normalizeCounts(
        spe.bad, size.factors = spe.bad$volume.sf, transform = "log", assay.type = "counts"
    )

#  ====== Dim reduction and intergration ====== 
set.seed(1054)

spe.good <- runPCA(spe.good, ncomponents=50, scale=TRUE, name="PCA", assay.type="volume_normcounts")
spe.bad <- runPCA(spe.bad, ncomponents=50, scale=TRUE, name="PCA", assay.type="volume_normcounts")

# run Harmony
spe.good <- RunHarmony(spe.good, group.by.vars="brnum", reduction.save="HARMONY_subset")
spe.bad <- RunHarmony(spe.bad, group.by.vars="brnum",  reduction.save="HARMONY_subset")

# runUMAP
spe.good <- runUMAP(spe.good, dimred="HARMONY_subset")
spe.bad <- runUMAP(spe.bad, dimred="HARMONY_subset")


# ======== Visualization ========



# UMAP
pdf(file.path(plot_dir, "proseg_UMAP_good_samples.pdf"), width = 5, height = 4)
scater::plotReducedDim(spe.good, dimred = "UMAP", colour_by = "brnum", scattermore=TRUE) +
    ggtitle("UMAP of good samples (Br9192 and Br9280)")
dev.off()

# UMAP
pdf(file.path(plot_dir, "proseg_UMAP_bad_samples.pdf"), width = 5, height = 4)
scater::plotReducedDim(spe.bad, dimred = "UMAP", colour_by = "brnum", scattermore=TRUE) +
    ggtitle("UMAP ofbad samples Br9017 and Br9206")
dev.off()

# save
saveRDS(spe.good, here("processed-data", "Xenium","04_dim_reduction", "spe_proseg_good.rds"))
saveRDS(spe.bad, here("processed-data", "Xenium","04_dim_reduction", "spe_proseg_bad.rds"))