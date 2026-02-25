library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)

# set directories
plots_dir <- here("plots", "VisiumHD", "04_clustering")
processed_dir <- here("processed-data", "VisiumHD", "04_clustering")

# load
sfe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "sfe_qc.rds"))


sfe <- logNormCounts(sfe)

# load
res <- readRDS(here(processed_dir, "rctd_results.rds"))

ws <- assay(res, "weights")
table(colSums(ws) == 0)
#  FALSE   TRUE 
# 193678  73814 

table(res$first_type)



# ========= Plottig first predicted cell type ==========
imgData(sfe) <- NULL # subsetting is throwing an error due to imgData

# subset to common cells (colData) between sfe and res
common_cells <- intersect(colnames(sfe), rownames(colData(res)))
sfe <- sfe[, common_cells]
res <- res[, common_cells]

sfe$first_type <- colData(res)$first_type

pal <- scCustomize::DiscretePalette_scCustomize(num_colors = length(unique(res$first_type)), palette = "polychrome")
pdf(here(plots_dir, "rctd_first_type.pdf"), width=12, height=6)
plotCoords(sfe, sample_id="sample_id", annotate="first_type", point_size=0.1, in_tissue=NULL) +
    scale_color_manual(values=pal) +
    ggtitle("RCTD first predicted cell type")
dev.off()

# subset sfe to drop Astrocyte cell types 
sfe.subset <- sfe[, !sfe$first_type %in% c("Astrocyte")]

# replot
pdf(here(plots_dir, "rctd_first_type_no_astrocytes.pdf"), width=12, height=6)
plotCoords(sfe.subset, sample_id="sample_id", annotate="first_type", point_size=0.1, in_tissue=NULL) +
    scale_color_manual(values=pal) +
    ggtitle("RCTD first predicted cell type (no Astrocytes)")
dev.off()


# now plot each cell type one per page using ggplot2. will need to get colData and facet around sample_id

df <- as.data.frame(colData(sfe))
df$cell_barcode <- rownames(df)
df$x <- spatialCoords(sfe)[,1]
df$y <- spatialCoords(sfe)[,2]
df <- df %>% select(cell_barcode, sample_id, first_type, x, y)



celltypes <- unique(df$first_type)

pdf(here(plots_dir, "rctd_first_type_all.pdf"), width = 12, height = 6)
for (ct in celltypes) {
  p <- ggplot(df, aes(x = x, y = y, color = first_type == ct)) +
    geom_point(size = 0.1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightgrey")) +
    facet_wrap(~sample_id) +
    ggtitle(paste0("RCTD first predicted cell type: ", ct)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  print(p)
}
dev.off()


# ======= comparing to max weights =========
# add proportion estimates as metadata
ws <- data.frame(t(as.matrix(ws)))
colData(sfe)[names(ws)] <- ws[colnames(sfe), ]

pdf(here(plots_dir, "rctd_celltype_proportions.pdf"), width=10, height=7)
    plotSpatialFeature(sfe, features="Oligodendrocyte", point_size=0.15, point_shape=15) +
    scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket")))
dev.off()



# ======== Plotting cell types ========

ids <- names(ws)[apply(ws, 1, which.max)]
ids <- gsub("\\.([A-z])", " \\1", ids)
idx <- match(colnames(sfe), rownames(ws))
table(sfe$RCTD_decon1<- factor(ids[idx]))

pdf(here(plots_dir, "rctd_celltype_max.pdf"), width=12, height=9)
    plotSpatialFeature(sfe, features="RCTD_decon1")
dev.off()



celltypes <- colnames(ws)
# plot one cell type per page
pdf(here(plots_dir, "celltype_weights_all.pdf"), width=12, height=6)
    for (ct in celltypes) {
        print(plotCoords(sfe, sample_id="sample_id", annotate=ct, point_size=0.1, in_tissue=NULL) +
        scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket"))) +
        ggtitle(ct))
    }
dev.off()