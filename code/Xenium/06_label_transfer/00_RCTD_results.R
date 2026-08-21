library(here)
library(patchwork)
library(scater)
library(scran)
library(ggspavis)
library(SpaceTrooper)
library(ggplot2)
library(dplyr)

# set directories
plots_dir <- here("plots", "Xenium", "06_label_transfer")
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")

# load xenium data
spe <- readRDS(file.path("processed-data", "Xenium", "04_clustering", "Banksy", "Banksy_integrated_res2.0_collapsed_v6.rds"))
spe

# load
res <- readRDS(here(processed_dir, "rctd_results.rds"))

ws <- assay(res, "weights")
table(colSums(ws) == 0)
#  FALSE   TRUE 
# 421643  60876 

table(res$first_type)



# ========= Adding first predicted cell type ==========
imgData(spe) <- NULL # subsetting is throwing an error due to imgData

# subset to common cells (colData) between spe and res
common_cells <- intersect(colnames(spe), rownames(colData(res)))
spe <- spe[, common_cells]
res <- res[, common_cells]

spe$first_type <- colData(res)$first_type


# ========== Adding broad cell type labels =============

sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce

unique(sce$broad_celltype)

table(sce$broad_celltype, sce$fine_celltype)

# add broad cell type labels to spe, based on broad - fine mapping in the sce
broad_fine_mapping <- colData(sce) %>% 
  as.data.frame() %>% 
  select(broad_celltype, fine_celltype) %>% 
  distinct()

spe$first_type_broad <- NA
for (i in 1:nrow(broad_fine_mapping)) {
  broad <- broad_fine_mapping$broad_celltype[i]
  fine <- broad_fine_mapping$fine_celltype[i]
  spe$first_type_broad[spe$first_type == fine] <- broad
}

table(spe$first_type_broad, spe$first_type)

# save
saveRDS(spe, file = here(processed_dir, "Banksy_integrated_res2.0_collapsed_v6_wRCTD_labels.rds"))



# ======= plotting =======
pal <- scCustomize::DiscretePalette_scCustomize(num_colors = length(unique(res$first_type)), palette = "polychrome")
pdf(here(plots_dir, "rctd_first_type.pdf"), width=12, height=6)
plotSpots(spe, sample_id="brnum", annotate="first_type", point_size=0.1, in_tissue=NULL) +
    scale_color_manual(values=pal) +
    ggtitle("RCTD first predicted cell type")
dev.off()

# make first_type a factor
spe$first_type <- factor(spe$first_type, levels=unique(res$first_type))
pdf(file = here(plots_dir, "rctd_first_type.pdf"), width = 12, height = 12)
plot_list_qc <- list()
for (sample in unique(spe$sample_id)) {
  spe_subset <- spe[, spe$sample_id == sample]
  p <- SpaceTrooper::plotCentroids(spe_subset, colourBy="first_type", size=0.1) + 
  ggtitle(sample) + 
  scale_color_manual(values=pal) 
  plot_list_qc[[sample]] <- p
}
wrap_plots(plot_list_qc, ncol=2)
dev.off()


# subset spe to drop Astrocyte cell types 
spe.subset <- spe[, !spe$first_type %in% c("Astrocyte")]

# replot
pdf(here(plots_dir, "rctd_first_type_no_astrocytes.pdf"), width=12, height=6)
plotCoords(spe.subset, sample_id="sample_id", annotate="first_type", point_size=0.1, in_tissue=NULL) +
    scale_color_manual(values=pal) +

# now with spacetrooper
pdf(file = here(plots_dir, "rctd_first_type_no_astrocytes_spacetrooper.pdf"), width = 12, height = 12)
plot_list_qc <- list()
for (sample in unique(spe.subset$sample_id)) {
  spe_subset2 <- spe.subset[, spe.subset$sample_id == sample]
  p <- SpaceTrooper::plotCentroids(spe_subset2, colourBy="first_type", size=0.1) + 
  ggtitle(sample) + 
  scale_color_manual(values=pal) 
  plot_list_qc[[sample]] <- p
}
wrap_plots(plot_list_qc, ncol=2)
dev.off()


# now plot each cell type one per page using ggplot2. will need to get colData and facet around sample_id

df <- as.data.frame(colData(spe))
df$cell_barcode <- rownames(df)
df$x <- spatialCoords(spe)[,1]
df$y <- spatialCoords(spe)[,2]
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

