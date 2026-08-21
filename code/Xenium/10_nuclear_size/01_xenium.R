library(crumblr)
library(variancePartition)
library(dreamlet)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(scattermore)


# save directories
processed_dir <- here("processed-data", "Xenium", "09_crumblr")
plot_dir <- here("plots", "Xenium", "09_crumblr")

# load xenium data
spe.old <- readRDS(here("processed-data","Xenium","04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe.new <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering", "Banksy_domains_v1.0.rds"))

# load label predictions (aligned to spe.old)
pred.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype.csv"))

# add predictions to spe.old first, while IDs still match
stopifnot(nrow(pred.broad) == ncol(spe.old))
spe.old$pred_broad_celltype <- pred.broad$pruned.labels
spe.old$pred_fine_celltype <- pred.fine$pruned.labels

# subset to common cells
common_cells <- intersect(colnames(spe.old), colnames(spe.new))
spe.old <- spe.old[, common_cells]
spe.new <- spe.new[, common_cells]

# copy predictions from spe.old to spe.new
spe.new$pred_broad_celltype <- spe.old$pred_broad_celltype
spe.new$pred_fine_celltype <- spe.old$pred_fine_celltype

# rename to spe
spe <- spe.new
rm(spe.old, spe.new)

colnames(spe) <- make.unique(colnames(spe), sep = "-")
rownames(spatialCoords(spe)) <- colnames(spe)


colnames(colData(spe))



# Plot nucleus area across broad and fine cell types
pdf(file.path(plot_dir, "nucleus_area_by_fine_celltype.pdf"), width = 12, height = 4)
scater::plotColData(spe, y="nucleus_area", x = "pred_fine_celltype", colour_by="pred_fine_celltype") + 
  ggtitle("Nucleus Area by Fine Cell Type") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()