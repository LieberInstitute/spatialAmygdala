library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("harmony")
library("Seurat")

# save directories
processed_dir <- here("processed-data", "Xenium", "04_dim_reduction")
plot_dir <- here("plots", "Xenium", "04_dim_reduction")

spe <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_proseg_5um_harmonized_singlecell.rds"))
spe


# convert to seurat object
colnames(spe) <- make.unique(colnames(spe))

obj <- as.Seurat(spe, counts = "counts", data = "cell_normcounts")
obj
# An object of class Seurat 
# 366 features across 1039305 samples within 1 assay 
# Active assay: originalexp (366 features, 0 variable features)
#  3 dimensional reductions calculated: PCA, HARMONY, UMAP



# ======= CCA integration ========

seurat_list <- SplitObject(obj, split.by = "brnum")


anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method= "LogNormalize",
  dims = 1:30
)

integrated_obj <- IntegrateData(
  anchorset = anchors
)

DefaultAssay(integrated_obj) <- "integrated"

integrated_obj <- ScaleData(integrated_obj)
integrated_obj <- RunPCA(integrated_obj)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30)

# plot umap
png(file.path(plot_dir, "seurat_cca_umap.png"), width = 5, height = 5, units = "in", res = 300)
DimPlot(integrated_obj, reduction = "umap", group.by = "brnum") +
  ggtitle("Seurat CCA UMAP") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

saveRDS(integrated_obj, file = file.path(processed_dir, "proseg_5um_cca_integrated_obj.rds"))
