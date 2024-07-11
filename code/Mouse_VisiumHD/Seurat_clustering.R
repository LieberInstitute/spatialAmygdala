# packages required for Visium HD
#install.packages("hdf5r")
#install.packages("arrow")
#
#remotes::install_github(repo = "satijalab/seurat-object", ref = "visium-hd")
#remotes::install_github(repo = "satijalab/seurat", ref = "visium-hd")
#
# NOTE: ^^ done

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)

# Load the data
# alter localdir to the path where you have saved the outs folder of the demo data
localdir <- here("raw-data","mouse_brain")
object <- Load10X_Spatial(data.dir = localdir)
object

colnames(object@meta.data)
#[1] "orig.ident"       "nCount_Spatial"   "nFeature_Spatial"


vln.plot <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial", pt.size.factor = 1.2) +
    theme(legend.position = "right")

png(here("plots","QC.png"), width = 10, height = 5, units = "in", res = 300)
vln.plot | count.plot
dev.off()

# normalize
DefaultAssay(object) <- "Spatial"
object <- NormalizeData(object)

# cluster
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, reduction.name = "pca")
object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
object <- FindClusters(object, resolution = 0.6, cluster.name = "seurat_cluster")

object <- RunUMAP(object, reduction = "pca", reduction.name = "umap", dims = 1:30)

# plotting
dim.plot <- DimPlot(object, reduction = "umap", group.by = "seurat_cluster", label = TRUE,
                    repel = T) + NoLegend()
cluster.plot <- SpatialDimPlot(object, group.by = "seurat_cluster", label = FALSE, pt.size.factor = 1.2) +
    theme(legend.position = "right")

png(here("plots","Mouse_VisiumHD", "UMAP_seurat_clusters.png"), width = 10, height = 10, units = "in", res = 300)
dim.plot
dev.off()

png(here("plots","Mouse_VisiumHD", "Spatial_seurat_clusters.png"), width = 10, height = 10, units = "in", res = 300)
cluster.plot
dev.off()

# save the object
saveRDS(object, here("processed-data","mouse_brain","mouse_brain_seurat.rds"))

