library(here)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

set.seed(12345)

# ===== Load =====
rat.amy <- readRDS(here("processed-data","snRNAseq","GSE212415_seurat.rds"))
sce.yu  <- readRDS(here("processed-data","snRNAseq","Yu_et_al","sce_amy_mouse.rds"))
sce.yu$orig_anno <- gsub("Mouse_", "", sce.yu$orig_anno)
sce.yu$space_anno <- gsub("Mouse_", "", sce.yu$space_anno)
# Sanity check the label columns
table(sce.yu$orig_anno)
table(sce.yu$space_anno)

# ===== Shared gene symbols =====
shared <- intersect(rownames(sce.yu), rownames(rat.amy))
message("Shared symbols: ", length(shared), " / Yu: ", nrow(sce.yu),
        " / Zhou: ", nrow(rat.amy))

# ===== Yu reference (Seurat) =====
yu.ref <- CreateSeuratObject(
  counts    = counts(sce.yu)[shared, ],
  meta.data = as.data.frame(colData(sce.yu))
)
yu.ref <- NormalizeData(yu.ref, verbose = FALSE) |>
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = 50, verbose = FALSE)

# ===== Zhou query — use RNA assay, restrict to shared features =====
DefaultAssay(rat.amy) <- "RNA"
rat.query <- subset(rat.amy, features = shared)
rat.query <- NormalizeData(rat.query, verbose = FALSE) |>
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = 50, verbose = FALSE)

# ===== Anchors (rPCA) =====
anchors <- FindTransferAnchors(
  reference = yu.ref,
  query     = rat.query,
  reduction = "rpca",
  dims      = 1:30,
  features  = intersect(VariableFeatures(yu.ref), VariableFeatures(rat.query))
)

# ===== Transfer both label sets =====
pred_cell <- TransferData(anchorset = anchors,
                          refdata = yu.ref$orig_anno,  dims = 1:30)
pred_space <- TransferData(anchorset = anchors,
                           refdata = yu.ref$space_anno, dims = 1:30)
colnames(pred_cell)  <- paste0("yu_cell_",  colnames(pred_cell))
colnames(pred_space) <- paste0("yu_space_", colnames(pred_space))

rat.amy <- AddMetaData(rat.amy, pred_cell[colnames(rat.amy), ])
rat.amy <- AddMetaData(rat.amy, pred_space[colnames(rat.amy), ])

saveRDS(rat.amy,
        here("processed-data","snRNAseq","zhou_with_yu_labels.rds"))

# ===== Quick QC =====
hist(rat.amy$yu_cell_prediction.score.max,
     breaks = 50, main = "Yu cell-type transfer confidence",
     xlab = "max prediction score")
table(rat.amy$yu_cell_predicted.id, rat.amy$cellType_broad)  # vs. Zhou's labels




library(ggplot2)
library(patchwork)
# ===== UMAP plotting =====
# Use Zhou's existing UMAP if present, otherwise compute one
if (!"umap" %in% Reductions(rat.amy)) {
  rat.amy <- FindVariableFeatures(rat.amy, nfeatures = 3000, verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 50, verbose = FALSE) |>
    RunUMAP(dims = 1:30, verbose = FALSE)
}
 
# Yu cell-type labels
p_cell <- DimPlot(
  rat.amy,
  group.by   = "yu_cell_predicted.id",
  reduction  = "umap",
  label      = TRUE,
  repel      = TRUE,
  raster     = TRUE,
  raster.dpi = c(1024, 1024),
  pt.size    = 0.2
) +
  ggtitle("Yu cell-type labels (transferred)") +
  theme(legend.position = "right",
        legend.text = element_text(size = 7)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
 
# Yu spatial labels
p_space <- DimPlot(
  rat.amy,
  group.by   = "yu_space_predicted.id",
  reduction  = "umap",
  label      = TRUE,
  repel      = TRUE,
  raster     = TRUE,
  raster.dpi = c(1024, 1024),
  pt.size    = 0.2
) +
  ggtitle("Yu spatial labels (transferred)") +
  theme(legend.position = "right",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
 
# Transfer confidence — diagnostic
p_conf <- FeaturePlot(
  rat.amy,
  features   = "yu_cell_prediction.score.max",
  reduction  = "umap",
  raster     = TRUE,
  raster.dpi = c(1024, 1024),
  pt.size    = 0.2
) +
  scale_color_viridis_c(option = "magma", limits = c(0, 1)) +
  ggtitle("Cell-type transfer confidence")
 
# Save
out_dir <- here("plots","snRNAseq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
 
pdf(file.path(out_dir, "zhou_umap_yu_labels.pdf"),
    width = 14, height = 10)
print(p_cell)
print(p_space)
print(p_conf)
print((p_cell | p_space) / p_conf +
        plot_annotation(title = "Zhou rat amygdala — Yu labels transferred"))
dev.off()
 
message("Saved: ", file.path(out_dir, "zhou_umap_yu_labels.pdf"))
 