library(here)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(scater)
library(harmony)
library(Matrix)

set.seed(12345)

# ===== Load =====
sce.q   <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_sce.rds"))
obj.ref <- readRDS(here("processed-data","snRNAseq","GSE195445_Human_obj.rds"))
sce.ref <- as.SingleCellExperiment(obj.ref)

# ===== Subset BOTH to inhibitory neurons =====
# Yu uses InN (its own naming)
ref.neu <- sce.ref[, sce.ref$celltype %in% c("InN")]
# Girgenti uses INH
q.neu   <- sce.q[,   sce.q$celltype  %in% c("INH")]

message("Yu neurons: ", ncol(ref.neu), " | Girgenti neurons: ", ncol(q.neu))

# ===== Shared genes =====
shared <- intersect(rownames(ref.neu), rownames(q.neu))
message("Shared genes: ", length(shared))

# ===== Build Seurat objects from lognorm, recompute HVGs in INH-only space =====
seu_from_lognorm <- function(mat, meta) {
  seu <- CreateSeuratObject(counts = mat, meta.data = meta)
  seu <- SetAssayData(seu, layer = "data", new.data = mat)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE) |>
         ScaleData(verbose = FALSE) |>
         RunPCA(npcs = 30, verbose = FALSE)
  seu
}

yu.ref <- seu_from_lognorm(logcounts(ref.neu)[shared, ],
                           as.data.frame(colData(ref.neu)))
girg.q <- seu_from_lognorm(assay(q.neu, "logcounts")[shared, ],
                           as.data.frame(colData(q.neu)))

# ===== Re-anchor + transfer in INH-only subspace =====
anchors <- FindTransferAnchors(
  reference = yu.ref,
  query     = girg.q,
  reduction = "rpca",
  dims      = 1:30,
  features  = intersect(VariableFeatures(yu.ref), VariableFeatures(girg.q))
)

pred <- TransferData(anchors, refdata = yu.ref$orig_anno, dims = 1:30)
colnames(pred) <- paste0("yu_inh_", colnames(pred))

# Attach refined labels to the INH SCE
colData(q.neu) <- cbind(colData(q.neu), DataFrame(pred[colnames(q.neu), ]))

# ===== Recompute dim reduction on INH-only (for visualization) =====
dec <- modelGeneVar(q.neu, block = q.neu$Sample.x)
hvg <- getTopHVGs(dec, n = 2000)
q.neu <- runPCA(q.neu, subset_row = hvg, ncomponents = 30,
                exprs_values = "logcounts")
q.neu <- RunHarmony(q.neu, group.by.vars = "Sample.x",
                    reduction = "PCA", reduction.save = "HARMONY")
q.neu <- runUMAP(q.neu, dimred = "HARMONY", n_dimred = 30, name = "UMAP_inh_0.2", min_dist=0.2)

saveRDS(q.neu, here("processed-data","snRNAseq","Lee_at_al","girgenti_INH.srds"))


# ===== UMAP =====
library(ggplot2); library(scattermore)
umap <- reducedDim(q.neu, "UMAP_inh_0.2")
yu_filt  <- gsub("^Human_", "", q.neu$yu_inh_predicted.id)

df <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2],
                 subtype = q.neu$subtype, yu_fine = yu_filt)
plot_umap <- function(col, title)
  ggplot(df, aes(UMAP1, UMAP2, color = .data[[col]])) +
    geom_scattermore(pointsize = 2) +
    labs(title = title, color = NULL) + theme_bw() +
    guides(color = guide_legend(override.aes = list(size = 3)))

out_plots <- here("plots","snRNAseq","girgenti")
pdf(file.path(out_plots, "umaps_INH_retransfer_0.2_2.pdf"), width = 9, height = 7)
print(plot_umap("subtype", "Girgenti INH subtypes"))
print(plot_umap("yu_fine", "Yu transferred labels (INH-only, >2%)"))
dev.off()








# Publication UMAP theme: clean, no grid/axes ticks, prominent legend
theme_umap <- function(base_size = 14) {
  theme_void(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 3,
                                     hjust = 0, margin = margin(b = 8)),
      legend.title    = element_text(face = "bold", size = base_size),
      legend.text     = element_text(size = base_size - 1),
      legend.key.size = unit(0.9, "lines"),
      legend.position = "right",
      plot.margin     = margin(10, 10, 10, 10)
    )
}
 
# Small UMAP axis arrows in the corner (replaces full axes in void theme)
umap_axis_arrows <- function(df, frac = 0.18) {
  x0 <- min(df$UMAP1); y0 <- min(df$UMAP2)
  xr <- diff(range(df$UMAP1)); yr <- diff(range(df$UMAP2))
  list(
    annotate("segment", x = x0, y = y0, xend = x0 + xr*frac, yend = y0,
             arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 2),
    annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + yr*frac,
             arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 2),
    annotate("text", x = x0 + xr*frac/2, y = y0 - yr*0.03,
             label = "UMAP1", size = 5, vjust = 1),
    annotate("text", x = x0 - xr*0.03, y = y0 + yr*frac/2,
             label = "UMAP2", size = 5, angle = 90, vjust = 0)
  )
}

cat_palette <- function(labels) {
  labs <- sort(unique(as.character(labels)))
  n <- length(labs)
  # 24-color set with maximal perceptual separation (Polychrome 'kelly'/'alphabet' style)
  big <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFD300",
           "#A65628","#F781BF","#1B9E77","#D95F02","#7570B3","#E7298A",
           "#66A61E","#E6AB02","#A6761D","#666666","#1F78B4","#33A02C",
           "#FB9A99","#FDBF6F","#CAB2D6","#6A3D9A","#B15928","#00CED1")
  if (n > length(big)) {
    pal <- grDevices::colorRampPalette(big)(n)        # interpolate if very many
  } else {
    pal <- big[seq_len(n)]
  }
  # keep "other" grey if present
  setNames(ifelse(labs == "other", "#CCCCCC", pal), labs)
}
 
plot_umap <- function(df, col, title, pal = NULL) {
  if (is.null(pal)) pal <- cat_palette(df[[col]])
  p <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[col]])) +
    geom_scattermore(pointsize = 5, pixels = c(2048, 2048)) +
    umap_axis_arrows(df) +
    labs(title = title, color = NULL) +
    theme_umap() +
    scale_color_manual(values = pal) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))
  p
}
 
# ---- usage ----
out_plots <- here::here("plots","snRNAseq","girgenti")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
 
# Save as PDF (vector frame + raster points) AND high-dpi pdf option
pdf(file.path(out_plots, "umaps_INH_retransfer_pretty.pdf"), width = 9, height = 7)
print(plot_umap(df, "subtype", "Lee INH subtypes"))
print(plot_umap(df, "yu_fine", "Yu transferred labels"))
dev.off()
 
# pdf at high dpi if you prefer raster for slides
ggsave(file.path(out_plots, "umap_subtype_pretty.pdf"),
       plot_umap(df, "subtype", "Lee INH subtypes"),
       width = 9, height = 7, dpi = 600)