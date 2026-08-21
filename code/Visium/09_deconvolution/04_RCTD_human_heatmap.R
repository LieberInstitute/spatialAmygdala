suppressPackageStartupMessages({
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(spacexr)
})

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load spe + RCTD
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
myRCTD <- readRDS(here(processed_dir, "rctd_results_human.rds"))

# load snRNA-seq reference (for broad class annotations)
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))

# ---------------------------
# Extract weights from new RCTD output
# ---------------------------
results = myRCTD@results
weights = lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
norm_weights <- normalize_weights(weights_df)
coords = myRCTD@spatialRNA@coords

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights_df)

# ---------------------------
# Drop ADARB2 from sce, CHAT domain from spe
# ---------------------------
sce <- sce[, !grepl("ADARB2", sce$fine_celltype)]
spe <- spe[, !grepl("CHAT", spe$BS_k16_Semisupervised_wAI)]

sce$fine_celltype <- factor(sce$fine_celltype)
sce$broad_celltype <- factor(sce$broad_celltype)
spe$BS_k16_Semisupervised_wAI <- factor(spe$BS_k16_Semisupervised_wAI)

# ---------------------------
# Create CGE/MGE/LGE broad classes
# ---------------------------
sce$RCTD_broad_celltype <- as.character(sce$broad_celltype)
sce$RCTD_broad_celltype[grepl("SST|PVALB", sce$fine_celltype)] <- "MGE"
sce$RCTD_broad_celltype[grepl("LAMP5|VIP|CCK|CALCR", sce$fine_celltype)] <- "CGE"
sce$RCTD_broad_celltype[
  grepl("Inhibitory", sce$broad_celltype) &
    !grepl("SST|PVALB|LAMP5|VIP|CCK|CALCR", sce$fine_celltype)
] <- "LGE"
sce$RCTD_broad_celltype <- factor(sce$RCTD_broad_celltype)

# ---------------------------
# Colors
# ---------------------------
anno_colors <- list(
  BroadClass = c(
    "CGE" = "#E15759",
    "MGE" = "#4E79A7",
    "LGE" = "#59A14F",
    "Excitatory" = "#F28E2B",
    "Non-neuronal" = "#B07AA1"
  )
)

pal <- c(
  AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6",
  PL = "#f1e438ff", BL = "#035185ff", LA = "#F4B400",
  CoA = "#5DA5DA", CeA = "#197d43ff", MeA = "#baf739ff",
  HPC = "#d6a8f8ff", CHAT = "#A0522D", Endothelial = "#444444",
  WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)




# ---------------------------
# Average weights per domain, then z-score
# ---------------------------
domains <- spe$BS_k16_Semisupervised_wAI
unique_domains <- sort(unique(domains))

ws_sum <- matrix(0, nrow = length(unique_domains), ncol = ncol(ws_mat))
rownames(ws_sum) <- unique_domains
colnames(ws_sum) <- colnames(ws_mat)

for (dom in unique_domains) {
  idx <- which(domains == dom)
  ws_sum[dom, ] <- colSums(ws_mat[idx, , drop = FALSE], na.rm = TRUE)
}

domain_spot_counts <- as.numeric(table(domains)[rownames(ws_sum)])
ws_sum_norm <- sweep(ws_sum, 1, domain_spot_counts, FUN = "/")
ws_sum_scaled <- base::scale(ws_sum_norm)

# ---------------------------
# Map fine cell types -> broad classes
# ---------------------------
map_df <- unique(data.frame(
  fine  = as.character(sce$fine_celltype),
  broad = as.character(sce$RCTD_broad_celltype),
  stringsAsFactors = FALSE
))

rctd_celltypes <- colnames(ws_sum_scaled)
broad_labels <- map_df$broad[match(rctd_celltypes, map_df$fine)]

keep_cols <- !is.na(broad_labels)
if (!any(keep_cols)) stop("No cell types mapped. Check naming consistency.")

ws_sum_scaled <- ws_sum_scaled[, keep_cols, drop = FALSE]
rctd_celltypes <- colnames(ws_sum_scaled)
broad_labels <- broad_labels[keep_cols]

annotation_col <- data.frame(BroadClass = broad_labels, stringsAsFactors = FALSE)
rownames(annotation_col) <- rctd_celltypes

# Reorder columns by broad class
desired_order <- c("Excitatory", "CGE", "MGE", "LGE", "Non-neuronal")
ordered_celltypes <- unlist(lapply(desired_order, function(cls) {
  sort(rownames(annotation_col)[annotation_col$BroadClass == cls])
}))

ws_sum_scaled_ordered <- ws_sum_scaled[, ordered_celltypes, drop = FALSE]
annotation_col_ordered <- annotation_col[ordered_celltypes, , drop = FALSE]

# Reorder rows (domains)
domain_order <- c("CoA", "BM", "PL", "BLD", "BL", "LA",
                   "HPC", "CLA", "MeA", "CeA", "AI", "Endothelial", "WM.1", "WM.2")
domain_order_present <- intersect(domain_order, rownames(ws_sum_scaled_ordered))
ws_sum_scaled_ordered <- ws_sum_scaled_ordered[domain_order_present, , drop = FALSE]

annotation_row <- data.frame(
  Domain = factor(domain_order_present, levels = domain_order_present),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- domain_order_present

# ---------------------------
# Build annotations + heatmap
# ---------------------------
col_ha <- HeatmapAnnotation(
  BroadClass = factor(annotation_col_ordered$BroadClass, levels = desired_order),
  which = "column",
  col = list(BroadClass = anno_colors$BroadClass),
  show_annotation_name = FALSE
)

row_ha <- rowAnnotation(
  Domain = annotation_row$Domain,
  col = list(Domain = pal),
  show_annotation_name = FALSE
)

rng <- range(ws_sum_scaled_ordered, na.rm = TRUE, finite = TRUE)
heatmap_colors <- colorRamp2(
  breaks = c(rng[1], 0, rng[2]),
  colors = c("purple", "white", "darkgreen")
)

ht <- Heatmap(
  ws_sum_scaled_ordered,
  name = "Z-score",
  col = heatmap_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  bottom_annotation = col_ha,
  right_annotation = row_ha,
  column_split = factor(annotation_col_ordered$BroadClass, levels = desired_order),
  heatmap_legend_param = list(title = "Z-score", legend_direction = "vertical"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- ws_sum_scaled_ordered[i, j]
    if (is.finite(v) && v > 2) {
      grid::grid.text("X", x = x, y = y,
                       gp = grid::gpar(fontsize = 10, fontface = "bold", col = "white"))
    }
  }
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
pdf(here(plots_dir, "rctd_heatmap_celltypes.pdf"), width = 12, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
