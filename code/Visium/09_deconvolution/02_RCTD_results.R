library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)
library(escheR)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))


# load
res <- readRDS(here(processed_dir, "rctd_results.rds"))
colnames(colData(res))
# [1] "cell_type_list" "conf_list"      "conv_all"       "conv_sub"      
# [5] "min_score"      "sample_id" 


ws <- assay(res, "weights")
table(colSums(ws) == 0)
#  FALSE   TRUE 
# 202038  21358 

# add proportion estimates as metadata
ws <- data.frame(t(as.matrix(ws)))
colData(spe)[names(ws)] <- ws[colnames(spe), ]

celltypes <- colnames(ws)

pdf(here(plots_dir, "rctd_celltype_proportions_Oligo.pdf"), width=15, height=10)
    plotCoords(spe, sample_id="sample_id", annotate="Oligodendrocyte", point_size=0.1) +
    scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket")))
dev.off()


# plot one cell type per page
pdf(here(plots_dir, "rctd_celltype_proportions_all.pdf"), width=15, height=10)
    for (ct in celltypes) {
        print(plotCoords(spe, sample_id="sample_id", annotate=ct, point_size=0.1) +
        scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket"))) +
        ggtitle(ct))
    }
dev.off()

type2plot <- c("DLK1_ZFHX3", "TAC1_PPP1R1B", "PENK_DRD2")
pdf(here(plots_dir, "rctd_celltype_proportions_CeA.pdf"), width=15, height=10)
    for (ct in celltypes) {
        print(plotCoords(spe, sample_id="sample_id", annotate=ct, point_size=0.1) +
        scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket"))) +
        ggtitle(ct))
    }
dev.off()




# =============== plotting cell type weights in spatial domains with escheR ===============



plot_feature_overlay <- function(spe, feature, pal, outdir = here("plots", "Visium", "08_marker_genes")) {
  # Check if feature exists in the object
  if (!feature %in% rownames(spe) && !feature %in% colnames(colData(spe))) {
    stop(paste("Feature", feature, "not found in assay or colData."))
  }

  # Create a working copy
  spe.sub <- spe

  # If it's a gene, compute log1p(counts)
  if (feature %in% rownames(spe.sub)) {
    spe.sub[[feature]] <- log1p(counts(spe.sub)[feature, ])
    fill_name <- "logcounts"
  } else {
    # Assume it's a cell type weight (continuous already)
    fill_name <- "weight"
  }

  # Output file
  outfile <- file.path(outdir, paste0("Spatial_", feature, ".pdf"))

  # Plot
  pdf(outfile, width = 8, height = 7)
  print(
    make_escheR(spe.sub) |>
      add_fill(var = feature, point_size = 1) |>  # continuous fill (gene or weight)
      add_ground(var = "BS_k16_relabel_ITC_smooth",
                 point_size = 0.5, stroke = 0.25, color_aes = TRUE) +
      scale_color_manual(values = pal, guide = "none") +
      scale_fill_gradient(low = "white", high = "black", name = fill_name) +
      ggtitle(bquote(italic(.(feature)))) +
      guides(
        fill = guide_colorbar(
          barwidth = 1.5,
          barheight = 15,
          title.position = "top"
        )
      ) +
      theme(
        plot.title = element_text(size = 22, margin = margin(b = 10)),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.key = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.spacing.y = unit(0.4, "lines"),
        legend.position = "right",
        legend.justification = "top",
        legend.background = element_blank(),
        plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
      )
  )
  dev.off()
}

# plot PENK, TAC1, and DLK1 neurons 
type2plot <- c("DLK1_ZFHX3", "TAC1_PPP1R1B", "PENK_DRD2")
for (ct in type2plot) plot_feature_overlay(spe, ct, pal)





# Desired sample order
sample_order <- c("Br2743","Br6423","Br6471","Br6660","Br8325","Br9192","Br9280")
spe$sample_id <- factor(as.character(spe$sample_id),
                        levels = sample_order,
                        ordered = TRUE)

# Features to plot (genes or cell type weights)
type2plot <- c("DLK1_ZFHX3", "TAC1_PPP1R1B", "PENK_DRD2")

# Output directory
outdir <- here("plots", "Visium", "08_marker_genes")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Loop over each feature (gene or cell type weight)
for (feature in type2plot) {
  plot_list <- list()

  # Loop through samples in the correct order
  for (sid in sample_order) {
    spe.sub <- spe[, spe$sample_id == sid]

    if (ncol(spe.sub) == 0) next  # skip missing samples

    # Determine if feature is a gene or weight
    if (feature %in% rownames(spe.sub)) {
      spe.sub[[feature]] <- log1p(counts(spe.sub)[feature, ])
      fill_name <- "logcounts"
    } else {
      fill_name <- "weight"
    }

    # Build plot
    p <- make_escheR(spe.sub) |>
      add_fill(var = feature, point_size = 1) |>
      add_ground(var = "BS_k16_relabel_ITC_smooth",
                 point_size = 0.5, stroke = 0.25, color_aes = TRUE) +
      scale_color_manual(values = pal, guide = "none") +
      scale_fill_gradient(low = "white", high = "black", name = fill_name) +
      ggtitle(paste0(sid, " - ", feature)) +
      guides(fill = guide_colorbar(
        barwidth = 1.5,
        barheight = 15,
        title.position = "top"
      )) +
      theme(
        plot.title = element_text(size = 18, margin = margin(b = 8)),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.position = "none",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
      )

    plot_list[[sid]] <- p
  }

  # Combine all sample plots (4 columns per row)
  combined_plot <- wrap_plots(plot_list, ncol = 4)

  # Save to PDF
  outfile <- file.path(outdir, paste0("Spatial_", feature, "_all_samples.pdf"))
  pdf(outfile, width = 8 * 4, height = 7 * ceiling(length(plot_list) / 4))
  print(combined_plot)
  dev.off()
}






# ======== Plotting cell types ========

ids <- names(ws)[apply(ws, 1, which.max)]
ids <- gsub("\\.([A-z])", " \\1", ids)
idx <- match(colnames(spe), rownames(ws))
table(spe$RCTD_decon1<- factor(ids[idx]))

# drop astrocytes (not in brain)
ws_no_astro <- ws[, !grepl("Astrocyte", colnames(ws))]
ids_no_astro <- names(ws_no_astro)[apply(ws_no_astro, 1, which.max)]
ids_no_astro <- gsub("\\.([A-z])", " \\1", ids_no_astro)
table(spe$RCTD_decon1_no_astro <- factor(ids_no_astro[idx]))

pal <- scCustomize::DiscretePalette_scCustomize(num_colors = length(unique(celltypes)), palette = "polychrome")

pdf(here(plots_dir, "rctd_celltype_max.pdf"), width=15, height=10)
    plotCoords(spe, annotate="RCTD_decon1_no_astro", sample_id="sample_id", point_size=0.05) +
    scale_color_manual(values=pal)
dev.off()














# ======================================
# Heatmap of average cell type proportions by spatial domain
# ======================================

# ==== Prepping snRNA-seq labels for heatmap ====

suppressPackageStartupMessages({
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(scCustomize)
})

# ---------------------------
# Load snRNA-seq + (assumes) spe + ws already in memory
# ---------------------------
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))

# Drop ADARB2 from sce, and CHAT domain from spe
sce <- sce[, !grepl("ADARB2", sce$fine_celltype)]
spe <- spe[, !grepl("CHAT", spe$BS_k16_Semisupervised_wAI)]

# Refactor
sce$fine_celltype <- factor(sce$fine_celltype)
sce$broad_celltype <- factor(sce$broad_celltype)
spe$BS_k16_Semisupervised_wAI <- factor(spe$BS_k16_Semisupervised_wAI)

# ---------------------------
# Create modified broad classes (CGE/MGE/LGE) for inhibitory
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
# Define colors for plotting
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


# Align barcodes
common_barcodes <- intersect(rownames(ws), colnames(spe))
message("Common barcodes: ", length(common_barcodes))

if (length(common_barcodes) == 0) {
  stop("No common barcodes between rownames(ws) and colnames(spe). Check orientation/names.")
}

spe_subset <- spe[, common_barcodes]

# Ensure matrix format and consistent spot names
ws_mat <- as.matrix(ws)

# IMPORTANT: ensure same ordering
ws_mat <- ws_mat[colnames(spe_subset), , drop = FALSE]
stopifnot(all(rownames(ws_mat) == colnames(spe_subset)))

# Get domain labels for each spot
domains <- spe_subset$BS_k16_Semisupervised_wAI
unique_domains <- sort(unique(domains))

# Initialize result matrix: rows = domains, columns = cell types
ws_sum <- matrix(0, nrow = length(unique_domains), ncol = ncol(ws_mat))
rownames(ws_sum) <- unique_domains
colnames(ws_sum) <- colnames(ws_mat)

# Sum weights per cell type for each domain
for (dom in unique_domains) {
  idx <- which(domains == dom)
  ws_sum[dom, ] <- colSums(ws_mat[idx, , drop = FALSE], na.rm = TRUE)
}

# Normalize by number of spots per domain (per-spot average)
domain_spot_counts <- as.numeric(table(domains)[rownames(ws_sum)])
ws_sum_norm <- sweep(ws_sum, 1, domain_spot_counts, FUN = "/")

# Z-score by cell type (across domains)
ws_sum_scaled <- base::scale(ws_sum_norm)

# ---------------------------
# Column annotations: map fine cell types -> broad classes
#   (DROP anything that does not map; no "Unknown")
# ---------------------------

map_df <- unique(data.frame(
  fine  = as.character(sce$fine_celltype),
  broad = as.character(sce$RCTD_broad_celltype),
  stringsAsFactors = FALSE
))

rctd_celltypes <- colnames(ws_sum_scaled)
broad_labels <- map_df$broad[match(rctd_celltypes, map_df$fine)]

# Drop unmapped celltypes entirely (this removes "Unknown" rather than showing it)
keep_cols <- !is.na(broad_labels)
if (!any(keep_cols)) {
  stop("All ws column names failed to map to sce$fine_celltype. Check naming consistency between ws and sce.")
}

ws_sum_scaled <- ws_sum_scaled[, keep_cols, drop = FALSE]
rctd_celltypes <- colnames(ws_sum_scaled)
broad_labels <- broad_labels[keep_cols]

annotation_col <- data.frame(BroadClass = broad_labels, stringsAsFactors = FALSE)
rownames(annotation_col) <- rctd_celltypes

# Reorder columns by broad class and alphabetically within
desired_order <- c("Excitatory", "CGE", "MGE", "LGE", "Non-neuronal")

ordered_celltypes <- unlist(lapply(desired_order, function(cls) {
  fine_cts <- rownames(annotation_col)[annotation_col$BroadClass == cls]
  sort(fine_cts)
}))

if (length(ordered_celltypes) == 0) {
  stop("No cell types retained after ordering. Check BroadClass labels or mapping.")
}

ws_sum_scaled_ordered <- ws_sum_scaled[, ordered_celltypes, drop = FALSE]
annotation_col_ordered <- annotation_col[ordered_celltypes, , drop = FALSE]

# ---------------------------
# Reorder and relabel rows (spatial domains) safely
# ---------------------------

domain_order <- c(
  "CoA", "BM", "PL", "BLD", "BL", "LA",
  "HPC", "CLA", "MeA", "CeA", "AI", "Endothelial", "WM.1", "WM.2"
)

domain_order_present <- intersect(domain_order, rownames(ws_sum_scaled_ordered))
if (length(domain_order_present) == 0) {
  stop("None of domain_order entries are present in rownames(ws_sum_scaled_ordered). Check domain labels.")
}

ws_sum_scaled_ordered <- ws_sum_scaled_ordered[domain_order_present, , drop = FALSE]

# Row annotations (Domain)
annotation_row <- data.frame(
  Domain = factor(domain_order_present, levels = domain_order_present),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- domain_order_present

# Domain palette (named!)
pal <- c( 
  AI        = "#D62728",  # strong red standout
  BM        = "#E67E22",  # burnt orange
  BLD      = "#9B59B6",  # violet
  PL     = "#f1e438ff",  # steel blue
  BL     = "#035185ff",  # sky blue
  LA         = "#F4B400",  # goldenrod
  CoA        = "#5DA5DA",  # blue-gray
  CeA        = "#197d43ff",  # emerald green
  MeA        = "#baf739ff",  # chartreuse
  HPC        = "#d6a8f8ff",  # deep royal purple (distinct from BA)
  CHAT       = "#A0522D",  # chestnut brown
  Endothelial   = "#444444",  # charcoal gray
  WM.1       = "#BBBBBB",  # light slate gray
  WM.2       = "#DDDDDD",   # mist gray
  CLA        = "#FF69B4"   # hot pink
)

# ---------------------------
# Build ComplexHeatmap annotations
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

# ---------------------------
# Heatmap colors (finite range!)
# ---------------------------
rng <- range(ws_sum_scaled_ordered, na.rm = TRUE, finite = TRUE)
if (!all(is.finite(rng))) {
  stop("Heatmap matrix has no finite values (all NA/Inf). Check ws_sum_scaled_ordered.")
}

heatmap_colors <- colorRamp2(
  breaks = c(rng[1], 0, rng[2]),
  colors = c("purple", "white", "darkgreen")
)

# ---------------------------
# Create heatmap
# ---------------------------
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
  heatmap_legend_param = list(
    title = "Z-score",
    legend_direction = "vertical"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- ws_sum_scaled_ordered[i, j]
    if (is.finite(v) && v > 2) {
      grid::grid.text(
        "X", x = x, y = y,
        gp = grid::gpar(fontsize = 10, fontface = "bold", col = "white")
      )
    }
  }
)


# ---------------------------
# Save to PDF
# ---------------------------
pdf(here(plots_dir, "rctd_heatmap_celltypes.pdf"), width = 12, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()










# ======= Box plots of cell type proportions by domain =======

library(dplyr)
library(ggplot2)
library(tidyr)

mge_celltypes <- rownames(annotation_col)[annotation_col$BroadClass == "MGE"]
mge_celltypes <- setdiff(mge_celltypes, "SST_TAC1")

common_barcodes <- intersect(rownames(ws), colnames(spe))

ws_mat <- as.matrix(ws[common_barcodes, , drop = FALSE])
spe_sub <- spe[, common_barcodes]

domains <- spe_sub$BS_k16_Semisupervised_wAI
samples <- spe_sub$sample_id

domains_of_interest <- c("LA", "BLD", "PL", "BM")
idx_keep <- domains %in% domains_of_interest

ws_mge <- ws_mat[idx_keep, mge_celltypes, drop = FALSE]
domains <- domains[idx_keep]
samples <- samples[idx_keep]

df <- as.data.frame(ws_mge)
df$Domain <- domains
df$Sample <- samples

df_long <- df %>%
  pivot_longer(cols = mge_celltypes, names_to = "CellType", values_to = "Weight")

df_norm <- df_long %>%
  group_by(Sample, Domain) %>%
  mutate(WeightNorm = Weight / sum(Weight, na.rm = TRUE)) %>%
  ungroup()

df_avg <- df_norm %>%
  group_by(Sample, Domain, CellType) %>%
  summarise(AvgWeight = mean(WeightNorm, na.rm = TRUE), .groups = "drop")

df_avg$Domain <- factor(df_avg$Domain, levels = c("LA", "BLD", "PL", "BM"))
domain_colors <- c("LA" = "blue", "BLD" = "red", "PL" = "green", "BM" = "purple")

# ---- 10. Plot ----
pdf(here(plots_dir, "rctd_mge_celltypes_boxplot_by_domain.pdf"), width = 8, height = 6)
ggplot(df_avg, aes(x = Domain, y = AvgWeight, fill = Domain)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = .3, size = 1, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free_y") +
  scale_fill_manual(values = domain_colors) +
  labs(
    title = "Average normalized MGE cell type weights",
    y = "Average normalized weight (per donor)",
    x = "Domain"
  )
dev.off()



# ======= CGE cell types =======

library(dplyr)
library(ggplot2)
library(tidyr)

cge_celltypes <- rownames(annotation_col)[annotation_col$BroadClass == "CGE"]

common_barcodes <- intersect(rownames(ws), colnames(spe))

ws_mat <- as.matrix(ws[common_barcodes, , drop = FALSE])
spe_sub <- spe[, common_barcodes]

domains <- spe_sub$BS_k16_Semisupervised_wAI
samples <- spe_sub$sample_id

domains_of_interest <- c("LA", "BLD", "PL", "BM")
idx_keep <- domains %in% domains_of_interest

ws_cge <- ws_mat[idx_keep, cge_celltypes, drop = FALSE]
domains <- domains[idx_keep]
samples <- samples[idx_keep]

df <- as.data.frame(ws_cge)
df$Domain <- domains
df$Sample <- samples

df_long <- df %>%
  pivot_longer(cols = cge_celltypes, names_to = "CellType", values_to = "Weight")

df_norm <- df_long %>%
  group_by(Sample, Domain) %>%
  mutate(WeightNorm = Weight / sum(Weight, na.rm = TRUE)) %>%
  ungroup()

df_avg <- df_norm %>%
  group_by(Sample, Domain, CellType) %>%
  summarise(AvgWeight = mean(WeightNorm, na.rm = TRUE), .groups = "drop")

df_avg$Domain <- factor(df_avg$Domain, levels = c("LA", "BLD", "PL", "BM"))
domain_colors <- c("LA" = "blue", "BLD" = "red", "PL" = "green", "BM" = "purple")

# ---- Plot ----
pdf(here(plots_dir, "rctd_cge_celltypes_boxplot_by_domain.pdf"), width = 8, height = 6)
ggplot(df_avg, aes(x = Domain, y = AvgWeight, fill = Domain)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = .3, size = 1, alpha = 0.8) +
  facet_wrap(~ CellType, scales = "free_y") +
  scale_fill_manual(values = domain_colors) +
  labs(
    title = "Average normalized CGE cell type weights",
    y = "Average normalized weight (per donor)",
    x = "Domain"
  )
dev.off()
