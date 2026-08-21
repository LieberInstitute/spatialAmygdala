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

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution","Allen_ITCs")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

dir <- "your_directory"
for (s in unique(spe$sample_id)) magick::image_write(magick::image_read(imgRaster(spe, sample_id = s, image_id = "lowres")), file.path(dir, paste0(s, "_lowres.tif")), format = "tiff")

# load
res <- readRDS(here(processed_dir, "rctd_Allen_ITCs_results.rds"))
colnames(colData(res))
# [1] "cell_type_list" "conf_list"      "conv_all"       "conv_sub"      
# [5] "min_score"      "sample_id" 

# plot first_type predicted from RCTD


# # plot
# pdf(here(plots_dir, "rctd_first_type.pdf"), width=15, height=10)
#     plotCoords(spe, annotate="first_type", sample_id="sample_id", point_size=0.05) 
# dev.off()


ws <- assay(res, "weights")
table(colSums(ws) == 0)

# add proportion estimates as metadata
ws <- data.frame(t(as.matrix(ws)))
colData(spe)[names(ws)] <- ws[colnames(spe), ]

celltypes <- colnames(ws)

pdf(here(plots_dir, "Allen_ITC_rctd_celltype_proportions.pdf"), width=15, height=10)
    plotCoords(spe, sample_id="sample_id", annotate="Oligodendrocyte", point_size=0.1) +
    scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket")))
dev.off()


# plot one cell type per page
pdf(here(plots_dir, "Allen_ITC_rctd_celltype_proportions_all.pdf"), width=15, height=10)
    for (ct in celltypes) {
        print(plotCoords(spe, sample_id="sample_id", annotate=ct, point_size=0.1) +
        scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket"))) +
        ggtitle(ct))
    }
dev.off()

type2plot <- c("DLK1_ZFHX3", "TAC1_PPP1R1B", "PENK_DRD2")
pdf(here(plots_dir, "Allen_ITC_rctd_celltype_proportions_CeA.pdf"), width=15, height=10)
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
# load snRNA-seq data
sce <- readRDS(here("processed-data", "snRNAseq", "sce_amy_macaque.rds"))
sce

# create a new broad cell type where "inhibitory" is broken into CGE, MGE, LGE
# All SST and PVALB are MGE
# LAMP5, VIP, and CCK are CGE
# inhibitory that are not in these categories are LGE
sce$fine_celltype <- factor(sce$fine_celltype)
sce$RCTD_broad_celltype <- as.character(sce$broad_celltype)

sce$RCTD_broad_celltype[grepl("SST|PVALB", sce$fine_celltype)] <- "MGE"
sce$RCTD_broad_celltype[grepl("LAMP5|VIP|CCK|CALCR", sce$fine_celltype)] <- "CGE"
sce$RCTD_broad_celltype[
  grepl("Inhibitory", sce$broad_celltype) & 
    !grepl("SST|PVALB|LAMP5|VIP|CCK|CALCR", sce$fine_celltype)
] <- "LGE"

sce$RCTD_broad_celltype <- factor(sce$RCTD_broad_celltype)
table(sce$RCTD_broad_celltype)

# ==== Define colors for plotting ====
anno_colors <- list(
  BroadClass = c(
    "CGE" = "#E15759",         # reddish
    "MGE" = "#4E79A7",         # blue
    "LGE" = "#59A14F",         # green
    "Excitatory" = "#F28E2B",  # orange
    "Non-neuronal" = "#B07AA1" # purple
  )
)

# === Load libraries ===
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(scCustomize)
library(here)

# ==== Summing & Normalizing cell type weights by spatial domain ====

# Align barcodes
common_barcodes <- intersect(rownames(ws), colnames(spe))
length(common_barcodes)  # Should be ~223,420

spe_subset <- spe[, common_barcodes]

# Ensure matrix format and consistent spot names
ws_mat <- as.matrix(ws)
stopifnot(all(rownames(ws_mat) == colnames(spe_subset)))

# Get domain labels for each spot
domains <- spe_subset$BS_k16_relabel_ITC_smooth
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

# === Normalize by number of spots per domain ===
domain_spot_counts <- table(domains)[rownames(ws_sum)]
ws_sum_norm <- sweep(ws_sum, 1, domain_spot_counts, FUN = "/")  # per-spot average

# Z-score by cell type (across domains)
ws_sum_scaled <- base::scale(ws_sum_norm)

# ==== Make column annotations ====

# Extract cell types in the same order as columns of the heatmap
rctd_celltypes <- colnames(ws_sum_scaled)

# Create a named vector mapping fine → broad
fine_to_broad <- setNames(as.character(sce$RCTD_broad_celltype), sce$fine_celltype)

# Map RCTD cell types to broad classes (e.g., CGE, MGE, etc.)
broad_labels <- unname(fine_to_broad[rctd_celltypes])

# Create annotation data frame
annotation_col <- data.frame(BroadClass = broad_labels)
rownames(annotation_col) <- rctd_celltypes

# ==== Reorder columns by broad class AND alphabetically within ====

desired_order <- c("Excitatory", "CGE", "MGE", "LGE", "Non-neuronal")

# For each broad class, get the fine celltypes sorted alphabetically
ordered_celltypes <- unlist(
  lapply(desired_order, function(cls) {
    fine_cts <- rownames(annotation_col)[annotation_col$BroadClass == cls]
    sort(fine_cts)
  })
)

# Reorder matrix and annotations (columns)
ws_sum_scaled_ordered <- ws_sum_scaled[, ordered_celltypes]
annotation_col_ordered <- annotation_col[ordered_celltypes, , drop = FALSE]

# ==== Reorder and relabel rows (spatial domains) ====

# Define desired domain order
domain_order <- c("aBA", "CoA", "BA", "BLVM.1", "BLVM.2", "LA",
                  "HPC", "ITC", "CeA", "MeA", "Other", "Meninges", "WM.1", "WM.2")

# Reorder rows of matrix
ws_sum_scaled_ordered <- ws_sum_scaled_ordered[domain_order, , drop = FALSE]

# ==== Row Annotations for Domains ====

annotation_row <- data.frame(Domain = factor(domain_order, levels = domain_order))
rownames(annotation_row) <- domain_order

# Color palette for domains
pal <- scCustomize::DiscretePalette_scCustomize(num_colors = length(domain_order), palette = "polychrome")

# Add to anno_colors list
anno_colors$Domain <- setNames(pal, domain_order)

# === Build ComplexHeatmap annotations ===

# Column annotation (BroadClass, bottom, no label)
col_ha <- HeatmapAnnotation(
  BroadClass = annotation_col_ordered$BroadClass,
  which = "column",
  col = list(BroadClass = anno_colors$BroadClass),
  show_annotation_name = FALSE
)

# Row annotation (Domain, right, no label)
row_ha <- rowAnnotation(
  Domain = annotation_row$Domain,
  col = list(Domain = anno_colors$Domain),
  show_annotation_name = FALSE
)

# ==== Define heatmap colors (symmetric, red/blue, white at 0) ====

heatmap_colors <- colorRamp2(
  breaks = c(min(ws_sum_scaled_ordered), 0, max(ws_sum_scaled_ordered)),
  colors = c("purple", "white", "darkgreen")
)

# reorder BroadClass factor levels for split order
annotation_col_ordered$BroadClass <- factor(
  annotation_col_ordered$BroadClass,
  levels = c("Excitatory", "CGE", "MGE", "LGE", "Non-neuronal")
)

# ==== Create heatmap ====
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
  column_split = annotation_col_ordered$BroadClass,
  heatmap_legend_param = list(
    title = "Z-score",
    legend_direction = "vertical"
  )
)

# ==== Save to PDF ====
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

domains <- spe_sub$BS_k16_relabel_ITC_smooth
samples <- spe_sub$sample_id

domains_of_interest <- c("LA", "BA", "aBA")
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

df_avg$Domain <- factor(df_avg$Domain, levels = c("LA", "BA", "aBA"))
domain_colors <- c("LA" = "blue", "BA" = "red", "aBA" = "purple")

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

domains <- spe_sub$BS_k16_relabel_ITC_smooth
samples <- spe_sub$sample_id

domains_of_interest <- c("LA", "BA", "aBA")
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

df_avg$Domain <- factor(df_avg$Domain, levels = c("LA", "BA", "aBA"))
domain_colors <- c("LA" = "blue", "BA" = "red", "aBA" = "purple")

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
