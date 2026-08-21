suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(patchwork)
  library(ggplot2)
  library(scater)
  library(pheatmap)
})

output_dir <- here("plots", "Xenium", "08_gex_correlations")
processed_dir_out <- here("processed-data", "Xenium", "08_gex_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir_out, recursive = TRUE, showWarnings = FALSE)

# ===== Load data =====
spe.xenium <- readRDS(file.path("processed-data", "Xenium", "05_spatial_clustering", "Banksy_domains_v1.0.rds"))

spe.visium <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

# ===== Subset both to common genes (Xenium panel ~366 genes) =====
common_genes <- intersect(rownames(spe.xenium), rownames(spe.visium))
message(paste("Number of common genes:", length(common_genes)))

spe.visium <- spe.visium[common_genes, ]
spe.xenium <- spe.xenium[common_genes, ]

# set rownames to gene_id
rownames(spe.visium) <- rowData(spe.visium)$gene_id
rownames(spe.xenium) <- rowData(spe.xenium)$ID

# ===== Set up registration variables =====
spe.xenium$layer <- spe.xenium$Banksy_domains
spe.visium$layer <- spe.visium$BS_k16_Semisupervised_wAI

# ===== Pseudobulk: Visium =====
spe_pseudo_vis <- registration_pseudobulk(
  spe.visium,
  var_registration = "layer",
  var_sample_id = "sample_id",
  min_ncells = 10
)

# ===== Pseudobulk: Xenium =====
spe_pseudo_xen <- registration_pseudobulk(
  spe.xenium,
  var_registration = "layer",
  var_sample_id = "brnum",
  min_ncells = 10
)

# ===== Registration model + stats: Visium (reference) =====
registration_mod_vis <- registration_model(spe_pseudo_vis, covars = NULL)
block_cor_vis <- registration_block_cor(spe_pseudo_vis, registration_model = registration_mod_vis)

results_enrichment_vis <- registration_stats_enrichment(
  spe_pseudo_vis,
  block_cor = block_cor_vis,
  covars = NULL,
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

# ===== Registration model + stats: Xenium (query) =====
registration_mod_xen <- registration_model(spe_pseudo_xen, covars = NULL)
block_cor_xen <- registration_block_cor(spe_pseudo_xen, registration_model = registration_mod_xen)

results_enrichment_xen <- registration_stats_enrichment(
  spe_pseudo_xen,
  block_cor = block_cor_xen,
  covars = NULL,
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

# ===== Build Visium modeling_results (reference) =====
results_pairwise_vis <- registration_stats_pairwise(
  spe_pseudo_vis,
  registration_model = registration_mod_vis,
  block_cor = block_cor_vis,
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

results_anova_vis <- registration_stats_anova(
  spe_pseudo_vis,
  block_cor = block_cor_vis,
  covars = NULL,
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

modeling_results_vis <- list(
  "enrichment" = results_enrichment_vis,
  "pairwise" = results_pairwise_vis,
  "anova" = results_anova_vis
)

modeling_results_xen <- list(
  "enrichment" = results_enrichment_xen,
  "pairwise" = registration_stats_pairwise(
    spe_pseudo_xen,
    registration_model = registration_mod_xen,
    block_cor = block_cor_xen,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  ),
  "anova" = registration_stats_anova(
    spe_pseudo_xen,
    block_cor = block_cor_xen,
    covars = NULL,
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
  )
)

# ===== Build Xenium query t-stat matrix for layer_stat_cor() =====
# layer_stat_cor() expects: rows = Ensembl gene IDs, cols = domain labels, values = t-stats
t_cols_xen <- grep("^t_", colnames(results_enrichment_xen), value = TRUE)
xenium_tstats <- as.matrix(results_enrichment_xen[, t_cols_xen])

# ===== Cross-platform correlation using layer_stat_cor() =====
cor_stats <- layer_stat_cor(
  stats = xenium_tstats,
  modeling_results = modeling_results_vis,
  model_type = "enrichment"
)

# ===== Plot with layer_stat_cor_plot() (clustered) =====
pdf(here(output_dir, "cross_platform_correlation_heatmap.pdf"), width = 10, height = 8)
layer_stat_cor_plot(cor_stats)
dev.off()

# ===== Plot with manual ordering (unclustered) =====
# Reorder rows (Xenium) and columns (Visium)
xenium_order <- c("Ctx", "EC", "Sub",  "LA", "PL", "BL", "BM_CoA", "MeA_AI", "CeA", "Endo", "Ventricle", "WM")
visium_order <- c("CLA", "HPC", "BLD", "LA", "PL", "BL", "BM", "CoA", "AI", "MeA", "CeA", "WM.2", "WM.1", "Endothelial")

# Subset to names that exist in cor_stats
xenium_order <- xenium_order[xenium_order %in% rownames(cor_stats)]
visium_order <- visium_order[visium_order %in% colnames(cor_stats)]

cor_stats_ordered <- cor_stats[xenium_order, visium_order, drop = FALSE]

pdf(here(output_dir, "cross_platform_correlation_heatmap_ordered.pdf"), width = 10, height = 8)
layer_stat_cor_plot(cor_stats_ordered, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


# From Xenium, drop Ctx, EC, Sub, Ventricle, Endo, WM
xenium_order <- setdiff(xenium_order, c("Ctx", "EC", "Sub", "Ventricle", "Endo", "WM"))

# visium drop CLA, HPC, BLD, WM.1, WM.2, Endothelial
visium_order <- setdiff(visium_order, c("CLA", "HPC", "BLD", "WM.1", "WM.2", "Endothelial"))

# Subset to names that exist in cor_stats
xenium_order <- xenium_order[xenium_order %in% rownames(cor_stats)]
visium_order <- visium_order[visium_order %in% colnames(cor_stats)]

cor_stats_ordered <- cor_stats[xenium_order, visium_order, drop = FALSE]

pdf(here(output_dir, "cross_platform_correlation_heatmap_ordered_dropped.pdf"), width = 3, height = 2.2)
layer_stat_cor_plot(cor_stats_ordered, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()


# ===== Scatterplots: LA and PL t-stat correlations =====
library(ggrepel)

# Map Ensembl IDs to gene symbols
# Xenium: rowData has $Symbol, Visium: rowData has $gene_name
xen_id2sym <- setNames(rowData(spe_pseudo_xen)$Symbol, rowData(spe_pseudo_xen)$gene_id)
vis_id2sym <- setNames(rowData(spe_pseudo_vis)$gene_name, rowData(spe_pseudo_vis)$gene_id)

# Build data frames with gene symbols
df_xen <- data.frame(
  ensembl = rownames(results_enrichment_xen),
  gene = xen_id2sym,
  t_LA_xen = results_enrichment_xen$t_stat_LA,
  t_PL_xen = results_enrichment_xen$t_stat_PL
)

df_vis <- data.frame(
  ensembl = rownames(results_enrichment_vis),
  gene = vis_id2sym,
  t_LA_vis = results_enrichment_vis$t_stat_LA,
  t_PL_vis = results_enrichment_vis$t_stat_PL
)

# Merge by gene symbol
df_merged <- merge(df_xen[, c("gene", "t_LA_xen", "t_PL_xen")],
                   df_vis[, c("gene", "t_LA_vis", "t_PL_vis")],
                   by = "gene")
df_merged <- df_merged[complete.cases(df_merged), ]
message(paste("Genes for scatterplot:", nrow(df_merged)))

# Helper
plot_domain_scatter <- function(df, t_xen_col, t_vis_col, domain_label) {
  r_val <- cor(df[[t_xen_col]], df[[t_vis_col]])

  df$avg_t <- (df[[t_xen_col]] + df[[t_vis_col]]) / 2
  top5 <- head(df[order(-df$avg_t), ], 5)

  ggplot(df, aes(x = .data[[t_vis_col]], y = .data[[t_xen_col]])) +
    geom_point(alpha = 0.6, size = 2, color = "grey40") +
    geom_point(data = top5, color = "red", size = 3) +
    geom_text_repel(data = top5, aes(label = gene),
                    color = "red", size = 4, max.overlaps = 20,
                    fontface = "italic") +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("R = ", round(r_val, 3)),
             hjust = 1.1, vjust = -0.5, size = 5, fontface = "bold") +
    labs(
      x = paste0("Visium ", domain_label, " t-statistic"),
      y = paste0("Xenium ", domain_label, " t-statistic"),
      title = paste0(domain_label, " enrichment t-statistics")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 13)
    )
}

p_PL <- plot_domain_scatter(df_merged, "t_PL_xen", "t_PL_vis", "PL")
p_LA <- plot_domain_scatter(df_merged, "t_LA_xen", "t_LA_vis", "LA")

# PL on top, LA on bottom
pdf(here(output_dir, "LA_PL_tstat_scatterplots.pdf"), width = 3, height = 6)
print(p_PL / p_LA)
dev.off()
