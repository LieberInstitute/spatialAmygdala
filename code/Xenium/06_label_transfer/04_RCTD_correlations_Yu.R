suppressPackageStartupMessages({
    library(here)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(SpatialFeatureExperiment)
    library(spacexr)
    library(ComplexHeatmap)
    library(circlize)
})

plot_dir <- here("plots", "Xenium", "06_label_transfer")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# Load reference and build ident -> broad class mapping
# =========================================================================
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce <- sce.amy

ident_to_celltype <- data.frame(
    ident = as.character(sce$ident),
    celltype = as.character(sce$celltype),
    stringsAsFactors = FALSE
) %>%
    distinct()

rm(sce.amy, sce)

## Strip "Human_" prefix
ident_to_celltype$ident <- gsub("^Human_", "", ident_to_celltype$ident)

## Build broad class: Excitatory / CGE / MGE / LGE / Non-neuronal
ident_to_celltype$RCTD_broad <- ident_to_celltype$celltype
ident_to_celltype$RCTD_broad[
    grepl("SST|PVALB", ident_to_celltype$ident)
] <- "MGE"
ident_to_celltype$RCTD_broad[
    grepl("LAMP5|VIP|CCK|CALCR", ident_to_celltype$ident)
] <- "CGE"
ident_to_celltype$RCTD_broad[
    ident_to_celltype$celltype == "InN" &
        !grepl("SST|PVALB|LAMP5|VIP|CCK|LHX8", ident_to_celltype$ident)
] <- "LGE"
ident_to_celltype$RCTD_broad[ident_to_celltype$celltype == "ExN"] <- "Excitatory"
ident_to_celltype$RCTD_broad[
    !ident_to_celltype$RCTD_broad %in% c("Excitatory", "CGE", "MGE", "LGE")
] <- "Non-neuronal"

# =========================================================================
# STEP 1: Build Xenium enrichment matrix (domain x fine cell type, z-scored)
# =========================================================================
spe_xen <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering",
                     "Banksy_domains_v1.0.rds"))
colnames(spe_xen) <- make.unique(colnames(spe_xen))

sample_ids <- c("Br9017", "Br9192", "Br9206", "Br9280")
rctd_xen_list <- lapply(sample_ids, function(s) {
    readRDS(here("processed-data", "Xenium", "06_label_transfer",
                  paste0("rctd_results_", s, ".rds")))
})
names(rctd_xen_list) <- sample_ids

## Match RCTD results to spe cells per sample
spe_xen$rctd_ident <- NA_character_
for (s in sample_ids) {
    spe_idx <- which(spe_xen$brnum == s)
    res_df <- rctd_xen_list[[s]]@results$results_df
    rctd_barcodes <- rownames(res_df)
    original_barcodes <- colnames(spe_xen)[spe_idx]
    m <- match(original_barcodes, rctd_barcodes)
    matched <- !is.na(m)
    message(s, ": matched ", sum(matched), " / ", length(spe_idx), " cells")
    spe_xen$rctd_ident[spe_idx[matched]] <- as.character(res_df$first_type[m[matched]])
}
spe_xen$rctd_ident <- gsub("^Human_", "", spe_xen$rctd_ident)

## Build enrichment matrix
has_rctd <- !is.na(spe_xen$rctd_ident)
df_xen <- data.frame(
    domain = as.character(spe_xen$Banksy_domains[has_rctd]),
    celltype = spe_xen$rctd_ident[has_rctd],
    stringsAsFactors = FALSE
)
obs_xen <- table(df_xen$domain, df_xen$celltype)
exp_xen <- outer(rowSums(obs_xen), colSums(obs_xen)) / sum(obs_xen)
log2_xen <- log2((obs_xen + 1) / (exp_xen + 1))
xen_mat <- base::scale(as.matrix(log2_xen))   # z-score per cell type

# =========================================================================
# STEP 2: Build Visium enrichment matrix from RCTD weights
# =========================================================================
spe_vis <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                         "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
rctd_vis <- readRDS(here("processed-data", "Visium", "09_deconvolution",
                          "rctd_results_human.rds"))

## Extract normalized weights
results_vis <- rctd_vis@results
weights_vis <- lapply(results_vis, function(x) x$all_weights)
weights_df <- as.data.frame(do.call(rbind, weights_vis))
norm_weights <- normalize_weights(weights_df)
coords_vis <- rctd_vis@spatialRNA@coords

spe_vis <- spe_vis[, colnames(spe_vis) %in% rownames(coords_vis)]
spe_vis <- spe_vis[, match(rownames(coords_vis), colnames(spe_vis))]

## Strip Human_ prefix from Visium cell type names to match Xenium
ws_vis <- as.matrix(norm_weights)
colnames(ws_vis) <- gsub("^Human_", "", colnames(ws_vis))

## Drop CHAT domain (as in your Visium script)
keep_vis <- !grepl("CHAT", spe_vis$BS_k16_Semisupervised_wAI)
spe_vis <- spe_vis[, keep_vis]
ws_vis <- ws_vis[keep_vis, , drop = FALSE]

## Sum weights per Visium domain, normalize by spot count, z-score per cell type
vis_domains <- as.character(spe_vis$BS_k16_Semisupervised_wAI)
unique_vis <- sort(unique(vis_domains))

vis_sum <- matrix(0, nrow = length(unique_vis), ncol = ncol(ws_vis),
                   dimnames = list(unique_vis, colnames(ws_vis)))
for (d in unique_vis) {
    vis_sum[d, ] <- colSums(ws_vis[vis_domains == d, , drop = FALSE], na.rm = TRUE)
}
vis_counts <- as.numeric(table(vis_domains)[rownames(vis_sum)])
vis_mat_raw <- sweep(vis_sum, 1, vis_counts, FUN = "/")
vis_mat <- base::scale(vis_mat_raw)   # z-score per cell type

# =========================================================================
# STEP 3: Harmonize domains — collapse Visium domains to Xenium naming
# =========================================================================
domain_map <- list(
    "BM_CoA"      = c("BM", "CoA"),
    "PL"          = "PL",
    "BL"          = c("BL", "BLD"),
    "LA"          = "LA",
    "MeA_AI"      = c("MeA", "AI"),
    "CeA"         = "CeA",
    "Endo"        = "Endothelial",
    "WM"          = c("WM.1", "WM.2")
)

vis_collapsed <- matrix(NA_real_,
                         nrow = length(domain_map), ncol = ncol(vis_mat),
                         dimnames = list(names(domain_map), colnames(vis_mat)))
for (xen_name in names(domain_map)) {
    vis_subs <- intersect(domain_map[[xen_name]], rownames(vis_mat))
    if (length(vis_subs) == 0) next
    vis_collapsed[xen_name, ] <- colMeans(vis_mat[vis_subs, , drop = FALSE],
                                           na.rm = TRUE)
}
vis_collapsed <- vis_collapsed[!apply(is.na(vis_collapsed), 1, all), , drop = FALSE]

# =========================================================================
# STEP 4: Align cell types and domains
# =========================================================================
shared_celltypes <- intersect(colnames(xen_mat), colnames(vis_collapsed))
shared_domains <- intersect(rownames(xen_mat), rownames(vis_collapsed))

message("Shared cell types: ", length(shared_celltypes))
message("Shared domains: ", length(shared_domains))

xen_aligned <- xen_mat[shared_domains, shared_celltypes, drop = FALSE]
vis_aligned <- vis_collapsed[shared_domains, shared_celltypes, drop = FALSE]

# =========================================================================
# STEP 5: Scatter plot — Visium z-score vs Xenium z-score
# =========================================================================
broad_map <- ident_to_celltype$RCTD_broad[
    match(shared_celltypes, ident_to_celltype$ident)
]

scatter_df <- expand.grid(
    domain = shared_domains,
    celltype = shared_celltypes,
    stringsAsFactors = FALSE
) %>%
    mutate(
        visium = as.numeric(vis_aligned[cbind(domain, celltype)]),
        xenium = as.numeric(xen_aligned[cbind(domain, celltype)]),
        broad = broad_map[match(celltype, shared_celltypes)]
    ) %>%
    filter(is.finite(visium), is.finite(xenium))

cor_val <- cor(scatter_df$visium, scatter_df$xenium, method = "pearson")
cor_spearman <- cor(scatter_df$visium, scatter_df$xenium, method = "spearman")

broad_colors <- c(
    "CGE" = "#E15759", "MGE" = "#4E79A7", "LGE" = "#59A14F",
    "Excitatory" = "#F28E2B", "Non-neuronal" = "#B07AA1"
)

p_scatter <- ggplot(scatter_df, aes(x = visium, y = xenium, color = broad)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = broad_colors, name = "Broad class") +
    labs(
        x = "Visium enrichment (z-score)",
        y = "Xenium enrichment (z-score)",
        title = "Cross-platform consistency: cell-type enrichment by domain",
        subtitle = sprintf("Pearson r = %.3f   |   Spearman ρ = %.3f   |   n = %d (domain × cell type)",
                           cor_val, cor_spearman, nrow(scatter_df))
    ) +
    theme_bw(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank()
    )

# =========================================================================
# STEP 6: Domain correlation heatmap (Visium domain x Xenium domain)
# =========================================================================
shared_ct_full <- intersect(colnames(xen_mat), colnames(vis_mat))
xen_full <- xen_mat[, shared_ct_full, drop = FALSE]
vis_full <- vis_mat[, shared_ct_full, drop = FALSE]

cor_mat <- matrix(NA_real_,
                   nrow = nrow(vis_full), ncol = nrow(xen_full),
                   dimnames = list(rownames(vis_full), rownames(xen_full)))
for (v in rownames(vis_full)) {
    for (x in rownames(xen_full)) {
        cor_mat[v, x] <- cor(vis_full[v, ], xen_full[x, ],
                              method = "pearson", use = "complete.obs")
    }
}

xen_order <- c("LA", "BL", "PL", "EC", "Sub", "Ctx",
                   "MeA_AI", "BM_CoA", "CeA", "Endo", "Ventricle", "WM")
vis_order <- c("CoA", "BM", "PL", "BLD", "BL", "LA",
                "HPC", "CLA", "MeA", "CeA", "AI", "Endothelial", "WM.1", "WM.2")

xen_order_present <- intersect(xen_order, colnames(cor_mat))
vis_order_present <- intersect(vis_order, rownames(cor_mat))
cor_mat <- cor_mat[vis_order_present, xen_order_present, drop = FALSE]

heat_col <- colorRamp2(
    breaks = c(-1, 0, 1),
    colors = c("#2166AC", "white", "#B2182B")
)

ht_cor <- Heatmap(
    cor_mat,
    name = "Pearson r",
    col = heat_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = "Visium domain",
    column_title = "Xenium domain",
    row_title_side = "left",
    column_title_side = "bottom",
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 45,
    heatmap_legend_param = list(title = "Pearson r"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        v <- cor_mat[i, j]
        if (is.finite(v)) {
            grid::grid.text(sprintf("%.2f", v), x = x, y = y,
                             gp = grid::gpar(fontsize = 7,
                                              col = ifelse(abs(v) > 0.5, "white", "black")))
        }
    }
)

# =========================================================================
# STEP 7: Save outputs
# =========================================================================
ggsave(here(plot_dir, "xenium_visium_scatter.pdf"), p_scatter,
        width = 9, height = 7)

pdf(here(plot_dir, "xenium_visium_domain_correlation.pdf"), width = 9, height = 8)
draw(ht_cor, heatmap_legend_side = "right")
dev.off()

## Combined panel
pdf(here(plot_dir, "xenium_visium_consistency_combined.pdf"), width = 16, height = 8)
grid::grid.newpage()
pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(1, 1))))
pushViewport(grid::viewport(layout.pos.col = 1))
print(p_scatter, newpage = FALSE)
popViewport()
pushViewport(grid::viewport(layout.pos.col = 2))
draw(ht_cor, newpage = FALSE, heatmap_legend_side = "right")
popViewport(2)
dev.off()

message("Plots saved to: ", plot_dir)
message("Overall Pearson correlation: ", round(cor_val, 3))
message("Overall Spearman correlation: ", round(cor_spearman, 3))