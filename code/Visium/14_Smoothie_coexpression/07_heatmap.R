library(here)
library(SpatialExperiment)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(matrixStats)

# ==============================================================================
# Configuration
# ==============================================================================
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

modules_df <- read.csv(here("processed-data", "Visium", "14_smoothie_coexpression",
                             "results", "modules_df_0.8_9.csv"))

plots_dir <- here("plots", "Visium", "14_smoothie_coexpression")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Cluster column — UPDATE THIS
clust_col <- "BS_k16_Semisupervised_wAI"  # <-- UPDATE

# Remove NA domains
keep_spots <- !is.na(spe[[clust_col]])
spe <- spe[, keep_spots]
domains <- factor(spe[[clust_col]])

# ==============================================================================
# Same modules and domain assignments as network plot
# ==============================================================================
# Keep all modules with >= 10 genes, except 12 and 29
drop_modules <- c(12, 29) # 12 = sex linked, 29 = low quality
MIN_GENES <- 10

keep_modules <- modules_df %>%
    group_by(module_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n >= MIN_GENES, !module_label %in% drop_modules) %>%
    pull(module_label)

cat("Keeping", length(keep_modules), "modules\n")

# Domain assignments for the 10 highlighted modules (others get "Other")
domain_lookup <- c(
    "7"  = "CLA",
    "5"  = "PL",
    "10" = "CoA",
    "16" = "HPC",
    "27" = "BLD",
    "6"  = "LA",
    "13" = "CeA",
    "19" = "MeA",
    "20" = "AI",
    "15" = "CHAT"
)

# Same domain color palette as network plot
pal <- c(
    AI          = "#D62728",
    BM          = "#E67E22",
    BLD         = "#9B59B6",
    PL          = "#f1e438ff",
    BL          = "#035185ff",
    LA          = "#F4B400",
    CoA         = "#5DA5DA",
    CeA         = "#197d43ff",
    MeA         = "#baf739ff",
    HPC         = "#d6a8f8ff",
    CHAT        = "#A0522D",
    Endothelial = "#444444",
    WM.1        = "#BBBBBB",
    WM.2        = "#DDDDDD",
    CLA         = "#FF69B4"
)

# ==============================================================================
# Filter to kept modules
# ==============================================================================
modules_filtered <- modules_df %>%
    filter(module_label %in% keep_modules)

mod_sizes <- modules_filtered %>%
    group_by(module_label) %>%
    summarise(n_genes = n(), .groups = "drop")

cat("Genes per module:\n")
print(mod_sizes)

# ==============================================================================
# Compute module scores (Smoothie method) using filtered genes
# ==============================================================================
lc <- logcounts(spe)
all_genes <- rownames(spe)

module_labels <- paste0("M", mod_sizes$module_label)

# Track which modules are domain-assigned (for bold formatting)
is_domain_module <- as.character(mod_sizes$module_label) %in% names(domain_lookup)

module_scores <- matrix(0, nrow = ncol(spe), ncol = nrow(mod_sizes))
colnames(module_scores) <- module_labels

for (i in seq_len(nrow(mod_sizes))) {
    mod <- mod_sizes$module_label[i]
    mod_genes <- modules_filtered$name[modules_filtered$module_label == mod]
    present <- mod_genes[mod_genes %in% all_genes]

    if (length(present) == 0) next

    mat <- as.matrix(lc[present, , drop = FALSE])
    row_maxs <- matrixStats::rowMaxs(mat)
    row_maxs[row_maxs == 0] <- 1
    module_scores[, i] <- colSums(mat / row_maxs)
}

# ==============================================================================
# Aggregate: mean module score per domain
# ==============================================================================
domain_module_mat <- matrix(NA, nrow = nlevels(domains), ncol = ncol(module_scores))
rownames(domain_module_mat) <- levels(domains)
colnames(domain_module_mat) <- colnames(module_scores)

for (d in levels(domains)) {
    idx <- which(domains == d)
    domain_module_mat[d, ] <- colMeans(module_scores[idx, , drop = FALSE], na.rm = TRUE)
}

# Z-score each module (column) across domains
scaled_mat <- scale(domain_module_mat)

# ==============================================================================
# Transpose: domains = columns, modules = rows
# ==============================================================================
scaled_mat <- t(scaled_mat)

# ==============================================================================
# Annotations
# ==============================================================================
# Domain colors for column annotation
domain_cols <- pal[colnames(scaled_mat)]
# Handle any domains not in pal
domain_cols[is.na(domain_cols)] <- "grey50"

col_ha <- HeatmapAnnotation(
    Domain = anno_simple(
        colnames(scaled_mat),
        col = setNames(pal[colnames(scaled_mat)], colnames(scaled_mat)),
        na_col = "grey50"
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
)
# Module domain color for row annotation (grey for unassigned)
mod_domain_cols <- sapply(as.character(mod_sizes$module_label), function(m) {
    if (m %in% names(domain_lookup)) pal[domain_lookup[m]] else "grey70"
})
names(mod_domain_cols) <- module_labels


# ==============================================================================
# Color scale
# ==============================================================================
col_fun <- colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")
)

# ==============================================================================
# Plot
# ==============================================================================
pdf(file.path(plots_dir, "module_domain_heatmap_selected.pdf"),
    width = 4.5,
    height= 7)

ht <- Heatmap(
    scaled_mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10,
                         fontface = ifelse(is_domain_module[match(rownames(scaled_mat), module_labels)],
                                           "bold", "plain")),
    column_names_gp = gpar(fontsize = 10),
    bottom_annotation = col_ha,
    column_title = "Spatial Domains",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 12),
    row_title = "Co-expression Modules",
    row_title_gp = gpar(fontsize = 12),
    heatmap_legend_param = list(
        title = "Scaled\nmodule\nscore",
        legend_height = unit(4, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
        val <- scaled_mat[i, j]
        if (!is.na(val) && val > 1.5) {
            grid.text("*", x, y, gp = gpar(fontsize = 14))
        }
    }
)

draw(ht, merge_legend = TRUE)
dev.off()

cat("Saved: module_domain_heatmap_selected.pdf\n")