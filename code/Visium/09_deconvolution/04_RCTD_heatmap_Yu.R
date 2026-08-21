suppressPackageStartupMessages({
    library(here)
    library(ComplexHeatmap)
    library(circlize)
    library(spacexr)
    library(dplyr)
})

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# load spe + RCTD
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
myRCTD <- readRDS(here(processed_dir, "rctd_results_human.rds"))

# load snRNA-seq reference
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce <- sce.amy

# ---------------------------
# Build ident -> broad class mapping (same as Xenium script)
# ---------------------------
ident_to_celltype <- data.frame(
    ident = as.character(sce$ident),
    celltype = as.character(sce$celltype),
    stringsAsFactors = FALSE
) %>%
    distinct()

## Strip "Human_" prefix
ident_to_celltype$ident <- gsub("^Human_", "", ident_to_celltype$ident)

## Build RCTD_broad: Excitatory / CGE / MGE / LGE / Non-neuronal
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

rm(sce.amy, sce)

# ---------------------------
# Extract weights from RCTD output
# ---------------------------
results <- myRCTD@results
weights <- lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
norm_weights <- normalize_weights(weights_df)
coords <- myRCTD@spatialRNA@coords

spe <- spe[, colnames(spe) %in% rownames(coords)]
spe <- spe[, match(rownames(coords), colnames(spe))]

## Strip "Human_" prefix from weight column names
ws_mat <- as.matrix(norm_weights)
colnames(ws_mat) <- gsub("^Human_", "", colnames(ws_mat))

## Normalize: dots -> spaces to match ident_to_celltype
colnames(ws_mat) <- gsub("\\.", " ", colnames(ws_mat))

# ---------------------------
# Drop CHAT domain
# ---------------------------
# keep <- !grepl("CHAT", spe$BS_k16_Semisupervised_wAI)
# spe <- spe[, keep]
# ws_mat <- ws_mat[keep, , drop = FALSE]
# spe$BS_k16_Semisupervised_wAI <- factor(as.character(spe$BS_k16_Semisupervised_wAI))

# ---------------------------
# Average weights per domain, then z-score
# ---------------------------
domains <- as.character(spe$BS_k16_Semisupervised_wAI)
unique_domains <- sort(unique(domains))

ws_sum <- matrix(0, nrow = length(unique_domains), ncol = ncol(ws_mat),
                  dimnames = list(unique_domains, colnames(ws_mat)))

for (dom in unique_domains) {
    idx <- which(domains == dom)
    ws_sum[dom, ] <- colSums(ws_mat[idx, , drop = FALSE], na.rm = TRUE)
}

domain_spot_counts <- as.numeric(table(domains)[rownames(ws_sum)])
ws_sum_norm <- sweep(ws_sum, 1, domain_spot_counts, FUN = "/")
ws_sum_scaled <- base::scale(ws_sum_norm)

# ---------------------------
# Map cell types -> broad classes using ident_to_celltype
# ---------------------------
rctd_celltypes <- colnames(ws_sum_scaled)
broad_labels <- ident_to_celltype$RCTD_broad[
    match(rctd_celltypes, ident_to_celltype$ident)
]

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
domain_order <- c( "LA","BLD", "BL", "PL", "HPC", "CLA", "MeA", "BM", "CoA", "CHAT", 
                     "CeA", "AI",  "Endothelial", "WM.1", "WM.2")
domain_order_present <- intersect(domain_order, rownames(ws_sum_scaled_ordered))
ws_sum_scaled_ordered <- ws_sum_scaled_ordered[domain_order_present, , drop = FALSE]

annotation_row <- data.frame(
    Domain = factor(domain_order_present, levels = domain_order_present),
    stringsAsFactors = FALSE
)
rownames(annotation_row) <- domain_order_present

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
    HPC = "#d6a8f8ff", Endothelial = "#444444",
    WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4", CHAT="black"
)

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
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(title = "Z-score", legend_direction = "vertical"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        v <- ws_sum_scaled_ordered[i, j]
        if (is.finite(v) && v > 2) {
            grid::grid.text("X", x = x, y = y,
                             gp = grid::gpar(fontsize = 10, fontface = "bold", col = "white"))
        }
    }
)

pdf(here(plots_dir, "rctd_heatmap_celltypes_Yu.pdf"), width = 12, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Visium heatmap saved to: ", plots_dir)