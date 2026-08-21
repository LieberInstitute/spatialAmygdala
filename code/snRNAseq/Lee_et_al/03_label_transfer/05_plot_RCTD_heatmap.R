suppressPackageStartupMessages({
    library(here)
    library(ComplexHeatmap)
    library(circlize)
    library(spacexr)
    library(dplyr)
})

plots_dir     <- here("plots", "snRNAseq", "girgenti")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# load spe + RCTD
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
myRCTD <- readRDS(here(processed_dir, "rctd_girgenti_subtype_subtype.rds"))

# ---------------------------
# Extract weights from RCTD output
# ---------------------------
results     <- myRCTD@results
weights     <- lapply(results, function(x) x$all_weights)
weights_df  <- data.frame(do.call(rbind, weights))
norm_weights <- normalize_weights(weights_df)
coords <- myRCTD@spatialRNA@coords

spe <- spe[, colnames(spe) %in% rownames(coords)]
spe <- spe[, match(rownames(coords), colnames(spe))]

ws_mat <- as.matrix(norm_weights)

# ---------------------------
# Average weights per domain, then z-score
# ---------------------------
domains        <- as.character(spe$BS_k16_Semisupervised_wAI)
unique_domains <- sort(unique(domains))

ws_sum <- matrix(0, nrow = length(unique_domains), ncol = ncol(ws_mat),
                 dimnames = list(unique_domains, colnames(ws_mat)))
for (dom in unique_domains) {
    idx <- which(domains == dom)
    ws_sum[dom, ] <- colSums(ws_mat[idx, , drop = FALSE], na.rm = TRUE)
}

domain_spot_counts <- as.numeric(table(domains)[rownames(ws_sum)])
ws_sum_norm   <- sweep(ws_sum, 1, domain_spot_counts, FUN = "/")
ws_sum_scaled <- base::scale(ws_sum_norm)

# ---------------------------
# Reorder rows (domains); keep all subtype columns
# ---------------------------
domain_order <- c("LA","BLD","BL","PL","HPC","CLA","MeA","BM","CoA","CHAT",
                  "CeA","AI","Endothelial","WM.1","WM.2")
domain_order_present <- intersect(domain_order, rownames(ws_sum_scaled))
ws_sum_scaled_ordered <- ws_sum_scaled[domain_order_present, , drop = FALSE]

annotation_row <- data.frame(
    Domain = factor(domain_order_present, levels = domain_order_present)
)
rownames(annotation_row) <- domain_order_present

# ---------------------------
# Colors
# ---------------------------
pal <- c(
    AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6",
    PL = "#f1e438ff", BL = "#035185ff", LA = "#F4B400",
    CoA = "#5DA5DA", CeA = "#197d43ff", MeA = "#baf739ff",
    HPC = "#d6a8f8ff", Endothelial = "#444444",
    WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4", CHAT = "black"
)

row_ha <- rowAnnotation(
    Domain = annotation_row$Domain,
    col = list(Domain = pal),
    show_annotation_name = FALSE
)

rng <- range(ws_sum_scaled_ordered, na.rm = TRUE, finite = TRUE)
heatmap_colors <- colorRamp2(c(rng[1], 0, rng[2]), c("purple", "white", "darkgreen"))

ht <- Heatmap(
    ws_sum_scaled_ordered,
    name            = "Z-score",
    col             = heatmap_colors,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,        
    show_column_names = TRUE,
    show_row_names    = TRUE,
    right_annotation  = row_ha,
    column_names_gp   = grid::gpar(fontsize = 9),
    row_names_gp      = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Z-score", legend_direction = "vertical"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        v <- ws_sum_scaled_ordered[i, j]
        if (is.finite(v) && v > 2)
            grid::grid.text("X", x, y,
                gp = grid::gpar(fontsize = 10, fontface = "bold", col = "white"))
    }
)

pdf(here(plots_dir, "rctd_heatmap_subtype_Lee.pdf"), width = 6, height = 4)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Saved: ", here(plots_dir, "rctd_heatmap_subtype_Lee.pdf"))