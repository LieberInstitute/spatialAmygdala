## =============================================================================
## CeA-domain neurons: CeM marker co-expression at single-cell resolution
## VisiumHD cell-level data (spe_cells), Br9280
##   1. subset to CeA spatial_domain
##   2. subset to neurons (drop glia/endothelial/microglia)
##   3. binarize CeM marker expression, plot co-expression heatmap + spatial map
## =============================================================================

library(here)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

processed_dir <- here("processed-data", "VisiumHD", "06_composition")
plot_dir      <- here("plots", "Visium", "16_CeA_analyses")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- parameters -------------------------------------------------------------
# CeM marker panel grounded in rodent + iScience (McCullough 2022) literature.
# Rodent Tac2 -> human TAC3 (neurokinin B); included both in case panel naming.
cem_markers <- c(
    "ISL1", "PRKCD", "PENK", "TAC3", "SST", "CRH", "CARTPT", "NTS"
)

# positivity threshold for binarization (logcounts)
pos_thresh <- 0  # any detected expression > 0; raise to be stricter

## ---- 1. load and subset to CeA --------------------------------------------
spe_cells <- readRDS(file.path(here("processed-data", "VisiumHD",
                                    "06_composition",
                                    "spe_cells_with_domain.rds")))

# subset to Br9280 + CeA spatial domain
spe_cea <- spe_cells[, !is.na(spe_cells$spatial_domain) &
                       spe_cells$spatial_domain == "CeA" &
                       spe_cells$sample_id == "Br9280_CeA"]

cat("Cells in CeA domain (Br9280):", ncol(spe_cea), "\n")

## ---- 2. subset to neurons --------------------------------------------------
# drop glia / endothelial / microglia by name prefix
# non_neuron_re <- "^(Astro|Oligo|OPC|Micro|Endo)"
# is_neuron <- !grepl(non_neuron_re, as.character(spe_cea$first_type))
# spe_neu   <- spe_cea[, is_neuron]
# spe_neu$first_type <- droplevels(factor(as.character(spe_neu$first_type)))

# cat("CeA neurons:", ncol(spe_neu), "\n")
# cat("Cell types present:\n"); print(sort(table(spe_neu$first_type), decreasing = TRUE))


## ---- 2. subset to PRKCD cells ----------------------------------------------
# is the "PRKCD" cell-type annotation one population, or is it lumping together
# multiple molecular types? Run the co-expression analysis within PRKCD-only.
is_prkcd  <- as.character(spe_cea$first_type) == "PRKCD"
spe_neu   <- spe_cea[, is_prkcd]
spe_neu$first_type <- droplevels(factor(as.character(spe_neu$first_type)))

cat("PRKCD-annotated CeA cells:", ncol(spe_neu), "\n")

## ---- 3. logcounts ----------------------------------------------------------
if (!"logcounts" %in% assayNames(spe_neu)) spe_neu <- scuttle::logNormCounts(spe_neu)

## ---- 4. marker expression matrix (genes x cells) ---------------------------
present <- cem_markers %in% rownames(spe_neu)
if (any(!present))
    warning("CeM markers not in rownames and skipped: ",
            paste(cem_markers[!present], collapse = ", "))
markers <- cem_markers[present]

mat <- as.matrix(logcounts(spe_neu)[markers, , drop = FALSE])    # genes x cells

# binarize: cell is positive for marker if logcounts > threshold
pos <- mat > pos_thresh                                          # logical, same shape

cat("\nFraction of CeA neurons positive per marker:\n")
print(round(rowMeans(pos), 3))

## ---- 5. Spearman correlation heatmap (continuous logcounts) ---------------
# symmetric: rank-based, handles zero-inflation better than Pearson on counts
spear <- cor(t(mat), method = "spearman")    # mat is genes x cells, transpose for cor()
diag(spear) <- NA                            # diagonal is uninformative; let color scale focus on off-diagonal
 
# tight color range centered at 0: sparse single-cell data rarely exceeds ±0.4
col_fun <- colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "white", "#B2182B"))
 
ht <- Heatmap(
    spear,
    name            = "Spearman",
    col             = col_fun,
    na_col          = "black",
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    show_row_names  = TRUE,
    show_column_names = TRUE,
    row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
    column_names_gp = gpar(fontsize = 9, fontface = "italic"),
    column_names_rot = 45,
    rect_gp         = gpar(col = "black", lwd = 0.4),
    heatmap_legend_param = list(at = c(-0.3, 0, 0.3))
)
 
pdf(file.path(plot_dir, "Br9280_CeA_neurons_CeM_spearman.pdf"),
    width = 4, height = 3)
draw(ht)
dev.off()

## ---- 6. PRKCD x PENK 4-category spatial map -------------------------------
# the headline pairwise comparison: are PRKCD+ and PENK+ the same cells?
if (all(c("PRKCD", "PENK") %in% markers)) {
    cat <- ifelse( pos["PRKCD", ] &  pos["PENK", ], "PRKCD+ PENK+",
            ifelse( pos["PRKCD", ] & !pos["PENK", ], "PRKCD+ only",
            ifelse(!pos["PRKCD", ] &  pos["PENK", ], "PENK+ only",
                                                     "double-neg")))
    cat <- factor(cat, levels = c("double-neg", "PRKCD+ only",
                                  "PENK+ only", "PRKCD+ PENK+"))

    cat_cols <- c("double-neg"   = "grey85",
                  "PRKCD+ only"  = "#009E73",
                  "PENK+ only"   = "#D55E00",
                  "PRKCD+ PENK+" = "#7B3294")

    coords <- spatialCoords(spe_neu)
    df <- data.frame(x = coords[, 1], y = coords[, 2], cat = cat)
    df <- df[order(df$cat), ]   # double-neg drawn first, positives on top

    # composition summary
    cat("\nPRKCD x PENK composition (CeA neurons):\n")
    print(round(prop.table(table(cat)), 3))

    p_xy <- ggplot(df, aes(x, y, color = cat)) +
        geom_point(size = 0.6, alpha = 0.9) +
        scale_color_manual(values = cat_cols, name = "Population") +
        coord_fixed() +
        theme_void(base_size = 11) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggtitle("Br9280 CeA neurons: PRKCD vs PENK")

    # composition bar inset
    bar_df <- as.data.frame(prop.table(table(cat)))
    colnames(bar_df) <- c("cat", "prop")
    p_bar <- ggplot(bar_df, aes(x = "", y = prop, fill = cat)) +
        geom_col(width = 0.5) +
        scale_fill_manual(values = cat_cols, guide = "none") +
        coord_flip() +
        theme_void() +
        labs(x = NULL, y = NULL)

    pdf(file.path(plot_dir, "Br9280_CeA_PRKCD_vs_PENK.pdf"),
        width = 7, height = 6)
    print(p_xy)
    print(p_bar)
    dev.off()
}

## ---- 7. save annotated object ----------------------------------------------
# attach positivity matrix back to colData for downstream
for (g in markers) spe_neu[[paste0(g, "_pos")]] <- pos[g, ]
saveRDS(spe_neu, file.path(processed_dir, "spe_cells_Br9280_CeA_neurons.rds"))

cat("\nDone.\n",
    "  coexpression heatmap : ", file.path(plot_dir, "Br9280_CeA_neurons_CeM_coexpression.pdf"), "\n",
    "  PRKCD x PENK map     : ", file.path(plot_dir, "Br9280_CeA_PRKCD_vs_PENK.pdf"), "\n",
    "  annotated cells      : ", file.path(processed_dir, "spe_cells_Br9280_CeA_neurons.rds"), "\n",
    sep = "")