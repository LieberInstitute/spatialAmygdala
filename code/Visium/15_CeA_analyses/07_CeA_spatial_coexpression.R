## =============================================================================
## Br9280 CeA single-cell figure panels
##   A. spatial domains, CeA highlighted
##   B. RCTD first_type, CeA subtypes highlighted
##   C. two-channel co-expression composite (red Pop1, green Pop2, yellow = both)
##   D. individual marker gene panels (vis_gene-style)
##
## Pop1 (red):   SST, TAC3, CRH       (CeL-overlap neuropeptide group)
## Pop2 (green): PENK, PRKCD, CARTPT  (CeM-PENK/CARTPT group; PRKCD as protein-level companion)
## =============================================================================

library(here)
library(SpatialExperiment)
library(scuttle)
library(ggplot2)
library(scales)
library(ggspavis)

processed_dir <- here("processed-data", "VisiumHD", "06_composition")
plot_dir      <- here("plots", "Visium", "16_CeA_analyses")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- parameters -------------------------------------------------------------
focal_sample <- "Br9280_CeA"

pop1_markers <- c("SST", "TAC3", "CRH")
pop2_markers <- c("PENK", "PRKCD", "CARTPT")
all_markers  <- unique(c(pop1_markers, pop2_markers))

# CeA RCTD subtypes (raw names, not anatomical divisions)
cea_cols <- c("DRD2 ISL1" = "#D55E00",
              "DRD2 PAX6" = "#0072B2",
              "PRKCD"     = "#009E73")

# spatial domain palette
domain_colors <- c(
    AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6",
    LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
    MeA = "#baf739ff", CHAT = "#A0522D",
    Ventricle = "#666666", Vascular = "#333333",
    WM = "#BBBBBB", Neuropil = "#DDDDDD"
)

## ---- 1. load and subset to Br9280 + CeA domain -----------------------------
spe_cells <- readRDS(file.path(processed_dir, "spe_cells_with_domain.rds"))
spe.subset <- spe_cells[, spe_cells$sample_id == focal_sample &
                           !is.na(spe_cells$spatial_domain)]

# rotate -90 degrees for plotting (matches earlier panels)
coords <- spatialCoords(spe.subset)
spatialCoords(spe.subset)[, 1] <- coords[, 2]
spatialCoords(spe.subset)[, 2] <- -coords[, 1]

## ---- Panel A: spatial domains, CeA opaque ----------------------------------
spe_A <- spe.subset
spe_A$is_cea_domain <- factor(
    as.character(spe_A$spatial_domain) == "CeA",
    levels = c(FALSE, TRUE)
)
spe_A <- spe_A[, order(spe_A$is_cea_domain)]

p_domain <- ggspavis::plotCoords(
    spe_A, annotate = "spatial_domain", in_tissue = NULL
) +
    #aes(alpha = is_cea_domain) +
    scale_color_manual(values = domain_colors, drop = FALSE, breaks = "CeA") +
    ggtitle("A. Spatial domains")

pdf(file.path(plot_dir, "Br9280_panelA_domains.pdf"), width = 4, height = 4)
print(p_domain)
dev.off()

## ---- Panel B: RCTD first_type, CeA subtypes opaque -------------------------
spe_B <- spe.subset
other_types <- setdiff(unique(as.character(spe_B$first_type)), names(cea_cols))
first_type_cols <- c(setNames(rep("grey70", length(other_types)), other_types),
                     cea_cols)
spe_B$first_type <- factor(as.character(spe_B$first_type),
                           levels = c(other_types, names(cea_cols)))
spe_B$is_cea <- factor(
    as.character(spe_B$first_type) %in% names(cea_cols),
    levels = c(FALSE, TRUE)
)
spe_B <- spe_B[, order(spe_B$is_cea)]

p_celltype <- ggspavis::plotCoords(
    spe_B, annotate = "first_type", in_tissue = NULL
) +
    aes(alpha = is_cea) +
    scale_color_manual(values = first_type_cols, drop = FALSE,
                       breaks = names(cea_cols)) +
    scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), guide = "none") +
    ggtitle("B. RCTD CeA subtypes")

pdf(file.path(plot_dir, "Br9280_panelB_first_type.pdf"), width = 4, height = 4)
print(p_celltype)
dev.off()

## ---- Panel C: two-channel population composite (CeA neurons only) ----------
# subset further: CeA domain + neurons (drop glia/endo/microglia)
spe_C <- spe.subset[, spe.subset$spatial_domain == "CeA"]
non_neuron_re <- "^(Astro|Oligo|OPC|Micro|Endo)"
spe_C <- spe_C[, !grepl(non_neuron_re, as.character(spe_C$first_type))]
if (!"logcounts" %in% assayNames(spe_C)) spe_C <- logNormCounts(spe_C)
 
present <- intersect(all_markers, rownames(spe_C))
miss    <- setdiff(all_markers, rownames(spe_C))
if (length(miss)) warning("Markers not found and skipped: ",
                          paste(miss, collapse = ", "))
 
# population score = sum of logcounts across each population's markers
p1 <- intersect(pop1_markers, present)
p2 <- intersect(pop2_markers, present)
stopifnot(length(p1) > 0, length(p2) > 0)
 
lc <- as.matrix(logcounts(spe_C))
score1 <- colSums(lc[p1, , drop = FALSE])    # Pop1: SST/TAC3/CRH
score2 <- colSums(lc[p2, , drop = FALSE])    # Pop2: PENK/PRKCD/CARTPT
 
# scale each score to [0,1] by its 95th percentile of nonzero values
norm01 <- function(x) {
    cap <- quantile(x[x > 0], 0.95, na.rm = TRUE)
    if (!is.finite(cap) || cap <= 0) return(rep(0, length(x)))
    pmin(x / cap, 1)
}
r <- norm01(score1)
g <- norm01(score2)
 
# floor: any expressing cell gets at least 0.3 intensity so it pops against black
floor_intensity <- function(x, floor = 0.3) {
    out <- numeric(length(x))
    nz  <- x > 0
    out[nz] <- floor + (1 - floor) * x[nz]
    out
}
r_f <- floor_intensity(r)
g_f <- floor_intensity(g)
 
# magenta + cyan = white (colorblind-safe two-channel mix):
#   Pop1 only -> magenta (R + B), Pop2 only -> cyan (G + B), both -> white
df <- data.frame(
    x   = spatialCoords(spe_C)[, 1],
    y   = -spatialCoords(spe_C)[, 2],
    r   = r_f,
    g   = g_f,
    any = r > 0 | g > 0
)
df$col <- rgb(
    red   = pmin(df$r, 1),                   # magenta contributes red
    green = pmin(df$g, 1),                   # cyan contributes green
    blue  = pmin(df$r + df$g, 1)             # both channels contribute blue
)
df$col[!df$any] <- "grey25"                  # non-expressing -> near-black on dark bg
 
# draw non-expressing cells first, expressing cells on top
df <- df[order(df$any), ]
 
# legend swatches: magenta / cyan / white
legend_df <- data.frame(
    x    = NA_real_, y = NA_real_,
    grp  = factor(c("SST/TAC3/CRH", "PENK/PRKCD/CARTPT", "co-expressing"),
                  levels = c("SST/TAC3/CRH", "PENK/PRKCD/CARTPT", "co-expressing"))
)
legend_cols <- c("SST/TAC3/CRH"      = rgb(1, 0, 1),    # magenta
                 "PENK/PRKCD/CARTPT" = rgb(0, 1, 1),    # cyan
                 "co-expressing"     = rgb(1, 1, 1))    # white
 
p_pop <- ggplot(df, aes(x, y)) +
    geom_point(color = df$col, size = 0.6) +
    geom_point(data = legend_df, aes(x, y, fill = grp), shape = 21,
               size = 4, colour = NA, na.rm = TRUE) +
    scale_fill_manual(values = legend_cols, name = "Population", drop = FALSE) +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(plot.background   = element_rect(fill = "black", color = NA),
          panel.background  = element_rect(fill = "black", color = NA),
          legend.background = element_rect(fill = "black", color = NA),
          legend.key        = element_rect(fill = "black", color = NA),
          legend.text       = element_text(color = "white"),
          legend.title      = element_text(color = "white"),
          plot.title        = element_text(color = "white")) +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    ggtitle("C. CeA neuron populations")
 
pdf(file.path(plot_dir, "Br9280_panelC_populations.pdf"), width = 5, height = 4)
print(p_pop)
dev.off()
 
# composition summary print
pop_class <- ifelse(r > 0 & g > 0, "double+",
              ifelse(r > 0,        "Pop1 (SST/TAC3/CRH)",
              ifelse(g > 0,        "Pop2 (PENK/PRKCD/CARTPT)",
                                   "double-")))
cat("\nPanel C population composition (CeA neurons):\n")
print(round(prop.table(table(pop_class)), 3))

### ---- Panel D: individual marker gene panels -------------------------------
# one panel per marker, expression on the CeA-neuron subset, viridis rocket
coords_C <- spatialCoords(spe_C)

pdf(file.path(plot_dir, "Br9280_panelD_marker_genes.pdf"),
    width = 4, height = 4)
for (g in present) {
    expr <- as.numeric(logcounts(spe_C)[g, ])
    df_g <- data.frame(x = coords_C[, 1], y = -coords_C[, 2], expr = expr)
    df_g <- df_g[order(df_g$expr), ]   # high expression drawn on top

    p_g <- ggplot(df_g, aes(x, y, color = expr)) +
        geom_point(size = 0.6, alpha = 0.9) +
        scale_color_gradientn(
            colours = viridisLite::rocket(10, direction = -1),
            name = "logcounts") +
        coord_fixed() +
        theme_void(base_size = 11) +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
              plot.title   = element_text(face = "italic")) +
        ggtitle(g)
    print(p_g)
}
dev.off()