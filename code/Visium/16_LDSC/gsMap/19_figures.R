suppressPackageStartupMessages({
    library("here")
    library("ComplexHeatmap")
    library("circlize")
    library("grid")
})

# ---------------------------
# Directories
# ---------------------------
# figdata/ sits next to 19_figures.py -- fix this if the gsMap code lives elsewhere
data_dir  <- here("code", "Visium", "16_LDSC", "gsMap", "figdata")
plots_dir <- here("processed-data", "Visium", "16_LDSC", "gsMap", "figures")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

arm    <- "test2"   # "test2" = functional-conditioned arm, "test1" = default baseline
p_file <- file.path(data_dir, paste0(arm, "_cauchy_matrix.csv"))

thr <- -log10(0.05)   # 1.301, the white point of the P scale

# ---------------------------
# Traits to keep, and their display names
# ---------------------------
# same 21 as the scDRS figure; keys are the gsMap spellings, which differ in
# places (BIP_PGC3 not BIP_2024, Alzheimer_v3 lowercase v, PTSD not PTSD_F3)
trait_labels <- c(
    GSCAN_SmkInit   = "Smoking initiation",
    Neuroticism     = "Neuroticism",
    BMI             = "BMI",
    EduYears        = "Education years",
    Anorexia        = "Anorexia",
    SCZ             = "Schizophrenia",
    ADHD            = "ADHD",
    Insomnia        = "Insomnia",
    Intelligence    = "Intelligence",
    MDD             = "Depression",
    GSCAN_DrnkWk    = "Drinks per week",
    BIP_PGC3        = "Bipolar disorder",
    Alzheimer_v3    = "Alzheimer's",
    Autism          = "Autism",
    PD              = "Parkinson's",
    PTSD            = "PTSD",
    Stroke_2022_Any = "Stroke",
    OUD_META        = "Opioid use disorder",
    Epilepsy_GGE    = "Epilepsy (GGE)",
    Height          = "Height",
    T2D             = "Type 2 diabetes"
)

# ---------------------------
# Color palette
# ---------------------------
pal <- c(
    IA          = "#D62728",
    BM          = "#E67E22",
    BLD         = "#9B59B6",
    PL          = "#f1e438",
    BL          = "#035185",
    LA          = "#F4B400",
    CoA         = "#5DA5DA",
    CeA         = "#197d43",
    MeA         = "#baf739",
    HPC         = "#d6a8f8",
    CHAT        = "#A0522D",
    Endothelial = "#444444",
    WM.1        = "#BBBBBB",
    WM.2        = "#DDDDDD",
    CLA         = "#FF69B4"
)

domain_order <- c(
    "CLA", "HPC", "BLD", "LA", "PL", "BL", "BM", "CoA",
    "IA", "CeA", "MeA", "CHAT", "Endothelial", "WM.1"
)

domain_drop <- "WM.2"

# ---------------------------
# Load Cauchy P matrix (traits x domains)
# ---------------------------
P <- read.csv(p_file, row.names = 1, check.names = FALSE)

colnames(P)[colnames(P) == "AI"] <- "IA"
P <- P[, !colnames(P) %in% domain_drop, drop = FALSE]

dim(P)
sort(rownames(P))

# match trait keys case-insensitively, so Alzheimer_v3 / Alzheimer_V3 both land
idx <- match(tolower(names(trait_labels)), tolower(rownames(P)))
if (anyNA(idx)) message("not in the cauchy matrix: ",
                        paste(names(trait_labels)[is.na(idx)], collapse = ", "))

keep_traits  <- names(trait_labels)[!is.na(idx)]
p_mat        <- as.matrix(P[idx[!is.na(idx)], , drop = FALSE])
rownames(p_mat) <- trait_labels[keep_traits]

# clip as the python does, then -log10; domains on rows
n_mat <- -log10(pmax(p_mat, 1e-300))
n_mat <- t(n_mat)

# ---------------------------
# Load Mantel-Haenszel OR (long) and mask saturated traits
# ---------------------------
or_long <- read.csv(file.path(data_dir, "domain_vs_rest_OR.csv"),
                    check.names = FALSE, stringsAsFactors = FALSE)
or_long <- or_long[or_long$arm == arm, ]
or_long$domain[or_long$domain == "AI"] <- "IA"

# a trait whose FDR call saturates (<2% or >98% of spots) has an undefined OR
# everywhere, not just in the cell that tripped it -- blank the whole trait
sat <- unique(or_long$trait[as.logical(or_long$saturated)])
if (length(sat)) message("OR undefined (saturated): ", paste(sat, collapse = ", "))
or_long$log2OR[or_long$trait %in% sat] <- NA_real_

o_mat <- tapply(or_long$log2OR, list(or_long$domain, or_long$trait), function(x) x[1])
f_mat <- tapply(or_long$fdr,    list(or_long$domain, or_long$trait), function(x) x[1])

# ---------------------------
# Shared row / column order
# ---------------------------
# columns ordered by mean -log10 P, strongest first (as in the python figure),
# and the SAME order is used for the OR panel so the two are comparable
domain_order <- intersect(domain_order, rownames(n_mat))
trait_order  <- names(sort(colMeans(n_mat, na.rm = TRUE), decreasing = TRUE))

n_mat <- n_mat[domain_order, trait_order]

# OR matrix is keyed by the raw gsMap trait names; relabel then align
colnames(o_mat) <- trait_labels[colnames(o_mat)]
colnames(f_mat) <- trait_labels[colnames(f_mat)]
o_mat <- o_mat[domain_order, trait_order]
f_mat <- f_mat[domain_order, trait_order]

stopifnot(identical(dimnames(n_mat), dimnames(o_mat)))

# FDR from the MH pooling is kept in the table but NOT plotted: pooled over ~1e5
# spots, its SE is tiny and nearly every cell clears 0.05, so the stars carried no
# information. Donor-level inference lives in gsmap_donor_meta.R.
s_mat <- ifelse(!is.na(f_mat) & f_mat < 0.05 & !is.na(o_mat), "*", "")
mean(s_mat == "*", na.rm = TRUE)

# ---------------------------
# Row annotation
# ---------------------------
dom_cols <- pal[domain_order]

row_ha <- rowAnnotation(
    Domain = factor(domain_order, levels = domain_order),
    col = list(Domain = dom_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(4, "mm")
)

# ---------------------------
# Color scales
# ---------------------------
# Both scales are capped at a percentile, not the max: one runaway cell
# (-log10 P reaches ~26 in the default arm) otherwise compresses the whole grid
# into the pale end. colorRamp2 clamps out-of-range values to the end colour, so
# nothing is lost -- the top tick is labelled with a >= to say so.
p_cap <- unname(ceiling(quantile(n_mat, 0.99, na.rm = TRUE)))
o_lim <- unname(ceiling(quantile(abs(o_mat), 0.95, na.rm = TRUE)))

sum(n_mat > p_cap, na.rm = TRUE)          # cells clamped at the top of the P scale
sum(abs(o_mat) > o_lim, na.rm = TRUE)     # cells clamped at either end of the OR scale

# P: grey below the 0.05 threshold, white AT it, orange-red above -- so the
# grey/orange boundary is significance, not the middle of the data range
p_cols <- colorRamp2(
    c(0, thr / 2, thr,
      thr + (p_cap - thr) / 3, thr + 2 * (p_cap - thr) / 3, p_cap),
    c("#D2D2D2", "#E4E4E4", "#FFFFFF", "#FDBE85", "#F16913", "#B30000")
)

o_cols <- colorRamp2(c(-o_lim, 0, o_lim), c("#762A83", "#F7F7F7", "#1B7837"))

# ticks: always mark 0 and the significance threshold, then round numbers up to the cap
p_ticks  <- c(0, thr, setdiff(pretty(c(thr, p_cap), 3), 0))
p_ticks  <- sort(unique(p_ticks[p_ticks <= p_cap]))
p_labels <- format(round(p_ticks, 1), trim = TRUE)
p_labels[p_ticks == thr] <- "1.3"
if (any(n_mat > p_cap, na.rm = TRUE)) {
    p_labels[length(p_labels)] <- paste0("\u2265", p_labels[length(p_labels)])
}

o_ticks  <- seq(-o_lim, o_lim, length.out = 5)
o_labels <- format(o_ticks, trim = TRUE)
if (any(abs(o_mat) > o_lim, na.rm = TRUE)) {
    o_labels[1] <- paste0("\u2264", o_labels[1])
    o_labels[length(o_labels)] <- paste0("\u2265", o_labels[length(o_labels)])
}

# ---------------------------
# Heatmaps
# ---------------------------
common <- list(
    left_annotation = row_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    rect_gp = gpar(col = "white", lwd = 0.5),
    width  = unit(7 * ncol(n_mat), "mm"),
    height = unit(7 * nrow(n_mat), "mm")
)

ht_p <- do.call(Heatmap, c(list(
    n_mat,
    name = "P",
    col = p_cols,
    na_col = "white",
    heatmap_legend_param = list(
        title = expression(-log[10] * italic(P)),
        at = p_ticks,
        labels = p_labels,
        legend_height = unit(3.5, "cm")
    )), common))

ht_o <- do.call(Heatmap, c(list(
    o_mat,
    name = "OR",
    col = o_cols,
    na_col = "white",     # undefined, not zero
    heatmap_legend_param = list(
        title = expression(log[2] * "OR"),
        at = o_ticks,
        labels = o_labels,
        legend_height = unit(3.5, "cm")
    )), common))

# ---------------------------
# Draw
# ---------------------------
w <- 6.5 + 0.29 * ncol(n_mat)
h <- 2.5 + 0.29 * nrow(n_mat)

pdf(file.path(plots_dir, paste0("gsmap_", arm, "_P.pdf")), width = w, height = h)
draw(ht_p,
     column_title = "gsMap domain enrichment (Cauchy-combined P, 7 donors)",
     column_title_gp = gpar(fontsize = 13),
     heatmap_legend_side = "right")
dev.off()

pdf(file.path(plots_dir, paste0("gsmap_", arm, "_log2OR.pdf")), width = w, height = h)
draw(ht_o,
     column_title = "gsMap spatial specificity, domain vs rest of amygdala",
     column_title_gp = gpar(fontsize = 13),
     heatmap_legend_side = "right")
dev.off()

message("Saved to: ", plots_dir)