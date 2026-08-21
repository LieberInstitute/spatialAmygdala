suppressPackageStartupMessages({
    library("ComplexHeatmap")
    library("circlize")
    library("grid")
})

# ---------------------------------------------------------------------------
# Supplementary: do psychiatric enrichments replicate across donors?
#
# Every other figure pools the 7 capture areas (cauchy_across_samples). This one
# deliberately does not: each panel is one domain, rows are psychiatric traits and
# columns are the 7 donors, so a row that is red in one column and grey in six is
# visibly a single-section result rather than a replicated one.
#
# Input:  figdata/donor_cauchy_long_test2.csv   (25_build_donor_cauchy_long.py)
# Output: gsmap_donor_replication_test2.pdf
# ---------------------------------------------------------------------------

proj      <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
data_dir  <- file.path(proj, "code/Visium/16_LDSC/gsMap/figdata")
plots_dir <- file.path(proj, "processed-data/Visium/16_LDSC/gsMap/figures")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

arm <- "test2"   # functional-conditioned arm
thr <- -log10(0.05)

long <- read.csv(file.path(data_dir, paste0("donor_cauchy_long_", arm, ".csv")),
                 check.names = FALSE, stringsAsFactors = FALSE)

# AI -> IA, as in the other figures
long$domain[long$domain == "AI"] <- "IA"

# ---------------------------
# Traits: the psychiatric set, best-powered GWAS per phenotype
# ---------------------------
# where a phenotype has several freezes, keep the one used in the main figures --
# otherwise near-duplicate rows (three MDDs, two SCZs) dominate the panel
trait_labels <- c(
    SCZ          = "Schizophrenia",
    BIP_PGC3     = "Bipolar disorder",
    MDD          = "Depression",
    ADHD         = "ADHD",
    PTSD_F3      = "PTSD",
    Autism       = "Autism",
    Anorexia     = "Anorexia"
)

n_psych_total <- length(unique(long$trait))     # before subsetting
long <- long[long$trait %in% names(trait_labels), ]
stopifnot(nrow(long) > 0)
n_psych_shown <- length(unique(long$trait))
message("psychiatric GWAS: showing ", n_psych_shown, " of ", n_psych_total,
        " (one freeze per phenotype)")

# ---------------------------
# Domains and donors
# ---------------------------
domain_order <- c("CLA", "HPC", "BLD", "LA", "PL", "BL", "BM", "CoA",
                  "IA", "CeA", "MeA", "CHAT", "Endothelial", "WM.1", "WM.2")
# intersect() would silently drop any domain missing from domain_order, so say so
dropped <- setdiff(unique(long$domain), domain_order)
if (length(dropped)) stop("domain_order omits: ", paste(dropped, collapse = ", "))
domain_order <- intersect(domain_order, unique(long$domain))
message("domains plotted: ", length(domain_order), " of ",
        length(unique(long$domain)), " in the table")

donor_order <- sort(unique(long$donor))

pal <- c(
    IA = "#D62728", BM = "#E67E22", BLD = "#9B59B6", PL = "#f1e438",
    BL = "#035185", LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43",
    MeA = "#baf739", HPC = "#d6a8f8", CHAT = "#A0522D",
    Endothelial = "#444444", WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)

# ---------------------------
# One trait x donor matrix per domain
# ---------------------------
trait_order <- names(trait_labels)

mk <- function(dom) {
    d <- long[long$domain == dom, ]
    m <- matrix(NA_real_, nrow = length(trait_order), ncol = length(donor_order),
                dimnames = list(trait_labels[trait_order], donor_order))
    idx <- cbind(match(d$trait, trait_order), match(d$donor, donor_order))
    m[idx] <- -log10(pmax(d$p_cauchy, 1e-300))
    m
}
mats <- lapply(domain_order, mk)
names(mats) <- domain_order

# shared colour scale across ALL panels, so a cell's colour means the same thing
# everywhere; capped at the 99th percentile rather than the max, since one runaway
# cell would otherwise flatten every other panel
all_vals <- unlist(mats, use.names = FALSE)
p_cap <- unname(ceiling(quantile(all_vals, 0.99, na.rm = TRUE)))
message("cells clamped at cap: ", sum(all_vals > p_cap, na.rm = TRUE),
        " of ", sum(!is.na(all_vals)))

# grey below 0.05, white AT it, orange-red above -- the grey/orange boundary is
# significance, not the middle of the data range
p_cols <- colorRamp2(
    c(0, thr / 2, thr,
      thr + (p_cap - thr) / 3, thr + 2 * (p_cap - thr) / 3, p_cap),
    c("#D2D2D2", "#E4E4E4", "#FFFFFF", "#FDBE85", "#F16913", "#B30000")
)

p_ticks  <- c(0, thr, setdiff(pretty(c(thr, p_cap), 3), 0))
p_ticks  <- sort(unique(p_ticks[p_ticks <= p_cap]))
p_labels <- format(round(p_ticks, 1), trim = TRUE)
p_labels[p_ticks == thr] <- "1.3"
if (any(all_vals > p_cap, na.rm = TRUE)) {
    p_labels[length(p_labels)] <- paste0("\u2265", p_labels[length(p_labels)])
}

# ---------------------------
# Assemble: one Heatmap per domain, wrapped over two rows of panels
# ---------------------------
# 14 domains in a single horizontal strip is ~10in wide and unreadable; splitting
# them over two stacked rows keeps cells square-ish and the figure page-shaped.
n_absent <- sum(is.na(unlist(mats, use.names = FALSE)))
message("structurally absent cells (domain not in that donor): ", n_absent)

panel <- function(dom, first_in_row, show_leg) {
    dom <- dom   # captured by cell_fun below
    Heatmap(
        mats[[dom]],
        name = paste0("P_", dom),
        col = p_cols,
        # absent must NOT look like p = 0.05 (which is white on this scale): draw a
        # diagonal slash so "domain not sampled in this donor" reads as a gap
        na_col = "#FFFFFF",
        cell_fun = function(j, i, x, y, w_, h_, fill) {
            if (is.na(mats[[dom]][i, j])) {
                grid.lines(x = c(x - w_ * 0.5, x + w_ * 0.5),
                           y = c(y - h_ * 0.5, y + h_ * 0.5),
                           gp = gpar(col = "#B0B0B0", lwd = 0.6))
            }
        },
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = first_in_row,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 6.5),
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 0.5),
        column_title = dom,
        column_title_gp = gpar(
            fontsize = 8, fill = pal[dom], border = NA,
            col = ifelse(dom %in% c("BL", "Endothelial", "CHAT", "CeA"), "white", "black")),
        width  = unit(3.0 * ncol(mats[[dom]]), "mm"),
        height = unit(5.5 * nrow(mats[[dom]]), "mm"),
        show_heatmap_legend = show_leg,
        heatmap_legend_param = list(
            title = expression(-log[10] * italic(P)),
            at = p_ticks, labels = p_labels,
            legend_height = unit(2.6, "cm"),
            title_gp = gpar(fontsize = 8),
            labels_gp = gpar(fontsize = 7)
        )
    )
}

split_at <- ceiling(length(domain_order) / 2)
rows_of  <- list(domain_order[seq_len(split_at)],
                 domain_order[(split_at + 1):length(domain_order)])

build_row <- function(doms, show_leg) {
    hl <- NULL
    for (j in seq_along(doms)) {
        h <- panel(doms[j], first_in_row = (j == 1),
                   show_leg = (show_leg && j == 1))
        hl <- if (is.null(hl)) h else hl + h
    }
    hl
}

# ComplexHeatmap refuses to nest a horizontal list inside a vertical one, so each
# row of panels is drawn into its own viewport rather than combined with %v%.
w <- 2.0 + 0.135 * max(lengths(rows_of)) * ncol(mats[[1]])
h <- 0.80 + 2 * (0.205 * length(trait_order) + 0.60)

pdf(file.path(plots_dir, paste0("gsmap_donor_replication_", arm, ".pdf")),
    width = w, height = h)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3,
             heights = unit.c(unit(1, "null"), unit(1, "null"), unit(4.5, "mm")))))

pushViewport(viewport(layout.pos.row = 1))
draw(build_row(rows_of[[1]], TRUE), newpage = FALSE,
     column_title = sprintf(
         "Psychiatric enrichment by donor (functional baseline, unpooled) \u2014 %d of %d GWAS",
         n_psych_shown, n_psych_total),
     column_title_gp = gpar(fontsize = 11),
     heatmap_legend_side = "right",
     ht_gap = unit(1.5, "mm"),
     padding = unit(c(1, 5, 4, 1), "mm"))   # b, l, t, r -- top pad clears the title
popViewport()

pushViewport(viewport(layout.pos.row = 2))
draw(build_row(rows_of[[2]], FALSE), newpage = FALSE,
     heatmap_legend_side = "right",
     ht_gap = unit(1.5, "mm"),
     padding = unit(c(1, 5, 1, 1), "mm"))
popViewport()

pushViewport(viewport(layout.pos.row = 3))
grid.text(paste0("One panel per domain (", length(domain_order), " of ",
                 length(unique(read.csv(file.path(data_dir, paste0(
                     "donor_cauchy_long_", arm, ".csv")))$domain)),
                 "); rows are traits, columns the 7 donors. ",
                 "Where a phenotype has several GWAS freezes only the best-powered is shown, ",
                 "so ", n_psych_shown, " of ", n_psych_total, " psychiatric GWAS appear. ",
                 "Slashed cells = domain absent from that donor (CLA 1/7, LA and IA 6/7). ",
                 "White = P = 0.05; grey below, orange-red above."),
          x = unit(3, "mm"), y = 0.5, just = "left",
          gp = gpar(fontsize = 5.6, col = "#555555"))
popViewport()

popViewport()
dev.off()

message("PDF size: ", round(w, 2), " x ", round(h, 2), " in")
message("Saved to: ", plots_dir)
