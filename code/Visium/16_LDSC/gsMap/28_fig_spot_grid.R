suppressPackageStartupMessages({
    library("ggplot2")
    library("data.table")
    library("ggrastr")
})

# ---------------------------------------------------------------------------
# Supplementary: per-spot gsMap enrichment -log10(P) for every psychiatric trait (rows) in every
# donor (columns), functional-conditioned arm.
#
# The pooled heatmaps say WHICH domain is enriched; this says whether the spatial
# pattern behind that number looks the same section to section. A trait whose
# signal is real should show the same anatomical shape in all 7 columns.
#
# Input:  figdata/spot_grid_test2_Psychiatric.csv.gz  (27_build_spot_grid_data.py)
#         figdata/spot_grid_test2_Psychiatric_meta.csv
# Output: gsmap_spot_grid_test2_Psychiatric.pdf
# ---------------------------------------------------------------------------

proj      <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
data_dir  <- file.path(proj, "code/Visium/16_LDSC/gsMap/figdata")
plots_dir <- file.path(proj, "processed-data/Visium/16_LDSC/gsMap/figures")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

arm <- "test2"
grp <- "Psychiatric"
thr <- -log10(0.05)

dt   <- fread(file.path(data_dir, sprintf("spot_grid_%s_%s.csv.gz", arm, grp)))
meta <- fread(file.path(data_dir, sprintf("spot_grid_%s_%s_meta.csv", arm, grp)))

vlim <- meta$vlim[1]
trait_order <- strsplit(meta$trait_order[1], "|", fixed = TRUE)[[1]]
donor_order <- strsplit(meta$donor_order[1], "|", fixed = TRUE)[[1]]

trait_labels <- c(
    SCZ = "Schizophrenia", BIP_PGC3 = "Bipolar", MDD = "Depression",
    ADHD = "ADHD", PTSD_F3 = "PTSD", Autism = "Autism", Anorexia = "Anorexia"
)

stopifnot(all(trait_order %in% names(trait_labels)))
dt[, trait := factor(trait_labels[trait], levels = trait_labels[trait_order])]
dt[, donor := factor(donor, levels = donor_order)]

# clamp rather than drop, so an extreme spot keeps its position and saturates in
# colour instead of leaving a hole in the tissue
dt[, nlp := pmin(nlp, vlim)]

# draw the most significant spots last so signal is not hidden under null spots
setorder(dt, nlp)

# Point size must match the SPOT PITCH, not be minimised: a Visium section spans
# ~130 spot columns across a ~1in panel, so each spot owns ~1/130 in. Drawing points
# smaller than that leaves gaps between them that alias into stripes and moire under
# rasterisation. pt_size below is that pitch converted to ggplot's mm-based size unit.
spot_pitch_in <- (1.02 * 2) / 130          # panel width in inches / spot columns
pt_size <- spot_pitch_in * 25.4 / .pt * 1.5 # -> ggplot size units, 1.5x for overlap

# 1.6M points as vectors makes a ~76 MB PDF that no viewer will open comfortably;
# rasterise the spot layer only, leaving all text and guides as vectors
p <- ggplot(dt, aes(x, y, colour = nlp)) +
    rasterise(geom_point(size = pt_size, stroke = 0, shape = 16), dpi = 600) +
    # grey below 0.05, white AT it, orange-red above -- the grey/white boundary IS
    # the significance threshold, so a reader sees significance without a legend lookup
    scale_colour_gradientn(
        colours = c("#D2D2D2", "#E4E4E4", "#FFFFFF", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"),
        values  = scales::rescale(c(0, thr * 0.75, thr,
                                    thr + (vlim - thr) * 0.25,
                                    thr + (vlim - thr) * 0.5,
                                    thr + (vlim - thr) * 0.75, vlim),
                                  from = c(0, vlim)),
        limits = c(0, vlim),
        # the top break must be vlim ITSELF, not round(vlim): rounding up puts it
        # outside limits and ggplot drops it silently, leaving the cap unlabelled
        breaks = c(0, thr, round(vlim / 2), vlim),
        labels = function(b) ifelse(abs(b - thr) < 1e-6, "1.3",
                             # 1 dp, matching the caption -- %.0f would print 12 for
                             # a cap of 11.554 and contradict "capped at 11.6"
                             ifelse(b >= vlim - 1e-6, sprintf("\u2265%.1f", vlim),
                                    sprintf("%g", round(b, 1)))),
        guide = guide_colourbar(barheight = unit(2.4, "cm"), barwidth = unit(0.25, "cm"),
                               title.position = "top")
    ) +
    facet_grid(trait ~ donor, switch = "y") +
    scale_y_reverse() +           # Visium pixel y increases downward
    coord_fixed() +
    labs(colour = expression(-log[10](italic(P))),
         title = sprintf("Psychiatric spatial enrichment by donor (functional baseline) \u2014 %d of 12 GWAS",
                         length(trait_order))) +
    theme_void(base_size = 9) +
    theme(
        strip.text.x = element_text(size = 7, margin = margin(b = 1.5)),
        strip.text.y.left = element_text(size = 7.5, angle = 0, hjust = 1,
                                        margin = margin(r = 2.5)),
        panel.spacing = unit(0.06, "lines"),
        plot.title = element_text(size = 10, hjust = 0, margin = margin(b = 4)),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        plot.margin = margin(4, 4, 4, 4)
    )

cap <- paste0(
    "Per-spot gsMap enrichment P under the full functional baseline; grey = above the ",
    "0.05 threshold, white at it, orange-red below. ",
    "One-sided scale shared across all ", length(trait_order) * length(donor_order),
    " panels (capped at ", sprintf("%.1f", vlim), ", the 99.5th percentile; ",
    sprintf("%.2f", meta$pct_clamped[1]), "% of spots clamped).\n",
    "Where a phenotype has several GWAS freezes only the best-powered is shown, so ",
    length(trait_order), " of 12 psychiatric GWAS appear. Coordinates centred and scaled ",
    "per donor, so panels are comparable in size but not in absolute dimensions; ",
    "y reversed for anatomical orientation."
)
# hard-wrap rather than trusting the device: an over-long caption line is silently
# clipped at the page edge instead of wrapping
cap <- paste(strwrap(gsub("\n", " ", cap), width = 165), collapse = "\n")

p <- p + labs(caption = cap) +
    theme(plot.caption = element_text(size = 5.4, colour = "#555555", hjust = 0,
                                      margin = margin(t = 5)))

w <- 1.35 + 1.02 * length(donor_order)
h <- 1.05 + 0.11 * (1 + lengths(regmatches(cap, gregexpr("\n", cap)))) +
     1.02 * length(trait_order)          # caption height scales with its line count
ggsave(file.path(plots_dir, sprintf("gsmap_spot_grid_%s_%s.pdf", arm, grp)),
       p, width = w, height = h, dpi = 450)

message("panels: ", length(trait_order), " x ", length(donor_order),
        "  spots: ", nrow(dt), "  vlim: ", vlim)
message("PDF size: ", round(w, 2), " x ", round(h, 2), " in")
message("Saved to: ", plots_dir)
