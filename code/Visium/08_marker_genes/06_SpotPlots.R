suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library("ComplexHeatmap")
    library("patchwork")
    library("ggplot2")
    library("escheR")
    library("BiocNeighbors")
    
})

# ============ Load data ============
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

# ============ Color palette ============
pal <- c(
  AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6", PL = "#f1e438ff",
  BL = "#035185ff", LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
  MeA = "#baf739ff", HPC = "#d6a8f8ff", CHAT = "#A0522D",
  Endothelial = "#444444", WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)

# ============ Subset to Br8325 ============
spe.subset <- spe[, spe$sample_id == "Br8325"]
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping]

# ============ Rescale spatialCoords to mm (ONCE) ============
one_area <- unique(spe.subset$capture_area)[1]
idx <- which(spe.subset$capture_area == one_area)
nn <- findKNN(spatialCoords(spe.subset)[idx, ], k = 1)
mm_per_unit <- 0.1 / median(nn$distance)  # 100 µm = 0.1 mm
spatialCoords(spe.subset) <- spatialCoords(spe.subset) * mm_per_unit

# ============ Scale bar parameters ============
xrange <- range(spatialCoords(spe.subset)[, 1])
yrange <- range(spatialCoords(spe.subset)[, 2])

bar_x   <- xrange[1] + diff(xrange) * 0.05
bar_y   <- yrange[2] - diff(yrange) * 0.02
bar_len <- 3
text_y  <- bar_y + diff(yrange) * 0.04

scalebar_layers <- list(
  annotate("segment", x = bar_x, xend = bar_x + bar_len,
           y = bar_y, yend = bar_y, linewidth = 1.5, color = "black"),
  annotate("text", x = bar_x + bar_len / 2, y = text_y,
           label = "3 mm", size = 3.5, color = "black")
)

# ============ Common themes ============
base_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 13),
  legend.text  = element_text(size = 8)
)

common_theme <- theme(
  plot.title = element_text(size = 22, margin = margin(b = 8)),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  legend.position = "bottom",
  legend.justification = "center",
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 11),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),
  plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
)

fill_scale_fullwidth <- function() {
  scale_fill_gradient(
    low = "white", high = "black", name = "logcounts",
    guide = guide_colorbar(
      direction = "horizontal", title.position = "left", title.hjust = 0,
      label.position = "bottom", barwidth = unit(.12, "npc"), barheight = unit(0.5, "cm")
    )
  )
}

# ============ 1. Histology with scale bar ============
samples <- unique(spe$sample_id)

pdf(here("plots", "Visium", "08_marker_genes", "Histology_AllSamples.pdf"), width = 6, height = 6)
for (s in samples) {
  print(
    ggspavis::plotVisium(spe, spots = FALSE, annotate = NULL, sample_id = s) +
      ggtitle(s) 
  )
}
dev.off()

# ============ 2. Individual marker gene plots (escheR) ============
plot_marker_gene <- function(spe, gene, pal, outdir = here("plots", "Visium", "08_marker_genes")) {
  if (!gene %in% rownames(spe)) stop(paste("Gene", gene, "not found."))
  spe[[gene]] <- log1p(counts(spe)[gene, ])

  pdf(file.path(outdir, paste0("Clusters_Br8325_", gene, ".pdf")), width = 6, height = 5)
  print(
    make_escheR(spe) |>
      add_fill(var = gene, point_size = 1) |>
      add_ground(var = "BS_k16_Semisupervised_wAI", point_size = 0.5, stroke = 0.25, color_aes = TRUE) +
      scale_color_manual(values = pal, guide = "none") +
      scale_fill_gradient(low = "white", high = "black", name = "logcounts") +
      ggtitle(bquote(italic(.(gene)))) +
      scalebar_layers +
      theme(
        plot.title = element_text(size = 22, margin = margin(b = 10)),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.title = element_text(size = 18), legend.text = element_text(size = 15),
        legend.key = element_blank(), legend.key.size = unit(1.2, "lines"),
        legend.spacing.y = unit(0.4, "lines"), legend.position = "right",
        legend.justification = "top", legend.background = element_blank(),
        plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
      ) +
      base_theme
  )
  dev.off()
}

markers <- c("AGBL1", "CNR1", "PEX5L", "UNCX", "GJB5","GJB3", "GLP2R",
             "PENK", "GPR88", "PDYN", "GRP", "TSHZ1", "FOXP2",
             "SATB1", "CYP26B1", "OTP", "SIM1", "LAMP5", "CRHR1", "CRHR2", "GABRG1", "GAL", 
             "BHMT", "GLYATL1", "MYRIP", "DMRT3","GLP2R", "ETV1", "ESR1", "MOXD1")

for (gene in markers) {
  plot_marker_gene(spe.subset, "CABP7", pal)
}



# ============ 3. Combined highlighted plot with scale bars ============
genes <- c("RGS4", "GJB3", "CYP26B1", "PDYN", "CNR1")
for (g in genes) {
  spe.subset[[g]] <- log1p(counts(spe.subset)[g, ])
}

mk_panel <- function(var, title_text, highlight_domain, show_scalebar = FALSE) {
  pal_highlight <- setNames(rep("white", length(pal)), names(pal))
  pal_highlight[highlight_domain] <- pal[highlight_domain]

  p <- make_escheR(spe.subset) |>
    add_fill(var = var, point_size = .75) |>
    add_ground(var = "BS_k16_Semisupervised_wAI",
               point_size = 0.75, stroke = 0.05, color_aes = TRUE) +
    scale_color_manual(values = pal_highlight, guide = "none") +
    fill_scale_fullwidth() +
    ggtitle(bquote(italic(.(title_text)))) +
    common_theme + base_theme + theme(legend.key.width = unit(.3, "npc"))

  if (show_scalebar) p <- p + scalebar_layers
  p
}

p1 <- mk_panel("RGS4",   "RGS4",   "BLD", show_scalebar = TRUE)
p2 <- mk_panel("GJB3",   "GJB3",   "PL", show_scalebar = TRUE)
p3 <- mk_panel("CYP26B1", "CYP26B1", "LA", show_scalebar = TRUE)
p4 <- mk_panel("PDYN",    "PDYN",    "CoA", show_scalebar = TRUE)
p5 <- mk_panel("CNR1",    "CNR1",    "BM", show_scalebar = TRUE)

pdf(here("plots", "Visium", "08_marker_genes", "MarkerGenes_Br8325_Combined_highlighted.pdf"),
    width = 20, height = 5)
(p1 | p2 | p3 | p4 | p5)
dev.off()