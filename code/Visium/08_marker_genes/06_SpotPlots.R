suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('ComplexHeatmap')
    library("patchwork")
    library("ggplot2")
    library("escheR")
})


spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe 


# ====== plot histology =======
pdf(here("plots","Visium", "08_marker_genes", "Histology_Br8325.pdf"), width = 6, height = 6)
ggspavis::plotVisium(spe, spots=FALSE, annotate=NULL, sample_id="Br8325")
dev.off()


pal <- c( 
  AI        = "#D62728",  # strong red standout
  BM        = "#E67E22",  # burnt orange
  BLD         = "#9B59B6",  # violet
  PL     = "#f1e438ff",  # steel blue
  BL     = "#035185ff",  # sky blue
  LA         = "#F4B400",  # goldenrod
  CoA        = "#5DA5DA",  # blue-gray
  CeA        = "#197d43ff",  # emerald green
  MeA        = "#baf739ff",  # chartreuse
  HPC        = "#d6a8f8ff",  # deep royal purple (distinct from BA)
  CHAT       = "#A0522D",  # chestnut brown
  Endothelial   = "#444444",  # charcoal gray
  WM.1       = "#BBBBBB",  # light slate gray
  WM.2       = "#DDDDDD",   # mist gray
  CLA        = "#FF69B4"   # hot pink
)


# get your domains (already a factor)
domains <- spe$BS_k16_Semisupervised_wAI

# ensure they’re a factor with the original order preserved
spe$spatial_domains <- factor(domains)

# reorder the palette to match the factor levels actually present in the data
pal_matched <- pal[names(pal) %in% levels(spe$spatial_domains)]
pal_matched <- pal_matched[match(levels(spe$spatial_domains), names(pal_matched))]

# now assign colors in that matching order
spe$spatial_domain_colors <- pal_matched[as.character(spe$spatial_domains)]

spe$spatial_domain_colors <- factor(spe$spatial_domain_colors,
                                    levels = pal_matched)

unique(spe$spatial_domain_colors)


# ============ Plotting broad marker genes without escheR ==============

spe.subset <- spe[, spe$sample_id == "Br8325"]
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping]

pdf(here("plots","Visium", "08_marker_genes", "LibrarySize_and_SNAP25_Br8325.pdf"), width = 12, height = 4)

base_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 13),
  legend.text  = element_text(size = 8)
)

p1 <- ggspavis::plotCoords(spe.subset, annotate = "sum_umi", pal="viridis", point_size=.2) +
  ggtitle("Library Size") +
  base_theme
p2 <- ggspavis::plotCoords(spe.subset, annotate = "SLC17A7", pal="viridis", point_size=.2) + 
  ggtitle("SLC17A7") +
  base_theme
p3 <- ggspavis::plotCoords(spe.subset, annotate = "PRKCD", pal="viridis", point_size=.2) + 
  ggtitle("") +
  base_theme

# display panels using patchwork
p1 | p2 | p3
dev.off()



# ============ Plotting marker genes with escheR ==============
# subset to Br8325
spe.subset <- spe[, spe$sample_id == "Br8325"]

# drop $exclude_overlapping   
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping]


plot_marker_gene <- function(spe, gene, pal, outdir = here("plots", "Visium", "08_marker_genes")) {
  # Ensure gene exists
  if (!gene %in% rownames(spe)) {
    stop(paste("Gene", gene, "not found in object."))
  }
  
  # Add log1p-transformed expression to colData
  spe[[gene]] <- log1p(counts(spe)[gene, ])
  
  # Construct output filename
  outfile <- file.path(outdir, paste0("Clusters_Br8325_", gene, ".pdf"))
  
  # Create PDF
  pdf(outfile, width = 6, height = 5)
  print(
  make_escheR(spe) |>
    add_fill(var = gene, point_size = 1) |>
    add_ground(var = "BS_k16_Semisupervised_wAI", point_size = 0.5, stroke = 0.25, color_aes = TRUE) +
    
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_gradient(low = "white", high = "black", name = "logcounts") +
    
    ggtitle(bquote(italic(.(gene)))) +  # italicize dynamically
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 10)),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.key = element_blank(),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "lines"),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_blank(),
      
      plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
    ) +
    base_theme
  )
  
  dev.off()

}


# genes to plot
markers <- c("FIBCD1", "CNR1", "PEX5L", "SLIT1", "GJB5", "GLP2R",
             "PENK", "GPR88","PDYN", "GRP",
              "TSHZ1", "FOXP2", "SATB1", "CYP26B1",
             "OTP", "SIM1", "LAMP5", "PDYN", "CRHR1", "CRHR2", "GABRG1")

for (gene in markers) {
  plot_marker_gene(spe.subset, gene, pal)
}

plot_marker_gene(spe.subset, "PDYN", pal)



# ============ Combined plot of selected marker genes with escheR ==============



pdf(here("plots","Visium", "08_marker_genes", "MarkerGenes_Br8325_Combined.pdf"),
    width = 20, height=5)

# add genes to colData (log1p counts)
genes <- c("PEX5L","LAMP5","CYP26B1","PDYN","CNR1")
for (g in genes) {
  spe.subset[[g]] <- log1p(counts(spe.subset)[g, ])
}

# Make each legend a full-width horizontal colorbar under its panel
fill_scale_fullwidth <- function() {
  scale_fill_gradient(
    low = "white", high = "black", name = "logcounts",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0,
      label.position = "bottom",
      # KEY: stretch to full available width of the panel
      barwidth = unit(.12, "npc"),
      barheight = unit(0.5, "cm")
    )
  )
}

common_theme <- theme(
  plot.title = element_text(size = 22, margin = margin(b = 8)),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),

  legend.position = "bottom",
  legend.justification = "center",
  legend.title = element_text(size = 14),
  legend.text  = element_text(size = 11),

  # remove extra padding around legend so it can span
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),

  # tight gutters, but keep a hair of bottom space
  plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
)

mk_panel <- function(var, title_text) {
  make_escheR(spe.subset) |>
    add_fill(var = var, point_size = .75) |>
    add_ground(
      var = "BS_k16_Semisupervised_wAI",
      point_size = 0.75, stroke = 0.05, color_aes = TRUE
    ) +
    scale_color_manual(values = pal, guide = "none") +
    fill_scale_fullwidth() +
    ggtitle(bquote(italic(.(title_text)))) +
    common_theme +
    base_theme + theme(legend.key.width = unit(.3, "npc"))
}

p1 <- mk_panel("PEX5L",   "PEX5L")
p2 <- mk_panel("LAMP5",   "LAMP5")
p3 <- mk_panel("CYP26B1", "CYP26B1")
p4 <- mk_panel("PDYN",    "PDYN")
p5 <- mk_panel("CNR1",    "CNR1")

(p1 | p2 | p3 | p4 | p5)

dev.off()
