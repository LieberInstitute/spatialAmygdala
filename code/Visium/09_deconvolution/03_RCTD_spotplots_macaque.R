library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SpatialFeatureExperiment)
library(scater)
library(scran)
library(Voyager)
library(ggspavis)
library(harmony)
library(escheR)

# set directories
plots_dir <- here("plots", "Visium", "09_deconvolution")
processed_dir <- here("processed-data", "Visium", "09_deconvolution")

# load
spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))


# load
res <- readRDS(here(processed_dir, "rctd_results_human.rds"))
colnames(colData(res))
# [1] "cell_type_list" "conf_list"      "conv_all"       "conv_sub"      
# [5] "min_score"      "sample_id" 


ws <- assay(res, "weights")
table(colSums(ws) == 0)
#  FALSE   TRUE 
# 202038  21358 

# add proportion estimates as metadata
ws <- data.frame(t(as.matrix(ws)))
colData(spe)[names(ws)] <- ws[colnames(spe), ]

celltypes <- colnames(ws)



# ======= Define color palette for spatial domains =======
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




# ============ Plotting marker genes with escheR ==============
# subset to Br8325
spe.subset <- spe[, spe$sample_id == "Br8325"]

# drop $exclude_overlapping   
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping]

base_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 13),
  legend.text  = element_text(size = 8)
)

plot_marker_celltype <- function(spe, celltype, pal, outdir = here("plots", "Visium", "09_deconvolution")) {
  
  # Construct output filename
  outfile <- file.path(outdir, paste0("Celltypes_Br8325_", celltype, ".pdf"))
  
  # Create PDF
  pdf(outfile, width = 5, height = 5)
  print(
  make_escheR(spe) |>
    add_fill(var = celltype, point_size = 1) |>
    add_ground(var = "BS_k16_Semisupervised_wAI", point_size = 0.5, stroke = 0.25, color_aes = TRUE) +
    
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_gradient(low = "white", high = "black", name = "logcounts") +
    
    ggtitle(bquote(italic(.(celltype)))) +  # italicize dynamically
    theme(
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
        )+
    base_theme
  )
  
  dev.off()

}


for (celltype in celltypes) {
  plot_marker_celltype(spe.subset, celltype, pal)
}








# ============ Plotting deconvolution cell types with escheR (combined) ==============

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(grid)      # unit()
  library(patchwork) # | operator
  # library(escheR)  # assuming already loaded in your session
})

# subset to Br8325
spe.subset <- spe[, spe$sample_id == "Br8325"]

# drop overlapping spots if that column exists / is used
spe.subset <- spe.subset[, !spe.subset$exclude_overlapping]

# cell types you want
celltypes <- c(
  "BA_MOXD1",
  "BA_MYRIP",
  "LA_RORB",
  "LA_ZBTB20",
  "aBA_ESR1",
  "DLK1_ZFHX3"
)

# --- safety check: make sure these exist in colData ---
missing_ct <- setdiff(celltypes, colnames(colData(spe.subset)))
if (length(missing_ct) > 0) {
  stop("These cell type columns are missing from spe.subset colData: ",
       paste(missing_ct, collapse = ", "))
}

# base theme (your same idea)
base_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(size = 13),
  legend.text  = element_text(size = 8)
)

# Make each legend a full-width horizontal colorbar under its panel
fill_scale_fullwidth <- function(name = "logcounts") {
  scale_fill_gradient(
    low = "white", high = "black", name = name,
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      title.hjust = 0,
      label.position = "bottom",
      # stretch across available panel width
      barwidth  = unit(1.2, "npc"),
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

  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(0, 0, 0, 0),

  plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
)

mk_celltype_panel <- function(var) {
  make_escheR(spe.subset) |>
    add_fill(var = var, point_size = 0.75) |>
    add_ground(
      var = "BS_k16_Semisupervised_wAI",
      point_size = 0.75,
      stroke = 0.05,
      color_aes = TRUE
    ) +
    scale_color_manual(values = pal, guide = "none") +
    fill_scale_fullwidth(name = "logcounts") +
    ggtitle(bquote(italic(.(var)))) +
    common_theme +
    base_theme +
    theme(legend.key.width = unit(0.3, "npc"))
}

# build plots (kept explicit so it's easy to reorder)
p1 <- mk_celltype_panel("BA_MOXD1")
p2 <- mk_celltype_panel("BA_MYRIP")
p3 <- mk_celltype_panel("LA_RORB")
p4 <- mk_celltype_panel("LA_ZBTB20")
p5 <- mk_celltype_panel("aBA_ESR1")
p6 <- mk_celltype_panel("DLK1_ZFHX3")

# output
outdir <- here("plots", "Visium", "09_deconvolution")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(outdir, "Celltypes_Br8325_Combined.pdf"),
    width = 24, height = 5)

print(p1 | p2 | p3 | p4 | p5 | p6)

dev.off()
