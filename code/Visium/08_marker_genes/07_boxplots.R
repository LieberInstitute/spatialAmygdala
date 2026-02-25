suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('EnhancedVolcano')
    library("patchwork")
})

outdir <- here("plots", "Visium", "08_marker_genes")

spe <- readRDS(here("processed-data", "Visium", "08_marker_genes","pseudo_BS_K16_final_labels_dropped.rds"))
spe <- spe[, spe$BS_k16_Semisupervised_wAI %in% c("BM","BA","PL", "BLVM","CeA","CoA","AI","LA","MeA")]

# Confirm matrix type
mat <- logcounts(spe)

# transpose to spots x genes
df <- as.data.frame(t(as.matrix(mat)))
df$domain <- spe$BS_k16_Semisupervised_wAI

# genes to plot
genes <- c("FIBCD1", "CNR1", "PEX5L", "SLIT1", "GJB5", "GLP2R",
             "PENK", "GPR88","PDYN", "GRP",
              "TSHZ1", "FOXP2", "SATB1", "CYP26B1",
             "OTP", "SIM1", "LAMP5", "CRHR1", "CRHR2")

# keep only valid genes
genes <- intersect(genes, colnames(df))

pal <- c( 
  AI        = "#D62728",  # strong red standout
  BM        = "#E67E22",  # burnt orange
  BLD      = "#9B59B6",  # violet
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


# confirm and fix duplicate gene names before transpose
gene_names <- rowData(spe)$gene_name
if (anyDuplicated(gene_names)) {
  message("Found duplicated gene names — making them unique.")
  gene_names <- make.unique(gene_names)
}
rownames(mat) <- gene_names

# transpose to spots x genes
df <- as.data.frame(t(as.matrix(mat)))

# attach domain info
df$domain <- spe$BS_k16_Semisupervised_wAI

# plotting boxplots
for (gene in genes) {

  p <- ggplot(df, aes(x = domain, y = .data[[gene]], fill = domain)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      color = "black",
      size = 0.4
    ) +
    ggbeeswarm::geom_quasirandom(
      size = 1.8,           # bigger points
      alpha = 0.4,          # more transparent
      shape = 21,           # hollow circle with fill
      color = "black",
      width = 0.25,
      dodge.width = 0.6
    ) +
    theme_bw() +
    scale_fill_manual(values = pal) +
    ylab("logcounts") +
    ggtitle(gene) +
    xlab("") +
    #ylim(c(0, 7.5)) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        color = "black"
      ),
      legend.position = "none",
    )

  ggsave(
    filename = file.path(outdir, paste0("boxplot_", gene, ".pdf")),
    plot = p,
    width = 3.5,
    height = 4,
    useDingbats = FALSE
  )
}
