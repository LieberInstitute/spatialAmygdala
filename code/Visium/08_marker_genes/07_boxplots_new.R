suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library("patchwork")
    library("ggplot2")
})

outdir <- here("plots", "Visium", "08_marker_genes")

# ============ Load data ============
spe <- readRDS(here("processed-data", "Visium", "08_marker_genes",
                     "pseudo_BS_K16_final_labels_dropped.rds"))
spe <- spe[, spe$BS_k16_Semisupervised_wAI %in%
             c("BM", "BL", "PL", "BLD", "CeA", "CoA", "AI", "LA", "MeA")]

# Load sig genes from enrichment model
sig_genes <- read.csv(here("processed-data", "Visium", "08_marker_genes",
                            "BS_k16_sig_genes_50_wBLVM.csv"))

# Build a lookup: for each gene, which domains are significantly enriched?
# Using FDR < 0.05 as threshold — adjust as needed
sig_lookup <- sig_genes %>%
  filter(fdr < 0.05) %>%
  group_by(gene) %>%
  summarise(sig_domains = list(unique(test)), .groups = "drop")

sig_map <- setNames(sig_lookup$sig_domains, sig_lookup$gene)

# ============ Palette ============
pal <- c(
  AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6", PL = "#f1e438ff",
  BL = "#035185ff", LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
  MeA = "#baf739ff", HPC = "#d6a8f8ff", CHAT = "#A0522D",
  Endothelial = "#444444", WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)

# ============ Prep expression matrix ============
mat <- logcounts(spe)
gene_names <- rowData(spe)$gene_name
if (anyDuplicated(gene_names)) {
  gene_names <- make.unique(gene_names)
}
rownames(mat) <- gene_names

df <- as.data.frame(t(as.matrix(mat)))
df$domain <- spe$BS_k16_Semisupervised_wAI

# ============ Genes to plot ============
genes <- c("FIBCD1", "CNR1", "PEX5L", "SLIT1", "GJB5", "GLP2R",
           "PENK", "GPR88", "PDYN", "GRP",
           "TSHZ1", "FOXP2", "SATB1", "CYP26B1",
           "OTP", "SIM1", "LAMP5", "CRHR1", "CRHR2")
genes <- intersect(genes, colnames(df))

# ============ Plot ============
for (gene in genes) {

  # Get significant domains for this gene (empty = none highlighted)
  sig_domains <- sig_map[[gene]] %||% character(0)

  # Build per-gene palette: significant domains get color, rest get grey
  pal_gene <- setNames(rep("grey80", length(pal)), names(pal))
  for (d in sig_domains) {
    if (d %in% names(pal)) pal_gene[d] <- pal[d]
  }

  p <- ggplot(df, aes(x = domain, y = .data[[gene]], fill = domain)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      color = "black",
      size = 0.4
    ) +
    ggbeeswarm::geom_quasirandom(
      size = 1.8,
      alpha = 0.4,
      shape = 21,
      color = "black",
      width = 0.25,
      dodge.width = 0.6
    ) +
    theme_bw() +
    scale_fill_manual(values = pal_gene) +
    ylab("logcounts") +
    ggtitle(gene) +
    xlab("") +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1, color = "black"
      ),
      legend.position = "none"
    )

  ggsave(
    filename = file.path(outdir, paste0("boxplot_", gene, "_highlighted.pdf")),
    plot = p,
    width = 3.5,
    height = 4,
    useDingbats = FALSE
  )
}



# ============ Combined boxplot: selected genes, amygdala domains only ============
amyg_domains <- c("AI", "BM", "BLD", "BL", "PL", "LA", "CoA", "CeA", "MeA")
combo_genes <- c("RGS4", "GJB3", "CYP26B1", "PDYN", "CNR1")
combo_genes <- intersect(combo_genes, colnames(df))

df_amyg <- df[df$domain %in% amyg_domains, ]
df_amyg$domain <- factor(df_amyg$domain, levels = amyg_domains)

plots <- list()
for (gene in combo_genes) {
  sig_domains <- sig_map[[gene]] %||% character(0)
  
  pal_gene <- setNames(rep("grey80", length(pal)), names(pal))
  for (d in sig_domains) {
    if (d %in% names(pal)) pal_gene[d] <- pal[d]
  }
  
  plots[[gene]] <- ggplot(df_amyg, aes(x = domain, y = .data[[gene]], fill = domain)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", size = 0.4) +
    ggbeeswarm::geom_quasirandom(
      size = 1.8, alpha = 0.4, shape = 21, color = "black",
      width = 0.25, dodge.width = 0.6
    ) +
    theme_bw() +
    scale_fill_manual(values = pal_gene) +
    ylab("logcounts") +
    ggtitle(gene) +
    xlab("") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
      legend.position = "none"
    )
}

p_combined <- plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]] | plots[[5]]

ggsave(
  filename = file.path(outdir, "boxplot_combined_amygdala_highlighted.pdf"),
  plot = p_combined,
  width = 16,
  height = 4,
  useDingbats = FALSE
)