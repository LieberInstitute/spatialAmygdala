suppressPackageStartupMessages({
  library("here")
  library("sessioninfo")
  library("SpatialExperiment")
  library("scater")
  library("spatialLIBD")
  library("dplyr")
  library("ComplexHeatmap")
  library("patchwork")
  library("RColorBrewer")
  library("SummarizedExperiment")
  library("circlize")
  library("grid")
  library("dreamlet")
})

spe <- readRDS(here(
  "processed-data","Visium","07_clustering","BayesSpace","MarkerGenes",
  "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"
))
spe

# ========= Pseudobulk dendrogram =========
pb <- aggregateToPseudoBulk(
  spe,
  assay = "counts",
  cluster_id = "BS_k16_Semisupervised_wAI",
  sample_id = "sample_id",
  verbose = FALSE
)

hcl  <- buildClusterTreeFromPB(pb, method = "ward.D")
dend <- as.dendrogram(hcl)  # domain dendrogram

# -------------------------------
# 0) Inputs
# -------------------------------
pal <- c(
  AI          = "#D62728",
  BM          = "#E67E22",
  BLD         = "#9B59B6",
  PL          = "#f1e438ff",
  BL          = "#035185ff",
  LA          = "#F4B400",
  CoA         = "#5DA5DA",
  CeA         = "#197d43ff",
  MeA         = "#baf739ff",
  HPC         = "#d6a8f8ff",
  CHAT        = "#A0522D",
  Endothelial = "#444444",
  WM.1        = "#BBBBBB",
  WM.2        = "#DDDDDD",
  CLA         = "#FF69B4"
)

markers <- unique(c(
  "NR4A2","GNB4",         # CLA
  "FEZF2","POU3F1",        # HPC
  "RGS4","PEX5L",         # BLD
  "CYP26B1","RORB",       # LA
  "CRHR1","MOXD1",        # PL
  "SLC1A2","SLC1A3",      # BL
  "NNAT","CNR1",          # BM
  "PDYN","CHRM3",         # CoA
  "TSHZ1","FOXP2",        # AI
  "PENK","PPP1R1B",       # CeA
  "OTP","CARTPT",         # MeA
  "CHAT","SLC5A7",        # CHAT
  "MT-CO1","MT-ND2",      # WM.2
  "HBB","HBA2",           # Endothelial
  "MBP","GFAP"            # WM.1
))
markers <- markers[markers %in% rownames(spe)]

# -------------------------------
# 1) Pseudobulk mean logcounts per gene × domain
# -------------------------------
celltype <- colData(spe)$BS_k16_Semisupervised_wAI
mat <- assay(spe, "logcounts")

mat_avg <- sapply(split(seq_len(ncol(mat)), celltype), function(idx) {
  rowMeans(mat[, idx, drop = FALSE])
})
mat_avg <- mat_avg[intersect(markers, rownames(mat_avg)), , drop = FALSE]
stopifnot(nrow(mat_avg) > 0, ncol(mat_avg) > 0)

# -------------------------------
# 2) Scale per gene across domains (genes x domains)
# -------------------------------
mat_scaled <- t(scale(t(mat_avg)))  # genes x domains
mat_scaled[is.na(mat_scaled)] <- 0

genes_hm <- rownames(mat_scaled)
domains  <- colnames(mat_scaled)

# ensure every domain has a color
missing_dom <- setdiff(domains, names(pal))
if (length(missing_dom) > 0) {
  pal <- c(pal, setNames(rep("grey60", length(missing_dom)), missing_dom))
}
pal_domain <- pal[domains]

# -------------------------------
# 3) Gene grouping (for gene-domain color stripe)
# -------------------------------
gene_groups <- list(
  CLA         = c("NR4A2","GNB4"),
  HPC         = c("FEZF2","POU3F1"),
  BLD         = c("RGS4","PEX5L"),
  LA          = c("CYP26B1","RORB"),
  PL          = c("CRHR1","MOXD1"),
  BL          = c("SLC1A2","SLC1A3"),
  BM          = c("NNAT","CNR1"),
  CoA         = c("PDYN","CHRM3"),
  AI          = c("TSHZ1","FOXP2"),
  CeA         = c("PENK","PPP1R1B"),
  MeA         = c("OTP","CARTPT"),
  CHAT        = c("CHAT","SLC5A7"),
  WM.2        = c("MT-CO1","MT-ND2"),
  Endothelial = c("HBB","HBA2"),
  WM.1        = c("MBP","GFAP")
)

gene_to_domain <- unlist(lapply(names(gene_groups), function(ct) {
  setNames(rep(ct, length(gene_groups[[ct]])), gene_groups[[ct]])
}))

gene_domain <- gene_to_domain[genes_hm]
gene_domain[is.na(gene_domain)] <- "Other"
gene_domain <- factor(gene_domain)

# add "Other" to pal for gene-domain coloring if needed
if (!("Other" %in% names(pal))) pal <- c(pal, Other = "grey80")

# -------------------------------
# 4) FLIP (transpose) + annotations (NO AMY/NON-AMY GROUPS)
# -------------------------------
mat_flip <- t(mat_scaled)  # rows = domains, cols = genes

domains_flip <- rownames(mat_flip)
genes_flip   <- colnames(mat_flip)

# domain annotation (rows): ONLY domain colors
row_anno_flip <- rowAnnotation(
  domain = factor(domains_flip, levels = domains_flip),
  col = list(domain = pal[domains_flip]),
  show_legend = TRUE
)

# gene annotation (columns): ONLY gene_domain colors
col_anno_flip <- HeatmapAnnotation(
  gene_domain = gene_domain,
  col = list(gene_domain = pal),
  show_legend = FALSE,
  which = "column"
)



# -------------------------------
# 5) Heatmap (flipped; dendrogram orders domains)
# -------------------------------
ht_flip <- Heatmap(
  mat_flip,
  name = "Z-score",

  left_annotation = row_anno_flip,  # domains
  top_annotation  = col_anno_flip,  # genes

  row_gap    = unit(1.25, "mm"),
  column_gap = unit(1.25, "mm"),

  cluster_rows    = dend,   # dendrogram on domains (rows)
  cluster_columns = FALSE,

  show_row_names = TRUE,
  show_column_names = TRUE,

  row_title    = "Spatial Domain",
  column_title = "Marker genes"
)

pdf(here("plots", "Visium", "08_marker_genes",
         "pseudobulk_heatmap_FLIPPED_no_AMY_groups.pdf"),
    width = 9, height = 5)
draw(ht_flip, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
