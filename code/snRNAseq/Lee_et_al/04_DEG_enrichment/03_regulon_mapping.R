suppressPackageStartupMessages({
  library(here)
  library(spatialLIBD)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(scuttle)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

out_plots <- here("plots","snRNAseq","girgenti")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)

# domain palette (shared)
pal <- c(
  AI="#D62728", BM="#E67E22", BLD="#9B59B6", PL="#f1e438ff",
  BL="#035185ff", LA="#F4B400", CoA="#5DA5DA", CeA="#197d43ff",
  MeA="#baf739ff", HPC="#d6a8f8ff", CHAT="#A0522D",
  Endothelial="#444444", WM.1="#BBBBBB", WM.2="#DDDDDD", CLA="#FF69B4"
)

# ==============================================================================
# SHARED INPUTS
# ==============================================================================
# Spatial modeling results (BS k16) + domain relabel
load(here("processed-data","Visium","08_marker_genes","BS_k16_modeling_results.Rdata"))
relabel <- c("aBA"="BM","ITC"="AI","Meninges"="Endothelial","Chat"="CHAT","BLVM"="BL")
rename_cols <- function(mat) {
  cn <- colnames(mat)
  for (old in names(relabel))
    cn <- gsub(paste0("(_)", old, "$"), paste0("\\1", relabel[old]), cn)
  colnames(mat) <- cn; mat
}
for (nm in names(modeling_results))
  modeling_results[[nm]] <- rename_cols(modeling_results[[nm]])

# Girgenti SCE: symbol <-> ensembl
sce <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_sce.rds"))
sym2ens <- setNames(rowData(sce)$featureid, rownames(sce))
ens2sym <- setNames(rownames(sce), rowData(sce)$featureid)

# INH GRN + disease TFs
supp_dir <- here("processed-data","snRNAseq","Lee_at_al","Supplementary Data")
inh <- read.csv(file.path(supp_dir, "Supplementary Data 24.csv"))
tf_rank <- read.csv(file.path(supp_dir, "Supplementary Data 31.csv"))
deg_col <- grep("DEG", colnames(tf_rank), value=TRUE, ignore.case=TRUE)[1]
tf_col  <- grep("TF|factor|motif", colnames(tf_rank), value=TRUE, ignore.case=TRUE)[1]
disease_tfs <- tf_rank |> arrange(desc(.data[[deg_col]])) |>
  slice_head(n=20) |> pull(.data[[tf_col]])
disease_tfs <- union(disease_tfs, c("KLF6","KLF7","KLF16"))

# ==============================================================================
# ANALYSIS 1 — Regulon AUD-target enrichment vs spatial domain markers
# ==============================================================================
inh_aud <- inh |> filter(TF %in% disease_tfs, !is.na(DEG) & DEG != "")
gene_list <- lapply(split(inh_aud$geneName, inh_aud$TF), function(g)
  unique(na.omit(sym2ens[unique(g)])))
gene_list <- gene_list[lengths(gene_list) >= 10]
cat("AUD-target set sizes:\n"); print(sort(lengths(gene_list), decreasing=TRUE))

enrich <- gene_set_enrichment(
  gene_list = gene_list,
  modeling_results = modeling_results,
  model_type = "enrichment"
)
 
domain_order <- c("LA","BLD","BL","PL","HPC","CLA","MeA","BM","CoA","CHAT",
                  "CeA","AI","Endothelial","WM.1","WM.2")
 
# spatialLIBD reverses the order, so use rev below
present <- intersect(rev(domain_order), unique(as.character(enrich$test)))
enrich$test <- factor(as.character(enrich$test), levels = rev(present))
 
pdf(file.path(out_plots,"regulon_AUDtarget_marker_enrichment.pdf"), width=10, height=8)
gene_set_enrichment_plot(enrich, PThresh=12, ORcut=2,
    plot_SetSize_bar = TRUE,
    model_colors = pal)
dev.off()
 
# ==============================================================================
# ANALYSIS 2 — Disease-TF expression across spatial domains (approach B)
# ==============================================================================
spe <- readRDS(here("processed-data","Visium","07_clustering","BayesSpace",
                    "MarkerGenes","spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
domains <- as.character(colData(spe)$BS_k16_Semisupervised_wAI)
keep <- !is.na(domains); spe <- spe[, keep]; domains <- domains[keep]

pb_se <- summarizeAssayByGroup(logcounts(spe), ids = domains, statistics = "mean")
pb    <- as.matrix(assay(pb_se, "mean"))
pb_z  <- t(scale(t(pb))); pb_z[!is.finite(pb_z)] <- 0

tf_present <- intersect(disease_tfs, rownames(pb_z))
if (length(tf_present) < 2)
  stop("Too few disease TFs found. head rownames: ", paste(head(rownames(pb_z)), collapse=", "))

Bs <- pb_z[tf_present, , drop = FALSE]
Bs <- Bs[matrixStats::rowVars(Bs) > 1e-8, , drop = FALSE]
stopifnot(nrow(Bs) >= 2, ncol(Bs) >= 2)

# ---- flip: domains as ROWS, TFs as COLUMNS ----
domain_order <- c("LA","BLD","BL","PL","HPC","CLA","MeA","BM","CoA","CHAT",
                  "CeA","AI","Endothelial","WM.1","WM.2")
M <- t(Bs)                                              # now domains x TFs
domain_order <- intersect(domain_order, rownames(M))    # keep present, in order
M <- M[domain_order, , drop = FALSE]

# domain row annotation (right side), colors from pal
missing_pal <- setdiff(rownames(M), names(pal))
if (length(missing_pal)) stop("pal missing domains: ", paste(missing_pal, collapse=", "))
dom_cols <- setNames(as.character(pal[rownames(M)]), rownames(M))
right_anno <- rowAnnotation(Domain = rownames(M), col = list(Domain = dom_cols),
                            show_annotation_name = FALSE, show_legend = FALSE)

safe_cluster <- function(m, dim = c("row","col")) {
  dim <- match.arg(dim)
  x <- if (dim == "col") t(m) else m
  d <- dist(x)
  if (!all(is.finite(d)) || sum(d) == 0) return(FALSE)
  hclust(d, method = "complete")
}

ht <- Heatmap(M, name = "scaled\nexpr",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  right_annotation = right_anno,
  cluster_rows    = FALSE,                 # keep domain_order
  cluster_columns = safe_cluster(M, "col"),# cluster the TFs
  row_order = domain_order,
  row_names_side = "right",                # domain labels on the right
  rect_gp = gpar(col = "grey85", lwd = 0.4),
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  column_title = "AUD disease-TF expression across domains")

pdf(file.path(out_plots, "regulon_TF_expression.pdf"), width = 8, height = 8)
draw(ht); dev.off()
# ==============================================================================
# ANALYSIS 3 — Network of one CeA regulon, highlighting CeA-marker targets
# ==============================================================================
FOCAL_TF    <- "KLF16"
FOCAL_DOMAIN <- "CeA"
TOP_N_MARKERS <- 100       # top CeA markers by t-stat

# CeA markers: top-N enrichment genes for the focal domain, ranked by t-stat
enr <- modeling_results$enrichment
t_col <- paste0("t_stat_", FOCAL_DOMAIN)
marker_ord <- order(enr[[t_col]], decreasing = TRUE)
ceA_markers_ens <- enr$ensembl[ head(marker_ord, TOP_N_MARKERS) ]
ceA_markers_sym <- unique(na.omit(ens2sym[ceA_markers_ens]))

# Focal TF's AUD-DEG targets (with edge weight = Correlation)
focal_edges <- inh |>
  filter(TF == FOCAL_TF, !is.na(DEG) & DEG != "") |>
  group_by(geneName) |> summarise(weight = max(Correlation), .groups="drop")

focal_edges$is_marker <- focal_edges$geneName %in% ceA_markers_sym
cat(sprintf("\n%s: %d AUD-DEG targets, %d are %s markers\n",
            FOCAL_TF, nrow(focal_edges), sum(focal_edges$is_marker), FOCAL_DOMAIN))

# Build graph: TF -> targets
edges_df <- data.frame(from = FOCAL_TF, to = focal_edges$geneName, weight = focal_edges$weight)
nodes_df <- data.frame(
  name = c(FOCAL_TF, focal_edges$geneName),
  type = c("TF", ifelse(focal_edges$is_marker,
                        paste0(FOCAL_DOMAIN, " marker target"), "other target"))
)
g <- graph_from_data_frame(edges_df, vertices = nodes_df, directed = TRUE)

node_cols <- c("TF" = "black")
node_cols[paste0(FOCAL_DOMAIN, " marker target")] <- unname(pal[FOCAL_DOMAIN])
node_cols["other target"] <- "grey80"

set.seed(1)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "grey85", alpha = 0.5, show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1)) +
  geom_node_point(aes(color = type, size = type)) +
  scale_color_manual(values = node_cols, name = NULL) +
  scale_size_manual(values = c("TF"=8, setNames(3, paste0(FOCAL_DOMAIN," marker target")),
                               "other target"=1.5), guide = "none") +
  # label TF + marker targets only (avoid clutter)
  geom_node_text(aes(label = ifelse(type != "other target", name, "")),
                 repel = TRUE, size = 3, max.overlaps = 30) +
  labs(title = sprintf("%s regulon (AUD-DEG targets) — %s markers highlighted",
                       FOCAL_TF, FOCAL_DOMAIN)) +
  theme_void() + theme(legend.position = "bottom",
                       plot.title = element_text(face="bold", hjust=0.5))

ggsave(file.path(out_plots, sprintf("regulon_network_%s_%s.pdf", FOCAL_TF, FOCAL_DOMAIN)),
       p, width = 9, height = 9)



       suppressPackageStartupMessages({
  library(here); library(spatialLIBD); library(SingleCellExperiment)
  library(dplyr); library(ComplexHeatmap); library(circlize)
})

out_plots <- here("plots","snRNAseq","girgenti")
pal <- c(
  AI="#D62728", BM="#E67E22", BLD="#9B59B6", PL="#f1e438ff",
  BL="#035185ff", LA="#F4B400", CoA="#5DA5DA", CeA="#197d43ff",
  MeA="#baf739ff", HPC="#d6a8f8ff", CHAT="#A0522D",
  Endothelial="#444444", WM.1="#BBBBBB", WM.2="#DDDDDD", CLA="#FF69B4"
)

# ---- modeling results + relabel ----
load(here("processed-data","Visium","08_marker_genes","BS_k16_modeling_results.Rdata"))
relabel <- c("aBA"="BM","ITC"="AI","Meninges"="Endothelial","Chat"="CHAT","BLVM"="BL")
rename_cols <- function(mat) {
  cn <- colnames(mat)
  for (old in names(relabel))
    cn <- gsub(paste0("(_)", old, "$"), paste0("\\1", relabel[old]), cn)
  colnames(mat) <- cn; mat
}
for (nm in names(modeling_results))
  modeling_results[[nm]] <- rename_cols(modeling_results[[nm]])

# ---- symbol -> ensembl ----
sce <- readRDS(here("processed-data","snRNAseq","Lee_at_al","girgenti_sce.rds"))
sym2ens <- setNames(rowData(sce)$featureid, rownames(sce))









# ==============================================================================
# EXC GRN (Supp Data 25)
# ==============================================================================
supp_dir <- here("processed-data","snRNAseq","Lee_at_al","Supplementary Data")
exc <- read.csv(file.path(supp_dir, "Supplementary Data 25.csv"))

# Disease TFs for EXC: data-driven = TFs with the most AUD-DEG targets in EXC GRN
exc_aud <- exc |> filter(!is.na(DEG) & DEG != "")
tf_aud_counts <- exc_aud |> dplyr::count(TF, sort = TRUE)
cat("Top EXC TFs by # AUD-DEG targets:\n"); print(head(tf_aud_counts, 25))

N_TF <- 20
exc_disease_tfs <- head(tf_aud_counts$TF,50)

# (optional) also force-include the INH KLF players for cross-class comparison
# exc_disease_tfs <- union(exc_disease_tfs, c("KLF6","KLF7","KLF16"))

# ==============================================================================
# Gene sets = each EXC disease TF's AUD-DEG targets
# ==============================================================================
exc_sets <- exc_aud |> filter(TF %in% exc_disease_tfs)
gene_list <- lapply(split(exc_sets$geneName, exc_sets$TF), function(g)
  unique(na.omit(sym2ens[unique(g)])))
gene_list <- gene_list[lengths(gene_list) >= 10]
cat("\nEXC AUD-target set sizes:\n"); print(sort(lengths(gene_list), decreasing=TRUE))

# ==============================================================================
# Fisher enrichment vs spatial domain markers
# ==============================================================================
enrich <- gene_set_enrichment(
  gene_list = gene_list,
  modeling_results = modeling_results,
  model_type = "enrichment"
)

# fixed domain order (match other panels)
domain_order <- c("LA","BLD","BL","PL","HPC","CLA","MeA","BM","CoA","CHAT",
                  "CeA","AI","Endothelial","WM.1","WM.2")
present <- intersect(rev(domain_order), unique(as.character(enrich$test)))
enrich$test <- factor(as.character(enrich$test), levels = rev(present))

pdf(file.path(out_plots,"regulon_EXC_AUDtarget_marker_enrichment.pdf"), width=4, height=6)
gene_set_enrichment_plot(enrich, PThresh=12, ORcut=2,
    plot_SetSize_bar = TRUE, model_colors = pal)
dev.off()

write.csv(enrich, here("processed-data","snRNAseq","dreamlet",
                       "regulon_EXC_AUDtarget_marker_enrichment.csv"), row.names=FALSE)



# ==============================================================================
# EXC disease-TF expression across spatial domains (approach B)
# ==============================================================================

domains <- as.character(colData(spe)$BS_k16_Semisupervised_wAI)
keep <- !is.na(domains); spe <- spe[, keep]; domains <- domains[keep]

pb_se <- summarizeAssayByGroup(logcounts(spe), ids = domains, statistics = "mean")
pb    <- as.matrix(assay(pb_se, "mean"))
pb_z  <- t(scale(t(pb))); pb_z[!is.finite(pb_z)] <- 0

tf_present <- intersect(exc_disease_tfs, rownames(pb_z))
if (length(tf_present) < 2)
  stop("Too few EXC TFs found. head rownames: ", paste(head(rownames(pb_z)), collapse=", "))

Bs <- pb_z[tf_present, , drop = FALSE]
Bs <- Bs[matrixStats::rowVars(Bs) > 1e-8, , drop = FALSE]
stopifnot(nrow(Bs) >= 2, ncol(Bs) >= 2)

# flip: domains as ROWS in fixed order, TFs as COLUMNS
domain_order <- c("LA","BLD","BL","PL","HPC","CLA","MeA","BM","CoA","CHAT",
                  "CeA","AI","Endothelial","WM.1","WM.2")
M <- t(Bs)
domain_order <- intersect(domain_order, rownames(M))
M <- M[domain_order, , drop = FALSE]

missing_pal <- setdiff(rownames(M), names(pal))
if (length(missing_pal)) stop("pal missing domains: ", paste(missing_pal, collapse=", "))
dom_cols <- setNames(as.character(pal[rownames(M)]), rownames(M))
right_anno <- rowAnnotation(Domain = rownames(M), col = list(Domain = dom_cols),
                            show_annotation_name = FALSE, show_legend = FALSE)

safe_cluster <- function(m, dim = c("row","col")) {
  dim <- match.arg(dim)
  x <- if (dim == "col") t(m) else m
  d <- dist(x)
  if (!all(is.finite(d)) || sum(d) == 0) return(FALSE)
  hclust(d, method = "complete")
}

ht <- Heatmap(M, name = "scaled\nexpr",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  right_annotation = right_anno,
  cluster_rows = FALSE, row_order = domain_order,
  cluster_columns = safe_cluster(M, "col"),
  row_names_side = "right",
  rect_gp = gpar(col = "grey85", lwd = 0.4),
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  column_title = "EXC AUD disease-TF expression across domains")

pdf(file.path(out_plots, "regulon_EXC_TF_expression.pdf"), width = 6, height = 7)
draw(ht); dev.off()