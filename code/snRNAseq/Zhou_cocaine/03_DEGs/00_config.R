#!/usr/bin/env Rscript


# ── Packages ─────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(spatialLIBD)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(biomaRt)
  library(here)
})


outdir <- here("plots", "snRNAseq", "Zhou_cocaine")
spe <- readRDS(here("processed-data", "Visium", "08_marker_genes","pseudo_BS_K16_final_labels_dropped.rds"))

load(here("processed-data", "Visium", "08_marker_genes", "BS_k16_modeling_results.Rdata"))
modeling_results

degs <- fread(here("processed-data", "snRNAseq", "Zhou_cocaine",
                   "41593_2023_1452_MOESM7_ESM.txt"))
head(degs)
#     gene        p_val  avg_log2FC pct.1 pct.2    p_val_adj   celltype
# 1   Pld5 6.375528e-66  0.08877314 0.212 0.129 1.102775e-61 Astrocytes
# 2 Zbtb16 7.468368e-27  0.14479552 0.285 0.200 1.291804e-22 Astrocytes
# 3 Rabep1 1.959923e-21 -0.15079978 0.509 0.593 3.390079e-17 Astrocytes
# 4   Lix1 1.052483e-18 -0.15627988 0.625 0.690 1.820480e-14 Astrocytes
# 5  Mtss2 1.154150e-18 -0.07517816 0.299 0.353 1.996332e-14 Astrocytes
# 6   Grm3 4.468725e-17 -0.11065874 0.831 0.863 7.729554e-13 Astrocytes
#          q_val     zscore
# 1 1.000193e-61   6.939796
# 2 5.858188e-23  11.224914
# 3 1.024909e-17 -11.384990
# 4 3.621259e-15 -11.804159
# 5 3.621259e-15  -5.600739
# 6 1.168423e-13  -8.314626

# ── Filter DEGs ──────────────────────────────────────────────────────────────
 
fdr_thresh   <- 0.10
log2fc_large <- 0.1
 
deg_sig <- degs[p_val_adj < fdr_thresh & abs(avg_log2FC) >= log2fc_large]
deg_sig[, direction := fifelse(avg_log2FC > 0, "up", "down")]
 
# ── Ortholog mapping: rat symbol → human Ensembl ID ──────────────────────────
 

library(orthogene)

gene_map <- orthogene::map_orthologs(
  genes          = rat_symbols,
  input_species  = "rat",
  output_species = "human",
  method         = "homologene"
)
ortho <- as.data.table(gene_map)
setnames(ortho, c("input_gene", "ortholog_gene"), c("rat_symbol", "human_symbol"))
ortho <- unique(ortho[!is.na(human_symbol) & human_symbol != ""])

# Map human symbols → Ensembl IDs using the atlas annotation
rd <- as.data.table(rowData(spe))
# Check your column names with: names(rd)
# Adjust "gene_name" and "gene_id" to match your rowData
symbol_to_ens <- unique(rd[, .(human_symbol = gene_name, ensembl_id = gene_id)])
ortho <- merge(ortho, symbol_to_ens, by = "human_symbol")
ortho <- unique(ortho[ensembl_id != ""])
 
# ── Build gene sets ──────────────────────────────────────────────────────────
 
deg_human <- merge(deg_sig, ortho, by.x = "gene", by.y = "rat_symbol",
                   allow.cartesian = TRUE)
 
make_set <- function(dt, ...) {
  ids <- unique(dt[..., ensembl_id])
  if (length(ids) >= 5) ids else NULL
}
 
gene_list <- Filter(Negate(is.null), list(
  Up_in_highAI   = make_set(deg_human, direction == "up"),
  Down_in_highAI = make_set(deg_human, direction == "down"),
  ExNeuron_DEGs  = make_set(deg_human, celltype == "ExNeuron"),
  InhNeuron_DEGs = make_set(deg_human, celltype == "InhNeuron"),
  ExNeuron_down  = make_set(deg_human, celltype == "ExNeuron" & direction == "down"),
  ExNeuron_up    = make_set(deg_human, celltype == "ExNeuron" & direction == "up")
))
 
# ── Enrichment (one domain vs rest) ──────────────────────────────────────────
 
enrich <- gene_set_enrichment(
  gene_list        = gene_list,
  modeling_results = modeling_results,
  model_type       = "enrichment",
  fdr_cut          = 0.1
)
fwrite(enrich, file.path(outdir, "enrichment_results.csv"))
 
pdf(file.path(outdir, "enrichment_plot.pdf"), width = 5, height = 6)
gene_set_enrichment_plot(
  enrichment       = enrich,
  xlabs            = names(gene_list),
  plot_SetSize_bar = TRUE,
  gene_list_length = lapply(gene_list, length)
)
dev.off()
 
# ── Enrichment reverse (genes depleted per domain) ───────────────────────────
 
enrich_rev <- gene_set_enrichment(
  gene_list        = gene_list,
  modeling_results = modeling_results,
  model_type       = "enrichment",
  fdr_cut          = 0.1,
  reverse          = TRUE
)
fwrite(enrich_rev, file.path(outdir, "enrichment_reverse_results.csv"))
 
pdf(file.path(outdir, "enrichment_reverse_plot.pdf"), width = 5, height = 6)
gene_set_enrichment_plot(
  enrichment       = enrich_rev,
  xlabs            = names(gene_list),
  plot_SetSize_bar = TRUE,
  gene_list_length = lapply(gene_list, length)
)
dev.off()
 
# ── Pairwise (domain A vs domain B) ─────────────────────────────────────────
 
if ("pairwise" %in% names(modeling_results)) {
  pairwise <- gene_set_enrichment(
    gene_list        = gene_list,
    modeling_results = modeling_results,
    model_type       = "pairwise",
    fdr_cut          = 0.1
  )
  fwrite(pairwise, file.path(outdir, "pairwise_results.csv"))
 

  p <- gene_set_enrichment_plot(
    enrichment       = pairwise,
    xlabs            = names(gene_list),
    plot_SetSize_bar = TRUE,
    gene_list_length = lapply(gene_list, length)
  )

}
 
pdf(file.path(outdir, "pairwise_plot.pdf"), width = 7, height = 12)
p
dev.off()

message("Done. Results in ", outdir)
 