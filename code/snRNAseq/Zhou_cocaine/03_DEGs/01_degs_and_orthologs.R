#!/usr/bin/env Rscript
# =============================================================================
# 01_degs_and_orthologs.R
# =============================================================================
# Read rat DEGs, filter to large-effect significant hits, map to human
# orthologs via biomaRt, and build the core gene sets.
#
# Inputs:  deg_file (from 00_config.R)
# Outputs: results/ortholog_map.csv
#          results/gene_sets.rds        (named list of human gene vectors)
#          results/gene_sets_summary.csv
#          results/deg_filtered.csv     (filtered DEG table with human orthologs)
# =============================================================================

source("00_config.R")


# ── 1. Read and filter ───────────────────────────────────────────────────────

message("═══ 01: Processing DEGs and mapping orthologs ═══\n")

deg <- fread(deg_file)
message(sprintf("Raw table: %d rows, %d genes, %d cell types",
                nrow(deg), uniqueN(deg$gene), uniqueN(deg$celltype)))

deg_sig <- deg[p_val_adj < fdr_thresh & abs(avg_log2FC) >= log2fc_large]
deg_sig[, direction := fifelse(avg_log2FC > 0, "up", "down")]

message(sprintf("After filtering (FDR < %g, |log2FC| >= %g): %d rows, %d unique genes\n",
                fdr_thresh, log2fc_large, nrow(deg_sig), uniqueN(deg_sig$gene)))

cat("Per cell type:\n")
print(deg_sig[, .(n_up = sum(direction == "up"),
                  n_down = sum(direction == "down")), by = celltype])


# ── 2. Ortholog mapping ─────────────────────────────────────────────────────

rat_genes <- unique(deg_sig$gene)
mappable  <- !grepl("^AABR07|^LOC[0-9]|^NEWGENE|^RGD[0-9]", rat_genes)
rat_symbols <- rat_genes[mappable]
message(sprintf("\n%d standard symbols to map, %d unmappable IDs dropped",
                sum(mappable), sum(!mappable)))

ortho <- tryCatch({
  rat_mart   <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  raw <- getLDS(
    attributes  = "rgd_symbol",
    filters     = "rgd_symbol",
    values      = rat_symbols,
    mart        = rat_mart,
    attributesL = "hgnc_symbol",
    martL       = human_mart,
    uniqueRows  = TRUE
  )
  out <- as.data.table(raw)
  setnames(out, c("rat_gene", "human_gene"))
  out <- unique(out[human_gene != ""])
  message(sprintf("biomaRt: %d rat → %d human symbols",
                  uniqueN(out$rat_gene), uniqueN(out$human_gene)))
  out
}, error = function(e) {
  message("biomaRt unavailable — using toupper() fallback")
  message("  Error: ", e$message)
  data.table(rat_gene = rat_symbols, human_gene = toupper(rat_symbols))
})

fwrite(ortho, file.path(out_dir, "ortholog_map.csv"))


# ── 3. Build gene sets ──────────────────────────────────────────────────────

deg_human <- merge(deg_sig, ortho, by.x = "gene", by.y = "rat_gene",
                   allow.cartesian = TRUE)

fwrite(deg_human, file.path(out_dir, "deg_filtered.csv"))

gene_sets <- list()

# A: Up in high AI — all cell types pooled
gs <- unique(deg_human[direction == "up", human_gene])
if (length(gs) >= 5) gene_sets[["up_in_highAI"]] <- gs

# B: Down in high AI — all cell types pooled
gs <- unique(deg_human[direction == "down", human_gene])
if (length(gs) >= 5) gene_sets[["down_in_highAI"]] <- gs

# C: Excitatory neuron DEGs (both directions)
gs <- unique(deg_human[celltype == "ExNeuron", human_gene])
if (length(gs) >= 5) gene_sets[["ExNeuron_DEGs"]] <- gs

# D: Inhibitory neuron DEGs (both directions)
gs <- unique(deg_human[celltype == "InhNeuron", human_gene])
if (length(gs) >= 5) gene_sets[["InhNeuron_DEGs"]] <- gs

# E: Excitatory — down only (the dominant direction, ~256 genes)
gs <- unique(deg_human[celltype == "ExNeuron" & direction == "down", human_gene])
if (length(gs) >= 5) gene_sets[["ExNeuron_down"]] <- gs

# F: Excitatory — up only
gs <- unique(deg_human[celltype == "ExNeuron" & direction == "up", human_gene])
if (length(gs) >= 5) gene_sets[["ExNeuron_up"]] <- gs

saveRDS(gene_sets, file.path(out_dir, "gene_sets.rds"))

set_info <- data.table(
  set_name = names(gene_sets),
  n_genes  = sapply(gene_sets, length)
)
fwrite(set_info, file.path(out_dir, "gene_sets_summary.csv"))

message("\nGene sets created:")
print(set_info)
message("\n✓ Done. Outputs in ", out_dir)
