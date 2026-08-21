#!/usr/bin/env Rscript
# =============================================================================
# 04_approach2_fisher.R
# =============================================================================
# Identify marker genes for each spatial domain, then test for over-
# representation of DEG gene sets among domain markers using Fisher's
# exact test.
#
# Inputs:  results/spe_scored.rds       (from 02)
#          results/gene_sets_filtered.rds (from 02)
#          results/ortholog_map.csv      (from 01)
# Outputs: results/domain_markers.rds
#          results/approach2_fisher_results.csv
#          results/approach2_overlap_genes.csv
#          results/approach2_heatmap.pdf
# =============================================================================

source("00_config.R")

message("═══ 04: Approach 2 — Fisher's exact test ═══\n")


# ── 1. Load data ────────────────────────────────────────────────────────────

spe       <- readRDS(file.path(out_dir, "spe_scored.rds"))
gene_sets <- readRDS(file.path(out_dir, "gene_sets_filtered.rds"))
ortho     <- fread(file.path(out_dir, "ortholog_map.csv"))

domains     <- factor(colData(spe)[[domain_col]])
atlas_genes <- rownames(spe)


# ── 2. Identify domain marker genes ─────────────────────────────────────────

message("Identifying domain markers (findMarkers, t-test, up only) ...")

spe$domain_factor <- domains

markers <- findMarkers(spe, groups = spe$domain_factor,
                       test.type = "t", pval.type = "any",
                       direction = "up", assay.type = "logcounts")

domain_markers <- list()
for (d in names(markers)) {
  mk <- as.data.frame(markers[[d]])
  top <- rownames(mk)[mk$FDR < 0.05 & mk$summary.logFC > 0]
  domain_markers[[d]] <- head(top, 500)
  message(sprintf("  %s: %d markers", d, length(domain_markers[[d]])))
}

saveRDS(domain_markers, file.path(out_dir, "domain_markers.rds"))


# ── 3. Fisher's exact tests ─────────────────────────────────────────────────

message("\nRunning Fisher's exact tests ...")

# Universe: genes present in atlas AND in the ortholog map
universe <- intersect(atlas_genes, unique(ortho$human_gene))
message(sprintf("  Universe: %d genes\n", length(universe)))

fisher_list   <- list()
overlap_list  <- list()

for (set_name in names(gene_sets)) {
  deg_in_u <- intersect(gene_sets[[set_name]], universe)

  for (d in names(domain_markers)) {
    mk_in_u <- intersect(domain_markers[[d]], universe)

    ol   <- intersect(deg_in_u, mk_in_u)
    a    <- length(ol)
    b    <- length(setdiff(deg_in_u, mk_in_u))
    c_   <- length(setdiff(mk_in_u, deg_in_u))
    d_   <- length(universe) - a - b - c_

    ft <- fisher.test(matrix(c(a, b, c_, d_), 2), alternative = "greater")

    fisher_list[[paste0(set_name, "||", d)]] <- data.table(
      gene_set   = set_name,
      domain     = d,
      overlap    = a,
      n_deg      = length(deg_in_u),
      n_markers  = length(mk_in_u),
      universe_n = length(universe),
      odds_ratio = ft$estimate,
      pval       = ft$p.value
    )

    if (a > 0) {
      overlap_list[[paste0(set_name, "||", d)]] <- data.table(
        gene_set = set_name,
        domain   = d,
        genes    = paste(sort(ol), collapse = ", "),
        n        = a
      )
    }
  }
}

fisher_all <- rbindlist(fisher_list)
fisher_all[, padj := p.adjust(pval, method = "BH")]
fisher_all <- fisher_all[order(pval)]
fwrite(fisher_all, file.path(out_dir, "approach2_fisher_results.csv"))

if (length(overlap_list) > 0) {
  fwrite(rbindlist(overlap_list),
         file.path(out_dir, "approach2_overlap_genes.csv"))
}

message("Significant results (padj < 0.05):")
sig <- fisher_all[padj < 0.05]
if (nrow(sig) > 0) {
  print(sig)
} else {
  message("  (none)")
}


# ── 4. Heatmap ──────────────────────────────────────────────────────────────

message("\nGenerating heatmap ...")

pmat <- dcast(fisher_all, gene_set ~ domain, value.var = "pval")
pm   <- as.matrix(pmat[, -1])
rownames(pm) <- pmat$gene_set

log10p <- pmin(-log10(pm), 10)
log10p[is.na(log10p)] <- 0

sig_stars <- ifelse(pm < 0.001, "***",
             ifelse(pm < 0.01,  "**",
             ifelse(pm < 0.05,  "*", "")))

pdf(file.path(out_dir, "approach2_heatmap.pdf"), width = 12, height = 6)
ht <- Heatmap(
  log10p,
  name = "-log10(p)",
  column_title = "Spatial domain",
  row_title = "DEG gene set",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  col = colorRamp2(c(0, 1.3, 5, 10),
                    c("grey95", "lightyellow", "orange", "red3")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sig_stars[i, j], x, y, gp = gpar(fontsize = 7))
  }
)
draw(ht)
dev.off()

message("\n✓ Done. Results in ", out_dir)
