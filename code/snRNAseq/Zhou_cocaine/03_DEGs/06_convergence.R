#!/usr/bin/env Rscript
# =============================================================================
# 06_convergence.R
# =============================================================================
# Identify spatial domains that emerge as significant from BOTH approaches:
#   - Top-ranked by composition-adjusted mean (Approach 1)
#   - Significant Fisher overlap (Approach 2)
#
# Domains that converge across both methods are the strongest candidates
# for harboring addiction-relevant transcriptional programs.
#
# Inputs:  results/approach1_emmeans.csv
#          results/approach1_lm_results.csv
#          results/approach2_fisher_results.csv
# Outputs: results/convergence_summary.csv
#          results/convergence_heatmap.pdf
# =============================================================================

source("00_config.R")

message("═══ 06: Convergence analysis ═══\n")


# ── 1. Load results from both approaches ────────────────────────────────────

lm_res    <- fread(file.path(out_dir, "approach1_lm_results.csv"))
emmeans   <- fread(file.path(out_dir, "approach1_emmeans.csv"))
fisher    <- fread(file.path(out_dir, "approach2_fisher_results.csv"))


# ── 2. Approach 1: rank domains by adjusted mean per gene set ────────────────

# For each gene set, rank domains from highest to lowest adjusted mean
emmeans[, rank := frank(-emmean), by = gene_set]

# Flag top 3 domains per gene set
a1_top <- emmeans[rank <= 3, .(gene_set, domain, emmean, SE, rank)]

message("Approach 1 — top 3 domains per gene set (by adjusted mean):")
print(a1_top[order(gene_set, rank)])


# ── 3. Approach 2: significant Fisher results ───────────────────────────────

a2_sig <- fisher[padj < 0.05, .(gene_set, domain,
                                 fisher_padj  = padj,
                                 odds_ratio,
                                 overlap)]

message(sprintf("\nApproach 2 — %d significant domain×set pairs (Fisher padj < 0.05)",
                nrow(a2_sig)))
if (nrow(a2_sig) > 0) print(a2_sig[order(fisher_padj)])


# ── 4. Find convergent hits ─────────────────────────────────────────────────

convergent <- merge(a1_top, a2_sig, by = c("gene_set", "domain"))

message(sprintf("\n═══ CONVERGENT DOMAINS: %d ═══", nrow(convergent)))

if (nrow(convergent) > 0) {
  convergent <- convergent[order(fisher_padj)]
  print(convergent)
  fwrite(convergent, file.path(out_dir, "convergence_summary.csv"))

  # Also report: which domains appear convergent for multiple gene sets?
  domain_counts <- convergent[, .N, by = domain][order(-N)]
  message("\nDomains convergent across multiple gene sets:")
  print(domain_counts)

} else {
  message("  No domains significant in both approaches.")
  message("  Consider:")
  message("    - Relaxing the emmeans rank threshold (top 5 instead of top 3)")
  message("    - Checking whether composition correction removed the signal")
  message("    - Looking at Approach 1 and 2 results independently")
  
  # Fallback: look for near-misses (top 5 + Fisher p < 0.1)
  a1_top5 <- emmeans[rank <= 5, .(gene_set, domain, emmean, rank)]
  a2_loose <- fisher[padj < 0.10, .(gene_set, domain,
                                      fisher_padj = padj, odds_ratio)]
  near_miss <- merge(a1_top5, a2_loose, by = c("gene_set", "domain"))
  if (nrow(near_miss) > 0) {
    message("\n  Near-misses (top 5 emmean + Fisher padj < 0.1):")
    print(near_miss[order(fisher_padj)])
    fwrite(near_miss, file.path(out_dir, "convergence_near_miss.csv"))
  }
}


# ── 5. Convergence heatmap ──────────────────────────────────────────────────

# Create a combined score: for each domain×set, compute
#   combined = -log10(fisher_padj) × (emmean rank <= 3)
# This highlights cells that are both highly ranked and Fisher-significant.

emm_ranked <- emmeans[, .(gene_set, domain, emmean,
                          rank = frank(-emmean)), by = gene_set]

combo <- merge(
  dcast(emm_ranked, domain ~ gene_set, value.var = "rank"),
  dcast(fisher, domain ~ gene_set, value.var = "padj"),
  by = "domain", suffixes = c("_rank", "_fisher")
)

# If we have results, make a simple combined visual
if (nrow(fisher) > 0 && nrow(emmeans) > 0) {
  message("\nGenerating convergence heatmap ...")
  
  # Build a matrix: -log10(fisher padj) for sets where domain is top 3
  gene_set_names <- unique(emmeans$gene_set)
  domain_names   <- unique(emmeans$domain)
  
  conv_mat <- matrix(0, nrow = length(gene_set_names),
                     ncol = length(domain_names),
                     dimnames = list(gene_set_names, domain_names))
  
  for (i in seq_len(nrow(emmeans))) {
    gs <- emmeans$gene_set[i]
    d  <- emmeans$domain[i]
    r  <- emmeans[gene_set == gs & domain == d, frank(-emmean), by = gene_set]$V1
    fp <- fisher[gene_set == gs & domain == d, padj]
    
    if (length(fp) > 0 && length(r) > 0) {
      # Score: -log10(p) weighted by rank (top domains get full weight)
      rank_weight <- ifelse(r <= 3, 1.0, ifelse(r <= 5, 0.5, 0.2))
      conv_mat[gs, d] <- min(-log10(fp), 10) * rank_weight
    }
  }
  
  pdf(file.path(out_dir, "convergence_heatmap.pdf"), width = 12, height = 6)
  ht <- Heatmap(
    conv_mat,
    name = "Convergence\nscore",
    column_title = "Spatial domain",
    row_title = "DEG gene set",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9),
    col = colorRamp2(c(0, 2, 5, 10),
                      c("grey95", "#FAEEDA", "#EF9F27", "#B2182B"))
  )
  draw(ht)
  dev.off()
}

message("\n✓ Done. Convergence results in ", out_dir)
