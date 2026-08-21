#!/usr/bin/env Rscript
# =============================================================================
# 03_approach1_linear_model.R
# =============================================================================
# Fit composition-corrected linear models to test whether DEG set scores
# differ across spatial domains after accounting for cell type composition.
#
#   Model:  UCell_score ~ domain + excitatory + inhibitory + astrocyte + ...
#   Test:   ANOVA comparing full model vs composition-only null model
#   Output: estimated marginal means per domain (adjusted for composition)
#
# Inputs:  results/composition_scores.csv  (from 02)
#          results/deg_scores.csv          (from 02)
# Outputs: results/approach1_lm_results.csv
#          results/approach1_emmeans.csv
#          results/approach1_pairwise.csv
#          results/approach1_heatmap.pdf
#          results/approach1_boxplots.pdf
# =============================================================================

source("00_config.R")

message("═══ 03: Approach 1 — Composition-corrected linear models ═══\n")


# ── 1. Load scores ──────────────────────────────────────────────────────────

comp <- fread(file.path(out_dir, "composition_scores.csv"))
degs <- fread(file.path(out_dir, "deg_scores.csv"))

# Merge into single model data frame
model_df <- merge(comp, degs, by = c("spot_id", "domain"))

# Identify columns
comp_vars <- setdiff(names(comp), c("spot_id", "domain"))
deg_vars  <- setdiff(names(degs), c("spot_id", "domain"))

model_df$domain <- factor(model_df$domain)

message(sprintf("  %d spots, %d domains", nrow(model_df), nlevels(model_df$domain)))
message(sprintf("  Composition covariates: %s", paste(comp_vars, collapse = ", ")))
message(sprintf("  Gene sets: %s\n", paste(deg_vars, collapse = ", ")))


# ── 2. Fit linear models ────────────────────────────────────────────────────

comp_rhs <- paste0("`", comp_vars, "`", collapse = " + ")

lm_results_list    <- list()
emmeans_list       <- list()
pairwise_list      <- list()

for (s in deg_vars) {

  f_full <- as.formula(paste0("`", s, "` ~ domain + ", comp_rhs))
  f_null <- as.formula(paste0("`", s, "` ~ ", comp_rhs))

  fit_full <- lm(f_full, data = model_df)
  fit_null <- lm(f_null, data = model_df)

  # F-test: does domain improve the model beyond composition?
  av <- anova(fit_null, fit_full)

  r2_full <- summary(fit_full)$r.squared
  r2_null <- summary(fit_null)$r.squared

  lm_results_list[[s]] <- data.table(
    gene_set  = s,
    f_stat    = av$F[2],
    df1       = av$Df[2],
    df2       = av$Res.Df[2],
    pval      = av$`Pr(>F)`[2],
    r2_full   = r2_full,
    r2_null   = r2_null,
    r2_domain = r2_full - r2_null
  )

  # Estimated marginal means (composition-adjusted domain means)
  emm <- emmeans(fit_full, ~ domain)
  emm_dt <- as.data.table(summary(emm))
  emm_dt$gene_set <- s
  emmeans_list[[s]] <- emm_dt

  # Pairwise domain contrasts
  pw <- as.data.table(pairs(emm, adjust = "BH"))
  pw$gene_set <- s
  pairwise_list[[s]] <- pw

  message(sprintf("  %s: F=%.1f, p=%.2e, R²(domain)=%.4f",
                  s, av$F[2], av$`Pr(>F)`[2], r2_full - r2_null))
}


# ── 3. Compile and save ─────────────────────────────────────────────────────

lm_results <- rbindlist(lm_results_list)
lm_results[, padj := p.adjust(pval, method = "BH")]
lm_results <- lm_results[order(pval)]
fwrite(lm_results, file.path(out_dir, "approach1_lm_results.csv"))

emmeans_all <- rbindlist(emmeans_list)
fwrite(emmeans_all, file.path(out_dir, "approach1_emmeans.csv"))

pairwise_all <- rbindlist(pairwise_list)
fwrite(pairwise_all, file.path(out_dir, "approach1_pairwise.csv"))

message("\nModel results:")
print(lm_results[, .(gene_set,
                      f_stat = round(f_stat, 1),
                      pval   = signif(pval, 3),
                      padj   = signif(padj, 3),
                      r2_domain = round(r2_domain, 4))])


# ── 4. Heatmap: adjusted means by domain ────────────────────────────────────

message("\nGenerating heatmap of adjusted means ...")

emm_wide <- dcast(emmeans_all, domain ~ gene_set, value.var = "emmean")
dm <- as.matrix(emm_wide[, -1])
rownames(dm) <- emm_wide$domain
dm_z <- scale(dm)

pdf(file.path(out_dir, "approach1_heatmap.pdf"), width = 10, height = 7)
ht <- Heatmap(
  t(dm_z),
  name = "Adjusted\nz-score",
  column_title = "Spatial domain",
  row_title = "DEG gene set",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  heatmap_legend_param = list(direction = "vertical")
)
draw(ht, heatmap_legend_side = "right")
dev.off()


# ── 5. Boxplots of raw scores ───────────────────────────────────────────────

message("Generating boxplots ...")

plot_df <- melt(as.data.table(model_df),
                id.vars = c("spot_id", "domain", comp_vars),
                measure.vars = deg_vars,
                variable.name = "gene_set", value.name = "score")

p <- ggplot(plot_df,
            aes(x = reorder(domain, score, FUN = median),
                y = score, fill = domain)) +
  geom_boxplot(outlier.size = 0.2, lwd = 0.3) +
  facet_wrap(~gene_set, scales = "free_y", ncol = 2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none",
        strip.text = element_text(size = 9, face = "bold")) +
  labs(x = "Spatial domain",
       y = "UCell score (raw)",
       title = "DEG set scores by domain (uncorrected for composition)",
       subtitle = "See heatmap / lm_results for composition-adjusted values")

ggsave(file.path(out_dir, "approach1_boxplots.pdf"), p,
       width = 12, height = 10)

message("\n✓ Done. Results in ", out_dir)
