#!/usr/bin/env Rscript
# =============================================================================
# 02_score_atlas.R
# =============================================================================
# Load the human spatial amygdala atlas, estimate cell type composition
# per spot using canonical markers, and score each spot for each DEG gene set.
#
# Inputs:  spe_file          (from 00_config.R)
#          results/gene_sets.rds  (from 01)
# Outputs: results/spe_scored.rds          (atlas with scores in colData)
#          results/composition_scores.csv   (per-spot composition estimates)
#          results/deg_scores.csv           (per-spot DEG set scores)
# =============================================================================

source("00_config.R")


# ── 1. Load atlas ────────────────────────────────────────────────────────────

message("═══ 02: Loading atlas and scoring spots ═══\n")

if (!file.exists(spe_file)) {
  stop("Atlas not found: ", spe_file,
       "\nUpdate spe_file in 00_config.R")
}

spe <- readRDS(spe_file)

if (inherits(spe, "Seurat")) {
  message("Converting Seurat → SingleCellExperiment ...")
  spe <- Seurat::as.SingleCellExperiment(spe)
}

stopifnot(domain_col %in% colnames(colData(spe)))
domains <- factor(colData(spe)[[domain_col]])
atlas_genes <- rownames(spe)

message(sprintf("Atlas: %d spots, %d genes, %d domains (%s)\n",
                ncol(spe), nrow(spe), nlevels(domains),
                paste(levels(domains), collapse = ", ")))


# ── 2. Load gene sets and filter to atlas genes ─────────────────────────────

gene_sets <- readRDS(file.path(out_dir, "gene_sets.rds"))

message("Filtering gene sets to atlas features:")
for (nm in names(gene_sets)) {
  n_before <- length(gene_sets[[nm]])
  gene_sets[[nm]] <- intersect(gene_sets[[nm]], atlas_genes)
  n_after <- length(gene_sets[[nm]])
  message(sprintf("  %s: %d → %d genes", nm, n_before, n_after))
  if (n_after < 5) {
    message(sprintf("    WARNING: <5 genes — dropping '%s'", nm))
    gene_sets[[nm]] <- NULL
  }
}

# Save the filtered version for downstream scripts
saveRDS(gene_sets, file.path(out_dir, "gene_sets_filtered.rds"))


# ── 3. Composition estimation ────────────────────────────────────────────────

message("\nEstimating cell type composition per spot ...")

mat <- as.matrix(logcounts(spe))

# Filter composition markers to those present
comp_filt <- lapply(composition_markers, function(g) intersect(g, atlas_genes))
comp_filt <- comp_filt[sapply(comp_filt, length) >= 2]

for (nm in names(comp_filt)) {
  message(sprintf("  %s: %d markers → %s",
                  nm, length(comp_filt[[nm]]),
                  paste(comp_filt[[nm]], collapse = ", ")))
}

comp_scores <- ScoreSignatures_UCell(mat, features = comp_filt, maxRank = 1500)
colnames(comp_scores) <- sub("_UCell$", "", colnames(comp_scores))

# Save
comp_df <- as.data.frame(comp_scores)
comp_df$spot_id <- colnames(spe)
comp_df$domain  <- as.character(domains)
fwrite(as.data.table(comp_df), file.path(out_dir, "composition_scores.csv"))

# Add to colData
for (s in colnames(comp_scores)) {
  colData(spe)[[paste0("comp_", s)]] <- comp_scores[, s]
}


# ── 4. DEG gene set scoring ─────────────────────────────────────────────────

message("\nScoring spots for DEG gene sets ...")

deg_scores <- ScoreSignatures_UCell(mat, features = gene_sets, maxRank = 1500)
colnames(deg_scores) <- sub("_UCell$", "", colnames(deg_scores))

# Save
deg_df <- as.data.frame(deg_scores)
deg_df$spot_id <- colnames(spe)
deg_df$domain  <- as.character(domains)
fwrite(as.data.table(deg_df), file.path(out_dir, "deg_scores.csv"))

# Add to colData
for (s in colnames(deg_scores)) {
  colData(spe)[[s]] <- deg_scores[, s]
}


# ── 5. Save scored atlas ────────────────────────────────────────────────────

message("\nSaving scored atlas ...")
saveRDS(spe, file.path(out_dir, "spe_scored.rds"))

message(sprintf("\n✓ Done. Scored %d gene sets × %d spots",
                ncol(deg_scores), nrow(deg_scores)))
message("  Composition covariates: ", paste(colnames(comp_scores), collapse = ", "))
message("  Outputs in ", out_dir)
