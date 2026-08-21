#!/usr/bin/env Rscript
# =============================================================================
# 05_spatial_plots.R
# =============================================================================
# Generate spatial maps of DEG set scores and composition estimates.
# Only meaningful for SpatialExperiment objects with spatial coordinates.
#
# Inputs:  results/spe_scored.rds   (from 02)
# Outputs: results/spatial_deg_scores.pdf
#          results/spatial_composition.pdf
#          results/spatial_domains.pdf
# =============================================================================

source("00_config.R")

message("═══ 05: Spatial visualization ═══\n")


# ── 1. Load scored atlas ────────────────────────────────────────────────────

spe <- readRDS(file.path(out_dir, "spe_scored.rds"))

if (!is(spe, "SpatialExperiment")) {
  message("Object is not a SpatialExperiment — skipping spatial plots.")
  message("(Boxplots and heatmaps from scripts 03/04 still apply.)")
  quit(save = "no")
}

coords  <- spatialCoords(spe)
domains <- colData(spe)[[domain_col]]

# Detect which score columns are present
all_cols   <- colnames(colData(spe))
deg_cols   <- intersect(all_cols,
                        c("up_in_highAI", "down_in_highAI",
                          "ExNeuron_DEGs", "InhNeuron_DEGs",
                          "ExNeuron_down", "ExNeuron_up"))
comp_cols  <- grep("^comp_", all_cols, value = TRUE)

message(sprintf("  %d spots, %d DEG score columns, %d composition columns",
                ncol(spe), length(deg_cols), length(comp_cols)))

# Handle multi-sample atlases: if there's a sample_id column, facet by it
has_samples <- "sample_id" %in% colnames(colData(spe))


# ── Helper: spatial plot ─────────────────────────────────────────────────────

make_spatial_plot <- function(spe, col_name, title, palette = "inferno") {
  df <- data.frame(
    x     = coords[, 1],
    y     = coords[, 2],
    value = colData(spe)[[col_name]]
  )
  if (has_samples) df$sample <- colData(spe)$sample_id

  p <- ggplot(df, aes(x, y, color = value)) +
    geom_point(size = 0.3, stroke = 0) +
    scale_color_viridis_c(option = palette) +
    coord_fixed() +
    theme_void(base_size = 9) +
    labs(title = title, color = col_name)

  if (has_samples) p <- p + facet_wrap(~sample)
  p
}


# ── 2. Spatial domain map ───────────────────────────────────────────────────

message("  Plotting spatial domains ...")

df_dom <- data.frame(
  x = coords[, 1], y = coords[, 2],
  domain = domains
)
if (has_samples) df_dom$sample <- colData(spe)$sample_id

p_dom <- ggplot(df_dom, aes(x, y, color = domain)) +
  geom_point(size = 0.3, stroke = 0) +
  coord_fixed() +
  theme_void(base_size = 10) +
  labs(title = "Spatial domains", color = domain_col)
if (has_samples) p_dom <- p_dom + facet_wrap(~sample)

ggsave(file.path(out_dir, "spatial_domains.pdf"), p_dom,
       width = 10, height = 8)


# ── 3. DEG score spatial maps ───────────────────────────────────────────────

if (length(deg_cols) > 0) {
  message("  Plotting DEG set scores ...")

  plots <- lapply(deg_cols, function(s) {
    make_spatial_plot(spe, s, s, palette = "inferno")
  })

  p_deg <- wrap_plots(plots, ncol = 2)
  ggsave(file.path(out_dir, "spatial_deg_scores.pdf"), p_deg,
         width = 14, height = 5 * ceiling(length(deg_cols) / 2))
}


# ── 4. Composition spatial maps ─────────────────────────────────────────────

if (length(comp_cols) > 0) {
  message("  Plotting composition scores ...")

  plots_c <- lapply(comp_cols, function(s) {
    make_spatial_plot(spe, s, sub("comp_", "", s), palette = "viridis")
  })

  p_comp <- wrap_plots(plots_c, ncol = 2)
  ggsave(file.path(out_dir, "spatial_composition.pdf"), p_comp,
         width = 14, height = 5 * ceiling(length(comp_cols) / 2))
}

message("\n✓ Done. Spatial plots in ", out_dir)
