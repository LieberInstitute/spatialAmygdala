## =============================================================================
## CeA division abundance across anterior-posterior (AP) levels
## abundance = per-sample fraction of ALL spots called CeC / CeL / CeM
## points = individual samples; lines = mean proportion per division across AP
## Expects spe$region (CeC/CeL/CeM/...) and spe$ap_level / spe$ap_short.
## =============================================================================

library(here)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)

processed_dir <- here("processed-data", "Visium", "09_deconvolution")
plot_dir      <- here("plots", "Visium", "09_deconvolution", "RCTD_human")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

region_levels <- c("CeC", "CeL", "CeM")
region_cols   <- c(CeC = "#D55E00", CeL = "#0072B2", CeM = "#009E73")  # Okabe-Ito
ap_order      <- c("Anterior-intermediate", "Intermediate",
                   "Posterior-intermediate", "Posterior")

## ---- 0. load labeled object & add AP level ---------------------------------
spe <- readRDS(file.path(processed_dir, "spe_with_region_labels.rds"))

ap_level <- c(Br9192 = "Anterior-intermediate", Br9280 = "Intermediate",
              Br2743 = "Intermediate",          Br6471 = "Intermediate",
              Br8325 = "Posterior-intermediate", Br6423 = "Posterior-intermediate",
              Br6660 = "Posterior")
ap_short <- c("Anterior-intermediate" = "AI", "Intermediate" = "I",
              "Posterior-intermediate" = "PI", "Posterior" = "P")
sid <- as.character(spe$sample_id)
stopifnot(all(sid %in% names(ap_level)))
spe$ap_level <- factor(ap_level[sid], levels = ap_order)
spe$ap_short <- factor(ap_short[as.character(spe$ap_level)],
                       levels = c("AI", "I", "PI", "P"))

## ---- 1. per-sample proportion of each division -----------------------------
cd <- as.data.frame(colData(spe))

# total spots per sample (denominator = ALL spots), and AP level per sample
sample_tot <- cd %>%
  group_by(sample_id) %>%
  summarise(n_total = n(),
            ap_level = first(ap_level),
            ap_short = first(ap_short), .groups = "drop")

# count of each CeA division per sample, completed with 0s where absent
div_counts <- cd %>%
  filter(region %in% region_levels) %>%
  count(sample_id, region, name = "n_region")

abund <- tidyr::expand_grid(sample_id = sample_tot$sample_id,
                            region    = factor(region_levels, levels = region_levels)) %>%
  left_join(div_counts, by = c("sample_id", "region")) %>%
  mutate(n_region = tidyr::replace_na(n_region, 0)) %>%
  left_join(sample_tot, by = "sample_id") %>%
  mutate(prop = n_region / n_total,
         region = factor(region, levels = region_levels),
         ap_level = factor(ap_level, levels = ap_order))

## ---- 2. mean proportion per AP level x division (for the lines) ------------
ap_mean <- abund %>%
  group_by(ap_level, region) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

## ---- 3. plot: abundance (x) vs AP level (y), lines per division ------------
# AP level on y running anterior (top) -> posterior (bottom)
abund$ap_level   <- factor(abund$ap_level,   levels = rev(ap_order))
ap_mean$ap_level <- factor(ap_mean$ap_level, levels = rev(ap_order))

p <- ggplot(abund, aes(prop, ap_level, color = region)) +
  # mean line per division across AP levels
  geom_path(data = ap_mean, aes(mean_prop, ap_level, group = region),
            linewidth = 1) +
  # mean marker
  geom_point(data = ap_mean, aes(mean_prop, ap_level),
             size = 5, shape = 18) +
  # individual samples (uniform size), jittered to separate overlaps
  geom_point(alpha = 0.6, size = 1.5,
             position = position_jitter(height = 0.12, width = 0)) +
  scale_color_manual(values = region_cols, name = "CeA division") +
  labs(x = "Proportion of spots in sample", y = "AP level",
       title = "CeA division abundance across anterior-posterior axis") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "CeA_abundance_by_AP.pdf"), p, width = 4.5, height = 5)
cat("Wrote:", file.path(plot_dir, "CeA_abundance_by_AP.pdf"), "\n")