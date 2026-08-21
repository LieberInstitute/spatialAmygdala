library(here)
library(ggplot2)
library(dplyr)
library(SpatialExperiment)
library(scater)
library(scran)
library(spacexr)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")
plot_dir <- here("plots", "VisiumHD", "05_label_transfer", "RCTD_human_Totty")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_cells.rds"))

myRCTD <- readRDS(file.path(processed_dir, "rctd_results_HDcells_Totty.rds"))

# ---- Align RCTD results to spe ----
# spe colnames repeat across samples (same cellid reused per sample), so plain
# match()/intersect() collapse everything into one sample's worth of cells.
# results_df disambiguated its rownames with make.unique() (.1/.2/.3/.4 suffixes),
# so we rebuild the identical key on spe before matching.
rd <- myRCTD@results$results_df

colnames(spe) <- make.unique(colnames(spe))
common <- intersect(colnames(spe), rownames(rd))
message("Cells matched: ", length(common), " of ", ncol(spe))

spe <- spe[, common]
rd  <- rd[common, ]

colData(spe)$first_type  <- rd$first_type
colData(spe)$second_type <- rd$second_type
colData(spe)$spot_class  <- rd$spot_class

# Sanity checks
print(table(spe$spot_class))
print(table(spe$sample_id, spe$spot_class))

# ---- Build plotting data frame ----
plot_df <- as.data.frame(spatialCoords(spe))
plot_df <- cbind(
  plot_df,
  as.data.frame(colData(spe)[, c("sample_id", "first_type", "spot_class")])
)

# Grey out low-confidence assignments so confident cell types read clearly.
# To keep doublet_uncertain visible, use only "reject" here.
# To color every cell, delete this line entirely.
plot_df$first_type[plot_df$spot_class %in% c("reject", "doublet_uncertain")] <- NA

# ---- Single PDF, single plot, faceted by sample ----
p <- ggplot(plot_df,
            aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = first_type)) +
  geom_point(size = 0.15, stroke = 0) +
  facet_wrap(~ sample_id, nrow = 1, scales = "free") +
  scale_color_discrete(na.value = "grey85") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = "Cell type assignments (RCTD first_type)", color = "Cell type") +
  theme_void() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "right",
    strip.text      = element_text(size = 12),
    aspect.ratio    = 1
  )

pdf(file.path(plot_dir, "celltypes_all_samples_Totty.pdf"), width = 28, height = 6)
print(p)
dev.off()













# ---- Build plotting data frame ----
plot_df <- as.data.frame(spatialCoords(spe))
plot_df <- cbind(
  plot_df,
  as.data.frame(colData(spe)[, c("sample_id", "first_type", "spot_class")])
)
 
# CeA cell types to highlight
cea_types <- c("SST_TAC1", "PENK_DRD2", "DLK1_ZFHX3")
 
# Everything that isn't a highlighted CeA type -> NA (drawn grey underneath)
plot_df$first_type <- as.character(plot_df$first_type)
plot_df$highlight <- ifelse(plot_df$first_type %in% cea_types, plot_df$first_type, NA)
plot_df$is_cea <- !is.na(plot_df$highlight)
 
# ---- Plot: grey background layer first, highlighted CeA types on top ----
p <- ggplot() +
  # background: all non-CeA cells, grey with alpha
  geom_point(
    data = plot_df[!plot_df$is_cea, ],
    aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres),
    color = "grey80", alpha = 0.5, size = 0.4, stroke = 0
  ) +
  # foreground: highlighted CeA cell types, colored and opaque
  geom_point(
    data = plot_df[plot_df$is_cea, ],
    aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = highlight),
    size = 0.8, stroke = 0
  ) +
  facet_wrap(~ sample_id, nrow = 1, scales = "free") +
  scale_color_manual(values = c(
    SST_TAC1   = "#E41A1C",
    PENK_DRD2  = "#377EB8",
    DLK1_ZFHX3 = "#4DAF4A"
  ), breaks = cea_types) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Putative CeA cell types (RCTD first_type)", color = "Cell type") +
  theme_void() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "right",
    strip.text      = element_text(size = 12),
    aspect.ratio    = 1
  )
 
pdf(file.path(plot_dir, "CeA_celltypes_all_samples_Totty.pdf"), width = 28, height = 6)
print(p)
dev.off()