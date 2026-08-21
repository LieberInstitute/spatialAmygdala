library(here)
library(ggplot2)
library(dplyr)
library(SpatialExperiment)
library(ggspavis)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")
plot_dir      <- here("plots", "Visium", "16_CeA_analyses")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# CeA cell-type colors (raw RCTD names; map to anatomical divisions:
#   DRD2 ISL1 -> CeC, DRD2 PAX6 -> CeL, PRKCD -> CeM)
cea_cols <- c("DRD2 ISL1" = "#D55E00",
              "DRD2 PAX6" = "#0072B2",
              "PRKCD"     = "#009E73")

# spatial domain palette
domain_colors <- c(
    AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6",
    LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
    MeA = "#baf739ff", CHAT = "#A0522D",
    Ventricle = "#666666", Vascular = "#333333",
    WM = "#BBBBBB", Neuropil = "#DDDDDD"
)

# --- Load and subset to Br9280 ---------------------------------------------
spe_cells <- readRDS(file.path(here("processed-data", "VisiumHD",
                                    "06_composition",
                                    "spe_cells_with_domain.rds")))
spe.subset <- spe_cells[, spe_cells$sample_id == "Br9280_CeA"]
spe.subset <- spe.subset[, !is.na(spe.subset$spatial_domain)]

# rotate -90 degrees for plotting
coords <- spatialCoords(spe.subset)
spatialCoords(spe.subset)[, 1] <- -coords[, 2]
spatialCoords(spe.subset)[, 2] <- -coords[, 1]

# === Plot 1: spatial domains, CeA domain full opacity ======================
spe.subset$is_cea_domain <- factor(
    as.character(spe.subset$spatial_domain) == "CeA",
    levels = c(FALSE, TRUE)
)
# reorder spots so non-CeA draw first, CeA on top (ggplot draws in row order)
spe.subset <- spe.subset[, order(spe.subset$is_cea_domain)]

p_domain <- ggspavis::plotCoords(
    spe.subset, annotate = "spatial_domain", in_tissue = NULL
) +
    aes(alpha = is_cea_domain) +
    scale_color_manual(values = domain_colors, drop = FALSE, breaks = "CeA") +
    scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), guide = "none")

pdf(here(plot_dir, "Br9280_domains.pdf"), width = 4, height = 4)
print(p_domain)
dev.off()

# === Plot 2: first_type cell labels, CeA subtypes full opacity =============
# full palette: CeA types get their colors, all other types grey
other_types <- setdiff(unique(as.character(spe.subset$first_type)), names(cea_cols))
first_type_cols <- c(
    setNames(rep("grey70", length(other_types)), other_types),
    cea_cols
)

# factor ordering: non-CeA first (drawn underneath in legend), CeA last
lev <- c(other_types, names(cea_cols))
spe.subset$first_type <- factor(as.character(spe.subset$first_type), levels = lev)

# alpha flag: CeA subtypes opaque, everything else faded
spe.subset$is_cea <- factor(
    as.character(spe.subset$first_type) %in% names(cea_cols),
    levels = c(FALSE, TRUE)
)

# reorder spots so non-CeA draw first and CeA draw on top
spe.subset <- spe.subset[, order(spe.subset$is_cea)]

p_celltype <- ggspavis::plotCoords(
    spe.subset, annotate = "first_type", in_tissue = NULL
) +
    aes(alpha = is_cea) +
    scale_color_manual(values = first_type_cols, drop = FALSE,
                       breaks = names(cea_cols)) +
    scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), guide = "none")

pdf(here(plot_dir, "Br9280_first_type_CeA.pdf"), width = 4, height = 4)
print(p_celltype)
dev.off()