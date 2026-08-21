library(here)
library(SpatialExperiment)
library(hoodscanR)
library(scico)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spdep)

plot_dir      <- here("plots", "VisiumHD", "05_label_transfer")
processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")

# load
spe <- readRDS(here(processed_dir, "spe_with_rctd_ITCmn.rds"))

# --- Filter: drop rejects AND low-confidence calls -------------------------
weight_cutoff <- 0.65

spe_nb <- spe[, spe$spot_class != "reject"]
message("Cells after dropping rejects: ", ncol(spe_nb))

# Capture the full cell-type universe BEFORE the confidence filter, so any
# subtype that's entirely low-confidence still gets a (zeroed) column in pm.
original_types <- sort(levels(droplevels(factor(spe_nb$first_type))))

hc_keep <- spe_nb$first_type_weight > weight_cutoff
hc_keep[is.na(hc_keep)] <- FALSE
message(sprintf("Keeping %d / %d cells with first_type_weight > %.2f (%.1f%%)",
                sum(hc_keep), length(hc_keep), weight_cutoff,
                100 * mean(hc_keep)))

# Per-type retention
before <- table(spe_nb$first_type)
spe_nb <- spe_nb[, hc_keep]
spe_nb$first_type <- droplevels(factor(spe_nb$first_type))
after  <- table(spe_nb$first_type)
print(data.frame(type     = names(before),
                 before   = as.integer(before),
                 after    = as.integer(after[names(before)]),
                 pct_kept = round(100 * as.integer(after[names(before)]) /
                                  as.integer(before), 1)))

# --- Parameters ------------------------------------------------------------
k_nn      <- 100
itc_cols  <- c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")
itc_types <- names(itc_cols)
all_types <- original_types
samples   <- sort(unique(as.character(spe_nb$sample_id)))

# --- Run hoodscanR per sample ----------------------------------------------
pm_list <- lapply(samples, function(s) {
    message("hoodscanR: ", s)
    sqe <- spe_nb[, spe_nb$sample_id == s]
    sqe <- readHoodData(sqe, anno_col = "first_type")
    nbs <- findNearCells(sqe, k = k_nn)
    pm  <- scanHoods(nbs$distance)
    pm  <- mergeByGroup(pm, nbs$cells)
    missing <- setdiff(all_types, colnames(pm))
    if (length(missing)) {
        pad <- matrix(0, nrow = nrow(pm), ncol = length(missing),
                      dimnames = list(rownames(pm), missing))
        pm <- cbind(pm, pad)
    }
    pm[, all_types, drop = FALSE]
})
names(pm_list) <- samples

pm_all <- do.call(rbind, pm_list)
message("Pooled probability matrix: ", nrow(pm_all), " cells × ",
        ncol(pm_all), " cell types")

# --- Full colocalization matrix -------------------------------------------
coloc <- cor(pm_all, use = "pairwise.complete.obs")

# --- Plot: full clustered colocalization heatmap --------------------------
cmax  <- max(abs(coloc[upper.tri(coloc)]))
col_fun <- colorRamp2(c(-cmax, 0, cmax), c("#3260a8", "white", "#a83232"))

type_anno <- ifelse(all_types %in% itc_types, all_types, "")
row_anno <- rowAnnotation(
    ITC = anno_simple(type_anno,
                      col = c(setNames(unname(itc_cols), itc_types),
                              setNames("white", "")),
                      border = FALSE),
    show_annotation_name = FALSE
)

h_full <- Heatmap(coloc,
    name = "colocalization\n(Pearson r)",
    col = col_fun,
    clustering_method_rows = "average",
    clustering_method_columns = "average",
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_names_side = "left", column_names_rot = 45,
    rect_gp = gpar(col = "grey85", lwd = 0.3),
    left_annotation = row_anno,
    column_title = sprintf("hoodscanR colocalization (k=%d, weight > %.2f, %d samples)",
                           k_nn, weight_cutoff, length(samples)))

pdf(here(plot_dir,
         sprintf("hoodscanR_colocalization_full_w%.2f.pdf", weight_cutoff)),
    width  = max(8, length(all_types) * 0.32 + 3),
    height = max(8, length(all_types) * 0.32 + 2))
draw(h_full); dev.off()

# --- Collapse Oligo subtypes ---------------------------------------------
print(all_types)
oligo_pattern <- "Oligo|^OL"
oligo_types <- grep(oligo_pattern, all_types, value = TRUE, ignore.case = TRUE)
stopifnot(length(oligo_types) > 0)
message("Collapsing into 'Oligo': ", paste(oligo_types, collapse = ", "))

pm_c <- pm_all
pm_c <- cbind(pm_c, Oligo = rowSums(pm_all[, oligo_types, drop = FALSE]))
pm_c <- pm_c[, !colnames(pm_c) %in% oligo_types]

collapsed_type <- as.character(spe_nb$first_type)
collapsed_type[collapsed_type %in% oligo_types] <- "Oligo"
spe_nb$collapsed_type <- factor(collapsed_type)

focal_types <- c("Oligo", itc_types)
focal_cols  <- c("Oligo" = "#7c5295", itc_cols)

# --- Focused colocalization: Oligo + ITC -----------------------------------
coloc_c     <- cor(pm_c, use = "pairwise.complete.obs")
coloc_focal <- coloc_c[focal_types, focal_types]

coloc_offdiag <- coloc_focal
diag(coloc_offdiag) <- NA
cmax <- max(abs(coloc_offdiag), na.rm = TRUE)
col_fun <- colorRamp2(c(-.1, 0, .1), c("#3260a8", "white", "#a83232"))

h_focal <- Heatmap(coloc_focal,
    name = "Pearson r",
    col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, w, h, fill) {
        if (i == j) {
            grid.rect(x, y, w, h, gp = gpar(fill = "black", col = "grey60"))
        }
    },
    row_names_side = "left", column_names_rot = 45,
    rect_gp = gpar(col = "grey60", lwd = 0.5),
    column_title = "Colocalization: ITC subtypes × Oligo")

pdf(here(plot_dir,
         sprintf("hoodscanR_colocalization_ITC_Oligo_w%.2f.pdf", weight_cutoff)),
    width = 5, height = 4)
draw(h_focal); dev.off()

# ==== Density analysis ====================================================
library(spatstat.geom)
library(spatstat.explore)

focal_types <- c("ITC_1", "ITC_2", "ITC_3")
focal_cols  <- c("ITC_1" = "#fd0d00",
                 "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")
samples     <- sort(unique(as.character(spe_nb$sample_id)))

# Build ppp from already-filtered spe_nb
build_ppp <- function(s) {
    keep <- spe_nb$sample_id == s &
            as.character(spe_nb$collapsed_type) %in%
                c(focal_types, "Oligo")
    cc   <- spatialCoords(spe_nb)[keep, , drop = FALSE]
    mks  <- factor(as.character(spe_nb$collapsed_type[keep]),
                   levels = c("Oligo", focal_types))
    ppp(x = cc[, 1], y = cc[, 2],
        window = ripras(cc[, 1], cc[, 2]),
        marks  = mks)
}
ppp_list <- lapply(samples, build_ppp); names(ppp_list) <- samples

r_unit  <- median(unlist(lapply(ppp_list, nndist)))
r_local <- 5 * r_unit

local_density <- do.call(rbind, lapply(samples, function(s) {
    p <- ppp_list[[s]]
    do.call(rbind, lapply(focal_types, function(t) {
        sub <- subset(p, marks == t)
        n   <- npoints(sub)
        if (n < 10) return(NULL)
        cp     <- closepairs(sub, rmax = r_local, what = "indices")
        counts <- tabulate(cp$i, nbins = n)
        data.frame(sample_id = s, cell_type = t,
                   density   = counts)
    }))
})) |>
    mutate(cell_type = factor(cell_type, levels = focal_types))

local_summary <- local_density |>
    group_by(sample_id, cell_type) |>
    summarise(median_density = median(density),
              n              = n(),
              .groups = "drop")

p <- ggplot(local_summary,
            aes(x = cell_type, y = median_density, fill = cell_type)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85,
                 colour = "grey20", linewidth = 0.3) +
    geom_jitter(width = 0.12, height = 0, size = 2.2, alpha = 0.85,
                shape = 21, colour = "grey20", stroke = 0.3) +
    scale_fill_manual(values = focal_cols, guide = "none") +
    labs(x = NULL,
         y = "Density (a.u.)",
         title = "Cell type density") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

pdf(here(plot_dir,
         sprintf("density_local_ITC_Oligo_boxplot_w%.2f_density.pdf",
                 weight_cutoff)),
    width = 3, height = 4)
print(p); dev.off()