library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggspavis)
library(spacexr)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")
plot_dir      <- here("plots", "VisiumHD", "05_label_transfer")

# --- Load spe --------------------------------------------------------------

spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control",
                    "spe_qc_cells.rds"))
spe

head(colnames(spe))
table(spe$sample_id)

# --- Load per-sample RCTD results ------------------------------------------

rctd_files <- list.files(
    processed_dir,
    pattern = "^rctd_results_HDcells_ITCmn_.+\\.rds$",
    full.names = TRUE
)
names(rctd_files) <- sub("^rctd_results_HDcells_ITCmn_(.+)\\.rds$", "\\1",
                         basename(rctd_files))

stopifnot(length(rctd_files) > 0)
message("Found ", length(rctd_files), " RCTD files: ",
        paste(names(rctd_files), collapse = ", "))

# RCTD script ran make.unique() on merged colnames before subsetting,
# leaving .1/.2 suffixes on per-sample barcodes. Strip them.
results_list <- lapply(names(rctd_files), function(sid) {
    rctd <- readRDS(rctd_files[[sid]])
    df   <- rctd@results$results_df

    bc <- sub("\\.\\d+$", "", rownames(df))
    stopifnot(!anyDuplicated(bc))

    df$sample_id <- sid
    df$barcode   <- bc

    w <- rctd@results$weights_doublet
    if (!is.null(w)) {
        df$first_type_weight  <- w[rownames(df), "first_type"]
        df$second_type_weight <- w[rownames(df), "second_type"]
    }
    rownames(df) <- NULL
    df
})
results_df <- do.call(rbind, results_list)

# Check cell type levels across samples
ct_levels <- lapply(rctd_files, function(f)
    levels(readRDS(f)@reference@cell_types))
if (!Reduce(identical, ct_levels)) {
    warning("Cell type levels differ across RCTD objects — check before pooling")
    cat("Union of cell types:\n");        print(Reduce(union,     ct_levels))
    cat("Intersection of cell types:\n"); print(Reduce(intersect, ct_levels))
}
celltypes <- Reduce(union, ct_levels)


# --- Subset spe and align results_df ---------------------------------------

print('Examining doublet mode results')

idx  <- match(spe_keys, results_df$key)
keep <- !is.na(idx)
spe         <- spe[, keep]
results_df  <- results_df[idx[keep], ]

colData(spe)$first_type        <- results_df$first_type
colData(spe)$second_type       <- results_df$second_type
colData(spe)$spot_class        <- results_df$spot_class
colData(spe)$first_type_weight <- results_df$first_type_weight

# --- Plot ITC subtypes per sample ------------------------------------------

itc_cols <- c("ITC_1" = "#fd0d00",
              "ITC_2" = "#f5b6b3",
              "ITC_3" = "#5d0500")

plot_df <- data.frame(
    cellid     = colnames(spe),
    sample_id  = spe$sample_id,
    x          = spatialCoords(spe)[, 1],
    y          = spatialCoords(spe)[, 2],
    first_type = results_df$first_type,
    spot_class = results_df$spot_class
) |>
    filter(spot_class != "reject",
           as.character(first_type) %in% names(itc_cols)) |>
    mutate(first_type = factor(as.character(first_type),
                               levels = names(itc_cols)))

samples <- sort(unique(plot_df$sample_id))
stopifnot(length(samples) > 0)

plots <- lapply(samples, function(s) {
    d <- plot_df[plot_df$sample_id == s, ]
    ggplot(d, aes(x = x, y = y, colour = first_type)) +
        geom_point(size = 0.6, alpha = 0.9) +
        scale_colour_manual(values = itc_cols, drop = FALSE) +
        coord_fixed() +
        scale_y_reverse() +
        labs(title = s,
             subtitle = paste0("n = ", nrow(d), " ITC spots"),
             colour = "ITC subtype") +
        theme_void(base_size = 11) +
        theme(plot.title    = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
        guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
})

n_col <- 3
n_row <- ceiling(length(plots) / n_col)

combined <- wrap_plots(plots, ncol = n_col) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")

pdf(here(plot_dir, "rctd_ITCmn_first_type_per_sample.pdf"),
    width = 4 * n_col, height = 4 * n_row)
print(combined)
dev.off()



# save spe
saveRDS(spe, here(processed_dir, "spe_with_rctd_ITCmn.rds"))

