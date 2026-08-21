library("SpatialExperiment")
library("here")
library("dplyr")
library("tidyr")
library("ggplot2")
library("BiocNeighbors")

out_dir  <- here("processed-data", "VisiumHD", "06_composition")
plot_dir <- here("plots",          "VisiumHD", "06_composition")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load both labelled objects --------------------------------------------
spe_dom   <- readRDS(here("processed-data", "VisiumHD", "04_clustering",
                          "spe_banksy_renamed.rds"))
spe_cells <- readRDS(here("processed-data", "VisiumHD", "05_label_transfer",
                          "spe_with_rctd_ITCmn.rds"))

# domain lookup keyed by sample_id::16um-bin-barcode
dom_by_key <- setNames(as.character(spe_dom$spatial_domain),
                       paste(spe_dom$sample_id, colnames(spe_dom), sep = "::"))
 
# =============================================================================
# Map each cell -> 16um bin -> spatial_domain
# =============================================================================

sids <- as.character(spe_cells$sample_id)
has_map <- "map" %in% colnames(colData(spe_cells)) &&
           sum(!vapply(spe_cells$map, is.null, logical(1))) > 0
 
if (has_map) {
    message("Using SpaceRanger `map` column (square_016um), area-weighted vote.")
 
    # --- unnest all per-cell maps into one long table (cell, sid, bin) ------
    mt <- rbindlist(
        lapply(seq_along(spe_cells$map), function(i) {
            m <- spe_cells$map[[i]]
            if (is.null(m) || !nrow(m) || !"square_016um" %in% colnames(m))
                return(NULL)
            data.table(cell    = i,
                       sid     = sids[i],
                       bin     = as.character(m$square_016um),
                       in_cell = if ("in_cell" %in% colnames(m)) m$in_cell else TRUE)
        }),
        use.names = TRUE, fill = TRUE
    )
 
    # free the heavy list-column now that it's unnested
    spe_cells$map <- NULL; gc()
 
    # --- restrict to in-cell 2um bins, attach domain -----------------------
    mt <- mt[in_cell == TRUE]
    mt[, dom := dom_by_key[paste(sid, bin, sep = "::")]]
    mt <- mt[!is.na(dom)]
 
    # --- area-weighted vote: count 2um bins per (cell, domain), take max ----
    votes <- mt[, .N, by = .(cell, dom)]
    setorder(votes, cell, -N)
 
    # flag ties (top two domains equal) for transparency
    tie_cells <- votes[, .(tie = (.N >= 2L) && (N[1] == N[2])), by = cell][tie == TRUE, cell]
    if (length(tie_cells))
        message(length(tie_cells), " cells had a tied domain vote (assigned by sort order).")
 
    winners <- votes[, .SD[1L], by = cell]   # top domain per cell
 
    cell_domain <- rep(NA_character_, ncol(spe_cells))
    cell_domain[winners$cell] <- winners$dom
    spe_cells$spatial_domain <- factor(cell_domain)
 
    rm(mt, votes, winners); gc()
 
} else {
    message("No usable `map` column -> per-sample nearest-bin spatial join.")
    dom_vec  <- rep(NA_character_, ncol(spe_cells))
    dist_vec <- rep(NA_real_,      ncol(spe_cells))
    for (s in shared) {
        ci <- which(spe_cells$sample_id == s); di <- which(spe_dom$sample_id == s)
        if (!length(di)) next
        nn <- queryKNN(spatialCoords(spe_dom)[di, , drop = FALSE],
                       spatialCoords(spe_cells)[ci, , drop = FALSE], k = 1)
        dom_vec[ci]  <- as.character(spe_dom$spatial_domain[di][nn$index[, 1]])
        dist_vec[ci] <- nn$distance[, 1]
    }
    spe_cells$spatial_domain <- factor(dom_vec)
    spe_cells$bin_dist       <- dist_vec
    message("Cell-to-bin distance (median per sample):")
    print(tapply(dist_vec, spe_cells$sample_id, median, na.rm = TRUE))
}
 
message("Unmapped cells: ", sum(is.na(spe_cells$spatial_domain)),
        " / ", ncol(spe_cells))
saveRDS(spe_cells, file.path(out_dir, "spe_cells_with_domain.rds"))
 
# =============================================================================
# Composition per domain
# =============================================================================
drop_domains <- c("DELETE", "Ventricle", "Vascular")
cell_meta <- as.data.frame(colData(spe_cells)[
        , c("sample_id", "spatial_domain", "first_type",
            "second_type", "spot_class")]) |>
    filter(!is.na(spatial_domain), !spatial_domain %in% drop_domains,
           spot_class != "reject") |>
    droplevels()
 
comp <- cell_meta |>
    count(spatial_domain, first_type, name = "n") |>
    group_by(spatial_domain) |>
    mutate(prop = n / sum(n), n_domain = sum(n)) |>
    ungroup()
comp_wide <- comp |>
    pivot_wider(id_cols = spatial_domain, names_from = first_type,
                values_from = prop, values_fill = 0)
 
write.csv(comp,      file.path(out_dir, "composition_long.csv"), row.names = FALSE)
write.csv(comp_wide, file.path(out_dir, "composition_wide.csv"), row.names = FALSE)
 
p_prop <- ggplot(comp, aes(spatial_domain, prop, fill = first_type)) +
    geom_col() +
    labs(x = NULL, y = "proportion of cells", fill = "cell type",
         title = "Cell-type composition per spatial domain") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
# ITC-subtype split by domain (interstitial vs WM question)
itc_cols <- c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")
itc_comp <- cell_meta |>
    filter(as.character(first_type) %in% names(itc_cols)) |>
    count(spatial_domain, first_type, name = "n") |>
    mutate(first_type = factor(as.character(first_type), levels = names(itc_cols)))
 
p_itc <- ggplot(itc_comp, aes(spatial_domain, n, fill = first_type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = itc_cols) +
    labs(x = NULL, y = "proportion of ITC cells", fill = "ITC subtype",
         title = "ITC subtype split by domain") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
pdf(file.path(plot_dir, "composition_per_domain.pdf"), width = 10, height = 6)
print(p_prop); print(p_itc); dev.off()
 
message("Done. Outputs in:\n  ", out_dir, "\n  ", plot_dir)
 