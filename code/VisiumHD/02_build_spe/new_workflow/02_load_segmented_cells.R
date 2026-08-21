library(here)
library(arrow)
library(dplyr)
library(VisiumIO)
library(SpatialExperiment)

.load_segmented <- function(outs_dir, sample_id) {
    message("Loading segmented cells: ", sample_id)

    spe <- TENxVisiumHD(
        segmented_outputs = file.path(outs_dir, "segmented_outputs"),
        format            = "h5",
        images            = "lowres"
    ) |> import()

    rownames(spe) <- make.unique(rowData(spe)$Symbol)
    if ("Type" %in% colnames(rowData(spe))) rowData(spe)$Type <- NULL
    spe$sample_id <- sample_id

    # ── Manually attach barcode mappings ─────────────────────────────────
    # barcode_mappings.parquet lives in outs/ root (not segmented_outputs/)
    map_path <- file.path(outs_dir, "barcode_mappings.parquet")

    if (!file.exists(map_path)) {
        warning("barcode_mappings.parquet not found for ", sample_id)
        return(spe)
    }

    map <- read_parquet(map_path) |>
        rename_with(tolower) |>
        filter(in_cell, cell_id %in% colnames(spe)) |>   # only assigned bins
        split(f = ~cell_id)

    spe$map <- map[match(colnames(spe), names(map))]

    message("  -> ", ncol(spe), " cells | map attached: ",
            sum(!sapply(spe$map, is.null)), " / ", ncol(spe))
    spe
}

sample_manifest <- list(
    Br9280_CeA = here("processed-data", "VisiumHD", "01_spaceranger",
                      "H1-W369TJK_A1", "outs"),
    Br8325_MeA = here("processed-data", "VisiumHD", "01_spaceranger",
                      "H1-937FVHX_A1", "outs"),
    Br8325_CeA = here("processed-data", "VisiumHD", "01_spaceranger",
                      "H1-937FVHX_D1", "outs"),
    Br9280_ITC = here("processed-data", "VisiumHD", "01_spaceranger",
                      "H1-HW9VGBW_A1", "outs"),
    Br9280_MeA = here("processed-data", "VisiumHD", "01_spaceranger",
                      "H1-HW9VGBW_D1", "outs")
)

spe_list <- mapply(
    .load_segmented,
    outs_dir  = sample_manifest,
    sample_id = names(sample_manifest),
    SIMPLIFY  = FALSE
)

spe_cells <- do.call(cbind, spe_list)

# convert to in-memory matrices (instead of DelayedArray-backed) for downstream compatibility
assays(spe_cells) <- endoapply(assays(spe_cells), as.matrix)

# Verify map is attached
message("colData cols: ", paste(names(colData(spe_cells)), collapse = ", "))
message("Map populated for: ",
        sum(!sapply(spe_cells$map, is.null)), " / ", ncol(spe_cells), " cells")
message("Map columns: ", paste(names(spe_cells$map[[1]]), collapse = ", "))

out_dir <- here("processed-data", "VisiumHD", "02_build_spe")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(spe_cells, file.path(out_dir, "spe_cells_combined.rds"))
message("Saved -> spe_cells_combined.rds")