suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("RColorBrewer")
    library("ggplot2")
    library("escheR")
})

## Load once
load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))

## All sample IDs present
sample_ids <- sort(unique(colData(spe)$sample_id))

## Loop one sample at a time
for (sid in sample_ids) {
    message("Processing sample: ", sid)

    # restrict to this sample
    spe.subset <- spe[, colData(spe)$sample_id == sid]

    # where the BayesSpace cluster csvs live for this sample
    bs_dir <- here::here("processed-data","Visium","07_clustering","BayesSpace","SVGs", sid)
    if (!dir.exists(bs_dir)) {
        message("  Skipping: directory not found")
        next
    }

    # subfolders named like BayesSpace_10, BayesSpace_12, etc.
    bs_folders <- list.files(bs_dir, full.names = TRUE)
    bs_folders <- bs_folders[grepl("BayesSpace_", basename(bs_folders))]
    if (length(bs_folders) == 0) {
        message("  Skipping: no BayesSpace_* folders in ", bs_dir)
        next
    }

    # ensure output dir exists
    out_dir <- here::here("plots","Visium","07_clustering","BayesSpace","SVGs", sid)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # loop over k’s for this sample
    for (bf in bs_folders) {
        k <- sub(".*BayesSpace_", "", basename(bf))
        bs_csv <- list.files(bf, pattern = "[.]csv$", full.names = TRUE)[1]
        if (is.na(bs_csv)) next

        bs_df <- read.csv(bs_csv)

        # add cluster column to colData(spe.subset)
        clustV <- paste0("BS_k", k)
        colData(spe.subset)[[clustV]] <- factor(bs_df$cluster)

        # plot PDF
        pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(
            length(unique(colData(spe.subset)[[clustV]]))
        )
        pdf(file = file.path(out_dir, paste0(clustV, "_stitched.pdf")),
            width = 10, height = 10)
        p <- make_escheR(spe.subset) |>
            add_fill(var = clustV, point_size = 1.75) +
            scale_fill_manual(values = pal) +
            ggtitle(paste0(sid, " • ", clustV))
        print(p)
        dev.off()
    }

}
