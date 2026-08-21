library(here)
library(ggplot2)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggspavis)
library(spacexr)

# ======= Load and prep data =======
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_cells.rds"))
spe <- logNormCounts(spe)

processed_dir <- here("processed-data", "VisiumHD", "05_label_transfer")

# Make sure gene names are in rownames
# (plotCoords uses feature_names to look up rownames)
genes <- c("TSHZ1", "DRD1", "FOXP2")
genes_present <- genes[genes %in% rownames(spe)]
genes_missing <- setdiff(genes, genes_present)
if (length(genes_missing) > 0) {
    message("Genes not found: ", paste(genes_missing, collapse = ", "))
}

# ======= Plot: one PDF per gene, one page per sample =======
plot_dir <- here("plots", "VisiumHD", "05_label_transfer", "RCTD_human_Yu")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

unique_samples <- sort(unique(spe$sample_id))

for (g in genes_present) {
    pdf_path <- file.path(plot_dir, paste0(g, "_expression.pdf"))
    pdf(file = pdf_path, width = 10, height = 10)

    for (s in unique_samples) {
        spe_sub <- spe[, spe$sample_id == s]

        p <- plotCoords(spe_sub,
                        annotate = g,
                        feature_names = "symbol",
                        assay_name = "logcounts",
                        pal = c("white", "red"),
                        point_size = 0.3,
                        y_reverse = TRUE) +
            ggtitle(paste0(s, " — ", g))

        print(p)
        cat("Plotted:", s, "-", g, "\n")
    }

    dev.off()
    cat("Saved:", pdf_path, "\n")
}

cat("\nDone!\n")