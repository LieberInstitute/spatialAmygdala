library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("harmony")

# Get input from command line
args <- commandArgs(trailingOnly = TRUE)
obj_name <- args[1]
message("Running for object: ", obj_name)

# Save directories
processed_dir <- here("processed-data", "Xenium", "04_dim_reduction")
plot_dir <- here("plots", "Xenium", "04_dim_reduction")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load the RData file
input_file <- here("processed-data","Xenium", "02_build_spe", paste0(obj_name, ".Rdata"))
load(input_file)  # Loads 'spe' object
spe

if ("Type" %in% colnames(rowData(spe))) {
    message("Detected Xenium object — subsetting to Gene Expression probes only.")
    gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
    spe <- spe[gene_expression_idx, ]
    rownames(spe) <- rowData(spe)$Symbol
    message("Subset dimensions: ", paste(dim(spe), collapse = " x "))
}

spe <- scuttle::addPerCellQC(spe)

# ======= Normalization ========
message("Normalizing: ", obj_name)

# ============================
# Normalization
# ============================

# Check available metrics
cd_names <- colnames(colData(spe))

if (all(c("cell_area", "nucleus_area") %in% cd_names)) {
    message("Using cell_area and nucleus_area for normalization (Xenium object).")

    # Compute scaling factors
    spe$cell_area.sf <- spe$cell_area / median(spe$cell_area, na.rm = TRUE)
    spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area, na.rm = TRUE)

    # QC plots
    png(file = file.path(plot_dir, paste0(obj_name, "_cell_area_scaling_factors.png")),
        width = 10, height = 5, units = "in", res = 300)
    hist(spe$cell_area.sf, breaks = 50, main = "Cell area scaling factor", xlab = "Scaling factor")
    dev.off()

    png(file = file.path(plot_dir, paste0(obj_name, "_nucleus_area_scaling_factors.png")),
        width = 10, height = 5, units = "in", res = 300)
    hist(spe$nucleus_area.sf, breaks = 50, main = "Nucleus area scaling factor", xlab = "Scaling factor")
    dev.off()

    # Apply normalization
    assay(spe, "nucleus_normcounts") <- scuttle::normalizeCounts(
        spe, size.factors = spe$nucleus_area.sf, transform = "log", assay.type = "counts"
    )
    assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(
        spe, size.factors = spe$cell_area.sf, transform = "log", assay.type = "counts"
    )

} else if ("volume" %in% cd_names) {
    message("Using volume for normalization (Proseg object).")

    # spe$volume.sf <- spe$volume / median(spe$volume, na.rm = TRUE)

    # # QC plot
    # png(file = file.path(plot_dir, paste0(obj_name, "_volume_scaling_factors.png")),
    #     width = 10, height = 5, units = "in", res = 300)
    # hist(spe$volume.sf, breaks = 50, main = "Volume scaling factor", xlab = "Scaling factor")
    # dev.off()

    # # Apply normalization (only one version)
    # assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(
    #     spe, size.factors = spe$volume.sf, transform = "log", assay.type = "counts"
    # )

    assay(spe, "cell_normcounts") <- log1p(counts(spe))

} else {
    stop("No valid normalization metric found in colData.")
}
spe

assay(spe, "cell_normcounts")[1:5, 1:5]
dim(assay(spe, "cell_normcounts"))




# ======== Dim reduction and Harmony batch correction ========
message("Starting Dim Reduction: ", obj_name)

set.seed(1000)

spe <- RunHarmony(spe, group.by.vars = "brnum", reduction.save="HARMONY_2")

spe <- runUMAP(spe, dimred = "HARMONY_2", name="UMAP_HARMONY_2")

spe <- runUMAP(spe, dimred = "PCA", name="UMAP_UNCORRECTED")
spe

# ======== Visualization ========
message("Starting Visualizations: ", obj_name)

# PCA plots
png(file.path(plot_dir, paste0(obj_name, "_Uncorrected_PCs_Brnum.png")), width = 8, height = 6, units = "in", res = 300)
plotReducedDim(spe, dimred = "PCA", ncomponents = 4, colour_by = "brnum", scattermore=TRUE)
dev.off()

png(file.path(plot_dir, paste0(obj_name, "_Corrected_PCs_Brnum.png")), width = 8, height = 6, units = "in", res = 300)
plotReducedDim(spe, dimred = "HARMONY", ncomponents = 4, colour_by = "brnum", scattermore=TRUE)
dev.off()

# UMAP plot
png(file.path(plot_dir, paste0(obj_name, "_Corrected_UMAP_Brnum.png")), width = 10, height = 10, units = "in", res = 300)
plotReducedDim(spe, dimred = "UMAP", ncomponents = 2, colour_by = "brnum", scattermore=TRUE)
dev.off()

# Marker genes
markers <- c("MOBP", "SLC17A7", "GAD1", "GULP1", "COL25A1", "PDYN", "TSHZ1", "SST")

plots <- lapply(markers, function(marker) {
  plotReducedDim(spe, dimred = "UMAP", colour_by = marker, by.assay.type = "cell_normcounts", scattermore=TRUE) +
    scale_color_gradient(low = "grey", high = "red") +
    ggtitle(marker)
})

png(file.path(plot_dir, paste0(obj_name, "_Corrected_UMAP_markers.png")), width = 20, height = 10, units = "in", res = 300)
cowplot::plot_grid(plotlist = plots, ncol = 4)
dev.off()
# Save
message("Saving: ", obj_name)
output_file <- file.path(processed_dir, paste0(obj_name, "_harmonized_singlecell.rds"))
saveRDS(spe, output_file)
