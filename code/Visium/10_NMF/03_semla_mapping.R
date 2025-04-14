suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("semla")
    library("visiumStitched")
    library("RcppML")
})

plot_dir <- here("plots", "Visium", "10_NMF")

spe <- readRDS(here("processed-data", "Visium", "09_marker_genes", "spe_Br8325_BS_manual.rds"))
spe
# dim: 36601 29885 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(1): symbol
# colnames(29885): AAACAAGTATCTCCCA-1_V13M06-387_A1
#   AAACACCAATAACTGC-1_V13M06-387_A1 ... TTGTTTCATTAGTCTA-1_V13M06-388_D1
#   TTGTTTCCATACAACT-1_V13M06-388_D1
# colData names(57): in_tissue array_row ... BS_k9 BS_manual
# reducedDimNames(4): PCA HARMONY UMAP-PCA UMAP-HARMONY
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

sce.excit <- readRDS(here("processed-data/snRNA-seq/sce_FINAL_all_celltypes.rds"))
sce.excit
# class: SingleCellExperiment 
# dim: 13842 166905 
# metadata(0):
# assays(2): counts logcounts
# rownames(13842): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
# colnames(166905): AAACCCACAGGTCCCA-1 AAACCCATCCATTGTT-1 ...
#   TTTGGTTGTGGCTAGA-1 TTTGTTGCACATGGTT-1
# colData names(34): orig.ident nCount_originalexp ... broad_celltype
#   fine_celltype
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
altExpNames(0):

unique(sce.excit$fine_celltype)
#  [1] "ESR1_ADRA1A"     "MEIS2_COL25A1"   "MEIS1_PARD3B"    "ZBTB20_SLC4A4"  
#  [5] "PEX5L_MYRIP"     "ST18_ABCA8"      "ADARB2_TRPS1"    "GULP1_TRHDE"    
#  [9] "RXFP1_KIAA1217"  "GRIK3_TNS3"      "SATB2_MPPED1"    "SLC17A8_ST8SIA2"

# drop RXFP1_KIAA1217 and SLC17A8_ST8SIA2 from the list of excitatory clusters
#sce.excit <- sce.excit[, sce.excit$fine_celltype != "RXFP1_KIAA1217" & sce.excit$fine_celltype != "SLC17A8_ST8SIA2"]

# ======= Converting to Seurat =======

# make rownames unique
rownames(spe) <- make.unique(rownames(spe))

as.Seurat <- function(spe,
    spatial_cols = c(
        "tissue" = "in_tissue",
        "row" = "array_row",
        "col" = "array_col",
        "imagerow" = "pxl_row_in_fullres_transformed",
        "imagecol" = "pxl_col_in_fullres_transformed"
    ),
    verbose = TRUE) {
    #   Seurat is only suggested
    BiocBaseUtils::checkInstalled("Seurat")

    SPOT_DIAMETER <- 55e-6

    #   Ensure all necessary columns are present in colData
    required_col_names <- c("tissue", "row", "col", "imagerow", "imagecol")
    if (!all(required_col_names %in% names(spatial_cols))) {
        missing_col_names <- required_col_names[!(required_col_names %in% names(spatial_cols))]
        stop(
            sprintf(
                "Expected the following named elements spatial_cols: '%s'",
                paste(missing_col_names, collapse = "', '")
            )
        )
    }

    required_cols <- spatial_cols[required_col_names]
    col_info <- cbind(colData(spe), SpatialExperiment::spatialCoords(spe))
    if (!all(required_cols %in% colnames(col_info))) {
        missing_cols <- required_cols[!(required_cols %in% colnames(col_info))]
        stop(
            sprintf(
                "Expected the following columns in colData(spe) or spatialCoords(spe): '%s'",
                paste(missing_cols, collapse = "', '")
            )
        )
    }

    #   Uniqueness of spot names
    if (any(duplicated(colnames(spe)))) {
        stop("Seurat requires colnames(spe) to be unique")
    }

    #   Low-res images must exist for each sample ID
    if (sum(imgData(spe)$image_id == "lowres") < length(unique(spe$sample_id))) {
        stop("Each sample ID must have a low-resolution image for conversion")
    }

    if (verbose) message("Running 'as.Seurat(spe)'...")
    seur <- Seurat::as.Seurat(spe)

    for (sample_id in unique(spe$sample_id)) {
        if (verbose) {
            message(
                sprintf(
                    "Adding spot coordinates and images for sample %s...",
                    sample_id
                )
            )
        }
        spe_small <- spe[, spe$sample_id == sample_id]

        coords <- col_info[, required_cols, drop = FALSE]
        colnames(coords) <- required_col_names
        coords$tissue <- as.integer(coords$tissue)
        coords <- as.data.frame(coords)

        this_img <- array(
            t(col2rgb(imgRaster(spe_small))),
            dim = c(dim(imgRaster(spe_small)), 3)
        ) / 256

        seur@images[[sample_id]] <- new(
            Class = "VisiumV1",
            image = this_img,
            scale.factors = Seurat::scalefactors(
                spot = NA, fiducial = NA, hires = NA,
                lowres = imgData(spe)[
                    imgData(spe_small)$image_id == "lowres", "scaleFactor"
                ]
            ),
            coordinates = coords,
            spot.radius = SPOT_DIAMETER / scaleFactors(spe_small),
            assay = "originalexp",
            key = paste0(sample_id, "_")
        )
    }

    if (verbose) message("Returning converted object...")
    return(seur)
}


#  === PR -> Seurat -> Semla ===
seur <- as.Seurat(spe)
# Use a different graphics device to avoid X11 display error
png(filename = here(plot_dir, "semla_update.png"))
semla <- UpdateSeuratForSemla(seur)
dev.off()
semla
# An object of class Seurat 
# 36601 features across 29885 samples within 1 assay 
# Active assay: originalexp (36601 features, 0 variable features)
#  2 layers present: counts, data
#  4 dimensional reductions calculated: PCA, HARMONY, UMAP.PCA, UMAP.HARMONY
#  1 image present: Br8325

# ==== convert SCE ====
# make colnames unique
colnames(sce.excit) <- make.unique(colnames(sce.excit))
seur.excit <- Seurat::as.Seurat(sce.excit)
seur.excit
# An object of class Seurat 
# 13842 features across 76519 samples within 1 assay 
# Active assay: originalexp (13842 features, 0 variable features)
#  2 layers present: counts, data
#  2 dimensional reductions calculated: PCA, UMAP


# ======== Semla Testing =========

png(here(plot_dir, "snRNA-seq_all_clusters_human.png"), width = 10, height = 10, units = "in", res = 300)
DimPlot(seur.excit, group.by = "fine_celltype", label = TRUE)
dev.off()

# plot spatial data
png(here(plot_dir, "Visium_all_clusters.png"), width = 10, height = 10, units = "in", res = 300)
MapLabels(semla, column_name="BS_manual")
dev.off()


# ======== NNLS ========

rm(spe)

# HVGs for NNLS
seur.excit <- seur.excit |> 
  FindVariableFeatures(nfeatures = 10000)

semla <- semla |>
  FindVariableFeatures(nfeatures = 10000)

# Run and time NNLS
ti <- Sys.time()
semla <- RunNNLS(object = semla, 
                            singlecell_object = seur.excit, 
                            groups = "fine_celltype",
                            singlecell_assay= "originalexp",
                            spatial_assay = "originalexp",)
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))


# Check available cell types
rownames(semla)
#  [1] "ADARB2-TRPS1"  "ESR1-ADRA1A"   "GRIK3-TNS3"    "GULP1-TRHDE"  
#  [5] "MEIS1-PARD3B"  "MEIS2-COL25A1" "PEX5L-MYRIP"   "SATB2-MPPED1" 
#  [9] "ST18-ABCA8"    "ZBTB20-SLC4A4"


# Plot selected cell types
DefaultAssay(semla) <- "celltypeprops"

selected_celltypes <- rownames(semla)
semla <- LoadImages(semla, image_height = 1e3)

# plot
pdf(here(plot_dir, "NNLS_celltypes_total_primate.pdf"), width = 10, height = 10)
plots <- lapply(seq_along(selected_celltypes), function(i) {
  MapFeatures(semla, pt_size = 1.3,
            features = selected_celltypes[i], image_use = NULL,
            arrange_features = "row", scale = "shared", 
            #override_plot_dims = TRUE,
            colors = RColorBrewer::brewer.pal(n = 10, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
}) |> setNames(nm = selected_celltypes)
plots
dev.off()



selcted_celltypes_for_multi <- c("ZFHX3-SCN5A", "TSHZ1.1", "SATB2-MPPED1", "Oligodendrocyte", "MEIS2-COL25A1", "ADARB2-TRPS1")

# plot all cell types
pdf(here(plot_dir, "NNLS_all_celltypes_total_primate.pdf"), width = 10, height = 10)
MapMultipleFeatures(semla, 
                    #image_use = "raw", 
                    pt_size = 1.5, max_cutoff = 0.99,
                    #override_plot_dims = TRUE, 
                    colors = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(12),
                    features = selcted_celltypes_for_multi) +
  plot_layout(guides = "collect")
dev.off()




cor_matrix <- FetchData(semla, selected_celltypes) |> 
  mutate_all(~ if_else(.x<0.1, 0, .x)) |>  # Filter lowest values (-> set as 0)
  cor()

diag(cor_matrix) <- NA
cor_matrix[is.na(cor_matrix)] <- 0
max_val <- max(cor_matrix, na.rm = T)
cols <- RColorBrewer::brewer.pal(7, "RdYlBu") |> rev(); cols[4] <- "white"
pdf(here(plot_dir, "NNLS_celltypes_correlation.pdf"), width = 10, height = 10)
pheatmap::pheatmap(cor_matrix, 
                   breaks = seq(-max_val, max_val, length.out = 100),
                   color=colorRampPalette(cols)(100),
                   cellwidth = 14, cellheight = 14, 
                   treeheight_col = 10, treeheight_row = 10, 
                   main = "Cell type correlation\nwithin spots")
dev.off()