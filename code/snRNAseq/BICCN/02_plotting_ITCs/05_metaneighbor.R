suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(MetaNeighbor)
    library(ggplot2)
    library(pheatmap)
    library(ComplexHeatmap)
    library(circlize)
})

# ---- Load data ----
sce  <- readRDS(here("processed-data", "snRNAseq", "BICCN", "sce_ITC_harmony.rds"))

# ---- 1. Unsupervised MetaNeighbor (broad types) ----

# MetaNeighbor expects study_id and cell_type in colData
sce$study_id  <- sce$dataset
sce$cell_type <- sce$celltype_broad

var_genes = variableGenes(dat = sce, exp_labels = sce$study_id)

aurocs <- MetaNeighborUS(
    var_genes   = var_genes,
    dat         = sce,
    study_id    = sce$study_id,
    cell_type   = sce$cell_type,
    fast_version = TRUE
)
 
# ---- ComplexHeatmap helper ----
# Parse MetaNeighbor labels (format: "study|celltype")
parse_mn_labels <- function(labels) {
    parts <- strsplit(labels, "\\|")
    data.frame(
        dataset  = sapply(parts, `[`, 1),
        celltype = sapply(parts, `[`, 2),
        stringsAsFactors = FALSE
    )
}
 
plot_auroc_heatmap <- function(aurocs, main = "",
                               cluster_method = "complete",
                                column_split=3,
                                row_split=3,
                               show_numbers = TRUE) {
    col_fun <- colorRamp2(c(0, 0.5, 1), c("#4575b4", "#f7f7f7", "#d73027"))

    # Parse labels BEFORE renaming
    lab_info <- parse_mn_labels(colnames(aurocs))

    # Now clean the display names
    rownames(aurocs) <- parse_mn_labels(rownames(aurocs))$celltype
    colnames(aurocs) <- parse_mn_labels(colnames(aurocs))$celltype
    # Dataset color palette
    datasets <- unique(lab_info$dataset)
    ds_cols <- setNames(
        c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")[seq_along(datasets)],
        datasets
    )
 
    # Annotations
    ha_col <- HeatmapAnnotation(
        Dataset  = lab_info$dataset,
        col = list(Dataset = ds_cols),
        annotation_name_side = "left",
        show_legend = TRUE
    )
    ha_row <- rowAnnotation(
        Dataset  = lab_info$dataset,
        col = list(Dataset = ds_cols),
        show_legend = FALSE
    )
 
    # Cell labels
    cell_fn <- if (show_numbers) {
        function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", aurocs[i, j]), x, y,
                            gp = grid::gpar(fontsize = 8,
                                            col = ifelse(aurocs[i, j] > 0.8 |
                                                         aurocs[i, j] < 0.2,
                                                         "white", "black")))
        }
    } else NULL
 
    # Clustering: group by celltype first, then by dataset within
    dist_mat <- as.dist(1 - aurocs)
    hc <- hclust(dist_mat, method = cluster_method)
 
    Heatmap(aurocs,
            name = "AUROC",
            col  = col_fun,
            column_title = main,
            column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
            cluster_rows    = hc,
            cluster_columns = hc,
            bottom_annotation  = ha_col,
            right_annotation = ha_row,
            column_split=column_split,
            row_split=row_split,
            #cell_fun = cell_fn,
            row_names_gp    = grid::gpar(fontsize = 9),
            column_names_gp = grid::gpar(fontsize = 9),
            border = TRUE,
            heatmap_legend_param = list(
                title = "AUROC",
                at = c(0, 0.25, 0.5, 0.75, 1)
            ))
}
 
# Plot broad heatmap
pdf(here("plots", "snRNAseq","ITCs", "metaneighbor_broad.pdf"),
    width = 7, height = 6)
draw(plot_auroc_heatmap(aurocs,
    main = "MetaNeighbor: ITC broad types across datasets",
    column_split=NULL,
    row_split=NULL,
    cluster_method = "complete"))
dev.off()
message("Broad-type heatmap saved.")
 
# ---- 2. Unsupervised MetaNeighbor (fine types) ----
# Test reproducibility at the fine-grained level
message("Running MetaNeighborUS on fine ITC types...")

# strip Humam_ from cell type fine
sce$celltype_fine <- gsub("^Human_", "", sce$celltype_fine) 

aurocs_fine <- MetaNeighborUS(
    var_genes   = var_genes,
    dat         = sce,
    study_id    = sce$study_id,
    cell_type   = sce$celltype_fine,
    fast_version = TRUE
)
 
message("AUROC matrix (fine types):")
print(round(aurocs_fine, 3))
 

pdf(here("plots", "snRNAseq","ITCs", "metaneighbor_fine.pdf"),
    width = 6, height = 5)
draw(plot_auroc_heatmap(aurocs_fine,
    main = "MetaNeighbor: ITC fine types across datasets",
    column_split=3,
    row_split=3,
    cluster_method = "complete"))
dev.off()
message("Fine-type heatmap saved.")
 









genes <- c("GRM1", "GABRA1", "DRD1", "TSHZ1", "FOXP2", "DRD3", "OPRM1", "PPP1R1B", "NOS1")
genes <- genes[genes %in% rownames(sce)]
 
plots <- lapply(genes, function(g) {
    plotReducedDim(sce, "UMAP_harmony", colour_by = g, point_size = 0.3) +
        ggtitle(g) +
        theme(plot.title = element_text(face = "italic"))
})
 
p <- wrap_plots(plots, ncol = 3)
 
pdf(here("plots", "snRNAseq", "ITCs", "ITC_gene_umaps.pdf"),
    width = 12, height = 10)
print(p)
dev.off()
message("Done.")
 