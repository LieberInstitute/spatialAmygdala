## 04_assign_ITC_broad.R
## Assign 3 broad ITC groups based on MetaNeighbor dendrogram

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scran)
    library(scuttle)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
})

sce <- readRDS(here("processed-data", "snRNAseq", "BICCN","sce_ITC_harmony.rds"))

# ---- Define groups ----
ITC_2 <- c("EMSN_232", "EMSN_230", "EMSN_231", "EMSN_233", "EMSN_234",
           "TSHZ1_PRKG1", "Human_TSHZ1 CALCRL")

ITC_1 <- c("EMSN_224", "EMSN_222", "EMSN_223",
           "EMSN_225", "EMSN_226")

ITC_3 <- c("EMSN_229", "EMSN_426", "EMSN_228", "EMSN_227",
           "TSHZ1_CPNE4", "Human_TSHZ1 SEMA3C")

# ---- Assign ----
sce$mn_cluster <- ifelse(sce$celltype_fine %in% ITC_1, "ITC_1",
                      ifelse(sce$celltype_fine %in% ITC_2, "ITC_2",
                      ifelse(sce$celltype_fine %in% ITC_3, "ITC_3", NA)))

# ---- 1. Find markers (1 vs all) ----
message("Running findMarkers...")
markers <- findMarkers(sce,
    groups = sce$mn_cluster,
    pval.type = "all",
    direction = "up"
)

# save top 50

top_n <- 50
out <- do.call(rbind, lapply(names(markers), function(g) {
    res <- as.data.frame(markers[[g]])
    res <- res[order(res$FDR), ]
    res <- head(res, top_n)
    # Keep only shared columns; rename pairwise logFCs generically
    lfc_cols <- grep("^logFC\\.", colnames(res), value = TRUE)
    res$mean_logFC <- rowMeans(res[, lfc_cols, drop = FALSE])
    res <- res[, c("p.value", "FDR", "summary.logFC", "mean_logFC")]
    res$gene <- rownames(res)
    res$group <- g
    res
}))
 
out_path <- here("processed-data", "snRNAseq", "ITC_broad_markers_top50.csv")
write.csv(out, out_path, row.names = FALSE)
message("Saved: ", out_path)
 
 
# ---- 4. UMAP colored by MetaNeighbor cluster (Siletti only) ----
sce_siletti <- sce[, sce$dataset == "BICCN"]

# Build a data frame from the UMAP coordinates
umap_df <- data.frame(
    UMAP1 = reducedDim(sce_siletti, "UMAP_harmony")[, 1],
    UMAP2 = reducedDim(sce_siletti, "UMAP_harmony")[, 2],
    cluster = sce_siletti$mn_cluster,  
    celltype = sce_siletti$celltype_fine
)

# Shuffle rows so no group is always on top
set.seed(42)
umap_df <- umap_df[sample(nrow(umap_df)), ]

# Cluster centroids for labels
centroids <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = umap_df, FUN = median)

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    scattermore::geom_scattermore(pointsize = 4, alpha = 1, pixels = c(2048, 2048)) +
    scale_colour_manual(values = c("ITC_1" = "#fd0d00",
                                 "ITC_2" = "#f5b6b3",
                                 "ITC_3" = "#5d0500")) +
    ggrepel::geom_label_repel(
        data = centroids,
        aes(label = cluster),
        size = 4, fontface = "bold",
        fill = "white", alpha = 0.8,
        label.size = 0.3,
        box.padding = 0.5,
        show.legend = FALSE
    ) +
    labs(title = "MetaNeighbor fine-type clusters (Siletti)",
         colour = "MN Cluster") +
    theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

pdf(here("plots", "snRNAseq", "ITCs", "umap_mn_clusters_siletti_scattermore.pdf"),
    width = 6, height = 6)
print(p)
dev.off()
message("UMAP saved.")



# ---- 6. Pseudobulk boxplots per cluster (Siletti) ----
library(patchwork)

genes_box <- c("TSHZ1", "FOXP2")

sce_siletti$mn_cluster <- factor(sce_siletti$mn_cluster, levels = c("ITC_1", "ITC_2", "ITC_3"))

# Sum raw counts per sample × cluster
pb_list <- lapply(levels(sce_siletti$mn_cluster), function(cl) {
    sub <- sce_siletti[, sce_siletti$mn_cluster == cl]
    samples <- unique(sub$sample_id)  # <-- change if needed

    do.call(rbind, lapply(samples, function(s) {
        cells <- sub[, sub$sample_id == s]
        raw_sums  <- rowSums(counts(cells)[genes_box, , drop = FALSE])
        lib_size  <- sum(counts(cells))
        log2cpm   <- log2(raw_sums / lib_size * 1e6 + 1)

        data.frame(
            sample  = s,
            cluster = cl,
            gene    = genes_box,
            log2cpm = log2cpm
        )
    }))
})

pb_df <- do.call(rbind, pb_list)
pb_df$cluster <- factor(pb_df$cluster)
pb_df$gene    <- factor(pb_df$gene, levels = genes_box)

# Drop samples that have 0 counts for any gene
samples_with_zero <- pb_df$sample[pb_df$log2cpm == 0]
pb_df <- pb_df[!pb_df$sample %in% samples_with_zero, ]

# Plot
p_box <- ggplot(pb_df, aes(x = cluster, y = log2cpm, fill = cluster)) +
geom_jitter(aes(colour = cluster), width = 0.15, size = 1.5, alpha = 0.8) +
geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("ITC_1" = "#fd0d00",
                                 "ITC_2" = "#f5b6b3",
                                 "ITC_3" = "#5d0500")) +
    scale_colour_manual(values = c("ITC_1" = "#fd0d00",
                                "ITC_2" = "#f5b6b3",
                                "ITC_3" = "#5d0500")) +
    labs(y = expression(log[2]~CPM), x = NULL) +
    theme_classic(base_size = 13) +
    theme(
        strip.text      = element_text(face = "bold.italic", size = 12),
        axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, NA)) +

pdf(here("plots", "snRNAseq", "ITCs", "pseudobulk_boxplots_siletti_donor.pdf"),
    width = 5, height = 4)
print(p_box)
dev.off()
message("Pseudobulk boxplots saved.")









library(ComplexHeatmap)
library(circlize)

# Marker genes per cluster
markers_itc <- list(
    ITC_1 = c("PPP1R1B", "PHACTR1", "PDE10A", "RARB", "CDH11"),
    ITC_2 = c("SHISA9", "VWC2L", "CADPS2", "GULP1", "PRKG1"),
    ITC_3 = c("DSCAM", "NRG3", "PREX2", "ADAMTSL1", "GRM3")
)
genes_hm <- unlist(markers_itc)

# Mean scaled logcounts per cluster
mean_expr <- sapply(levels(sce_siletti$mn_cluster), function(cl) {
    sub <- sce_siletti[, sce_siletti$mn_cluster == cl]
    rowMeans(logcounts(sub)[genes_hm, , drop = FALSE])
})

# Scale across clusters (rows)
mat_scaled <- t(scale(t(mean_expr)))

# Row annotation showing cluster membership
row_ann <- HeatmapAnnotation(
    Marker = rep(names(markers_itc), lengths(markers_itc)),
    col = list(Marker = c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")),
    which = "row",
    show_legend = TRUE,
    show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#1178df", "white", "#e70e27"))

ht <- Heatmap(mat_scaled,
    name = "Scaled\nlogcounts",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 10, fontface = "italic"),
    column_names_gp = gpar(fontsize = 11),
    column_title = "ITC subtype markers",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    left_annotation = row_ann,
    row_split = rep(names(markers_itc), lengths(markers_itc)),
    row_gap = unit(2, "mm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
)

pdf(here("plots", "snRNAseq", "ITCs", "itc_subtype_markers_heatmap.pdf"),
    width = 3, height = 5)
draw(ht)
dev.off()

message("ITC subtype markers heatmap saved.")










# Gene sets
genes_spiny  <- c("PPP1R1B", "OPRM1", "DRD1", "DRD2", "DRD3", "NOS1")
genes_aspiny <- c("DSCAM", "GRM1", "GRM3", "CALB1", "GABRA1", "ADCYAP1R1")

# Helper to build pseudobulk df for a given gene set
make_pb <- function(genes) {
    pb_list <- lapply(levels(sce_siletti$mn_cluster), function(cl) {
        sub <- sce_siletti[, sce_siletti$mn_cluster == cl]
        samples <- unique(sub$sample_id)
        do.call(rbind, lapply(samples, function(s) {
            cells <- sub[, sub$sample_id == s]
            raw_sums  <- rowSums(counts(cells)[genes, , drop = FALSE])
            lib_size  <- sum(counts(cells))
            log2cpm   <- log2(raw_sums / lib_size * 1e6 + 1)
            data.frame(sample = s, cluster = cl, gene = genes, log2cpm = log2cpm)
        }))
    })
    pb_df <- do.call(rbind, pb_list)
    samples_with_zero <- pb_df$sample[pb_df$log2cpm == 0]
    pb_df <- pb_df[!pb_df$sample %in% samples_with_zero, ]
    pb_df$cluster <- factor(pb_df$cluster, levels = c("ITC_1", "ITC_2", "ITC_3"))
    pb_df$gene    <- factor(pb_df$gene, levels = genes)
    pb_df
}

itc_cols <- c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")

make_boxplot <- function(pb_df, title) {
    ggplot(pb_df, aes(x = cluster, y = log2cpm, fill = cluster)) +
        geom_jitter(aes(colour = cluster), width = 0.15, size = 1.5, alpha = 0.8) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
        facet_wrap(~ gene, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = itc_cols) +
        scale_colour_manual(values = itc_cols, guide = "none") +
        coord_cartesian(ylim = c(0, NA)) +
        labs(y = expression(log[2]~CPM), x = NULL, title = title) +
        theme_classic(base_size = 13) +
        theme(
            plot.title       = element_text(face = "bold", size = 14),
            strip.text       = element_text(face = "bold.italic", size = 12),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            legend.position  = "none",
            panel.grid.minor = element_blank()
        )
}

pb_spiny  <- make_pb(genes_spiny)
pb_aspiny <- make_pb(genes_aspiny)

p_spiny  <- make_boxplot(pb_spiny,  "Spiny markers")
p_aspiny <- make_boxplot(pb_aspiny, "Aspiny markers")

# Save
pdf(here("plots", "snRNAseq", "ITCs", "pseudobulk_boxplots_spiny_markers.pdf"),
    width = 10, height = 3)
print(p_spiny)
dev.off()

pdf(here("plots", "snRNAseq", "ITCs", "pseudobulk_boxplots_aspiny_markers.pdf"),
    width = 10, height = 3)
print(p_aspiny)
dev.off()















# ============================================================
# RECEPTOR-ONLY GENE LISTS
# ============================================================

genes_glutamate <- list(
    AMPA      = c("GRIA1", "GRIA2", "GRIA3", "GRIA4"),
    NMDA      = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B"),
    Kainate   = c("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5"),
    Delta     = c("GRID1", "GRID2"),
    mGluR_I   = c("GRM1", "GRM5"),
    mGluR_II  = c("GRM2", "GRM3"),
    mGluR_III = c("GRM4", "GRM6", "GRM7", "GRM8")
)

genes_gaba <- list(
    GABAA_alpha = c("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6"),
    GABAA_beta  = c("GABRB1", "GABRB2", "GABRB3"),
    GABAA_gamma = c("GABRG1", "GABRG2", "GABRG3"),
    GABAA_delta = c("GABRD"),
    GABAA_epsilon = c("GABRE"),
    GABAA_theta = c("GABRQ"),
    GABAA_pi    = c("GABRP"),
    GABAA_rho   = c("GABRR1", "GABRR2", "GABRR3"),
    GABAB       = c("GABBR1", "GABBR2")
)

genes_catecholamine <- list(
    Dopamine_R = c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5"),
    Alpha1_AR  = c("ADRA1A", "ADRA1B", "ADRA1D"),
    Alpha2_AR  = c("ADRA2A", "ADRA2B", "ADRA2C"),
    Beta_AR    = c("ADRB1", "ADRB2", "ADRB3")
)

# ============================================================
# DOT PLOT FUNCTION
# ============================================================

make_dotplot <- function(sce, gene_list, title) {
    genes <- unlist(gene_list)
    genes <- genes[genes %in% rownames(sce)]
    subfam <- rep(names(gene_list), lengths(gene_list))
    subfam <- subfam[unlist(gene_list) %in% genes]

    # Compute per-cluster: pct expressed + mean expression
    dot_df <- do.call(rbind, lapply(levels(sce$mn_cluster), function(cl) {
        sub <- sce[, sce$mn_cluster == cl]
        mat <- logcounts(sub)[genes, , drop = FALSE]
        do.call(rbind, lapply(seq_along(genes), function(i) {
            vals <- mat[i, ]
            data.frame(
                gene    = genes[i],
                cluster = cl,
                pct_exp = mean(vals > 0) * 100,
                avg_exp = mean(vals)
            )
        }))
    }))


    dot_df$gene    <- factor(dot_df$gene, levels = rev(genes))
    dot_df$cluster <- factor(dot_df$cluster, levels = c("ITC_1", "ITC_2", "ITC_3"))
    dot_df$subfamily <- factor(subfam[match(dot_df$gene, genes)], levels = names(gene_list))

    ggplot(dot_df, aes(x = cluster, y = gene)) +
        geom_point(aes(size = pct_exp, colour = avg_exp)) +
        scale_size_continuous(range = c(0.5, 6), name = "% Expressing",
                              breaks = c(10, 25, 50, 75)) +
        scale_colour_gradient(low = "white", high = "#B2182B",
                      name = "Mean\nExpression") +
        facet_grid(subfamily ~ ., scales = "free_y", space = "free_y",
                   switch = "y") +
        labs(title = title, x = NULL, y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
            plot.title        = element_text(face = "bold", size = 14),
            axis.text.y       = element_text(face = "italic", size = 9),
            axis.text.x       = element_text(angle = 45, hjust = 1),
            strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
            strip.placement   = "outside",
            panel.grid.major  = element_line(colour = "grey90"),
            panel.grid.minor  = element_blank(),
            legend.position   = "right"
        )
}

# ============================================================
# GENERATE AND SAVE
# ============================================================

p_glut <- make_dotplot(sce_siletti, genes_glutamate, "Glutamate receptors")
p_gaba <- make_dotplot(sce_siletti, genes_gaba, "GABA receptors")
p_cat  <- make_dotplot(sce_siletti, genes_catecholamine, "Catecholamine receptors")

pdf(here("plots", "snRNAseq", "ITCs", "dotplot_glutamate_receptors.pdf"),
    width = 5, height = 10)
print(p_glut)
dev.off()

pdf(here("plots", "snRNAseq", "ITCs", "dotplot_gaba_receptors.pdf"),
    width = 5, height = 10)
print(p_gaba)
dev.off()

pdf(here("plots", "snRNAseq", "ITCs", "dotplot_catecholamine_receptors.pdf"),
    width = 5, height = 6)
print(p_cat)
dev.off()

message("Receptor dot plots saved.")