
suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scuttle)
    library(edgeR)
    library(variancePartition)
    library(ComplexHeatmap)
    library(circlize)
})

sce <- readRDS(here("processed-data", "snRNAseq", "BICCN", "sce_ITC_harmony.rds"))

# ---- Define groups ----
ITC_2 <- c("EMSN_232", "EMSN_230", "EMSN_231", "EMSN_233", "EMSN_234",
           "TSHZ1_PRKG1", "Human_TSHZ1 CALCRL")
ITC_1 <- c("EMSN_224", "EMSN_222", "EMSN_223", "EMSN_225", "EMSN_226")
ITC_3 <- c("EMSN_229", "EMSN_426", "EMSN_228", "EMSN_227",
           "TSHZ1_CPNE4", "Human_TSHZ1 SEMA3C")

sce$mn_cluster <- ifelse(sce$celltype_fine %in% ITC_1, "ITC_1",
                  ifelse(sce$celltype_fine %in% ITC_2, "ITC_2",
                  ifelse(sce$celltype_fine %in% ITC_3, "ITC_3", NA)))

sce_siletti <- sce[, sce$dataset == "BICCN"]
sce_siletti$mn_cluster <- factor(sce_siletti$mn_cluster,
                                  levels = c("ITC_1", "ITC_2", "ITC_3"))

# refactor $sample_id and $donor
sce_siletti$sample_id <- factor(sce_siletti$sample_id)
sce_siletti$donor     <- factor(sce_siletti$donor)

# Majority donor per sample
donor_majority <- sce_siletti |>
    colData() |>
    as.data.frame() |>
    dplyr::count(sample_id, donor) |>
    dplyr::slice_max(n, by = sample_id, n = 1) |>
    dplyr::select(sample_id, donor)

donor_lookup <- setNames(donor_majority$donor, donor_majority$sample_id)

# ---- Pseudobulk: one per sample × cluster ----
sce_siletti$pb_id <- paste0(sce_siletti$sample_id, "_", sce_siletti$mn_cluster)

pb <- aggregateAcrossCells(sce_siletti,
    ids        = sce_siletti$pb_id,
    use.assay.type = "counts"
)

# Build metadata
meta <- data.frame(
    pb_id   = pb$pb_id,
    cluster = factor(pb$mn_cluster),
    donor   = factor(donor_lookup[pb$sample_id]),
    sample  = pb$sample_id,
    row.names = pb$pb_id
)

table(meta$donor)
# H19.30.001 H19.30.002 H18.30.002 H18.30.001 
#        125        122        132          2 


# Drop pseudobulks with very few cells
keep_pb <- pb$ncells >= 10
message("Keeping ", sum(keep_pb), " / ", length(keep_pb), " pseudobulk samples")
pb   <- pb[, keep_pb]
meta <- meta[keep_pb, ]

message("Pseudobulk samples per cluster:")
print(table(meta$cluster))
message("Donors: ", length(unique(meta$donor)))

# Drop samples with missing donor
has_donor <- !is.na(meta$donor)
message("Dropping ", sum(!has_donor), " pseudobulks with NA donor")
pb   <- pb[, has_donor]
meta <- meta[has_donor, ]

# Then proceed with DGE
dge <- DGEList(counts = counts(pb))
dge <- dge[filterByExpr(dge, group = meta$cluster), ]
dge <- calcNormFactors(dge)

# ---- dream: mixed model ----
form <- ~ 0 + cluster + (1 | donor)

vobj <- voomWithDreamWeights(dge, form, meta)
fit  <- dream(vobj, form, meta)

# ---- Contrasts: each cluster vs rest ----
contrasts <- makeContrastsDream(form, meta,
    ITC_1_vs_rest = "(clusterITC_1 - (clusterITC_2 + clusterITC_3) / 2)",
    ITC_2_vs_rest = "(clusterITC_2 - (clusterITC_1 + clusterITC_3) / 2)",
    ITC_3_vs_rest = "(clusterITC_3 - (clusterITC_1 + clusterITC_2) / 2)"
)

fit2 <- dream(vobj, form, meta, L = contrasts)

# ---- Extract results ----
de_itc1 <- topTable(fit2, coef = "ITC_1_vs_rest", number = Inf)
de_itc2 <- topTable(fit2, coef = "ITC_2_vs_rest", number = Inf)
de_itc3 <- topTable(fit2, coef = "ITC_3_vs_rest", number = Inf)

message("ITC_1 vs rest: ", sum(de_itc1$adj.P.Val < 0.05), " DEGs (FDR < 0.05)")
message("ITC_2 vs rest: ", sum(de_itc2$adj.P.Val < 0.05), " DEGs (FDR < 0.05)")
message("ITC_3 vs rest: ", sum(de_itc3$adj.P.Val < 0.05), " DEGs (FDR < 0.05)")

# ---- Save ----
de_results <- list(
    ITC_1_vs_rest = de_itc1,
    ITC_2_vs_rest = de_itc2,
    ITC_3_vs_rest = de_itc3
)
saveRDS(de_results, here("processed-data", "snRNAseq", "BICCN", "dream_ITC_DE.rds"))
message("Done.")





plot_volcano <- function(de, title, n_label = 15,
                          fc_thresh = 1, fdr_thresh = 0.05) {
    de$gene <- rownames(de)
    de$neg_log10_fdr <- -log10(de$adj.P.Val)

    de$sig <- case_when(
        de$adj.P.Val < fdr_thresh & de$logFC >  fc_thresh ~ "Up",
        de$adj.P.Val < fdr_thresh & de$logFC < -fc_thresh ~ "Down",
        TRUE ~ "NS"
    )

    # Top genes to label
    top <- de[de$sig == "Up", ]
    top <- head(top[order(-top$logFC), ], n_label)

    ggplot(de, aes(x = logFC, y = neg_log10_fdr, colour = sig)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_label_repel(data = top, aes(label = gene),
                         size = 3, fontface = "italic",
                         max.overlaps = 20, label.size = 0.2,
                         fill = "white", alpha = 0.8,
                         show.legend = FALSE) +
        geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed",
                   colour = "grey40", linewidth = 0.3) +
        geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed",
                   colour = "grey40", linewidth = 0.3) +
        scale_colour_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey80"),
                            name = NULL) +
        labs(x = expression(log[2]~FC), y = expression(-log[10]~FDR),
             title = title,
             subtitle = paste0("Up: ", sum(de$sig == "Up"),
                               "  Down: ", sum(de$sig == "Down"))) +
        ggpubr::theme_pubr(base_size = 13) +
        theme(
            plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            legend.position = "right"
        )
}

p1 <- plot_volcano(de_results$ITC_1_vs_rest, "ITC_1 vs rest")
p2 <- plot_volcano(de_results$ITC_2_vs_rest, "ITC_2 vs rest")
p3 <- plot_volcano(de_results$ITC_3_vs_rest, "ITC_3 vs rest")

pdf(here("plots", "snRNAseq", "ITCs", "dream_volcanos.pdf"),
    width = 15, height = 5)
print(p1 | p2 | p3)
dev.off()
message("Volcano plots saved.")






library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


run_go_top <- function(de, n = 100) {
    de$gene <- rownames(de)
    
    # Top n upregulated genes by logFC, passing FDR
    sig <- de[de$adj.P.Val < 0.05, ]
    sig <- sig[order(-sig$logFC), ]
    top_genes <- head(sig$gene, n)
    
    message("  Using ", length(top_genes), " genes")
    
    enrichGO(
        gene          = top_genes,
        universe      = de$gene,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        minGSSize     = 15,
        maxGSSize     = 500,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.1
    )
}

message("ITC_1:")
go_itc1 <- run_go_top(de_results$ITC_1_vs_rest)
message("ITC_2:")
go_itc2 <- run_go_top(de_results$ITC_2_vs_rest)
message("ITC_3:")
go_itc3 <- run_go_top(de_results$ITC_3_vs_rest)

message("ITC_1: ", nrow(go_itc1), " terms")
message("ITC_2: ", nrow(go_itc2), " terms")
message("ITC_3: ", nrow(go_itc3), " terms")

# ---- Plot top terms ----
plot_go <- function(go_res, title, n = 15) {
    if (is.null(go_res) || nrow(go_res) == 0) {
        return(ggplot() + ggtitle(paste0(title, "\nNo significant terms")) +
                   theme_void())
    }

    dotplot(go_res, showCategory = n) +
        ggtitle(title) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p1 <- plot_go(go_itc1, "ITC_1 vs rest")
p2 <- plot_go(go_itc2, "ITC_2 vs rest")
p3 <- plot_go(go_itc3, "ITC_3 vs rest")

pdf(here("plots", "snRNAseq", "ITCs", "GO_enrichment.pdf"),
    width = 8, height = 10)
print(p1)
print(p2)
print(p3)
dev.off()

# ---- Save results ----
saveRDS(list(ITC_1 = go_itc1, ITC_2 = go_itc2, ITC_3 = go_itc3),
        here("processed-data", "snRNAseq", "BICCN", "GO_ITC_enrichment.rds"))
message("Done.")








library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)
library(enrichplot)

# ---- GSEA per cluster ----
run_gsea <- function(de) {
    de$gene <- rownames(de)

    # Ranked gene list by t-statistic (more stable than logFC)
    gene_list <- setNames(de$t, de$gene)
    gene_list <- sort(gene_list, decreasing = TRUE)

    gseGO(
        geneList     = gene_list,
        OrgDb        = org.Hs.eg.db,
        keyType      = "SYMBOL",
        ont          = "BP",
        minGSSize    = 15,
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        eps          = 1e-300
    )
}

message("Running GSEA...")
gsea_itc1 <- run_gsea(de_results$ITC_1_vs_rest)
gsea_itc2 <- run_gsea(de_results$ITC_2_vs_rest)
gsea_itc3 <- run_gsea(de_results$ITC_3_vs_rest)

message("ITC_1: ", nrow(gsea_itc1), " significant terms")
message("ITC_2: ", nrow(gsea_itc2), " significant terms")
message("ITC_3: ", nrow(gsea_itc3), " significant terms")

# ---- Dotplots: activated vs suppressed ----
plot_gsea_dot <- function(gsea_res, title, n = 15) {
    if (is.null(gsea_res) || nrow(gsea_res) == 0) {
        return(ggplot() + ggtitle(paste0(title, "\nNo significant terms")) +
                   theme_void())
    }
    dotplot(gsea_res, showCategory = n, split = ".sign") +
        facet_grid(~ .sign) +
        ggtitle(title) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

pdf(here("plots", "snRNAseq", "ITCs", "GSEA_dotplots.pdf"),
    width = 12, height = 10)
print(plot_gsea_dot(gsea_itc1, "ITC_1 vs rest"))
print(plot_gsea_dot(gsea_itc2, "ITC_2 vs rest"))
print(plot_gsea_dot(gsea_itc3, "ITC_3 vs rest"))
dev.off()

# ---- Comparative: merge results across clusters ----
go_itc1_df <- as.data.frame(go_itc1) %>% dplyr::mutate(cluster = "ITC_1")
go_itc2_df <- as.data.frame(go_itc2) %>% dplyr::mutate(cluster = "ITC_2")
go_itc3_df <- as.data.frame(go_itc3) %>% dplyr::mutate(cluster = "ITC_3")

all_go <- rbind(go_itc1_df, go_itc2_df, go_itc3_df)

parse_ratio <- function(x) {
    sapply(strsplit(x, "/"), function(r) as.numeric(r[1]) / as.numeric(r[2]))
}

# Fold enrichment
all_go$FoldEnrichment <- parse_ratio(all_go$GeneRatio) / parse_ratio(all_go$BgRatio)

top_terms <- all_go %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_min(p.adjust, n = 10) %>%
    dplyr::pull(Description) %>%
    unique()

fe_df <- all_go %>%
    dplyr::filter(Description %in% top_terms) %>%
    dplyr::select(Description, FoldEnrichment, cluster) %>%
    tidyr::pivot_wider(names_from = cluster, values_from = FoldEnrichment) %>%
    tibble::column_to_rownames("Description")

fe_mat <- as.matrix(fe_df)
fe_mat[is.na(fe_mat)] <- 0

ht <- Heatmap(fe_mat,
    name = "Fold\nEnrichment",
    col  = colorRamp2(c(0, max(fe_mat) / 2, max(fe_mat)),
                      c("#f7f7f7", "#fc8d59", "#d73027")),
    column_title = "GO enrichment: top 100 upregulated genes per cluster",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 11),
    border = TRUE,
    heatmap_legend_param = list(title = "Fold\nEnrichment")
)

pdf(here("plots", "snRNAseq", "ITCs", "GO_top100_heatmap.pdf"),
    width = 8, height = 10)
draw(ht)
dev.off()
message("GO heatmap saved.")








# ---- Extract core enrichment genes from key GO terms ----
library(ComplexHeatmap)
library(circlize)

gsea_results <- list(ITC_1 = gsea_itc1, ITC_2 = gsea_itc2, ITC_3 = gsea_itc3)

# Pick terms to interrogate
terms_of_interest <- c(
    "synapse organization",
    "cell projection morphogenesis",
    "neuron projection morphogenesis",
    "cell-cell adhesion",
    "oxidative phosphorylation",
    "cellular respiration",
    "regulation of nervous system development"
)

# Extract core enrichment genes per term per cluster
extract_core_genes <- function(gsea_res, term) {
    res_df <- as.data.frame(gsea_res)
    row <- res_df[res_df$Description == term, ]
    if (nrow(row) == 0) return(character(0))
    strsplit(row$core_enrichment, "/")[[1]]
}

# Build a list of genes per term
term_genes <- lapply(terms_of_interest, function(term) {
    genes_all <- unique(unlist(lapply(gsea_results, extract_core_genes, term = term)))
    genes_all
})
names(term_genes) <- terms_of_interest

# Print counts
lapply(names(term_genes), function(t) {
    message(t, ": ", length(term_genes[[t]]), " core genes")
})

# ---- Pseudobulk expression per cluster ----

make_term_heatmap <- function(genes, term_name, sce_sub) {
    genes <- genes[genes %in% rownames(sce_sub)]
    if (length(genes) < 3) return(NULL)

    pb <- sapply(levels(factor(sce_sub$mn_cluster)), function(cl) {
        Matrix::rowMeans(logcounts(sce_sub)[genes, sce_sub$mn_cluster == cl, drop = FALSE])
    })

    pb_scaled <- t(scale(t(pb)))
    max_val <- max(abs(pb_scaled), na.rm = TRUE)

    Heatmap(pb_scaled,
        name = "z-score",
        col  = colorRamp2(c(-max_val, 0, max_val), c("#4575b4", "#f7f7f7", "#d73027")),
        column_title = term_name,
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
        column_names_gp = gpar(fontsize = 11),
        border = TRUE,
        show_row_dend = TRUE,
        heatmap_legend_param = list(title = "Scaled\nexpression")
    )
}

# ---- Plot each term ----
pdf(here("plots", "snRNAseq", "ITCs", "GSEA_core_genes_heatmaps.pdf"),
    width = 5, height = 12)

for (term in terms_of_interest) {
    genes <- term_genes[[term]]
    ht <- make_term_heatmap(genes, term, sce_siletti)
    if (!is.null(ht)) {
        draw(ht)
    }
}

dev.off()
message("Core gene heatmaps saved.")

# ---- Also make a combined "greatest hits" heatmap ----
# Top 10 core genes per term by absolute logFC in the relevant cluster
de_results <- readRDS(here("processed-data", "snRNAseq", "BICCN", "dream_ITC_DE.rds"))

pick_top_core <- function(genes, de, n = 10) {
    genes <- genes[genes %in% rownames(de)]
    if (length(genes) == 0) return(character(0))
    sub <- de[genes, ]
    head(rownames(sub[order(-abs(sub$logFC)), ]), n)
}

# ITC_1 enriched terms: oxphos, respiration
itc1_genes <- unique(c(
    pick_top_core(term_genes[["oxidative phosphorylation"]], de_results$ITC_1_vs_rest),
    pick_top_core(term_genes[["cellular respiration"]], de_results$ITC_1_vs_rest)
))

# ITC_3 enriched terms: synapse, projection
itc3_genes <- unique(c(
    pick_top_core(term_genes[["synapse organization"]], de_results$ITC_3_vs_rest),
    pick_top_core(term_genes[["cell projection morphogenesis"]], de_results$ITC_3_vs_rest)
))

combined_genes <- c(itc1_genes, itc3_genes)
gene_group <- c(rep("OxPhos / Respiration\n(ITC_1 enriched)", length(itc1_genes)),
                rep("Synapse / Projection\n(ITC_3 enriched)", length(itc3_genes)))

pb_combo <- sapply(levels(factor(sce_siletti$mn_cluster)), function(cl) {
    rowMeans(logcounts(sce_siletti)[combined_genes, sce_siletti$mn_cluster == cl, drop = FALSE])
})
pb_combo <- t(scale(t(pb_combo)))
max_val <- max(abs(pb_combo), na.rm = TRUE)

ha <- rowAnnotation(
    Term = gene_group,
    col = list(Term = c("OxPhos / Respiration\n(ITC_1 enriched)" = "#d73027",
                         "Synapse / Projection\n(ITC_3 enriched)" = "#4575b4")),
    show_legend = TRUE
)

ht_combo <- Heatmap(pb_combo,
    name = "z-score",
    col  = colorRamp2(c(-max_val, 0, max_val), c("#4575b4", "#f7f7f7", "#d73027")),
    column_title = "Core enrichment genes: ITC_1 vs ITC_3 signatures",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    left_annotation = ha,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_split = factor(gene_group, levels = unique(gene_group)),
    row_gap = unit(3, "mm"),
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    column_names_gp = gpar(fontsize = 11),
    border = TRUE,
    heatmap_legend_param = list(title = "Scaled\nexpression")
)

pdf(here("plots", "snRNAseq", "ITCs", "GSEA_combined_core_heatmap.pdf"),
    width = 5, height = 10)
draw(ht_combo)
dev.off()
message("Combined heatmap saved.")