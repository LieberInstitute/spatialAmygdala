library(here)
library(SpatialExperiment)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# ==============================================================================
# Configuration
# ==============================================================================
spe <- readRDS(here("processed-data", "Visium", "07_clustering", "BayesSpace",
                     "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))

results_dir <- here("processed-data", "Visium", "14_smoothie_coexpression", "results")
plots_dir <- here("plots", "Visium", "14_smoothie_coexpression", "pathway_analysis")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

modules_df <- read.csv(file.path(results_dir, "modules_df_0.8_9.csv"))

clust_col <- "BS_k16_Semisupervised_wAI"
keep_spots <- !is.na(spe[[clust_col]])
spe <- spe[, keep_spots]
domains <- factor(spe[[clust_col]])

# After setting domains, drop unwanted ones
exclude <- c("Endothelial", "WM.1", "WM.2")
keep_dom <- !domains %in% exclude
spe <- spe[, keep_dom]
domains <- droplevels(domains[keep_dom])

# ==============================================================================
# Modules
# ==============================================================================
analyze_modules <- c(
    "6"  = "LA",
    "27" = "BLD",
    "5"  = "PL",
    "10" = "CoA",
    "13" = "CeA",
    "19" = "MeA",
    "20" = "AI",
    "16" = "HPC"
)

# ==============================================================================
# Run GO MF enrichment (or load from saved results)
# ==============================================================================
go_mf_csv <- file.path(plots_dir, "GO_MF_results.csv")

if (file.exists(go_mf_csv)) {
    cat("Loading saved GO MF results...\n")
    mf_df <- read.csv(go_mf_csv)
} else {
    cat("Running GO MF enrichment...\n")
    gene_list <- lapply(names(analyze_modules), function(mod) {
        modules_df$name[modules_df$module_label == as.integer(mod)]
    })
    names(gene_list) <- paste0("M", names(analyze_modules), " - ", analyze_modules)

    cc_go_mf <- compareCluster(
        gene_list,
        fun = "enrichGO",
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
    )
    cc_go_mf_simp <- simplify(cc_go_mf, cutoff = 0.7)
    mf_df <- as.data.frame(cc_go_mf_simp)
    write.csv(mf_df, go_mf_csv, row.names = FALSE)
}

# ==============================================================================
# Select terms — EDIT THESE
# ==============================================================================
selected_terms <- c(
    "gated channel activity",
    "channel regulator activity",
    "voltage-gated channel activity",
    "calmodulin-dependent protein kinase activity",
    "G-protein alpha-subunit binding",
    "G protein-coupled amine receptor activity",
    "G protein-coupled peptide receptor activity",
    "peptide receptor activity"
)

# ==============================================================================
# Extract genes per term
# ==============================================================================
term_genes <- mf_df %>%
    filter(Description %in% selected_terms) %>%
    dplyr::select(Description, geneID) %>%
    distinct() %>%
    mutate(gene = strsplit(geneID, "/")) %>%
    tidyr::unnest(gene) %>%
    dplyr::select(Description, gene) %>%
    distinct()

# Remove genes not in SPE
term_genes <- term_genes %>% filter(gene %in% rownames(spe))

cat("Terms:", length(unique(term_genes$Description)), "\n")
cat("Total genes:", length(unique(term_genes$gene)), "\n")
for (t in selected_terms) {
    n <- sum(term_genes$Description == t)
    cat("  ", t, ":", n, "genes\n")
}

# ==============================================================================
# Compute mean logcounts per domain
# ==============================================================================
lc <- logcounts(spe)
all_genes_unique <- unique(term_genes$gene)

domain_expr <- matrix(NA, nrow = length(all_genes_unique), ncol = nlevels(domains))
rownames(domain_expr) <- all_genes_unique
colnames(domain_expr) <- levels(domains)

for (d in levels(domains)) {
    idx <- which(domains == d)
    domain_expr[, d] <- rowMeans(as.matrix(lc[all_genes_unique, idx, drop = FALSE]))
}

# Z-score each gene across domains
scaled_expr <- t(scale(t(domain_expr)))

# ==============================================================================
# Top 5 genes per term by max absolute z-score
# ==============================================================================
gene_max_z <- data.frame(gene = rownames(scaled_expr),
                          max_z = apply(abs(scaled_expr), 1, max, na.rm = TRUE))
top_genes <- term_genes %>%
    left_join(gene_max_z, by = "gene") %>%
    group_by(Description) %>%
    slice_max(max_z, n = 5, with_ties = FALSE) %>%
    ungroup()
term_genes <- top_genes %>% dplyr::select(Description, gene)

# ==============================================================================
# Order genes by term
# ==============================================================================
gene_order <- term_genes %>%
    mutate(Description = factor(Description, levels = selected_terms)) %>%
    arrange(Description, gene) %>%
    pull(gene) %>%
    unique()

scaled_expr <- scaled_expr[gene_order, , drop = FALSE]

# Row split by term
gene_to_term <- data.frame(gene = gene_order) %>%
    left_join(
        term_genes %>%
            mutate(Description = factor(Description, levels = selected_terms)) %>%
            arrange(Description) %>%
            distinct(gene, .keep_all = TRUE),
        by = "gene"
    )

row_split <- factor(gene_to_term$Description, levels = selected_terms)

# ==============================================================================
# Annotations
# ==============================================================================
# Term colors for row annotation
term_colors <- c(
    "gated channel activity" = "#F4B400",
    "channel regulator activity" = "#F4B400",
    "voltage-gated channel activity" = "#f1e438ff",
    "calmodulin-dependent protein kinase activity" = "#f1e438ff",
    "G-protein alpha-subunit binding" = "#197d43ff",
    "G protein-coupled amine receptor activity" = "#197d43ff",
    "G protein-coupled peptide receptor activity" = "#baf739ff",
    "peptide receptor activity" = "#baf739ff"
)

row_ha <- rowAnnotation(
    Term = anno_simple(
        as.character(gene_to_term$Description),
        col = term_colors
    ),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
        title = "GO Term",
        title_gp = gpar(fontsize = 9, fontface = "bold"),
        labels_gp = gpar(fontsize = 7)
    )
)

# Domain colors for column annotation
full_pal <- c(
    AI = "#D62728", BM = "#E67E22", BLD = "#9B59B6", PL = "#f1e438ff",
    BL = "#035185ff", LA = "#F4B400", CoA = "#5DA5DA", CeA = "#197d43ff",
    MeA = "#baf739ff", HPC = "#d6a8f8ff", CHAT = "#A0522D",
    Endothelial = "#444444", WM.1 = "#BBBBBB", WM.2 = "#DDDDDD", CLA = "#FF69B4"
)

col_colors <- sapply(colnames(scaled_expr), function(d) {
    if (d %in% names(full_pal)) full_pal[d] else "grey70"
})

col_ha <- HeatmapAnnotation(
    Domain = anno_simple(
        colnames(scaled_expr),
        col = setNames(col_colors, colnames(scaled_expr))
    ),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

# ==============================================================================
# Color scale
# ==============================================================================
col_fun <- colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")
)

# ==============================================================================
# Plot — genes on Y, domains on X (tall format)
# ==============================================================================
pdf(file.path(plots_dir, "GO_MF_gene_heatmap_long.pdf"),
    width = max(6, ncol(scaled_expr) * 0.4 + 3),
    height = max(8, nrow(scaled_expr) * 0.25 + 4))

ht <- Heatmap(
    scaled_expr,
    name = "Z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    clustering_method_columns = "ward.D2",
    row_split = row_split,
    row_gap = unit(2, "mm"),
    row_title_gp = gpar(fontsize = 7, fontface = "bold"),
    row_title_rot = 0,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    left_annotation = row_ha,
    bottom_annotation = col_ha,
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
        title = "Scaled\nexpression",
        legend_height = unit(3, "cm")
    )
)

draw(ht, merge_legend = TRUE,
     column_title = "Gene Expression Underlying Selected GO MF Terms",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

cat("\nSaved: GO_MF_gene_heatmap_long.pdf\n")