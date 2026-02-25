library(NMFscape)
library(here)
library(SingleCellExperiment)
library(Seurat)

# set random seed for reproducibility
set.seed(12345)

# load zhou RDS
sce <- readRDS(here("processed-data","snRNAseq","Zhou_cocaine","Zhou_sce_NMF50.rds"))
sce
# class: SingleCellExperiment 
# dim: 17297 163003 
# metadata(1): NMF_basis
# assays(2): counts logcounts
# rownames(17297): AABR07000156.1 Lrp11 ... AABR07043200.1 Pomp
# rowData names(0):
# colnames(163003): AAACCCAAGAAACCCG-1_1 AAACCCACAAAGCACG-1_1 ...
#   TTTGTTGTCTTCGTAT-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(36): orig.ident nCount_RNA ... batch ident
# reducedDimNames(1): NMF
# mainExpName: RNA
# altExpNames(2): SCT integrated

# ======= Find DEPs =======

deps_results <- FindAllDEPs(
  sce,
  "ident",
  nmf_name = "NMF",
  test = "wilcox",
  pval.type = "some",
  log.p = FALSE
)
deps_results
# List of length 13
# names(13): Cck+/Vip+ Astrocytes Oligodendrocytes ... Reln+ Endothelial Pvalb+


# Volcano
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_DEPs_volcano.pdf"), width=10, height=8)
plotDEPsDots(
  deps_results,
  fold_threshold = 1.5,
  pval_threshold = 0.05
)
dev.off()



# ======= Find cocaine DEPs =======

cocaine_deps <- FindAllDEPs(
  sce,
  "addiction.index",
  nmf_name = "NMF",
  test = "wilcox",
  pval.type = "some",
  log.p = FALSE
)
cocaine_deps







library(dplyr)
library(tidyr)
library(rstatix)

# Get NMF scores
nmf_scores <- reducedDim(sce, "NMF")

# Build data frame
df <- data.frame(
  nmf_scores,
  cell_type = colData(sce)$ident,
  condition = colData(sce)$addiction.index
)

# Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("NMF_"),
    names_to = "factor",
    values_to = "score"
  )

# Check which groups have enough cells
group_counts <- df_long %>%
  group_by(cell_type, factor, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n, values_fill = 0)

# Look at it
print(group_counts)

# Filter to cell types with enough cells in each condition (e.g., >= 10)
min_cells <- 10

valid_groups <- df_long %>%
  group_by(cell_type, condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
  filter(if_all(-cell_type, ~ . >= min_cells)) %>%
  pull(cell_type)

df_filtered <- df_long %>%
  filter(cell_type %in% valid_groups)

# Now run the tests
results <- df_filtered %>%
  group_by(cell_type, factor) %>%
  wilcox_test(score ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

effect_sizes <- df_filtered %>%
  group_by(cell_type, factor) %>%
  wilcox_effsize(score ~ condition)

results_full <- results %>%
  left_join(effect_sizes, by = c("cell_type", "factor", "group1", "group2", ".y."))

# Filter for drug-responsive factors
drug_responsive <- results_full %>%
  filter(p.adj < 0.05, abs(effsize) > 0.3) %>%
  arrange(cell_type, p.adj)

drug_responsive



pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_All_NMF_Factors_Heatmap.pdf"), width=12, height=10)
results_full %>%
  filter(group2 == "none") %>%
  ggplot(aes(x = cell_type, y = factor, fill = effsize)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.adj < 0.05 & abs(effsize) > 0.3, "*", "")), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~group1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "NMF factor changes with cocaine exposure by cell type",
    subtitle = "* = FDR < 0.05 and |effect size| > 0.3",
    fill = "Effect size"
  )
dev.off()











pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_Cocaine_DEPs_volcano.pdf"), width=10, height=8)
plotDEPsDots(
  cocaine_deps,
  fold_threshold = 0.1,
  pval_threshold = 0.05
)
dev.off()


top_genes <- getTopFeatures(sce, name = "NMF", n = 100)
weights <- getBasis(sce, name = "NMF")
head(weights)
#                      NMF_41       NMF_42       NMF_43       NMF_44       NMF_45
# AABR07000156.1 3.705804e-06 6.222450e-06 0.000000e+00 1.365481e-07 0.000000e+00
# Lrp11          0.000000e+00 2.717238e-06 3.537897e-04 3.056219e-04 9.117222e-05
# Pcmt1          0.000000e+00 4.331683e-05 0.000000e+00 0.000000e+00 8.883517e-05
# Nup43          5.913128e-06 2.469178e-05 0.000000e+00 8.147241e-06 9.293978e-05
# Lats1          2.108658e-04 9.294333e-05 2.128215e-05 1.007896e-05 9.922039e-05
# Katna1         8.602765e-05 6.379303e-05 1.163818e-04 2.336913e-05 4.572465e-05
#                      NMF_46       NMF_47       NMF_48       NMF_49       NMF_50
# AABR07000156.1 0.000000e+00 0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00
# Lrp11          0.000000e+00 0.0006307106 4.311436e-04 0.000000e+00 0.000000e+00
# Pcmt1          1.237874e-04 0.0000000000 0.000000e+00 1.346548e-04 1.951010e-04
# Nup43          1.044455e-06 0.0000000000 5.020552e-05 5.089717e-06 0.000000e+00
# Lats1          1.265723e-04 0.0000000000 2.487062e-04 5.609086e-05 0.000000e+00
# Katna1         1.217979e-04 0.0000000000 2.691757e-05 1.697504e-05 5.290427e-05


library(clusterProfiler)
library(msigdbr)
library(enrichplot)

weights <- getBasis(sce, name = "NMF")

nmf3_weights <- weights[, "NMF_15"]
nmf3_ranks <- sort(nmf3_weights, decreasing = TRUE)

# GO BP pathways
msigdb_rat <- msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "GO:BP")
m_t2g <- msigdb_rat %>%
  dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(
  geneList = nmf3_ranks,
  TERM2GENE = m_t2g,
  pvalueCutoff = 0.05,
  scoreType = "pos",
  verbose = FALSE
)

# 1. Quick overview - what's enriched?
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_GO_BP_dotplot.pdf"), width=8, height=6)
dotplot(gsea_res, showCategory = 15)
dev.off()

# 2. Sanity check top hit with classic GSEA plot
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_GO_BP_gseaplot2.pdf"), width=8, height=6)
gseaplot2(gsea_res, geneSetID = 1, pvalue_table = TRUE)
dev.off()

# 3. If too many redundant terms, cluster them
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_GO_BP_treeplot.pdf"), width=8, height=6)
gsea_res_sim <- pairwise_termsim(gsea_res)
treeplot(gsea_res_sim)
dev.off()










library(clusterProfiler)
library(ReactomePA)

weights <- getBasis(sce, name = "NMF")

nmf3_weights <- weights[, "NMF_43"]
nmf3_ranks <- sort(nmf3_weights, decreasing = TRUE)

# ReactomePA needs Entrez IDs
library(org.Rn.eg.db)

gene_map <- bitr(
  names(nmf3_ranks),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)

# Rebuild ranked list with Entrez IDs
nmf3_entrez <- nmf3_ranks[names(nmf3_ranks) %in% gene_map$SYMBOL]
names(nmf3_entrez) <- gene_map$ENTREZID[match(names(nmf3_entrez), gene_map$SYMBOL)]

# Run GSEA with Reactome
gsea_reactome <- gsePathway(
  geneList = nmf3_entrez,
  organism = "rat",
  pvalueCutoff = 0.05,
  scoreType = "pos",
  verbose = FALSE
)

# Same workflow
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_dotplot.pdf"), width=8, height=6)
dotplot(gsea_reactome, showCategory = 15)
dev.off()

pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_gseaplot2.pdf"), width=8, height=6)
gseaplot2(gsea_reactome, geneSetID = 1, pvalue_table = TRUE)
dev.off()








library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(org.Rn.eg.db)

# =============================================================================
# 1. Look at top genes directly
# =============================================================================

weights <- getBasis(sce, name = "NMF")
nmf3_weights <- weights[, "NMF_43"]
nmf3_ranks <- sort(nmf3_weights, decreasing = TRUE)

# What are the actual top genes?
head(nmf3_ranks, 30)

# =============================================================================
# 2. Run GSEA
# =============================================================================

gene_map <- bitr(
  names(nmf3_ranks),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)

nmf3_entrez <- nmf3_ranks[names(nmf3_ranks) %in% gene_map$SYMBOL]
names(nmf3_entrez) <- gene_map$ENTREZID[match(names(nmf3_entrez), gene_map$SYMBOL)]

gsea_reactome <- gsePathway(
  geneList = nmf3_entrez,
  organism = "rat",
  pvalueCutoff = 0.05,
  scoreType = "pos",
  eps = 0,  # better p-value estimation
  verbose = FALSE
)

# =============================================================================
# 3. Handle pathway redundancy (Reactome version)
# =============================================================================

gsea_reactome_sim <- pairwise_termsim(gsea_reactome)

# Treeplot - clusters similar pathways
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_treeplot.pdf"), width=8, height=6)
treeplot(gsea_reactome_sim, showCategory = 20) +
  ggtitle("NMF_3: Reactome Pathway Clusters")
dev.off()

# Filtered dotplot - remove very specific/broad terms
gsea_filtered <- gsea_reactome
gsea_filtered@result <- gsea_reactome@result %>%
  filter(p.adjust < 0.05) %>%
  filter(setSize > 50 & setSize < 500)

pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_filtered_dotplot.pdf"), width=8, height=6)
dotplot(gsea_filtered, showCategory = 12) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  ggtitle("NMF_3: Reactome Pathways (filtered)")
dev.off()

# =============================================================================
# 4. Leading edge analysis - find hub genes
# =============================================================================

leading_edge_genes <- gsea_reactome@result %>%
  filter(p.adjust < 0.05) %>%
  pull(core_enrichment) %>%
  strsplit("/") %>%
  unlist()

leading_edge_counts <- table(leading_edge_genes) %>%
  sort(decreasing = TRUE) %>%
  as.data.frame() %>%
  setNames(c("ENTREZID", "n_pathways"))

leading_edge_counts$SYMBOL <- gene_map$SYMBOL[
  match(leading_edge_counts$ENTREZID, gene_map$ENTREZID)
]

# Top hub genes
head(leading_edge_counts, 20)

# Plot hub genes
top_hubs <- head(leading_edge_counts, 20)
pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_LeadingEdge_HubGenes.pdf"), width=8, height=6)
ggplot(top_hubs, aes(x = reorder(SYMBOL, n_pathways), y = n_pathways)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of enriched pathways",
    title = "NMF_3: Hub genes across Reactome pathways"
  ) +
  theme_minimal()
dev.off()

# =============================================================================
# 5. Gene-concept network
# =============================================================================

pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_cnetplot.pdf"), width=8, height=6)
cnetplot(
  gsea_filtered,
  showCategory = 5,
  cex_label_gene = 0.6,
  cex_label_category = 0.8
) +
  ggtitle("NMF_3: Gene-Pathway Network")
dev.off()

# =============================================================================
# 6. Heatmap of genes across pathways
# =============================================================================

pdf(here("plots","snRNAseq","Zhou_cocaine","NMF50_NMF3_GSEA_Reactome_heatplot.pdf"), width=8, height=6)
heatplot(gsea_filtered, showCategory = 10) +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))
dev.off()

# =============================================================================
# 7. Summary table
# =============================================================================

nmf3_summary <- gsea_reactome@result %>%
  filter(p.adjust < 0.05) %>%
  select(Description, setSize, enrichmentScore, NES, pvalue, p.adjust) %>%
  arrange(p.adjust)

nmf3_summary