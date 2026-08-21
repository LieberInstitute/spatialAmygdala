#!/usr/bin/env Rscript
# =============================================================================
# Compositional analysis of cell types across spatial domains using crumblr
# =============================================================================
# 
# Input: Xenium data with spatial domain assignments and predicted cell type labels
# Goal:  Identify cell types enriched in each spatial domain
# 
# Strategy:
#   1. Aggregate cell type counts per FOV (or spatial unit) within each domain
#   2. Apply crumblr CLR transform with precision weights
#   3. Omnibus F-test: does composition vary across domains at all?
#   4. 1-vs-all contrasts: which cell types are enriched in each specific domain?
#   5. Multivariate tree-based testing for extra power
# =============================================================================

library(crumblr)
library(variancePartition)
library(dreamlet)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)


# save directories
processed_dir <- here("processed-data", "Xenium", "09_crumblr")
plot_dir <- here("plots", "Xenium", "09_crumblr")

# load xenium data
spe.old <- readRDS(here("processed-data","Xenium","04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe.new <- readRDS(here("processed-data", "Xenium", "05_spatial_clustering", "Banksy_domains_v1.0.rds"))

# load label predictions (aligned to spe.old)
pred.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype.csv"))

# add predictions to spe.old first, while IDs still match
stopifnot(nrow(pred.broad) == ncol(spe.old))
spe.old$pred_broad_celltype <- pred.broad$pruned.labels
spe.old$pred_fine_celltype <- pred.fine$pruned.labels

# subset to common cells
common_cells <- intersect(colnames(spe.old), colnames(spe.new))
spe.old <- spe.old[, common_cells]
spe.new <- spe.new[, common_cells]

# copy predictions from spe.old to spe.new
spe.new$pred_broad_celltype <- spe.old$pred_broad_celltype
spe.new$pred_fine_celltype <- spe.old$pred_fine_celltype

# rename to spe
spe <- spe.new
rm(spe.old, spe.new)

colnames(spe) <- make.unique(colnames(spe), sep = "-")
rownames(spatialCoords(spe)) <- colnames(spe)


colnames(colData(spe))
#  [1] "cell_id"                          "transcript_counts"               
#  [3] "control_probe_counts"             "control_codeword_counts"         
#  [5] "unassigned_codeword_counts"       "deprecated_codeword_counts"      
#  [7] "total_counts"                     "nucleus_area"                    
#  [9] "sample_id"                        "AspectRatio"                     
# [11] "Area_um"                          "brnum"                           
# [13] "sum"                              "detected"                        
# [15] "total"                            "control_sum"                     
# [17] "control_detected"                 "target_sum"                      
# [19] "target_detected"                  "ctrl_total_ratio"                
# [21] "log2Ctrl_total_ratio"             "log2AspectRatio"                 
# [23] "SignalDensity"                    "log2SignalDensity"               
# [25] "Area_um_outlier_mc"               "Area_um_outlier_sc"              
# [27] "QC_score"                         "low_qcscore"                     
# [29] "cell_area.sf"                     "nucleus_area.sf"                 
# [31] "clust_HARMONY_M0_lam0.8_k50_res2" "Banksy_res2.0_collapsed_v6"      
# [33] "pred_broad_celltype"              "pred_fine_celltype" 

# celltype = "pred_fine_celltype" 
# domain = "Banksy_res2.0_collapsed_v6"
#!/usr/bin/env Rscript

# =============================================================================
# 2. Define variables
# =============================================================================

celltype_col <- "pred_fine_celltype"
domain_col   <- "Banksy_domains"

# remove cells with NA cell type labels (pruned by SingleR)
na_ct <- is.na(colData(spe)[[celltype_col]])
cat("Removing", sum(na_ct), "cells with NA cell type labels\n")
spe <- spe[, !na_ct]

# check what we have
cat("Cell types:\n")
print(sort(table(colData(spe)[[celltype_col]]), decreasing = TRUE))
cat("\nSpatial domains:\n")
print(sort(table(colData(spe)[[domain_col]]), decreasing = TRUE))
cat("\nBrain numbers:\n")
print(table(spe$brnum))


# =============================================================================
# 3. Aggregate to pseudobulk using dreamlet
# =============================================================================
#
# Create a composite sample_id = brnum:domain so that each
# brain × domain combination is a separate "sample".
# This gives biological replication (multiple brains) within each domain.
#
# cluster_id = cell type label (aggregateToPseudoBulk counts cells per type)

spe$sample_domain <- paste0(spe$brnum, "_", colData(spe)[[domain_col]])

pb <- aggregateToPseudoBulk(spe,
  assay      = "counts",
  cluster_id = celltype_col,
  sample_id  = "sample_domain",
  verbose    = FALSE
)

# cellCounts() extracts the cell count matrix: rows = samples, cols = cell types
counts_mat <- cellCounts(pb)
cat("\nPseudobulk cell count matrix:", nrow(counts_mat), "samples x",
    ncol(counts_mat), "cell types\n")

# Extract and augment sample metadata
sample_info <- as.data.frame(colData(pb))
sample_info$brnum  <- sub("_.*", "", rownames(sample_info))
sample_info$domain <- sub("^[^_]+_", "", rownames(sample_info))

cat("\nSamples per domain:\n")
print(table(sample_info$domain))
cat("\nSamples per brain:\n")
print(table(sample_info$brnum))


# =============================================================================
# 4. Apply crumblr transformation
# =============================================================================

# crumblr expects rows = cell types, columns = samples
cobj <- crumblr(counts_mat)

# cobj$E       = CLR-transformed values (cell types x samples)
# cobj$weights = precision weights (inverse variance)

# Build sample_info from the rownames of the count matrix
sample_info <- data.frame(
  sample_id = rownames(counts_mat),
  brnum     = sub("_.*", "", rownames(counts_mat)),
  domain    = sub("^[^_]+_", "", rownames(counts_mat)),
  row.names = rownames(counts_mat)
)

# Then for PCA, just cbind since row order matches
pca <- prcomp(t(standardize(cobj)))
df_pca <- data.frame(pca$x, sample_info[colnames(cobj), ])
# =============================================================================
# 5. PCA of composition
# =============================================================================

# standardize() gives approximately equal sampling variance per observation
# which improves PCA performance (as per crumblr vignette)
pca <- prcomp(t(standardize(cobj)))
df_pca <- data.frame(pca$x, sample_id = rownames(pca$x)) %>%
  left_join(sample_info %>% mutate(sample_id = rownames(sample_info)),
            by = "sample_id")

p_pca <- ggplot(df_pca, aes(PC1, PC2, color = domain, shape = brnum)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "PCA of cell type composition",
       color = "Domain", shape = "Brain")

ggsave(file.path(plot_dir, "composition_pca.pdf"), p_pca,
       width = 8, height = 6)
cat("Saved: composition_pca.pdf\n")


# =============================================================================
# 6. Variance partitioning
# =============================================================================

# How much compositional variation is explained by domain vs brain?
form_vp <- ~ brnum + domain
vp <- fitExtractVarPartModel(cobj, form_vp, sample_info)

p_vp <- plotPercentBars(vp) +
  ggtitle("Variance in cell type composition")

ggsave(file.path(plot_dir, "variance_partition_FixedEffects.pdf"), p_vp,
       width = 8, height = 6)


# =============================================================================
# 7. Differential composition testing
# =============================================================================

# Cell-means model: ~ 0 + domain + brnum
#   - domain as cell-means (no intercept) makes 1-vs-all contrasts clean
#   - brnum controls for brain-to-brain variation
# dream() handles the precision weights from crumblr

sample_info$domain <- factor(sample_info$domain)
sample_info$brnum  <- factor(sample_info$brnum)

fit <- dream(cobj, ~ 0 + domain + (1|brnum), sample_info)
fit <- eBayes(fit)

# --- 7a. Omnibus F-test: does each cell type vary across domains? -------------
domain_coefs <- grep("^domain", colnames(coef(fit)), value = TRUE)
cat("\nDomain coefficients in model:\n")
print(domain_coefs)

omnibus_results <- topTable(fit, coef = domain_coefs, number = Inf,
                            sort.by = "F")

cat("\n=== Omnibus F-test: which cell types vary across domains? ===\n")
print(omnibus_results)

write.csv(omnibus_results,
          file.path(processed_dir, "omnibus_Ftest_results.csv"))


# --- 7b. One-vs-all contrasts for each domain --------------------------------

# Make domain names syntactically valid before modeling
sample_info$domain <- make.names(sample_info$domain)
sample_info$domain <- factor(sample_info$domain)

colnames(fit$coefficients) <- make.names(colnames(fit$coefficients))
colnames(fit$design) <- make.names(colnames(fit$design))

domain_levels <- levels(sample_info$domain)
K <- length(domain_levels)

contrast_strings <- sapply(domain_levels, function(d) {
  others <- setdiff(domain_levels, d)
  d_coef <- paste0("domain", d)
  others_coefs <- paste0("domain", others)
  paste0(d_coef, " - (", paste(others_coefs, collapse = " + "), ")/",
         length(others))
})
names(contrast_strings) <- paste0(domain_levels, "_vs_rest")

contrast_matrix <- makeContrasts(contrasts = contrast_strings,
                                 levels = colnames(coef(fit)))

cat("\nContrast matrix (domain columns):\n")
print(contrast_matrix[domain_coefs, ])

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results for each domain
results_list <- lapply(seq_len(K), function(i) {
  tt <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
  tt$cell_type <- rownames(tt)
  tt$domain <- domain_levels[i]
  tt$contrast <- colnames(contrast_matrix)[i]
  tt
})

results_all <- bind_rows(results_list)

cat("\n=== 1-vs-All results (significant at FDR < 0.05) ===\n")
sig <- results_all %>% filter(adj.P.Val < 0.05) %>% arrange(domain, adj.P.Val)
if (nrow(sig) > 0) {
  print(sig %>% select(domain, cell_type, logFC, adj.P.Val), n = 50)
} else {
  cat("No significant results at FDR < 0.05\n")
  cat("Top results per domain:\n")
  top_per_domain <- results_all %>%
    group_by(domain) %>%
    slice_min(P.Value, n = 3) %>%
    select(domain, cell_type, logFC, P.Value, adj.P.Val)
  print(as.data.frame(top_per_domain))
}

write.csv(results_all,
          file.path(processed_dir, "one_vs_all_results.csv"),
          row.names = FALSE)


# =============================================================================
# 8. Heatmap of effect sizes
# =============================================================================

plot_df <- results_all %>%
  mutate(
    sig = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p_heatmap <- ggplot(plot_df, aes(x = domain, y = cell_type, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "logFC\n(CLR)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Spatial Domain",
    y = "Cell Type",
    title = "Cell type enrichment across spatial domains",
    subtitle = "crumblr 1-vs-all contrasts | * p<0.05, ** p<0.01, *** p<0.001"
  )

ggsave(file.path(plot_dir, "domain_celltype_enrichment_heatmap.pdf"),
       p_heatmap, width = 10, height = 8)
cat("Saved: domain_celltype_enrichment_heatmap.pdf\n")


# =============================================================================
# 9. Multivariate tree-based testing
# =============================================================================

# Build cell type hierarchy from pseudobulk gene expression
# following the dreamlet integration vignette:
# buildClusterTreeFromPB() computes hierarchical clustering
# from the pseudobulked expression profiles
hcl <- buildClusterTreeFromPB(pb)

# Plot the tree itself
pdf(file.path(plot_dir, "celltype_hierarchy.pdf"), width = 10, height = 6)
plot(hcl, main = "Cell type hierarchy from pseudobulk expression")
dev.off()

# Run treeTest for each 1-vs-all contrast
tree_results <- list()
for (i in seq_len(K)) {
  coef_name <- colnames(contrast_matrix)[i]
  
  tryCatch({
    res_tree <- treeTest(fit2, cobj, hcl, coef = coef_name)
    res_tree$contrast <- coef_name
    tree_results[[coef_name]] <- res_tree
    
    # Plot tree with p-values
    p_tree <- plotTreeTest(res_tree) +
      ggtitle(paste("Tree test:", domain_levels[i], "vs rest"))
    ggsave(file.path(plot_dir, paste0("treeTest_", domain_levels[i], ".pdf")),
           p_tree, width = 10, height = 6)
    
    # Plot tree with regression coefficients (logFC)
    p_beta <- plotTreeTestBeta(res_tree) +
      ggtitle(paste(domain_levels[i], "vs rest"))
    ggsave(file.path(plot_dir, paste0("treeTestBeta_", domain_levels[i], ".pdf")),
           p_beta, width = 10, height = 6)
    
    cat("Saved tree plots for", domain_levels[i], "\n")
  }, error = function(e) {
    cat("treeTest failed for", coef_name, ":", conditionMessage(e), "\n")
  })
}

saveRDS(tree_results,
        file.path(processed_dir, "treeTest_results.rds"))


# =============================================================================
# 10. Save all outputs
# =============================================================================

saveRDS(cobj, file.path(processed_dir, "crumblr_cobj.rds"))
saveRDS(fit, file.path(processed_dir, "crumblr_fit_cellmeans.rds"))
saveRDS(fit2, file.path(processed_dir, "crumblr_fit_contrasts.rds"))
saveRDS(pb, file.path(processed_dir, "pseudobulk_pb.rds"))

cat("\n=== Done! ===\n")
sessionInfo()