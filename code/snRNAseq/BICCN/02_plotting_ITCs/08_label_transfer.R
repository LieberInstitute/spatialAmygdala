## 05_symphony_label_transfer.R
## Transfer ITC subtype labels from BICCN reference to query using Symphony

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    library(scater)
    library(symphony)
    library(harmony)
    library(ggplot2)
    library(patchwork)
})

# ---- Load integrated object with labels ----
sce <- readRDS(here("processed-data", "snRNAseq", "BICCN", "sce_ITC_harmony.rds"))

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


ref_idx   <- sce$dataset == "BICCN"
sce_ref   <- sce[, ref_idx]
sce_query <- sce[, !ref_idx]

stopifnot(!is.null(sce_ref$mn_cluster))
message("Reference: ", ncol(sce_ref), " cells")
message("Query:     ", ncol(sce_query), " cells")
print(table(sce_ref$mn_cluster))

# ---- HVGs (same as integration script) ----
dec  <- modelGeneVar(sce, block = sce$dataset)
hvgs <- getTopHVGs(dec, n = 2000)

# ---- 1. Scale reference expression (center + scale per gene) ----
ref_exp <- as.matrix(logcounts(sce_ref)[hvgs, ])

vargenes_means_sds <- data.frame(
    symbol = hvgs,
    mean   = rowMeans(ref_exp),
    stddev = apply(ref_exp, 1, sd)
)

# Drop zero-variance genes
keep <- vargenes_means_sds$stddev > 1e-6
hvgs <- hvgs[keep]
vargenes_means_sds <- vargenes_means_sds[keep, ]
ref_exp <- ref_exp[keep, ]

ref_exp_scaled <- t(scale(t(ref_exp)))  # matches what Symphony expects

# ---- 2. PCA on scaled data ----
pca_res <- irlba::prcomp_irlba(t(ref_exp_scaled), n = 20, center = FALSE)
loadings <- pca_res$rotation   # genes x PCs
Z_pca    <- pca_res$x          # cells x PCs

# Store in SCE for consistency
reducedDim(sce_ref, "PCA") <- Z_pca

# ---- 3. Run Harmony on reference PCs ----
ref_metadata <- data.frame(
    cell_id    = colnames(sce_ref),
    sample_id  = sce_ref$sample_id,
    mn_cluster = sce_ref$mn_cluster
)

ref_harmObj <- harmony::HarmonyMatrix(
    data_mat         = reducedDim(sce_ref, "PCA"),
    meta_data        = ref_metadata,
    vars_use         = "sample_id",
    theta            = 2,
    nclust           = 100,
    max.iter.harmony = 20,
    return_object    = TRUE,
    do_pca           = FALSE
)

# ---- 4. Build Symphony reference ----
reference <- buildReferenceFromHarmonyObj(
    ref_harmObj,
    ref_metadata,
    vargenes_means_sds,
    loadings,
    verbose      = TRUE,
    do_umap      = FALSE  # we'll compute UMAP ourselves
)
message("Reference built.")

# ---- 5. Map query cells ----
query_exp <- as.matrix(logcounts(sce_query)[hvgs, ])

query_metadata <- data.frame(
    cell_id   = colnames(sce_query),
    sample_id = sce_query$sample_id
)

query <- mapQuery(
    query_exp,
    query_metadata,
    reference,
    vars         = "sample_id",
    do_normalize = FALSE,
    do_umap      = FALSE
)
message("Query mapped.")

# ---- 6. Predict labels with kNN ----
query <- knnPredict(
    query,
    reference,
    train_labels = reference$meta_data$mn_cluster,
    k = 25,
    confidence = TRUE
)

query_labels <- query$meta_data$cell_type_pred_knn
query_conf   <- query$meta_data$cell_type_pred_knn_prob

message("Label transfer complete.")
print(table(query_labels))

# ---- 7. Add predictions back to full SCE ----
sce$mn_cluster_transfer <- NA_character_
sce$transfer_confidence <- NA_real_

sce$mn_cluster_transfer[ref_idx]  <- as.character(sce_ref$mn_cluster)
sce$transfer_confidence[ref_idx]  <- 1.0
sce$mn_cluster_transfer[!ref_idx] <- query_labels
sce$transfer_confidence[!ref_idx] <- query_conf

# Remap the query labels
sce$mn_cluster_transfer <- dplyr::recode(sce$mn_cluster_transfer,
    "1" = "ITC_1", "2" = "ITC_2", "3" = "ITC_3"
)

table(sce$mn_cluster_transfer, sce$dataset)

# ---- 8. Combined UMAP in Symphony space ----
# Combine harmonized embeddings: reference Z_corr + query Z
ref_emb   <- t(reference$Z_corr)   # cells x dims
query_emb <- t(query$Z)            # cells x dims
combined_emb <- rbind(ref_emb, query_emb)
rownames(combined_emb) <- c(colnames(sce_ref), colnames(sce_query))
combined_emb <- combined_emb[colnames(sce), ]

# Store in SCE and run UMAP via scater
reducedDim(sce, "SYMPHONY") <- combined_emb
set.seed(49287)
sce <- runUMAP(sce, dimred = "SYMPHONY", name = "UMAP_symphony")


# ---- 9. Diagnostic plots ----
message("Generating plots...")

itc_cols <- c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")

umap_df <- data.frame(
    UMAP1      = reducedDim(sce, "UMAP_symphony")[, 1],
    UMAP2      = reducedDim(sce, "UMAP_symphony")[, 2],
    cluster    = sce$mn_cluster_transfer,
    dataset    = sce$dataset,
    confidence = sce$transfer_confidence,
    source     = ifelse(ref_idx, "Reference", "Query")
)
set.seed(42)
umap_df <- umap_df[sample(nrow(umap_df)), ]

p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_colour_manual(values = itc_cols) +
    labs(title = "Symphony label transfer", colour = "ITC subtype") +
    theme_minimal(base_size = 13) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p2 <- p1 + facet_wrap(~ source) + ggtitle("Reference vs Query")
p3 <- p1 + facet_wrap(~ dataset) + ggtitle("By dataset")

p4 <- ggplot(umap_df[umap_df$source == "Query", ],
             aes(x = UMAP1, y = UMAP2, colour = confidence)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_colour_viridis_c(option = "magma", direction = -1) +
    labs(title = "Transfer confidence (query)", colour = "Prob") +
    theme_minimal(base_size = 13) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank())

p5 <- ggplot(umap_df[umap_df$source == "Query", ],
             aes(x = confidence, fill = cluster)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = itc_cols) +
    labs(title = "Confidence distribution", x = "Max kNN probability", y = "Cells") +
    theme_minimal(base_size = 13)

pdf(here("plots", "snRNAseq", "ITCs", "symphony_label_transfer.pdf"),
    width = 12, height = 10)
print((p1 | p4) / p2 / p3)
print(p5)
dev.off()
message("Plots saved.")


# ---- Stacked bar plots: ITC subtype proportions per sample ----

prop_df <- as.data.frame(table(
    dataset = sce$dataset,
    cluster = sce$mn_cluster_transfer
))

prop_df <- do.call(rbind, lapply(split(prop_df, prop_df$dataset), function(d) {
    d$pct <- d$Freq / sum(d$Freq) * 100
    d
}))

prop_df$cluster <- factor(prop_df$cluster, levels = c("ITC_1", "ITC_2", "ITC_3"))

itc_cols <- c("ITC_1" = "#fd0d00", "ITC_2" = "#f5b6b3", "ITC_3" = "#5d0500")

p_bar <- ggplot(prop_df, aes(x = dataset, y = pct, fill = cluster)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = itc_cols) +
    labs(y = "% of cells", x = NULL, fill = "ITC subtype",
         title = "ITC subtype composition by dataset") +
    theme_classic(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", size = 14),
        axis.text.x        = element_text(size = 11),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank()
    )

pdf(here("plots", "snRNAseq", "ITCs", "stacked_barplot_itc_by_dataset.pdf"),
    width = 4, height = 5)
print(p_bar)
dev.off()

message("Stacked bar plot saved.")



# ---- 11. Save ----
saveRDS(sce, here("processed-data", "snRNAseq", "BICCN",
                  "sce_ITC_harmony_transferred.rds"))
saveRDS(reference, here("processed-data", "snRNAseq", "BICCN",
                         "symphony_reference.rds"))
message("Done.")
sessioninfo::session_info()