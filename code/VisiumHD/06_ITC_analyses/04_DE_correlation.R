library(here)
library(scran)
library(scater)
library(ggplot2)

# --------------------------------------------------------------------------
# 1. Load snRNA-seq reference and get logFC between the two TSHZ1 subtypes
# --------------------------------------------------------------------------
load(here("processed-data", "snRNAseq", "Yu_et_al", "yu_sce_gtf.rda"))
sce <- sce.amy

# drop "Humam_" from identity labels to match RCTD reference
sce$ident <- sub("^Human_", "", sce$ident)

## Subset to the two TSHZ1 populations
sce_tshz1 <- sce[, sce$ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C")]
sce_tshz1$ident <- factor(sce_tshz1$ident)

## findMarkers: TSHZ1_CALCRL vs TSHZ1_SEMA3C
markers_snrna <- findMarkers(sce_tshz1, groups = sce_tshz1$ident,
                              test.type = "t", direction = "any",
                              pval.type = "all", lfc = 0)

## Extract logFC — the result for TSHZ1_CALCRL gives logFC relative to TSHZ1_SEMA3C
res_snrna <- as.data.frame(markers_snrna[["TSHZ1 CALCRL"]])
res_snrna$gene <- rownames(res_snrna)
## The logFC column name will be "logFC.TSHZ1_SEMA3C"
snrna_lfc_col <- grep("logFC", colnames(res_snrna), value = TRUE)[1]

# --------------------------------------------------------------------------
# 2. Load VisiumHD + RCTD results, assign labels, then get logFC
# --------------------------------------------------------------------------
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_cells.rds"))
colnames(spe) <- make.unique(colnames(spe))

spe <- logNormCounts(spe)

rctd <- readRDS(here("processed-data", "VisiumHD", "05_label_transfer", "rctd_results_HDcells.rds"))
res_df <- rctd@results$results_df

## Match RCTD labels to spe
m <- match(colnames(spe), rownames(res_df))
spe$rctd_ident <- res_df$first_type[m]

# drop "Human_" from RCTD labels to match snRNA-seq
spe$rctd_ident <- sub("^Human_", "", spe$rctd_ident)

## Subset to the two TSHZ1 populations
spe_tshz1 <- spe[, which(spe$rctd_ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C"))]
spe_tshz1$rctd_ident <- factor(spe_tshz1$rctd_ident)

markers_hd <- findMarkers(spe_tshz1, groups = spe_tshz1$rctd_ident,
                           test.type = "t", direction = "any",
                           pval.type = "all", lfc = 0)

res_hd <- as.data.frame(markers_hd[["TSHZ1 CALCRL"]])
res_hd$gene <- rownames(res_hd)
hd_lfc_col <- grep("logFC", colnames(res_hd), value = TRUE)[1]

# --------------------------------------------------------------------------
# 3. Merge on shared genes and correlate
# --------------------------------------------------------------------------
merged <- merge(res_snrna[, c("gene", snrna_lfc_col)],
                res_hd[, c("gene", hd_lfc_col)],
                by = "gene")
colnames(merged) <- c("gene", "snRNAseq_logFC", "HD_logFC")

cor_val <- cor(merged$snRNAseq_logFC, merged$HD_logFC, method = "pearson")
message("Pearson r = ", round(cor_val, 3))

# --------------------------------------------------------------------------
# 4. Plot (panel U style)
# --------------------------------------------------------------------------
p <- ggplot(merged, aes(x = HD_logFC, y = snRNAseq_logFC)) +
    geom_point(size = 0.5, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("r = ", round(cor_val, 2)), size = 5) +
    labs(x = "HD: AI.1 vs AI.2\n(logFC)",
         y = "snRNAseq: AI.1 vs AI.2\n(logFC)",
         title = "AI subtype correlation") +
    theme_classic(base_size = 14)

ggsave(here("plots", "VisiumHD", "06_ITC_analyses", "TSHZ1_subtype_logFC_correlation.pdf"),
       p, width = 5, height = 5)





# --------------------------------------------------------------------------
# 1. Load snRNA-seq reference
# --------------------------------------------------------------------------

## Build AI flag: TSHZ1_CALCRL and TSHZ1_SEMA3C are the AI subtypes
sce$is_AI <- ifelse(sce$ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C"), "AI", "Other")

## --- Broad: AI vs all ---
markers_broad_snrna <- findMarkers(sce, groups = factor(sce$is_AI),
                                    test.type = "t", direction = "any",
                                    pval.type = "all", lfc = 0)
res_broad_snrna <- as.data.frame(markers_broad_snrna[["AI"]])
res_broad_snrna$gene <- rownames(res_broad_snrna)
broad_snrna_lfc_col <- grep("logFC", colnames(res_broad_snrna), value = TRUE)[1]

## --- Subtype: TSHZ1_CALCRL vs TSHZ1_SEMA3C ---
sce_tshz1 <- sce[, sce$ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C")]
sce_tshz1$ident <- factor(sce_tshz1$ident)

markers_sub_snrna <- findMarkers(sce_tshz1, groups = sce_tshz1$ident,
                                  test.type = "t", direction = "any",
                                  pval.type = "all", lfc = 0)
res_sub_snrna <- as.data.frame(markers_sub_snrna[["TSHZ1 CALCRL"]])
res_sub_snrna$gene <- rownames(res_sub_snrna)
sub_snrna_lfc_col <- grep("logFC", colnames(res_sub_snrna), value = TRUE)[1]

# --------------------------------------------------------------------------
# 2. Load VisiumHD + RCTD results
# --------------------------------------------------------------------------
## Drop cells with no RCTD assignment
spe <- spe[, !is.na(spe$rctd_ident)]

## AI flag
spe$is_AI <- ifelse(spe$rctd_ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C"), "AI", "Other")

## Make sure logcounts exists
if (!"logcounts" %in% assayNames(spe)) spe <- scater::logNormCounts(spe)

## --- Broad: AI vs all ---
markers_broad_hd <- findMarkers(spe, groups = factor(spe$is_AI),
                                 test.type = "t", direction = "any",
                                 pval.type = "all", lfc = 0)
res_broad_hd <- as.data.frame(markers_broad_hd[["AI"]])
res_broad_hd$gene <- rownames(res_broad_hd)
broad_hd_lfc_col <- grep("logFC", colnames(res_broad_hd), value = TRUE)[1]

## --- Subtype: TSHZ1_CALCRL vs TSHZ1_SEMA3C ---
spe_tshz1 <- spe[, spe$rctd_ident %in% c("TSHZ1 CALCRL", "TSHZ1 SEMA3C")]
spe_tshz1$rctd_ident <- factor(spe_tshz1$rctd_ident)

markers_sub_hd <- findMarkers(spe_tshz1, groups = spe_tshz1$rctd_ident,
                               test.type = "t", direction = "any",
                               pval.type = "all", lfc = 0)
res_sub_hd <- as.data.frame(markers_sub_hd[["TSHZ1 CALCRL"]])
res_sub_hd$gene <- rownames(res_sub_hd)
sub_hd_lfc_col <- grep("logFC", colnames(res_sub_hd), value = TRUE)[1]

# --------------------------------------------------------------------------
# 3. Merge and correlate
# --------------------------------------------------------------------------
## Broad
merged_broad <- merge(res_broad_snrna[, c("gene", broad_snrna_lfc_col)],
                       res_broad_hd[, c("gene", broad_hd_lfc_col)],
                       by = "gene")
colnames(merged_broad) <- c("gene", "snRNAseq_logFC", "HD_logFC")
cor_broad <- cor(merged_broad$snRNAseq_logFC, merged_broad$HD_logFC, method = "pearson")

## Subtype
merged_sub <- merge(res_sub_snrna[, c("gene", sub_snrna_lfc_col)],
                     res_sub_hd[, c("gene", sub_hd_lfc_col)],
                     by = "gene")
colnames(merged_sub) <- c("gene", "snRNAseq_logFC", "HD_logFC")
cor_sub <- cor(merged_sub$snRNAseq_logFC, merged_sub$HD_logFC, method = "pearson")

message("Broad (AI vs all) Pearson r = ", round(cor_broad, 3))
message("Subtype (CALCRL vs SEMA3C) Pearson r = ", round(cor_sub, 3))

# --------------------------------------------------------------------------
# 4. Plot
# --------------------------------------------------------------------------
p_broad <- ggplot(merged_broad, aes(x = HD_logFC, y = snRNAseq_logFC)) +
    geom_point(size = 0.5, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("r = ", round(cor_broad, 2)), size = 5) +
    labs(x = "HD: AI vs all\n(logFC)",
         y = "snRNAseq: AI vs all\n(logFC)",
         title = "Broad AI correlation") +
    theme_classic(base_size = 14)

p_sub <- ggplot(merged_sub, aes(x = HD_logFC, y = snRNAseq_logFC)) +
    geom_point(size = 0.5, alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("r = ", round(cor_sub, 2)), size = 5) +
    labs(x = "HD: AI.1 vs AI.2\n(logFC)",
         y = "snRNAseq: AI.1 vs AI.2\n(logFC)",
         title = "AI subtype correlation") +
    theme_classic(base_size = 14)

library(patchwork)
p_combined <- p_broad + p_sub + plot_annotation(tag_levels = "A")

plots_dir <- here("plots", "VisiumHD", "06_ITC_analyses")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(here(plots_dir, "AI_logFC_correlations.pdf"),
       p_combined, width = 10, height = 5)