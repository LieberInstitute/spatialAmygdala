#  ---------------------------------------------------------------------------
#  Zhou et al. cocaine snRNA-seq — DE reanalysis
#
#  - Aggregate to pseudobulk at Yu spatial-label resolution
#    (BLA / CEA / COA-MEA / IA / Cortical-IN-like / non-neuron)
#  - Identify and exclude outlier samples via sample-sample correlation
#  - Test cocaine vs naive per subregion with dreamlet (limma-voom on
#    pseudobulk), with batch as a fixed effect
#  ---------------------------------------------------------------------------

library(here)
library(Seurat)
library(SingleCellExperiment)
library(dreamlet)
library(variancePartition)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggrepel)
library(matrixStats)
library(pheatmap)

set.seed(12345)

# ===========================================================================
# Config
# ===========================================================================
LABEL_COL <- "yu_space_predicted.id"
SCORE_COL <- "yu_space_prediction.score.max"
RUN_TAG   <- "yu_space"

out_plots <- here("plots","snRNAseq","dreamlet")
out_data  <- here("processed-data","snRNAseq","dreamlet")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_data,  recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Load + filter cells by label-transfer confidence
# ===========================================================================
rat.amy <- readRDS(here("processed-data","snRNAseq","zhou_with_yu_labels.rds"))
sce <- as.SingleCellExperiment(rat.amy, assay = "RNA")
assays(sce) <- list(counts = counts(sce))

keep <- colData(sce)[[SCORE_COL]] > 0.5
message("Cells passing ", SCORE_COL, " > 0.5: ",
        sum(keep), " / ", ncol(sce),
        " (", round(100 * mean(keep), 1), "%)")
sce <- sce[, keep]

# Factor setup
sce$treatment <- factor(sce$treatment, levels = c("naive", "cocaine"))
sce$sample    <- factor(sce$sample)
sce$batch     <- factor(sce$batch)
colData(sce)[[LABEL_COL]] <- factor(colData(sce)[[LABEL_COL]])

print(table(colData(sce)[[LABEL_COL]]))

# ===========================================================================
# 2. Pseudobulk aggregation
# ===========================================================================
pb <- aggregateToPseudoBulk(
  sce,
  assay      = "counts",
  cluster_id = LABEL_COL,
  sample_id  = "sample",
  verbose    = TRUE
)
colData(pb)$sample <- rownames(colData(pb))

# Initial process (used only for QC — outliers identified before final fit)
res.proc_qc <- processAssays(
  pb,
  formula          = ~ treatment + batch,
  min.cells        = 20,
  normalize.method = "TMM"
)

cd_pb <- as.data.frame(colData(pb))
pdf(file.path(out_plots, "00_qc_pca.pdf"), width = 9, height = 7)
for (ct in assayNames(res.proc_qc)) {
  vobj <- assay(res.proc_qc, ct)$E
  vobj <- vobj[rowVars(vobj) > 0, , drop = FALSE]
  pc <- prcomp(t(vobj), scale. = TRUE)
  df <- as.data.frame(pc$x[, 1:2]) |>
    rownames_to_column("sample") |>
    mutate(treatment = cd_pb$treatment[match(sample, cd_pb$sample)],
           batch     = cd_pb$batch[match(sample, cd_pb$sample)])
  p <- ggplot(df, aes(PC1, PC2, color = treatment, shape = batch)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = sample), size = 2.5, max.overlaps = 20) +
    theme_bw() + ggtitle(paste0("PCA — ", ct))
  print(p)
}
dev.off()
 
# ===== Drop outliers =====
# Identified from PCA above — these three samples sit far from the main
# cluster across all subregions. Likely repurposed from other studies.
drop_samples <- c("Rat_Opioid_HS_1", "Rat_Opioid_HS_2", "Rat_Amygdala_787A_all_seq")
pb_clean <- pb[, !rownames(colData(pb)) %in% drop_samples]
 
message("Samples retained: ", ncol(pb_clean), " / ", ncol(pb))
print(table(colData(pb_clean)$treatment, colData(pb_clean)$batch))
 
# ===== Re-process + fit =====
res.proc <- processAssays(
  pb_clean,
  formula          = ~ treatment,
  min.cells        = 20,
  min.samples      = 6,
  min.count        = 5,
  normalize.method = "TMM"
)
 
res.vp <- fitVarPart(res.proc, formula = ~ treatment)
p_vp <- plotVarPart(sortCols(res.vp), label.angle = 60)
ggsave(file.path(out_plots, paste0("01_varPart_", RUN_TAG, ".pdf")),
       p_vp, width = 12, height = 6)
 
res.dl <- dreamlet(res.proc, formula = ~ treatment)
 
tt_cvn <- topTable(res.dl, coef = "treatmentnaive", number = Inf) |>
  as_tibble() |>
  mutate(logFC = -logFC, t = -t, contrast = "cocaine_v_naive")
 
# ===== Results =====
deg_summary <- tt_cvn |>
  group_by(assay) |>
  summarise(
    n_tested = n(),
    n_up     = sum(adj.P.Val < 0.05 & logFC > 0),
    n_down   = sum(adj.P.Val < 0.05 & logFC < 0),
    n_total  = n_up + n_down,
    .groups  = "drop"
  ) |>
  arrange(desc(n_total))
print(deg_summary)
 
# Volcanos
tt_vol <- tt_cvn |>
  group_by(assay) |>
  mutate(rank_p = rank(P.Value),
         is_sig = adj.P.Val < 0.05,
         lab    = ifelse(rank_p <= 10 & is_sig, ID, NA_character_)) |>
  ungroup()
 
p_vol <- ggplot(tt_vol, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = is_sig), size = 0.6, alpha = 0.6) +
  geom_text_repel(aes(label = lab), size = 2.5,
                  max.overlaps = Inf, segment.size = 0.2) +
  scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
  facet_wrap(~ assay, scales = "free") +
  labs(title = "Cocaine vs. Naive",
       subtitle = "Positive logFC = up in cocaine") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_plots, paste0("06_volcanos_", RUN_TAG, ".pdf")),
       p_vol, width = 14, height = 10)
 
# ===== Save =====
saveRDS(res.dl,   file.path(out_data, paste0("dreamlet_cvn_", RUN_TAG, ".rds")))
saveRDS(res.proc, file.path(out_data, paste0("processAssays_", RUN_TAG, ".rds")))
saveRDS(pb_clean, file.path(out_data, paste0("pseudobulk_clean_", RUN_TAG, ".rds")))
write.csv(tt_cvn,      file.path(out_data, paste0("DE_cvn_", RUN_TAG, ".csv")), row.names = FALSE)
write.csv(deg_summary, file.path(out_data, paste0("DEG_summary_cvn_", RUN_TAG, ".csv")), row.names = FALSE)
 