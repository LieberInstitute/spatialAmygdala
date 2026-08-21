#  ---------------------------------------------------------------------------
#  Maximum-power test: pseudobulk at sample level (collapse all cell types)
#  One pseudobulk profile per rat. Tests cocaine vs naive at bulk resolution.
#  If there's any signal in the data, this finds it.
#  ---------------------------------------------------------------------------

library(here)
library(Seurat)
library(SingleCellExperiment)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)

set.seed(12345)

out_plots <- here("plots","snRNAseq","dreamlet")
out_data  <- here("processed-data","snRNAseq","dreamlet")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_data,  recursive = TRUE, showWarnings = FALSE)

# ===== Load =====
rat.amy <- readRDS(here("processed-data","snRNAseq","zhou_with_yu_labels.rds"))
sce <- as.SingleCellExperiment(rat.amy, assay = "RNA")
assays(sce) <- list(counts = counts(sce))

# Drop outlier samples (identified from earlier PCA)
drop_samples <- c("Rat_Opioid_HS_1", "Rat_Opioid_HS_2", "Rat_Amygdala_787A_all_seq")
sce <- sce[, !sce$sample %in% drop_samples]
message("Cells retained: ", ncol(sce))

# ===== Aggregate counts per sample =====
samples <- unique(as.character(sce$sample))
bulk_counts <- sapply(samples, function(s) {
  rowSums(counts(sce)[, sce$sample == s, drop = FALSE])
})
colnames(bulk_counts) <- samples
message("Bulk matrix: ", nrow(bulk_counts), " genes x ", ncol(bulk_counts), " samples")

# Sample metadata
meta <- as.data.frame(colData(sce)) |>
  distinct(sample, treatment, batch, label) |>
  filter(sample %in% samples)
meta <- meta[match(samples, meta$sample), ]
meta$treatment <- factor(meta$treatment, levels = c("naive", "cocaine"))
meta$batch     <- factor(meta$batch)
print(meta)
print(table(meta$treatment, meta$batch))

# ===== edgeR/limma-voom DE =====
dge <- DGEList(counts = bulk_counts)

# Filter low-expression genes — standard edgeR filter
keep <- filterByExpr(dge, group = meta$treatment, min.count = 10, min.total.count = 15)
dge <- dge[keep, , keep.lib.sizes = FALSE]
message("Genes after filtering: ", nrow(dge))

dge <- calcNormFactors(dge, method = "TMM")

# Model 1: with batch
design <- model.matrix(~ treatment + batch, data = meta)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
tt_batch <- topTable(fit, coef = "treatmentcocaine", number = Inf) |>
  rownames_to_column("gene") |>
  as_tibble()

# Model 2: without batch (sanity check)
design_nb <- model.matrix(~ treatment, data = meta)
v_nb <- voom(dge, design_nb, plot = FALSE)
fit_nb <- lmFit(v_nb, design_nb)
fit_nb <- eBayes(fit_nb)
tt_nobatch <- topTable(fit_nb, coef = "treatmentcocaine", number = Inf) |>
  rownames_to_column("gene") |>
  as_tibble()

# ===== Summary =====
summarise_degs <- function(tt, name) {
  tt |>
    summarise(
      model       = name,
      n_tested    = n(),
      n_up_05     = sum(adj.P.Val < 0.05 & logFC > 0),
      n_down_05   = sum(adj.P.Val < 0.05 & logFC < 0),
      n_up_10     = sum(adj.P.Val < 0.10 & logFC > 0),
      n_down_10   = sum(adj.P.Val < 0.10 & logFC < 0),
      n_nominal   = sum(P.Value < 0.01)
    )
}

deg_summary <- bind_rows(
  summarise_degs(tt_batch,   "with_batch"),
  summarise_degs(tt_nobatch, "no_batch")
)
print(deg_summary)

# ===== P-value distribution =====
phist_df <- bind_rows(
  tt_batch   |> mutate(model = "with_batch"),
  tt_nobatch |> mutate(model = "no_batch")
)
p_phist <- ggplot(phist_df, aes(P.Value)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white") +
  facet_wrap(~ model) +
  theme_bw() +
  labs(title = "P-value distribution — sample-level pseudobulk DE",
       subtitle = "Peaked near 0 = real signal; flat = none")
ggsave(file.path(out_plots, "10_bulk_pvalue_hist.pdf"),
       p_phist, width = 10, height = 5)

# ===== Volcanos =====
make_volcano <- function(tt, title) {
  tt <- tt |>
    mutate(is_sig = adj.P.Val < 0.05,
           lab = ifelse(rank(P.Value) <= 20, gene, NA_character_))
  ggplot(tt, aes(logFC, -log10(P.Value))) +
    geom_point(aes(color = is_sig), size = 1, alpha = 0.6) +
    geom_text_repel(aes(label = lab), size = 3, max.overlaps = Inf) +
    scale_color_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey60")) +
    labs(title = title,
         subtitle = "Top 20 by raw p-value labeled") +
    theme_bw() + theme(legend.position = "none")
}

pdf(file.path(out_plots, "11_bulk_volcanos.pdf"), width = 10, height = 8)
print(make_volcano(tt_batch,   "Sample-level pseudobulk — with batch"))
print(make_volcano(tt_nobatch, "Sample-level pseudobulk — no batch"))
dev.off()

# ===== Save =====
write.csv(tt_batch,   file.path(out_data, "DE_bulk_with_batch.csv"),  row.names = FALSE)
write.csv(tt_nobatch, file.path(out_data, "DE_bulk_no_batch.csv"),    row.names = FALSE)
write.csv(deg_summary, file.path(out_data, "DEG_summary_bulk.csv"),   row.names = FALSE)

message("Done. Top genes (with batch):")
print(head(tt_batch, 20))