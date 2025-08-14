library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")

# save directories
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
plot_dir <- here("plots", "Xenium", "06_label_transfer")


# ===== Load predictions ======
xenium.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype.csv"))
xenium.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype.csv"))

proseg.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "proseg_pred_broad_celltype.csv"))
proseg.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "proseg_pred_fine_celltype.csv"))

xenium.broad_n100 <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype_n100.csv"))
xenium.fine_n100 <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype_n100.csv"))

# ======= Scores ggplot ======

library(tidyverse)

# Helper function: extract matching scores and label per dataset
prepare_score_df <- function(df, source_label) {
  df %>%
    as_tibble() %>%
    mutate(cell_id = row_number()) %>%
    pivot_longer(cols = starts_with("scores."), names_to = "celltype", values_to = "score") %>%
    mutate(celltype_clean = str_remove(celltype, "scores\\."),
           celltype_clean = str_replace_all(celltype_clean, "\\.", "-"),
           label_clean = str_replace_all(labels, "\\.", "-")) %>%
    filter(celltype_clean == label_clean) %>%
    mutate(source = source_label)
}


# ----- Broad cell types -----

# Prepare both datasets
df_xenium <- prepare_score_df(xenium.broad, "Xenium")
df_proseg <- prepare_score_df(proseg.broad, "Proseg")
df_xenium_n100 <- prepare_score_df(xenium.broad_n100, "Xenium (n=100)")

# Combine
df_combined <- bind_rows(df_xenium, df_proseg, df_xenium_n100)

# Plot
pdf(file = here(plot_dir, "xenium_vs_proseg_label_transfer_scores_broad.pdf"), width = 10, height = 6)
ggplot(df_combined, aes(x = label_clean, y = score, fill = source)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label", fill = "Source")
dev.off()


# ----- Fine cell types -----

# Prepare both datasets
df_xenium_fine <- prepare_score_df(xenium.fine, "Xenium")
df_proseg_fine <- prepare_score_df(proseg.fine, "Proseg")

# Combine
df_combined_fine <- bind_rows(df_xenium_fine, df_proseg_fine)

# Plot
pdf(file = here(plot_dir, "xenium_vs_proseg_label_transfer_scores_fine.pdf"), width = 10, height = 6)
ggplot(df_combined_fine, aes(x = label_clean, y = score, fill = source)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label", fill = "Source")
dev.off()




# ====== Heatmap of all predictions ======

library(pheatmap)


# 1. Prepare long-format data
df_long <- xenium.fine %>%
  as_tibble() %>%
  mutate(cell_id = row_number()) %>%
  pivot_longer(cols = starts_with("scores."), names_to = "predicted_celltype", values_to = "score") %>%
  mutate(predicted_celltype = str_remove(predicted_celltype, "scores\\."),
         predicted_celltype = str_replace_all(predicted_celltype, "\\.", "-"),
         true_label = str_replace_all(xenium.fine$labels[cell_id], "\\.", "-"))

# 2. Compute mean or median score for each (true_label x predicted_celltype)
confusion_mat <- df_long %>%
  group_by(true_label, predicted_celltype) %>%
  summarize(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = predicted_celltype, values_from = median_score) %>%
  column_to_rownames("true_label")

# 3. Convert to matrix
confusion_mat_matrix <- as.matrix(confusion_mat)

# Plot heatmap
pdf(file = here(plot_dir, "xenium_fine_off_diagonal_median_scores_heatmap.pdf"), width = 10, height = 8)
pheatmap(confusion_mat_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = viridis::viridis(100),
         main = "Off-diagonal Median Scores (Xenium Fine)",
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8)
dev.off()
