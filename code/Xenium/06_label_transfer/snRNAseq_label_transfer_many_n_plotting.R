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


# Load broad labels
p3  <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n3.csv"))
p5  <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n5.csv"))
p10 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n10.csv"))
p15 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n15.csv"))
p25 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n25.csv"))
p35 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n35.csv"))
p50 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n50.csv"))
p100 <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_broad_celltype_n100.csv"))

# load fine labels
p3_fine  <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n3.csv"))
p5_fine  <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n5.csv"))
p10_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n10.csv"))
p15_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n15.csv"))
p25_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n25.csv"))
p35_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n35.csv"))
p50_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n50.csv"))
p100_fine <- read.csv(here("processed-data/Xenium/06_label_transfer/xenium_5um_pred_fine_celltype_n100.csv"))


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
# Prepare individually
df3   <- prepare_score_df(p3,   "n=3")
df5   <- prepare_score_df(p5,   "n=5")
df10  <- prepare_score_df(p10,  "n=10")
df15  <- prepare_score_df(p15,  "n=15")
df25  <- prepare_score_df(p25,  "n=25")
df35  <- prepare_score_df(p35,  "n=35")
df50  <- prepare_score_df(p50,  "n=50")
df100 <- prepare_score_df(p100, "n=100")

# Combine all
df_combined <- bind_rows(df3, df5, df10, df15, df25, df35, df50, df100)

# Order by de.n
df_combined$source <- factor(df_combined$source, levels = paste0("n=", c(3, 5, 10, 15, 25, 35, 50, 100)))

# Plot
pdf(file = here(plot_dir, "xenium_label_scores_many_n_broad.pdf"), width = 20, height = 6)
ggplot(df_combined, aes(x = label_clean, y = score, fill = source)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label", fill = "Source")
dev.off()


# ---- Fine cell types -----
# Prepare individually
df3_fine   <- prepare_score_df(p3_fine,   "n=3")
df5_fine   <- prepare_score_df(p5_fine,   "n=5")
df10_fine  <- prepare_score_df(p10_fine,  "n=10")
df15_fine  <- prepare_score_df(p15_fine,  "n=15")
df25_fine  <- prepare_score_df(p25_fine,  "n=25")
df35_fine  <- prepare_score_df(p35_fine,  "n=35")
df50_fine  <- prepare_score_df(p50_fine,  "n=50")
df100_fine <- prepare_score_df(p100_fine, "n=100")

# Combine all
df_combined_fine <- bind_rows(df3_fine, df5_fine, df10_fine, df15_fine, df25_fine, df35_fine, df50_fine, df100_fine)

# Order by de.n
df_combined_fine$source <- factor(df_combined_fine$source, levels = paste0("n=", c(3, 5, 10, 15, 25, 35, 50, 100)))

# Plot
pdf(file = here(plot_dir, "xenium_label_scores_many_n_fine.pdf"), width = 20, height = 6)
ggplot(df_combined_fine, aes(x = label_clean, y = score, fill = source)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label", fill = "Source")
dev.off()





# Compute cell-type means for each source
celltype_means <- df_combined_fine %>%
  group_by(source, label_clean) %>%
  summarize(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

# Plot
pdf(file = here(plot_dir, "xenium_label_scores_celltype_mean_vs_n_fine.pdf"), width = 8, height = 5)
ggplot(celltype_means, aes(x = source, y = mean_score)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightgray") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  theme_bw() +
  labs(x = "Number of Marker Genes (de.n)", y = "Mean Score per Cell Type",
       title = "Mean Assigned Score per Fine Cell Type (by de.n)")
dev.off()
