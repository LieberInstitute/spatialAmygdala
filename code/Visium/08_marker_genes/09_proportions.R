suppressPackageStartupMessages({
  library("here")
  library("sessioninfo")
  library("SpatialExperiment")
  library("scater")
  library("spatialLIBD")
  library("dplyr")
  library("ComplexHeatmap")
  library("patchwork")
  library("RColorBrewer")
  library("SummarizedExperiment")
  library("circlize")
  library("grid")
  library("dreamlet")
})

spe <- readRDS(here("processed-data","Visium","07_clustering","BayesSpace","MarkerGenes","spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe

# reorder donors: Br9192, 9280, 2743, 6471, 8325, 6423, 6660
spe$sample_id <- factor(spe$sample_id, levels = c("Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660"))

# ========= Domain Proportions across Samples Barplot =========

# get counts per domain per sample
domain_counts <- table(spe$sample_id, spe$BS_k16_Semisupervised_wAI)
domain_counts_df <- as.data.frame(domain_counts)
colnames(domain_counts_df) <- c("Sample", "Domain", "Count")

# calculate proportions
domain_props_df <- domain_counts_df %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count))  
domain_props_df

# plot
library(ggplot2)
p <- ggplot(domain_props_df, aes(x = Sample, y = Proportion, fill = Domain)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Domain Proportions Across Samples",
       x = "Sample ID",
       y = "Proportion of Spots") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = pal)   

ggsave(
  filename = here("plots","Visium", "08_marker_genes", "Domain_Proportions_Across_Samples.pdf"),
  plot = p,
  width = 5,
  height = 6,
  useDingbats = FALSE
)