

suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(ggspavis)
  library(scCustomize)
  library(patchwork)
  library(ggplot2)
})

output_dir <- here("plots", "Xenium", "05_spatial_clustering")

# load
spe <- readRDS(file = here(processed_dir, "Banksy_integrated_res2.0_collapsed_v6_wRCTD_labels.rds"))



# NOTE: This object was made using the RCTD labels with ~450k cells. This is half the total amount. 
# RCTD likely excluded many when ran with default paramters that I think use UMI_min=100. Need to re run.

# ======== plotting domains =======

num_groups <- length(unique(spe$Banksy_res2.0_collapsed_v6))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )


# Plot with ggspavis
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6_RCTD.pdf"))
pdf(pdf_path, width = 15, height = 15)
ggspavis::plotSpots(spe,
    annotate = 'Banksy_res2.0_collapsed_v6',
    in_tissue = NULL,
    sample_id = "brnum"
  ) +
    scale_color_manual(values = colors) +
    ggtitle(paste0("res=2.0 collapsed")) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2, 
                    byrow = TRUE, 
                    override.aes = list(size=5)
                    )
            )
dev.off()

# ========= Stacked bars of cell cell type proportions per domain ========

# exclude non-neurons
spe.neurons <- spe[, spe$first_type_broad != "Non-neuronal"]

num_groups <- length(unique(spe.neurons$first_type))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )


# get proportion of cell types $first_type, per domain 
celltype_counts <- table(spe.neurons$Banksy_res2.0_collapsed_v6, spe.neurons$first_type)
celltype_props <- prop.table(celltype_counts, margin = 1)
celltype_props_df <- as.data.frame(celltype_props)
colnames(celltype_props_df) <- c("Domain", "CellType", "Proportion")

# Plot stacked bar plot
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6_RCTD_celltype_proportions.pdf"))
pdf(pdf_path, width = 10, height = 6)
ggplot(celltype_props_df, aes(x = Domain, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle("Proportion of Cell Types per Domain")
dev.off()


# ======= Proportions excitatory per domain =======


# get proportion of excitatory cell types per domain
spe.excitatory <- spe.neurons[, spe.neurons$first_type_broad == "Excitatory"]

num_groups <- length(unique(spe.excitatory$first_type))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )


celltype_counts_excitatory <- table(spe.excitatory$Banksy_res2.0_collapsed_v6, spe.excitatory$first_type)
celltype_props_excitatory <- prop.table(celltype_counts_excitatory, margin = 1)
celltype_props_excitatory_df <- as.data.frame(celltype_props_excitatory)
colnames(celltype_props_excitatory_df) <- c("Domain", "CellType", "Proportion")

# Plot stacked bar plot for excitatory cell types
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6_RCTD_excitatory_celltype_proportions.pdf"))
pdf(pdf_path, width = 10, height = 6)
ggplot(celltype_props_excitatory_df, aes(x = Domain, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle("Proportion of Excitatory Cell Types per Domain")
dev.off()


# ======= Proportions inhibitory per domain =======

# get proportion of inhibitory cell types per domain
spe.inhibitory <- spe.neurons[, spe.neurons$first_type_broad == "Inhibitory"]   

num_groups <- length(unique(spe.inhibitory$first_type))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )

celltype_counts_inhibitory <- table(spe.inhibitory$Banksy_res2.0_collapsed_v6, spe.inhibitory$first_type)
celltype_props_inhibitory <- prop.table(celltype_counts_inhibitory, margin = 1)
celltype_props_inhibitory_df <- as.data.frame(celltype_props_inhibitory)
colnames(celltype_props_inhibitory_df) <- c("Domain", "CellType", "Proportion")

# Plot stacked bar plot for inhibitory cell types
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6_RCTD_inhibitory_celltype_proportions.pdf"))
pdf(pdf_path, width = 10, height = 6)
ggplot(celltype_props_inhibitory_df, aes(x = Domain, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle("Proportion of Inhibitory Cell Types per Domain")
dev.off()

# ====== Propotion broad cell types per domain ======

num_groups <- length(unique(spe.neurons$first_type_broad))
colors <- scCustomize::scCustomize_Palette(
    num_groups = num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
  )

# get proportion of cell types $first_type, per domain
celltype_counts <- table(spe.neurons$Banksy_res2.0_collapsed_v6, spe.neurons$first_type_broad)
celltype_props <- prop.table(celltype_counts, margin = 1)
celltype_props_df <- as.data.frame(celltype_props)
colnames(celltype_props_df) <- c("Domain", "CellType", "Proportion")

# Plot stacked bar plot
pdf_path <- file.path(output_dir, paste0("Banksy_integrated_lambda_res2.0_collapsed_v6_RCTD_broad_celltype_proportions.pdf"))
pdf(pdf_path, width = 10, height = 6)
ggplot(celltype_props_df, aes(x = Domain, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle("Proportion of Broad Cell Types per Domain")
dev.off()
