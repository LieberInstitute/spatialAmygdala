suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(spatialLIBD)
  library(patchwork)
  library(ggplot2)
  library(scater)
  library(pheatmap)
})

output_dir <- here("plots", "Xenium", "08_gex_correlations")
processed_dir_out <- here("processed-data", "Xenium", "08_gex_correlations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir_out, recursive = TRUE, showWarnings = FALSE)

# ===== Load data =====
# load xenium data
spe.old <- readRDS(here("processed-data","Xenium","04_dim_reduction", "spe_xenium_5um_harmonized_singlecell.rds"))
spe.new <- readRDS(file.path("processed-data", "Xenium", "05_spatial_clustering", "Banksy_domains_v1.0.rds"))

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
# [33] "Banksy_domains" 





# ===== Stacked bar plot: subset to amygdala domains =====

# Define domains to keep
amygdala_domains <- c("BL", "LA", "PL", "BM_CoA", "MeA_AI", "CeA")

# ----- Build data frame -----
df_bar <- data.frame(
  domain = spe$Banksy_domains,
  celltype = spe$pred_broad_celltype
) %>%
  filter(!is.na(celltype), !is.na(domain), domain %in% amygdala_domains)

# Set domain order
df_bar$domain <- factor(df_bar$domain, levels = amygdala_domains)

# Compute proportions per domain
df_bar_summary <- df_bar %>%
  group_by(domain, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(domain) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()

# ----- Plot -----
p_stacked <- ggplot(df_bar_summary, aes(x = prop, y = domain, fill = celltype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Excitatory" = "#E76F51",
                               "Inhibitory" = "#264653",
                               "Non-neuronal" = "#A9A9A9"),
                    name = "Cell type") +
  labs(x = "Percent (%)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "domain_celltype_stacked_bar_amygdala.pdf"),
       p_stacked, width = 4, height = 5)





# ===== Dot plots: size = cell type proportion, color = raw mean expression =====
# Neuronal only, one PDF per donor, grey-to-blue color scale

# ----- Define marker-domain pairs -----
marker_domain_pairs <- data.frame(
  gene = c("PDYN", "COL25A1", "CYP26B1", "PENK"),
  domain = c("PL", "PL", "LA", "CeA"),
  stringsAsFactors = FALSE
)

# ----- Extract expression and metadata -----
genes <- marker_domain_pairs$gene
expr_mat <- as.matrix(assay(spe, "nucleus_normcounts")[genes, , drop = FALSE])

df <- data.frame(
  cell_id = colnames(spe),
  domain = spe$Banksy_domains,
  celltype = spe$pred_broad_celltype,
  donor = spe$brnum,
  t(expr_mat),
  check.names = FALSE
)

# Remove NAs and non-neuronal
df <- df %>%
  filter(!is.na(celltype), !is.na(domain), celltype != "Non-neuronal")

# Pivot to long
df_long <- df %>%
  pivot_longer(cols = all_of(genes),
               names_to = "gene",
               values_to = "logcounts")

# Filter to relevant domain per gene
df_plot <- df_long %>%
  inner_join(marker_domain_pairs, by = c("gene", "domain"))

# Add facet label
df_plot$facet_label <- paste0(df_plot$gene, " (", df_plot$domain, ")")

# ----- Helper function -----
make_donor_dotplot <- function(df_plot, donor_id) {
  
  df_donor <- df_plot %>%
    filter(donor == donor_id)
  
  # Cell type proportions per domain (neuronal only)
  domain_props <- df_donor %>%
    distinct(cell_id, .keep_all = TRUE) %>%
    group_by(domain, celltype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(domain) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup()
  
  # Mean expression per gene x domain x celltype
  df_summ <- df_donor %>%
    group_by(gene, domain, celltype, facet_label) %>%
    summarise(
      mean_expr = mean(logcounts, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Merge proportions
  df_summ <- df_summ %>%
    left_join(domain_props %>% select(domain, celltype, prop),
              by = c("domain", "celltype"))
  
  # Plot
  ggplot(df_summ, aes(x = celltype, y = facet_label)) +
    geom_point(aes(size = prop, color = mean_expr)) +
    scale_color_gradient(low = "grey85", high = "red", name = "Mean\nexpr") +
    scale_size_continuous(range = c(1, 10), name = "% of neurons\nin domain",
                          breaks = c(10, 25, 50, 75)) +
    labs(x = NULL, y = NULL, title = donor_id) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold")
    )
}

# ----- Generate and save -----
p_br9280 <- make_donor_dotplot(df_plot, "Br9280")
p_br9206 <- make_donor_dotplot(df_plot, "Br9206")

ggsave(file.path(output_dir, "marker_gene_dotplot_Br9280_neuronal.pdf"),
       p_br9280, width = 5, height = 4)

ggsave(file.path(output_dir, "marker_gene_dotplot_Br9206_neuronal.pdf"),
       p_br9206, width = 5, height = 4)

