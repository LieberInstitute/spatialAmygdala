suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("scater")
    library("spatialLIBD")
    library("dplyr")
    library('ComplexHeatmap')
    library("patchwork")
    library("ggplot2")
    library("escheR")
})


spe <- readRDS(here("processed-data","Visium", "07_clustering", "BayesSpace", "MarkerGenes", "spe_harmony_markers_BS_k16_Semisupervised_wAI.rds"))
spe 


# ======== Plotting domains =======

pal <- c( 
  AI        = "#D62728",  # strong red standout
  BM        = "#E67E22",  # burnt orange
  BLD      = "#9B59B6",  # violet
  PL     = "#f1e438ff",  # steel blue
  BL     = "#035185ff",  # sky blue
  LA         = "#F4B400",  # goldenrod
  CoA        = "#5DA5DA",  # blue-gray
  CeA        = "#197d43ff",  # emerald green
  MeA        = "#baf739ff",  # chartreuse
  HPC        = "#d6a8f8ff",  # deep royal purple (distinct from BA)
  CHAT       = "#A0522D",  # chestnut brown
  Endothelial   = "#444444",  # charcoal gray
  WM.1       = "#BBBBBB",  # light slate gray
  WM.2       = "#DDDDDD",   # mist gray
  CLA        = "#FF69B4"   # hot pink
)


# ========= Pseudobulk dendrogram =========
library(dreamlet)

pb <- aggregateToPseudoBulk(spe,
  assay = "counts",
  cluster_id = "BS_k16_Semisupervised_wAI",
  sample_id = "sample_id",
  verbose = FALSE
)

# Hierarchical clustering of cell types
hcl <- buildClusterTreeFromPB(pb, method="ward.D")

pdf(here("plots", "Visium", "08_marker_genes", "Dendrogram_domains_pseudobulk.pdf"),
    width = 6, height = 4)
plot(hcl, hang=-1)
dev.off()


# ========= Plotting with Image =========

pdf(here("plots", "Visium", "08_marker_genes", "Domains_with_image_all_samples.pdf"),
    width = 8, height = 7)
spatialLIBD::vis_clus(
    spe,
    sampleid = unique(spe$sample_id)[5],
    clustervar = 'BS_k16_Semisupervised_wAI',
    is_stitched = TRUE,
    colors = pal
) + guides(
      fill = guide_legend(
        override.aes = list(size = 6),
        ncol = 1,
        byrow = TRUE
      )
    ) +
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 10)),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.key = element_blank(),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "lines"),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_blank(),
      plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
    )
dev.off()


# Desired order
sample_order <- c("Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660")

# Reorder factor levels
spe$sample_id <- factor(spe$sample_id, levels = sample_order)

spe$sample_id <- factor(as.character(spe$sample_id),
                        levels = sample_order,
                        ordered = TRUE)

# Get all unique sample IDs
sample_ids <- unique(spe$sample_id)
sample_ids

ordered_ids <- levels(droplevels(spe$sample_id))


# now plot with image but loot through all samples
plot_list <- list()
pdf(here("plots", "Visium", "08_marker_genes", "Domains_with_image_loop_all_samples.pdf"),
     width = 10 * 3, height = 10)
for (sid in ordered_ids) {
  p <- spatialLIBD::vis_clus(
    spe,
    sampleid = sid,
    clustervar = 'BS_k16_Semisupervised_wAI',
    is_stitched = TRUE,
    colors = pal
  ) + guides(
      fill = guide_legend(
        override.aes = list(size = 6),
        ncol = 1,
        byrow = TRUE
      )
    ) +
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 10)),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.key = element_blank(),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "lines"),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_blank(),
      plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
    )
  plot_list[[sid]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 4) &
  theme(plot.margin = margin(1,1,1,1))

print(combined_plot)

dev.off()


# plot only 8325
pdf(here("plots", "Visium", "08_marker_genes", "Domains_with_image_Br8325.pdf"),
    width = 8, height = 7)
spatialLIBD::vis_clus(
    spe,
    sampleid = "Br8325",
    clustervar = 'BS_k16_Semisupervised_wAI',
    is_stitched = TRUE,
    colors = pal
) + guides(
      fill = guide_legend(
        override.aes = list(size = 6),
        ncol = 1,
        byrow = TRUE
      )
    ) +
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 10)),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.key = element_blank(),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "lines"),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_blank(),
      plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
    )
dev.off()



# ====== gross Anterior, Middle, Posterior grouping and plots =======

library(dplyr)

spe$AP_axis <- dplyr::case_when(
  spe$sample_id %in% c("Br9192") ~ "Anterior-intermediate",
  spe$sample_id %in% c("Br9280", "Br2743", "Br6471") ~ "Intermediate",
  spe$sample_id %in% c("Br8325", "Br6423") ~ "Posterior-intermediate",
  spe$sample_id %in% c("Br6660") ~ "Posterior",
  TRUE ~ NA_character_
)

spe$AP_axis <- factor(
  spe$AP_axis,
  levels = c("Anterior-intermediate", "Intermediate", "Posterior-intermediate", "Posterior")
)

# get domain proportions per AP group
domain_props_sample <- colData(spe) |>
  as.data.frame() |>
  group_by(sample_id, AP_axis, BS_k16_Semisupervised_wAI) |>
  summarise(counts = n(), .groups = "drop_last") |>
  mutate(proportion = counts / sum(counts)) |>
  ungroup()


domain_props_mean <- domain_props_sample |>
  group_by(AP_axis, BS_k16_Semisupervised_wAI) |>
  summarise(mean_prop = mean(proportion), .groups = "drop")

library(ggplot2)

pdf(here("plots", "Visium", "08_marker_genes", "Domain_composition_AP_axis_barplot.pdf"),
    width = 6, height = 4)
ggplot(domain_props_mean,
       aes(x = AP_axis, y = mean_prop,
           fill = BS_k16_Semisupervised_wAI)) +
  geom_col() +
  theme_bw() +
  ylab("Mean proportion of spots") +
  xlab("AP axis") +
  scale_y_continuous(labels = scales::percent_format()) +
  ggtitle("Average domain composition across AP axis") +
  scale_fill_manual(values = pal, name = "Spatial Domains")
dev.off()


itc_prop <- domain_props_sample |>
  filter(BS_k16_Semisupervised_wAI == "AI")

pdf(here("plots", "Visium", "08_marker_genes", "ITC_proportion_AP_axis_spot_jitter_boxplot.pdf"),
    width = 3.5, height = 3.5)
ggplot(itc_prop, aes(x = AP_axis, y = proportion)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.7, fill = "red") +
  geom_jitter(width = 0.12, height = 0, size = 2.5, alpha = 0.8) +  # each point = one sample
  labs(x = "AP axis", y = "ITC proportion (per sample)") +
  theme_bw() 
dev.off()

# # spotplot of example samples from each AP group
# example_samples <- c("Br9192", "Br8325", "Br6660")
# spe_example <- spe[, spe$sample_id %in% example_samples]

# # loop through and make plots
# plot_list <- list()
# for (sid in example_samples) {
#   spe.subset <- spe_example[, spe_example$sample_id == sid]

#   p <- make_escheR(spe.subset) |>
#     add_fill(var = "BS_k16_Semisupervised_wAI", point_size = 1.25) +
#     scale_fill_manual(values = pal, name = "Spatial Domains") +
#     ggtitle(sid) +
#     guides(
#       fill = guide_legend(
#         override.aes = list(size = 6),
#         ncol = 1,
#         byrow = TRUE
#       )
#     ) +
#     theme(
#       plot.title = element_text(size = 22, margin = margin(b = 10)),
#       panel.border = element_rect(colour = "black", fill = NA, size = 1),
#       legend.title = element_text(size = 18),
#       legend.text = element_text(size = 15),
#       legend.key = element_blank(),
#       legend.key.size = unit(1.2, "lines"),
#       legend.spacing.y = unit(0.4, "lines"),
#       legend.position = "right",
#       legend.justification = "top",
#       legend.background = element_blank(),
#       plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
#     )

#   plot_list[[sid]] <- p
# }

# # Combine all plots side by side
# combined_plot <- wrap_plots(plot_list, ncol = 3)

# # Save to PDF
# pdf(here("plots", "Visium", "08_marker_genes", "Domains_example_samples_AP_axis.pdf"),
#     width = 8 * 2.5, height = 7)
# print(combined_plot)
# dev.off()





# ====== plotting clusters ========


# Desired order
sample_order <- c("Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660")

# Reorder factor levels
spe$sample_id <- factor(spe$sample_id, levels = sample_order)

spe$sample_id <- factor(as.character(spe$sample_id),
                        levels = sample_order,
                        ordered = TRUE)


# Get all unique sample IDs
sample_ids <- unique(spe$sample_id)
sample_ids

ordered_ids <- levels(droplevels(spe$sample_id))

# Create a list to hold plots
plot_list <- list()

# Loop through each sample_id
for (sid in ordered_ids) {
  spe.subset <- spe[, spe$sample_id == sid]

  p <- make_escheR(spe.subset) |>
    add_fill(var = "BS_k16_Semisupervised_wAI", point_size = 1) +
    scale_fill_manual(values = pal, name = "Spatial Domains") +
    ggtitle(sid) +
    guides(
      fill = guide_legend(
        override.aes = list(size = 6),
        ncol = 1,
        byrow = TRUE
      )
    ) +
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 10)),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.key = element_blank(),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.y = unit(0.4, "lines"),
      legend.position = "right",
      legend.justification = "top",
      legend.background = element_blank(),
      plot.margin = margin(t = 10, r = 70, b = 10, l = 40)
    )

  plot_list[[sid]] <- p
}

# Combine all plots side by side
combined_plot <- wrap_plots(plot_list, ncol = 4)

# Save to PDF
pdf(here("plots", "Visium", "08_marker_genes", "Clusters_all_samples.pdf"),
    width = 8 * 4, height = 7*2)
print(combined_plot)
dev.off()







# -----------------------------
# 1) Desired order
# -----------------------------
sample_order <- c("Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660")

spe$sample_id <- factor(as.character(spe$sample_id),
                        levels = sample_order,
                        ordered = TRUE)

ordered_ids <- levels(droplevels(spe$sample_id))

# -----------------------------
# 2) Compute per-sample x/y ranges
#    (used to keep TRUE spatial scale consistent)
# -----------------------------
get_xy <- function(spe_obj) {
  xy <- as.data.frame(spatialCoords(spe_obj))
  xy <- xy[, 1:2, drop = FALSE]
  colnames(xy) <- c("x", "y")
  xy
}

range_df <- lapply(ordered_ids, function(sid) {
  spe_sub <- spe[, spe$sample_id == sid]
  xy <- get_xy(spe_sub)

  data.frame(
    sample_id = sid,
    x_rng = diff(range(xy$x, na.rm = TRUE)),
    y_rng = diff(range(xy$y, na.rm = TRUE))
  )
}) |>
  bind_rows()

# guard against any weird zero ranges
range_df$x_rng <- pmax(range_df$x_rng, 1)
range_df$y_rng <- pmax(range_df$y_rng, 1)

# -----------------------------
# 3) Make per-sample plots (NO legends; we collect once)
# -----------------------------
plot_list <- list()

for (sid in ordered_ids) {
  spe.subset <- spe[, spe$sample_id == sid]

  p <- make_escheR(spe.subset) |>
    add_fill(var = "BS_k16_Semisupervised_wAI", point_size = 1) +
    scale_fill_manual(values = pal, name = "Spatial Domains") +
    ggtitle(sid) +
    coord_equal(expand = FALSE) +   # IMPORTANT: preserves true x/y scale in each panel
    theme(
      plot.title = element_text(size = 22, margin = margin(b = 8)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      legend.position = "none",     # IMPORTANT: collect later
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

  plot_list[[sid]] <- p
}

# -----------------------------
# 4) Assemble with TRUE scale consistency across rows
#    Key idea:
#    - widths within a row proportional to x-range
#    - add a spacer to row2 so total row widths match row1
#    - set row heights proportional to (max y-range)/(total x-range of that row)
#      because coord_equal ties height to width
# -----------------------------
row1_ids <- sample_order[1:4]
row2_ids <- sample_order[5:7]

w1 <- range_df$x_rng[match(row1_ids, range_df$sample_id)]
w2 <- range_df$x_rng[match(row2_ids, range_df$sample_id)]

y1 <- range_df$y_rng[match(row1_ids, range_df$sample_id)]
y2 <- range_df$y_rng[match(row2_ids, range_df$sample_id)]

# Spacer width so row2 doesn't get "blown up" to fill the row width
sp_w <- sum(w1) - sum(w2)
sp_w <- ifelse(sp_w > 0, sp_w, 1)

row1 <- plot_list[[row1_ids[1]]] + plot_list[[row1_ids[2]]] +
        plot_list[[row1_ids[3]]] + plot_list[[row1_ids[4]]] +
        plot_layout(widths = w1)

row2 <- plot_list[[row2_ids[1]]] + plot_list[[row2_ids[2]]] +
        plot_list[[row2_ids[3]]] + plot_spacer() +
        plot_layout(widths = c(w2, sp_w))

# Height ratios that preserve 1:1 scale ACROSS rows
h_rel1 <- max(y1) / sum(w1)
h_rel2 <- max(y2) / (sum(w2) + sp_w)

combined_plot <-
  (row1 / row2) +
  plot_layout(heights = c(h_rel1, h_rel2), guides = "collect") &
  theme(
    legend.position = "right",
    legend.justification = "top",
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 15),
    legend.key = element_blank(),
    legend.key.size = unit(1.2, "lines"),
    legend.spacing.y = unit(0.4, "lines")
  ) &
  guides(
    fill = guide_legend(override.aes = list(size = 6), ncol = 1, byrow = TRUE)
  )

# -----------------------------
# 5) Save
# -----------------------------
pdf(here("plots", "Visium", "08_marker_genes", "Clusters_all_samples.pdf"),
    width = 32, height = 14)
print(combined_plot)
dev.off()




















# ======= Plotting immature marker genes ======

markers <- c("NR2F2", "SP8", "PROX1", "MKI67", "TBR1", "DCX", "NCAM1", "BCL2", "SOX11", "ROBO1")

pdf(here("plots", "Visium", "08_marker_genes", "Immature_markers_Heatmap.pdf"), width = 8, height = 6)
plotGroupedHeatmap(spe, features=markers,
    group="BS_k16_Semisupervised_wAI")
dev.off()


pdf(here("plots", "Visium", "08_marker_genes", "Immature_markers_Violins.pdf"), width = 8, height = 6)
plotExpression(
  spe,
  features = markers,
  x = "BS_k16_Semisupervised_wAI",
  exprs_values = "logcounts"
)
dev.off()