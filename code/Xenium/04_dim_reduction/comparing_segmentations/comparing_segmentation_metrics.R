library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("harmony")

# Save directories
processed_dir <- here("processed-data", "Xenium", "04_dim_reduction")
plot_dir <- here("plots", "Xenium", "04_dim_reduction")

# Load the RData file
ps.5um <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_proseg_5um_harmonized_singlecell.rds"))
ps.5um

ps.15um <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_proseg_15um_harmonized_singlecell.rds"))
ps.15um

xs.5um <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_xenium_5um_harmonized_singlecell.rds"))
xs.5um

xs.15um <- readRDS(here("processed-data/Xenium/04_dim_reduction/spe_xenium_15um_harmonized_singlecell.rds"))
xs.15um



# ==== Add QC metrics to each SPE ====
# NOTE: QC metrics already present from previous steps 

ps.5um$Segmentation <- "Proseg_5um"
ps.15um$Segmentation <- "Proseg_15um"
xs.5um$Segmentation <- "Xenium_5um"
xs.15um$Segmentation <- "Xenium_15um"



# ==== Combine colData into one DataFrame ====
qc_df <- rbind(
  as.data.frame(colData(ps.5um)[, c("Segmentation", "detected", "sum")]),
  as.data.frame(colData(ps.15um)[, c("Segmentation", "detected", "sum")]),
  as.data.frame(colData(xs.5um)[, c("Segmentation", "detected", "sum")]),
  as.data.frame(colData(xs.15um)[, c("Segmentation", "detected", "sum")])
)
colnames(qc_df) <- c("Segmentation", "DetectedGenes", "LibrarySize")


# ==== Summarize total counts and total cells ====
summary_df <- qc_df %>%
  group_by(Segmentation) %>%
  summarise(
    TotalTranscripts = sum(LibrarySize),
    TotalCells = n(),
    .groups = "drop"
  )


# ==== Plot 1: Mean number of detected genes per cell ====
p1 <- ggplot(qc_df, aes(x = Segmentation, y = DetectedGenes)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "red") +
  labs(y = "Detected Genes per Cell", x = NULL, title = "Mean Detected Genes per Cell") +
  theme_bw()

# ==== Plot 2: Mean library size per cell ====
p2 <- ggplot(qc_df, aes(x = Segmentation, y = LibrarySize)) +
  geom_violin(trim = FALSE, fill = "lightgreen") +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "red") +
  labs(y = "Library Size per Cell", x = NULL, title = "Mean Library Size per Cell") +
  theme_bw()

# ==== Plot 3: Total number of transcripts ====
p3 <- ggplot(summary_df, aes(x = Segmentation, y = TotalTranscripts, fill = Segmentation)) +
  geom_bar(stat = "identity") +
  labs(y = "Total Transcripts", x = NULL, title = "Total Transcripts per Dataset") +
  theme_bw() + theme(legend.position = "none")

# ==== Plot 4: Total number of detected cells ====
p4 <- ggplot(summary_df, aes(x = Segmentation, y = TotalCells, fill = Segmentation)) +
  geom_bar(stat = "identity") +
  labs(y = "Total Cells", x = NULL, title = "Total Cells per Dataset") +
  theme_bw() + theme(legend.position = "none")


# ==== Combine plots ====
png(file.path(plot_dir, "comparing_segmentation_metrics.png"), width = 10, height = 10, units = "in", res = 300)
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()





# =============== Visualizationg metrics per sample ================


# ==== Combine colData into one DataFrame with brnum ====
qc_df <- rbind(
  as.data.frame(colData(ps.5um)[, c("Segmentation", "brnum", "detected", "sum")]),
  as.data.frame(colData(ps.15um)[, c("Segmentation", "brnum", "detected", "sum")]),
  as.data.frame(colData(xs.5um)[, c("Segmentation", "brnum", "detected", "sum")]),
  as.data.frame(colData(xs.15um)[, c("Segmentation", "brnum", "detected", "sum")])
)
colnames(qc_df) <- c("Segmentation", "Sample", "DetectedGenes", "LibrarySize")

# ==== Summarize per-sample metrics ====
sample_summary <- qc_df %>%
  group_by(Segmentation, Sample) %>%
  summarise(
    MeanDetectedGenes = mean(DetectedGenes),
    MeanLibrarySize = mean(LibrarySize),
    TotalTranscripts = sum(LibrarySize),
    TotalCells = n(),
    .groups = "drop"
  )


# Map segmentation to box colors
box_colors <- c("Proseg_5um" = "grey70", "Proseg_15um" = "grey70",
                "Xenium_5um" = "grey30", "Xenium_15um" = "grey30")

# ==== Plot 1: Mean number of detected genes per cell ====
p1 <- ggplot(sample_summary, aes(x = Segmentation, y = MeanDetectedGenes)) +
  geom_boxplot(aes(fill = Segmentation), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = box_colors, guide = "none") +
  geom_jitter(aes(color = Sample, shape = Sample), width = 0.2, size = 3) +
  labs(y = "Mean Detected Genes per Cell", x = NULL, title = "Detected Genes per Cell") +
  theme_bw()

# ==== Plot 2: Mean library size per cell ====
p2 <- ggplot(sample_summary, aes(x = Segmentation, y = MeanLibrarySize)) +
  geom_boxplot(aes(fill = Segmentation), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = box_colors, guide = "none") +
  geom_jitter(aes(color = Sample, shape = Sample), width = 0.2, size = 3) +
  labs(y = "Mean Library Size per Cell", x = NULL, title = "Library Size per Cell") +
  theme_bw()

# ==== Plot 3: Total number of transcripts ====
p3 <- ggplot(sample_summary, aes(x = Segmentation, y = TotalTranscripts)) +
  geom_boxplot(aes(fill = Segmentation), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = box_colors, guide = "none") +
  geom_jitter(aes(color = Sample, shape = Sample), width = 0.2, size = 3) +
  labs(y = "Total Transcripts", x = NULL, title = "Total Transcripts per Dataset") +
  theme_bw()

# ==== Plot 4: Total number of detected cells ====
p4 <- ggplot(sample_summary, aes(x = Segmentation, y = TotalCells)) +
  geom_boxplot(aes(fill = Segmentation), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = box_colors, guide = "none") +
  geom_jitter(aes(color = Sample, shape = Sample), width = 0.2, size = 3) +
  labs(y = "Total Cells", x = NULL, title = "Total Cells per Dataset") +
  theme_bw()

# ==== Combine plots ====
png(file.path(plot_dir, "comparing_segmentation_metrics_samples_greys.png"),
    width = 12, height = 10, units = "in", res = 300)
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()




# =========== Comparing differences in cell size =========

colnames(colData(ps.5um))
#  [1] "cell"             "original_cell_id" "centroid_x"       "centroid_y"      
#  [5] "centroid_z"       "fov"              "cluster"          "volume"          
#  [9] "scale"            "population"       "sample_id"        "in_tissue"       
# [13] "brnum"            "sum"              "detected"         "total"           
# [17] "sum"              "detected"         "total"            "Segmentation"  

colnames(colData(xs.5um))
#  [1] "cell_id"                    "transcript_counts"         
#  [3] "control_probe_counts"       "control_codeword_counts"   
#  [5] "unassigned_codeword_counts" "deprecated_codeword_counts"
#  [7] "total_counts"               "cell_area"                 
#  [9] "nucleus_area"               "sample_id"                 
# [11] "brnum"                      "sum"                       
# [13] "detected"                   "total"                     
# [15] "Segmentation"   


# Assign cell size metric for each object
ps.5um$CellSize <- ps.5um$volume
ps.15um$CellSize <- ps.15um$volume
xs.5um$CellSize <- xs.5um$cell_area
xs.15um$CellSize <- xs.15um$cell_area

# Combine into single dataframe
cellsize_df <- rbind(
  as.data.frame(colData(ps.5um)[, c("Segmentation", "brnum", "sum", "CellSize")]),
  as.data.frame(colData(ps.15um)[, c("Segmentation", "brnum", "sum", "CellSize")]),
  as.data.frame(colData(xs.5um)[, c("Segmentation", "brnum", "sum", "CellSize")]),
  as.data.frame(colData(xs.15um)[, c("Segmentation", "brnum", "sum", "CellSize")])
)
colnames(cellsize_df) <- c("Segmentation", "Sample", "LibrarySize", "CellSize")

# Scatter plot function
make_scatter <- function(data, title) {
  ggplot(data, aes(x = CellSize, y = LibrarySize, color = Sample, shape = Sample)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title,
         x = "Cell Size (volume or area)",
         y = "Library Size (Total Counts)") +
    theme_bw()
}

# Subset datasets
ps_5_plot <- make_scatter(filter(cellsize_df, Segmentation == "Proseg_5um"), "Proseg 5µm")
ps_15_plot <- make_scatter(filter(cellsize_df, Segmentation == "Proseg_15um"), "Proseg 15µm")
xs_5_plot <- make_scatter(filter(cellsize_df, Segmentation == "Xenium_5um"), "Xenium 5µm")
xs_15_plot <- make_scatter(filter(cellsize_df, Segmentation == "Xenium_15um"), "Xenium 15µm")

# Combine plots
png(file.path(plot_dir, "cellsize_vs_librarysize_scatter.png"), width = 10, height = 10, units = "in", res = 300)
cowplot::plot_grid(ps_5_plot, ps_15_plot, xs_5_plot, xs_15_plot,
                   ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()



# ============== Plotting UMaps of each sample ==============

all_spe <- list(
  Proseg_5um = ps.5um,
  Proseg_15um = ps.15um,
  Xenium_5um = xs.5um,
  Xenium_15um = xs.15um
)

# ==== Create list to hold UMAP plots ====
umap_plots <- list()

# ==== Generate UMAPs for each dataset x sample ====
for (dataset_name in names(all_spe)) {
  spe_obj <- all_spe[[dataset_name]]
  
  for (sample_id in unique(spe_obj$brnum)) {
    # Subset to one sample
    spe_subset <- spe_obj[, spe_obj$brnum == sample_id]
    
    # Create UMAP plot
    p <- plotUMAP(spe_subset, colour_by = "brnum", scattermore=TRUE) +
      ggtitle(paste(dataset_name, "-", sample_id)) +
      theme(legend.position = "none",
            plot.title = element_text(size = 10, face = "bold"))
    
    # Store with unique name
    umap_plots[[paste(dataset_name, sample_id, sep = "_")]] <- p
  }
}



# ==== Save grid to file ====
png(file.path(plot_dir, "UMAPs_by_dataset_and_sample.png"),
    width = 16, height = 16, units = "in", res = 300)
cowplot::plot_grid(plotlist = umap_plots, ncol = 4)

dev.off()

png(file.path(plot_dir, "UMAPs_by_dataset_and_sample_legend.png"),
    width = 5, height = 5, units = "in", res = 300)
plotUMAP(ps.5um, colour_by = "brnum", point_size=0.3) +
  ggtitle("Proseg 5µm - Legend") +
  theme(legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"))
dev.off()

