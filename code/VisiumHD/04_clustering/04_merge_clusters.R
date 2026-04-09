library("SpatialFeatureExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("ggspavis")
library("scCustomize")
library("ggplot2")

# ======= Load base spe =======
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_016.rds"))

# drop out of tissue for all except Br8325_MeA
spe <- spe[, spe$in_tissue | spe$sample_id == "Br8325_MeA"]

# ======= Load Banksy output for resolution 0.6 =======
banksy_dir <- here("processed-data", "VisiumHD", "04_clustering", "016_wAI_markers")
banksy_file <- file.path(banksy_dir, "Banksy_integrated_lambda_0.8_res0.6.rds")

if (!file.exists(banksy_file)) {
    stop("Banksy file for resolution 0.6 not found: ", banksy_file)
}

cat("Loading Banksy output for resolution 0.6\n")
sfe.sub <- readRDS(banksy_file)

# Find the cluster column
clust_cols <- grep("clust_HARMONY_M1_lam0\\.8", colnames(colData(sfe.sub)), value = TRUE)
if (length(clust_cols) == 0) {
    stop("No cluster column found matching pattern clust_HARMONY_M1_lam0.8")
}
clust_col <- clust_cols[1]
cat("Cluster column:", clust_col, "\n")

# Transfer labels to full spe
spe[[clust_col]] <- factor(sfe.sub[[clust_col]])

# ======= Define cluster renaming =======
# Default labels for all samples
default_labels <- c(
    "1"  = "MeA",
    "2"  = "BM",
    "3"  = "Neuropil",
    "4"  = "BLD",
    "5"  = "BM",
    "6"  = "WM",
    "7"  = "Ventricle",
    "8"  = "BLD",
    "9"  = "CeA",
    "10" = "CoA",
    "11" = "Neuropil",
    "12" = "BM",
    "13" = "AI",
    "14" = "14",
    "15" = "DELETE",
    "16" = "BM",
    "17" = "BM",
    "18" = "Neuropil",
    "19" = "CHAT",
    "20" = "DELETE",
    "21" = "DELETE",
    "22" = "DELETE"
)

# Sample-specific overrides: list of sample_id -> list of cluster -> label
sample_overrides <- list(
    "Br9280_ITC" = c("2" = "LA"),
    "Br8325_CeA" = c("12" = "BLD")
)

# ======= Apply renaming per spot (respecting sample-specific overrides) =======
original_clusters <- as.character(spe[[clust_col]])
sample_ids <- as.character(spe$sample_id)
new_labels <- character(length(original_clusters))

for (i in seq_along(original_clusters)) {
    cl <- original_clusters[i]
    sid <- sample_ids[i]

    # Check for sample-specific override first
    if (sid %in% names(sample_overrides) && cl %in% names(sample_overrides[[sid]])) {
        new_labels[i] <- sample_overrides[[sid]][cl]
    } else {
        new_labels[i] <- default_labels[cl]
    }
}

# ======= Remove DELETE spots =======
keep <- new_labels != "DELETE"
cat("Removing", sum(!keep), "spots marked DELETE\n")
spe <- spe[, keep]
new_labels <- new_labels[keep]

# Store renamed labels
label_col <- "spatial_domain"
spe[[label_col]] <- factor(new_labels)

cat("Final spatial domains:\n")
print(table(spe[[label_col]]))

# ======= Determine colors =======
domain_levels <- sort(unique(new_labels))
num_groups <- length(domain_levels)
cat("Number of spatial domains:", num_groups, "\n")

colors <- scCustomize_Palette(
    num_groups,
    ggplot_default_colors = FALSE,
    color_seed = 123
)
names(colors) <- domain_levels

# ======= Plot (same style as original holistic plot) =======
plots_dir <- here("plots", "VisiumHD", "04_clustering", "016_wAI_markers")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(plots_dir, "Banksy_integrated_lambda_0.8_res0.6_renamed.pdf")
pdf(file = pdf_path, width = 40, height = 10)

p <- plotCoords(spe, annotate = label_col, in_tissue = NULL, sample_id = "sample_id") +
    scale_color_manual(values = colors) +
    ggtitle(paste0("Banksy Integrated Clustering (lambda=0.8, res=0.6, ", num_groups, " spatial domains)")) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

print(p)
dev.off()
cat("\nDone! Renamed plot saved to:", pdf_path, "\n")