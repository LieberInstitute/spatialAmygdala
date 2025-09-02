suppressPackageStartupMessages({
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("RColorBrewer")
    library("ggplot2")
    library("patchwork")
    library("DeconvoBuddies")
})

out_dir <- here("processed-data", "Visium", "08_marker_genes")

load(here("processed-data","Visium", "06_batch_correction", "spe_harmony.Rdata"))
dim(spe)

# drop low quality samples
spe <- spe[, !colData(spe)$sample_id %in% c("Br9469", "Br9017", "Br9206")]

# get folders in cluster_Csv
bs_folders <- list.files(here::here("processed-data","Visium", "07_clustering", "BayesSpace","SVGs","cluster_csv_new"), full.names = TRUE)
bs_folders
# [1] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/HVGs/cluster_csv/BayesSpace_10"
# [2] "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/08_clustering/BayesSpace/HVGs/cluster_csv/BayesSpace_12"

# get the number of clusters at the end of the folder name
bs_k <- gsub(".*BayesSpace_", "", bs_folders)

# loop through each folder, open the csv inside, and add clusters to spe colData
for (i in seq_along(bs_folders)) {
    bs_folder <- bs_folders[i]
    bs_csv <- list.files(bs_folder, full.names = TRUE)
    bs_csv <- bs_csv[grepl("csv", bs_csv)]
    bs_csv <- bs_csv[1]
    bs_df <- read.csv(bs_csv)

    colData(spe)[[paste0("BS_k", bs_k[i])]] <- factor(bs_df$cluster)
}

colnames(colData(spe))
# [43] "BS_k10"                      "BS_k11"                     
# [45] "BS_k12"                      "BS_k13"                     
# [47] "BS_k14"                      "BS_k15"                     
# [49] "BS_k16"                      "BS_k17"                     
# [51] "BS_k18"                      "BS_k19"                     
# [53] "BS_k2"                       "BS_k20"                     
# [55] "BS_k3"                       "BS_k4"                      
# [57] "BS_k5"                       "BS_k6"                      
# [59] "BS_k7"                       "BS_k8"                      
# [61] "BS_k9"   



# ======== MARKER GENES =========

# remove mito genes
mt_genes <- grepl("^MT-", rownames(spe), ignore.case = FALSE)
spe <- spe[!mt_genes, ]


# Detect BS_k* cluster columns 
bs_cols <- grep("^BS_k\\d+$", colnames(colData(spe)), value = TRUE)
if (length(bs_cols) == 0) stop("No BS_k* columns detected in colData(spe).")

# number of top genes to save
top_n <- 100

set.seed(123)
for (cl in bs_cols) {
  message("Running findMarkers on: ", cl)
  groups <- colData(spe)[[cl]]
  ok <- !is.na(groups)

  if (length(unique(groups[ok])) < 2L) {
    message("  Skipping ", cl, " (fewer than 2 groups after removing NAs).")
    next
  }

  markers <- scran::findMarkers(
    x         = spe[, ok],
    groups    = groups[ok],
    test.type = "t",
    pval.type = "any",    # faster, only need ranks
    full.stats = FALSE,
    sorted    = TRUE,
    direction = "up"
  )

  cl_names <- names(markers)
  if (all(suppressWarnings(!is.na(as.numeric(cl_names))))) {
    cl_names <- as.character(sort(as.numeric(cl_names)))
  }

  wide <- data.frame(row.names = seq_len(top_n))
  for (g in cl_names) {
    genes <- rownames(markers[[g]])
    if (length(genes) > top_n) genes <- genes[seq_len(top_n)]
    length(genes) <- top_n  # pad with NA if needed
    wide[[paste0("cluster", g)]] <- genes
  }

  fn <- file.path(out_dir, paste0(cl, "__top", top_n, "_genes.csv"))
  write.csv(wide, fn, row.names = FALSE)
  message("  Wrote: ", fn)
}

message("Done. Files in: ", normalizePath(out_dir))