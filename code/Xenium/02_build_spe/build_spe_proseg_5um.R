library(SpatialExperimentIO)
library(SpatialExperiment)
library(here)
library(Seurat)

xenium_data <- here("processed-data","Xenium", "01_resegment","proseg_5um")

# get folder names inside xenium data dir
xenium_names <- list.dirs(xenium_data, full.names = FALSE, recursive = FALSE)
xenium_names



# ======= Function to read into SPE =======
# modified from: https://github.com/dcjones/proseg/blob/main/extra/proseg-to-seurat.R

ProsegToSpatialExperiment <- function (
    proseg_output_path, 
    sample_id = "sample",
    expected_counts_basename = "expected-counts",
    cell_metadata_basename = "cell-metadata"
) {
  ## Load expected counts
  expected_counts_path <- file.path(proseg_output_path, paste0(expected_counts_basename, ".csv.gz"))
  if (!file.exists(expected_counts_path)) {
    expected_counts_path <- file.path(proseg_output_path, paste0(expected_counts_basename, ".parquet"))
    if (!file.exists(expected_counts_path)) {
      stop("Can't find expected-counts file.")
    }
    expected_counts <- arrow::read_parquet(expected_counts_path)
  } else {
    expected_counts <- read.csv(expected_counts_path, header = TRUE, sep = ",")
  }

  ## Load cell metadata
  cell_metadata_path <- file.path(proseg_output_path, paste0(cell_metadata_basename, ".csv.gz"))
  if (!file.exists(cell_metadata_path)) {
    cell_metadata_path <- file.path(proseg_output_path, paste0(cell_metadata_basename, ".parquet"))
    if (!file.exists(cell_metadata_path)) {
      stop("Can't find cell-metadata file.")
    }
    cell_metadata <- arrow::read_parquet(cell_metadata_path)
  } else {
    cell_metadata <- read.csv(cell_metadata_path, header = TRUE, sep = ",")
  }

  ## Exclude cells with undefined centroids
  mask <- is.finite(cell_metadata$centroid_x) & is.finite(cell_metadata$centroid_y)
  cell_metadata <- cell_metadata[mask, ]
  expected_counts <- expected_counts[mask, ]

  ## Format counts matrix
  counts <- Matrix::Matrix(t(as.matrix(expected_counts)), sparse = TRUE)

  ## Spatial coordinates
  spatialCoords <- as.matrix(cell_metadata[, c("centroid_x", "centroid_y")])
  colnames(spatialCoords) <- c("x", "y")
  rownames(spatialCoords) <- rownames(cell_metadata) <- colnames(counts)

  ## Build SpatialExperiment without imgData
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = cell_metadata,
    spatialCoords = spatialCoords,
    sample_id = sample_id
  )

  ## Add sample ID and in_tissue flag
  spe$sample_id <- sample_id
  spe$in_tissue <- 1

  return(spe)
}



# ======== Read into proseg to SPE ========
spe.1 <- ProsegToSpatialExperiment(
  proseg_output_path = file.path(xenium_data, xenium_names[1]),
  sample_id = xenium_names[1]
)
spe.1$brnum <- "Br9206"
spe.1

spe.2 <- ProsegToSpatialExperiment(
  proseg_output_path = file.path(xenium_data, xenium_names[2]),
  sample_id = xenium_names[2]
)
spe.2$brnum <- "Br9017"
spe.2

spe.3 <- ProsegToSpatialExperiment(
  proseg_output_path = file.path(xenium_data, xenium_names[3]),
  sample_id = xenium_names[3]
)
spe.3$brnum <- "Br9192"
spe.3

spe.4 <- ProsegToSpatialExperiment(
  proseg_output_path = file.path(xenium_data, xenium_names[4]),
  sample_id = xenium_names[4]
)
spe.4$brnum <- "Br9280"
spe.4


# combine all the spatial experiments
spe <- cbind(spe.1, spe.2, spe.3, spe.4)
spe

# save
save(spe, file = here("processed-data","Xenium", "02_build_spe", "spe_proseg_5um.Rdata"))