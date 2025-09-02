# bash: conda activate r-seurat-beta

# Install beta version of Seurat and SeuratObject for reading segmentations into new conda R environment
#devtools::install_github("satijalab/seurat-object", ref = "spaceranger-4.0", force = TRUE)
#devtools::install_github("satijalab/seurat", ref = "spaceranger-4.0", force = TRUE)

library(Seurat)
library(SeuratObject)
library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(S4Vectors)
library(SpatialExperiment)

data.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/H1-W369TJK_A1/outs/"

# load data
obj <- Load10X_Spatial(data.dir = data.dir, bin.size = c(8,16,"polygons"))

# save object
saveRDS(obj, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br9280_CeA.rds"))



# NOTE - I came back tot his and ran with "module load conda_R/4.4". Not the cost conda environment above
# load
obj <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "obj_Br9280_CeA.rds"))



## EDIT THESE TWO LINES ONLY
assays_to_convert <- c("Spatial.008um", "Spatial.016um", "Spatial.Polygons")  # order matters
images_for_assays <- c("slice1.008um", "slice1.016um", "slice1.polygons")  # same order

stopifnot(length(assays_to_convert) == length(images_for_assays))

mk_spimg <- function(obj, img_id, coords) {
  # try Seurat's lowres image; else make a tiny blank canvas
  img_try <- try(Seurat::GetImage(obj, image = img_id, mode = "lowres"), silent = TRUE)
  if (!inherits(img_try, "try-error") && !is.null(img_try)) {
    SpatialExperiment::SpatialImage(x = as.raster(img_try))
  } else {
    blank <- as.raster(matrix(NA, nrow = 10, ncol = 10))
    SpatialExperiment::SpatialImage(x = blank)
  }
}

spe_list <- vector("list", length(assays_to_convert))
names(spe_list) <- assays_to_convert

for (i in seq_along(assays_to_convert)) {
  assay  <- assays_to_convert[i]
  img_id <- images_for_assays[i]

  coords <- Seurat::GetTissueCoordinates(obj, image = img_id)
  # If your columns aren't "x","y", switch to c("imagecol","imagerow")
  spatialCoords <- as.matrix(coords[, c("x","y")])
  cells <- rownames(coords)
  rownames(spatialCoords) <- cells

  counts <- SeuratObject::LayerData(obj[[assay]], layer = "counts")[, cells, drop = FALSE]
  meta <- obj@meta.data[cells, , drop = FALSE]
  meta$sample_id <- img_id

  rd <- S4Vectors::DataFrame(gene_id = rownames(counts))

  spimg <- mk_spimg(obj, img_id, coords)
  sf <- try(obj@images[[img_id]]@scale.factors$lowres, silent = TRUE)
  if (inherits(sf, "try-error") || is.null(sf)) sf <- 1

  imgData <- S4Vectors::DataFrame(
    sample_id   = img_id,   # must match colData$sample_id
    image_id    = img_id,
    data        = I(list(spimg)),   # <- SpatialImage object
    scaleFactor = as.numeric(sf)
  )

  spe_list[[i]] <- SpatialExperiment::SpatialExperiment(
    assays        = S4Vectors::SimpleList(counts = counts),
    rowData       = rd,
    colData       = S4Vectors::DataFrame(meta),
    spatialCoords = spatialCoords,
    imgData       = imgData
  )
  spe_list[[i]]$in_tissue <- 1L
}

spe_list


# save each spe seperately. 
for (i in seq_along(spe_list)) {
  assay <- assays_to_convert[i]
  spe <- spe_list[[i]]
  saveRDS(spe, here("processed-data", "VisiumHD", "02_build_spe", paste0("Br9280_CeA_", assay, "_spe.rds")))
}