library(here)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sf)
library(Matrix)
library(S4Vectors)
library(SpatialFeatureExperiment)
library(SpatialExperiment)

# ======== Util for reading Visium HD cell segmentations ========
# source: https://github.com/estellad/EuroBioC2025OSTAWorkshop/blob/devel/vignettes/utils.R
# I assume you could also just install the EuroBioC2025OSTAWorkshop package

flip_sf_Y <- function(sf, type = "POINT", img_height, scalef){
  st_geometry(sf) <- st_sfc( 
    lapply(st_geometry(sf), function(geom) {
      coords <- st_coordinates(geom)
      coords[, 2] <- img_height/scalef - coords[, 2] 
      if(type == "POINT"){
        st_point(coords)
      }else{ # type == "POLYGON"
        st_polygon(list(matrix(coords[, 1:2], ncol = 2)))
      }
    }),
    crs = st_crs(sf)
  )
  
  return(sf)
}


shift_to_origin <- function(sf_object) {
  
  bbox <- st_bbox(sf_object)
  min_x <- bbox["xmin"]
  min_y <- bbox["ymin"]
  
  cat("Original bounding box:\n")
  print(bbox)
  cat("\nShifting by: x =", -min_x, ", y =", -min_y, "\n")
  
  # Subtract the minimums from all coordinates
  shifted_sf <- sf_object
  st_geometry(shifted_sf) <- st_geometry(sf_object) - c(min_x, min_y)
  
  cat("\nShifted bounding box:\n")
  print(st_bbox(shifted_sf))
  
  return(shifted_sf)
}


readVisiumHDCellSeg <- function(td, res){
  if(res == "hires"){
    png_name = "tissue_hires_image.png"
    scalef_name = "tissue_hires_scalef"
  }else{ # res = "lowres"
    png_name = "tissue_lowres_image.png"
    scalef_name = "tissue_lowres_scalef"
  }
  
  ## Coord ---------------------------------------------------------------
  # Read the GeoJSON file of cell segmentation
  geo_data <- st_read(file.path(td, "segmented_outputs", 
                                "cell_segmentations.geojson"))

  geo_data <- shift_to_origin(geo_data)
  scalef <- rjson::fromJSON(file = file.path(td, "segmented_outputs", "spatial", 
                                             "scalefactors_json.json"))
  image <- magick::image_read(file.path(td, "segmented_outputs", 
                                "spatial", png_name))
  info <- magick::image_info(image)
  
  # Get centroid
  st_crs(geo_data) <- NA
  rownames(geo_data) <- geo_data$cell_id
  centroids <- st_centroid(geo_data)
  centroids$cell_id <- as.character(centroids$cell_id)
  
  
  ## Countmat ------------------------------------------------------------
  countmat_file <- file.path(td, "segmented_outputs", 
                             "filtered_feature_cell_matrix.h5")
  vhdcellsce <- DropletUtils::read10xCounts(countmat_file, col.names = TRUE)
  
  colnames(vhdcellsce) <- sub("cellid_0*([0-9]+)-.*", "\\1", vhdcellsce$Barcode)
  
  rownames(vhdcellsce) <- rowData(vhdcellsce)$Symbol
  rownames(vhdcellsce) <- make.unique(rownames(vhdcellsce))
  
  
  ## Matching coordinates with count matrix ------------------------------
  # Subset to common cells
  common_cells <- intersect(centroids$cell_id, colnames(vhdcellsce))
  
  geo_data <- geo_data[geo_data$cell_id %in% common_cells, ]
  centroids <- centroids[centroids$cell_id %in% common_cells, ]
  vhdcellsce <- vhdcellsce[, colnames(vhdcellsce) %in% common_cells]
  
  # # Sanity check on ordering of cell ids
  # all(as.character(centroids$cell_id) == colnames(vhdcellsce))
  
  # Extract coordinates
  coords <- st_coordinates(centroids)
  colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

  # normalize coordinates by subtracting min
  coords <- as.data.frame(coords)
  coords$pxl_col_in_fullres <- coords$pxl_col_in_fullres - min(coords$pxl_col_in_fullres)
  coords$pxl_row_in_fullres <- coords$pxl_row_in_fullres - min(coords$pxl_row_in_fullres)
  print(min(coords$pxl_col_in_fullres))
  print(min(coords$pxl_row_in_fullres))


  ## Construct SPE -------------------------------------------------------
  vhd <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = as(counts(vhdcellsce), "dgCMatrix")),
    rowData = rowData(vhdcellsce),
    colData = cbind(colData(vhdcellsce), coords),
    spatialCoordsNames = colnames(coords),
    scaleFactors = scalef$tissue_hires_scalef, 
    imageSources = file.path(td, "segmented_outputs", "spatial", png_name), 
    loadImage = TRUE
  )
  imgData(vhd)$image_id <- res
  
  # vhd
  
  ## Coerce to SFE -------------------------------------------------------
  vhdsfe <- toSpatialFeatureExperiment(vhd)
  # vhdsfe
  
  ## Add polygons to colGeometries
  colGeometries(vhdsfe)$cellSeg <- geo_data
  
  
  # Flip Y in colGeometries -------------------------------------------------
  # Flip Centroid Y
  # centroids <- colGeometries(vhdsfe)$centroids
  # colGeometries(vhdsfe)$updatecentroids <- 
  #   flip_sf_Y(centroids, type = "POINT", img_height = info$height, 
  #             scalef = imgData(vhdsfe)$scaleFactor) 
  
  # # Flip Cellseg Y
  # cellSeg <- colGeometries(vhdsfe)$cellSeg
  # colGeometries(vhdsfe)$updatecellSeg <- 
  #   flip_sf_Y(cellSeg, type = "POLYGON", img_height = info$height, 
  #             scalef = imgData(vhdsfe)$scaleFactor) 
  
  # Flip Centroid Y in `spatialCoords()`
  # spatialCoords(vhdsfe)[, "pxl_row_in_fullres"] <- 
  #   info$height/imgData(vhdsfe)$scaleFactor - 
  #   spatialCoords(vhdsfe)[, "pxl_row_in_fullres"]
  
  # Add high res coords to `colData()`
  # vhdsfe[[paste0("pxl_col_in_", res)]] <- 
  #   spatialCoords(vhdsfe)[, "pxl_col_in_fullres"] * imgData(vhdsfe)$scaleFactor
  # vhdsfe[[paste0("pxl_row_in_", res)]] <- 
  #   spatialCoords(vhdsfe)[, "pxl_row_in_fullres"] * imgData(vhdsfe)$scaleFactor
  
  return(vhdsfe)
}



# ========== File to load ==========

# first sample
data1.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/H1-W369TJK_A1/outs"
data2.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_A1/outs"
data3.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_D1/outs"
data4.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_A1/outs"
data5.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_D1/outs"


sfe <- readVisiumHDCellSeg(data1.dir, res = "hires")
sfe$Sample <- NULL
sfe$Sample <- "Br9280_CeA"
sfe <- changeSampleIDs(sfe, c(sample01="H1-W369TJK_A1"))

sfe2 <- readVisiumHDCellSeg(data2.dir, res = "hires")   
sfe2$Sample <- NULL
sfe2$Sample <- "Br8325_MeA"
sfe2 <- changeSampleIDs(sfe2, c(sample01="H1-937FVHX_A1"))

sfe3 <- readVisiumHDCellSeg(data3.dir, res = "hires")
sfe3$Sample <- NULL
sfe3$Sample <- "Br8325_CeA"
sfe3 <- changeSampleIDs(sfe3, c(sample01="H1-937FVHX_D1"))

sfe4 <- readVisiumHDCellSeg(data4.dir, res = "hires")
sfe4$Sample <- NULL
sfe4$Sample <- "Br9280_ITC"
sfe4 <- changeSampleIDs(sfe4, c(sample01="H1-HW9VGBW_A1"))

sfe5 <- readVisiumHDCellSeg(data5.dir, res = "hires")
sfe5$Sample <- NULL
sfe5$Sample <- "Br9280_MeA"
sfe5 <- changeSampleIDs(sfe5, c(sample01="H1-HW9VGBW_D1"))




# save each sfe seperately.
saveRDS(sfe, here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_CeA.rds"))
saveRDS(sfe2, here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br8325_MeA.rds"))
saveRDS(sfe3, here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br8325_CeA.rds"))
saveRDS(sfe4, here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_ITC.rds"))
saveRDS(sfe5, here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_MeA.rds"))


# load
sfe <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_CeA.rds"))
sfe2 <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br8325_MeA.rds"))
sfe3 <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br8325_CeA.rds"))
sfe4 <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_ITC.rds"))
sfe5 <- readRDS(here("processed-data", "VisiumHD", "02_build_spe", "sfe_Br9280_MeA.rds"))

# merge using cbind
vhd_all <- SpatialFeatureExperiment::cbind(sfe, sfe2, sfe3, sfe4, sfe5)
vhd_all
# class: SpatialFeatureExperiment 
# dim: 18085 284986 
# metadata(0):
# assays(1): counts
# rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
# rowData names(3): ID Symbol Type
# colnames(284986): 4 7 ... 63230-2 63232-2
# colData names(5): Barcode sample_id pxl_col_in_hires pxl_row_in_hires
#   Sample
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor

# unit:
# Geometries:
# colGeometries: centroids (POINT), cellSeg (POLYGON), updatecentroids (POINT), updatecellSeg (POLYGON) 

# Graphs:
# H1-W369TJK_A1: 
# H1-937FVHX_A1: 
# H1-937FVHX_D1: 
# H1-HW9VGBW_A1: 
# H1-HW9VGBW_D1:

# save
saveRDS(vhd_all, here("processed-data", "VisiumHD", "02_build_spe", "sfe_combined.rds"))