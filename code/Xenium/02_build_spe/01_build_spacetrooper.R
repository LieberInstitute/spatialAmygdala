library(here)
library(SpatialExperiment)
library(SpaceTrooper)


processed_dir <- here("processed-data", "Xenium", "02_build_spe")
xenium_data <- here("processed-data", "Xenium", "01_resegment", "5um_segmentations")

# get folder names inside xenium data dir
xenium_names <- list.dirs(xenium_data, full.names = FALSE, recursive = FALSE)
xenium_names
# [1] "20240425__170523__0022862"                            
# [2] "20240425__170523__0023004"                            
# [3] "output-XETG00117__0023153__Region_1__20240329__171204"
# [4] "output-XETG00117__0023154__Region_1__20240329__171203"

# ======= Read in all files into seperate SFE objects =======



spe.1 <- readXeniumSPE(here("processed-data","Xenium", "01_resegment", "5um_segmentations", xenium_names[1], "outs"),
    type="HDF5",
    coordNames=c("x_centroid", "y_centroid"),
    boundariesType="parquet",
    computeMissingMetrics=TRUE,
    keepPolygons=TRUE,
    countsFilePattern="cell_feature_matrix",
    metadataFPattern="cells.csv.gz",
    polygonsFPattern="cell_boundaries",
    polygonsCol="polygons",
    txPattern="transcripts",
    )
spe.1 
spe.1$brnum <- "Br9206"
spe.1$sample_id <- xenium_names[1]


spe.2 <- readXeniumSPE(here("processed-data","Xenium", "01_resegment", "5um_segmentations", xenium_names[2], "outs"),
    type="HDF5",
    coordNames=c("x_centroid", "y_centroid"),
    boundariesType="parquet",
    computeMissingMetrics=TRUE,
    keepPolygons=TRUE,
    countsFilePattern="cell_feature_matrix",
    metadataFPattern="cells.csv.gz",
    polygonsFPattern="cell_boundaries",
    polygonsCol="polygons",
    txPattern="transcripts",
    )
spe.2
spe.2$brnum <- "Br9017"
spe.2$sample_id <- xenium_names[2]


spe.3 <- readXeniumSPE(here("processed-data","Xenium", "01_resegment", "5um_segmentations", xenium_names[3], "outs"),
    type="HDF5",
    coordNames=c("x_centroid", "y_centroid"),
    boundariesType="parquet",
    computeMissingMetrics=TRUE,
    keepPolygons=TRUE,
    countsFilePattern="cell_feature_matrix",
    metadataFPattern="cells.csv.gz",
    polygonsFPattern="cell_boundaries",
    polygonsCol="polygons",
    txPattern="transcripts",
    )
spe.3
spe.3$brnum <- "Br9192"
spe.3$sample_id <- xenium_names[3]


spe.4 <- readXeniumSPE(here("processed-data","Xenium", "01_resegment", "5um_segmentations", xenium_names[4], "outs"),
    type="HDF5",
    coordNames=c("x_centroid", "y_centroid"),
    boundariesType="parquet",
    computeMissingMetrics=TRUE,
    keepPolygons=TRUE,
    countsFilePattern="cell_feature_matrix",
    metadataFPattern="cells.csv.gz",
    polygonsFPattern="cell_boundaries",
    polygonsCol="polygons",
    txPattern="transcripts",
    )
spe.4
spe.4$brnum <- "Br9280"
spe.4$sample_id <- xenium_names[4]


# ======= Save SFE objects =======


# merge
spe <- cbind(spe.1, spe.2, spe.3, spe.4)
spe
# class: SpatialExperiment 
# dim: 541 1018068 
# metadata(8): polygons technology ... polygons technology
# assays(1): counts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0344
#   DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018068): aaaaeomf-1 aaaajkhp-1 ... oimbboka-1 oimbcgpk-1
# colData names(13): cell_id transcript_counts ... polygons brnum
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id


# saveRDS
save(spe, file = here("processed-data","Xenium", "02_build_spe", "spe_xenium_5um_spacetrooper.Rdata"))

