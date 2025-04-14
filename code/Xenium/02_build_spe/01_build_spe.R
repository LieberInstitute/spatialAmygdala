library(SpatialExperimentIO)
library(SpatialExperiment)
library(here)

processed_dir <- here("processed-data", "Xenium", "02_build_spe")
xenium_data <- here("processed-data", "Xenium", "01_resegment")

# get folder names inside xenium data dir
xenium_names <- list.dirs(xenium_data, full.names = FALSE, recursive = FALSE)
xenium_names
# [1] "20240425__170523__0022862"                            
# [2] "20240425__170523__0023004"                            
# [3] "output-XETG00117__0023153__Region_1__20240329__171204"
# [4] "output-XETG00117__0023154__Region_1__20240329__171203"


spe.1 <- readXeniumSXE(here("processed-data","Xenium", "01_resegment", xenium_names[1], "outs"))
spe.1 
spe.1$brnum <- "Br9206"
spe.1$sample_id <- xenium_names[1]
# class: SpatialExperiment 
# dim: 541 287576 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(287576): aaaadgkh-1 aaaadlfe-1 ... oihlfjab-1 oihliodj-1
# colData names(10): cell_id transcript_counts ... nucleus_area sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):

spe.2 <- readXeniumSXE(here("processed-data","Xenium", "01_resegment", xenium_names[2], "outs"))
spe.2
spe.2$brnum <- "Br9017"
spe.2$sample_id <- xenium_names[2]
# class: SpatialExperiment 
# dim: 541 252805 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(252805): aaaaecde-1 aaaaiook-1 ... oijepiab-1 oijfcgml-1
# colData names(10): cell_id transcript_counts ... nucleus_area sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):

spe.3 <- readXeniumSXE(here("processed-data","Xenium", "01_resegment", xenium_names[3], "outs"))
spe.3
spe.3$brnum <- "Br9192"
spe.3$sample_id <- xenium_names[3]
# class: SpatialExperiment 
# dim: 541 279464 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(279464): aaaajbpi-1 aaabched-1 ... oilldfpe-1 oillfmdl-1
# colData names(10): cell_id transcript_counts ... nucleus_area sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

spe.4 <- readXeniumSXE(here("processed-data","Xenium", "01_resegment", xenium_names[4], "outs"))
spe.4
spe.4$brnum <- "Br9280"
spe.4$sample_id <- xenium_names[4]
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):
# class: SpatialExperiment 
# dim: 541 198224 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(198224): aaaadokp-1 aaaafdph-1 ... oihobmbk-1 oihoeehh-1
# colData names(10): cell_id transcript_counts ... nucleus_area sample_id
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(0):


# combine all the spatial experiments
spe <- cbind(spe.1, spe.2, spe.3, spe.4)
spe
# class: SpatialExperiment 
# dim: 541 1018069 
# metadata(0):
# assays(1): counts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018069): aaaadgkh-1 aaaadlfe-1 ... oihobmbk-1 oihoeehh-1
# colData names(11): cell_id transcript_counts ... sample_id brnum
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# save
save(spe, file = here("processed-data","Xenium", "02_build_spe", "spe_combined.Rdata"))