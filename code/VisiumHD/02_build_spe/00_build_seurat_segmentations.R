# I installed the beta version of Seurat and SeuratObject for reading segmentations into distinct R library
#devtools::install_github("satijalab/seurat-object", ref = "spaceranger-4.0", force = TRUE)
#devtools::install_github("satijalab/seurat", ref = "spaceranger-4.0", force = TRUE)

# to reaccess run: 
# module load conda_R/4.4
# export SEURAT_BETA_LIB=$HOME/Rlibs-seurat-beta-4.4
# R

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






# first sample
data1.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/H1-W369TJK_A1/outs/"
obj1  <- Load10X_Spatial(data.dir = data1.dir, bin.size = c(8,16,"polygons"))
saveRDS(obj1, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br9280_CeA.rds"))
rm(obj1)

# last four samples
data2.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_A1/H1-937FVHX_A1/outs/"
obj2  <- Load10X_Spatial(data.dir = data2.dir, bin.size = c(8,16,"polygons"))
saveRDS(obj2, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br8325_MeA.rds"))
rm(obj2)

data3.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_D1/H1-937FVHX_D1/outs/"
obj3  <- Load10X_Spatial(data.dir = data3.dir, bin.size = c(8,16,"polygons"))
saveRDS(obj3, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br8325_CeA.rds"))
rm(obj3)

data4.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_A1/H1-HW9VGBW_A1/outs/"
obj4  <- Load10X_Spatial(data.dir = data4.dir, bin.size = c(8,16,"polygons"))
saveRDS(obj4, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br9280_ITC.rds"))
rm(obj4)

data5.dir <- "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_D1/H1-HW9VGBW_D1/outs/"
obj5  <- Load10X_Spatial(data.dir = data5.dir, bin.size = c(8,16,"polygons"))
saveRDS(obj5, file = here("processed-data", "VisiumHD", "02_build_spe", "obj_Br9280_MeA.rds"))
rm(obj5)



