
#library(remotes)
#remotes::install_github("pachterlab/SpatialFeatureExperiment", ref = "devel")

library(here)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(Voyager)

# define directories
processed_dir <- 
processed_dir <- here("processed-data", "VisiumHD")
plot_dir <- here("plots","VisiumHD", "mouse_brain")


hd_dir <-  c("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/H1-W369TJK_A1/outs/binned_outputs")



sfe <- read10xVisium(samples=hd_dir, sample_id="H1-W369TJK_A1",
    type = "HDF5", data = "filtered")


sfe <- readVisiumHD(data_dir=hd_dir,
                  #sample_id = "H1-W369TJK_A1",
                  bin_size = "16", # this defines which of 1:3 resolutions to load
                  type = "HDF5", # Note, "sparse" -> takes longer to load
                  data = "filtered", # spots under tissue
                  images = "lowres" 
                  )
