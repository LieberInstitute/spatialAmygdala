#   Create input images to ImageJ, and write sample info for use later when
#   reading in ImageJ outputs and building the SpatialExperiment
args <- commandArgs(trailingOnly = TRUE)
group <- args[1]
capture_area <- eval(parse(text = args[2]))
group
capture_area

setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
library(here)
library(tidyverse)
library(visiumStitched)
library(sessioninfo)

imagej_dir = here('processed-data', 'visium_stitching', 'NacUtils', 'imagej', group)
info_out_path = here('processed-data', 'visium_stitching', 'NacUtils', paste0(group,'.csv'))

sample_info = tibble(
    group = group,
    capture_area = capture_area,
    imagej_xml_path = file.path(imagej_dir, paste0(group, '.xml')),
    imagej_image_path = file.path(imagej_dir, paste0(group, '.png')),
    spaceranger_dir = here('processed-data', '01_spaceranger', capture_area, 'outs', 'spatial')
)

in_dir = file.path(imagej_dir, 'input')
dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)

rescale_imagej_inputs(sample_info, in_dir) |>
    write_csv(info_out_path)

session_info()