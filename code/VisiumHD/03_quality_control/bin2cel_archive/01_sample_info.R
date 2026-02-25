# This code was modified from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/code/04_bin2cell/
#   Gather together different sources of information to form a unified table
#   of sample information

library(tidyverse)
library(here)
library(sessioninfo)

out_path = here('raw-data', 'sample_info', 'VisiumHD_sample_info.csv')

tibble(
        sample_id = c('H1-W369TJK_A1'),
        donor = c('Br9280'),
        region = c('CeA'),
        raw_image = c(
            here('raw-data', 'images', 'vis-hd','40x_HE_AMY_s003.tif')
        ),
        spaceranger_out = c(
            here('processed-data', 'VisiumHD', '01_spaceranger', 'H1-W369TJK_A1')
        )
    ) |>
    write_csv(out_path)

session_info()
