setwd("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/")
suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
#suppressPackageStartupMessages(library("anndata"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
library(HDF5Array)
library(dplyr)


#spe <- readRDS(here("processed-data/rds/02_visium_qc","qc_spe_w_spg_N63.rds"))
#spe = readRDS(
#  "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/publication/spe_n24_example-samples.rds"
#  )
#old-object
#spe_path= "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/VisiumHD/Samui_app/spe_cells_with_domain.rds"
spe_path= "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/07_img_cropping/spe_cropped.rds"
spe =readRDS(spe_path)
coords <- DataFrame(spatialCoords(spe))
colData(spe) <- cbind(colData(spe), coords)
## Load SpD data ----
#finalized_spd <- readRDS(here("processed-data/rds/spatial_cluster", "PRECAST","test_clus_label_df_semi_inform_k_2-16.rds"))

## Attach SpD label to spe ----
#col_data_df <- colData(spe) |>
#  data.frame() |>
#  left_join(
#    finalized_spd,
#    by = c("key"),
#    relationship = "one-to-one"
#  )

#rownames(col_data_df) <- colnames(spe)
#colData(spe) <- DataFrame(col_data_df)

# Call SPG spots ----
#spe$pnn_pos <- ifelse(spe$spg_PWFA > 0.05, TRUE, FALSE)
## NOTE: neuropil spot are spots doesn't have DAPI staining
#spe$neuropil_pos <- ifelse(spe$spg_PDAPI > 0.05,FALSE, TRUE)
#spe$neun_pos <- ifelse(spe$spg_PNeuN > 0.05 & spe$spg_PNeuN < 0.3,TRUE, FALSE)
#spe$vasc_pos <- ifelse(spe$spg_PClaudin5 > 0.05 & spe$spg_PClaudin5 < 0.20,TRUE, FALSE)

spe_out <- here("processed-data", "samui", "visiumHD", "spe_visiumHD.h5ad")

write_anndata <- function(sce, out_path) {
  invisible(
    basiliskRun(
      fun = function(sce, filename) {
        library("zellkonverter")
        library("reticulate")

        # Convert SCE to AnnData:
        adata <- SCE2AnnData(sce)

        #  Write AnnData object to disk
        adata$write(filename = filename)

        return()
      },
      env = zellkonverterAnnDataEnv(),
      sce = sce,
      filename = out_path
    )
  )
}


#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe)$spatial <- spatialCoords(spe)

write_anndata(spe, spe_out)
#saveHDF5SummarizedExperiment(spe)
#writeHDF5Array()

#writeH5AD(spe, spe_out)

session_info()

#brnums <- unique(colData(spe)$brnum)
sample_ids <- unique(colData(spe)$sample_id)
sample_ids <- as.character(sample_ids)
#writeLines(brnums, con = '/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/11_samui/sample_list.txt')
writeLines(sample_ids, con = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/samui/visiumHD/sample_ids-visium.txt')

