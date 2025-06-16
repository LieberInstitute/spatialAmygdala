setwd('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')
suppressPackageStartupMessages({
    library("dplyr")
    library("purrr")
    library("Seurat")
    library("here")
    library("sessioninfo")
    library("SpatialExperiment")
    library("PRECAST")
})

#start from spe without batch correction
load(here("processed-data", "04_normalization", "spe_stitched_norm.Rdata"))

colnames(spe) <- spe$key

seuList <- unique(spe$sample_id) |>
    set_names(unique(spe$sample_id)) |>
    map(.f = function(id) {
        tmp_spe <- spe[, spe$sample_id == id]

        tmp_spe$row <- tmp_spe$array_row
        tmp_spe$col <- tmp_spe$array_col

        # browser()
        CreateSeuratObject(
            counts=as.matrix(counts(tmp_spe)),
            meta.data=data.frame(colData(tmp_spe)),
            project="AMY")
    })

set.seed(1)
preobj <- CreatePRECASTObject(seuList = seuList, gene.number=4000, selectGenesMethod='HVGs',
                              premin.spots = 1, premin.features=1, postmin.spots=1, postmin.features=1)
preobj@seulist

PRECASTObj <- AddAdjList(preobj, platform = "Visium")

PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE,  maxIter = 30, verbose = TRUE)

K <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

PRECASTObj <- PRECAST(PRECASTObj, K = K)

save(PRECASTObj, file = here("processed-data", "08_clustering", "PRECAST","HVGs","Rdata_objects", paste0("PRECASTObj_",K,".Rdata")))