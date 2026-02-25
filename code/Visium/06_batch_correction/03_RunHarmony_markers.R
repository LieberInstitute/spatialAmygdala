library("here")
library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("harmony")

processed_data <- here("processed-data", "Visium", "06_batch_correction")

# load data with PCA embeddings
load(here("processed-data", "Visium", "05_dim_reduction", "spe_stitched_pca_markers.Rdata"))
spe

# drop low quality samples
spe <- spe[, !colData(spe)$sample_id %in% c("Br9469", "Br9017", "Br9206")]

# ===== Run Harmony batch correction across sample_id and slide_id =====

set.seed(195)

message("running Harmony (sample) - ", Sys.time())
spe <- RunHarmony(spe, 
                group.by.vars = c("sample_id"), 
                reduction = "PCA", 
                assay.use = "logcounts", 
                reduction.save="PCA-HARMONY_sample"
                )

# save
saveRDS(spe, file = here(processed_data, "spe_harmony_markers.rds"))