library("here")
library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("harmony")

processed_data <- here("processed-data", "Visium", "06_batch_correction")

# load data with PCA embeddings
load(here("processed-data", "Visium", "05_dim_reduction", "spe_stitched_pca.Rdata"))
spe

# view colData
colnames(colData(spe))

unique(spe$capture_area)

# make new Slide column by dropping _X
spe$slide_id <- gsub("_.*", "", spe$capture_area)
unique(spe$slide_id)

# ===== Run Harmony batch correction across sample_id and slide_id =====

set.seed(195)
message("running Harmony (sample + slide) - ", Sys.time())
spe <- RunHarmony(spe, 
                group.by.vars = c("sample_id", "slide_id"), 
                reduction = "PCA", 
                assay.use = "logcounts", 
                reduction.save="PCA-HARMONY_sample_slide"
                )


message("running Harmony (sample) - ", Sys.time())
spe <- RunHarmony(spe, 
                group.by.vars = c("sample_id"), 
                reduction = "PCA", 
                assay.use = "logcounts", 
                reduction.save="PCA-HARMONY_sample"
                )

# save
saveRDS(spe, file = here(processed_data, "spe_harmony.rds"))