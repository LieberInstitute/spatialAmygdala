
suppressPackageStartupMessages({
    library("here")
    library("SpatialExperiment")
    library("spatialLIBD")
    library("rtracklayer")
    library("lobstr")
    library("sessioninfo")
})

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_AMY.rda"))

sample_info <- data.frame(slide = as.factor(REDCap_AMY$slide))
sample_info$array <- as.factor(REDCap_AMY$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_AMY$sample, "-"), `[`, 1))
sample_info$species <- as.factor(REDCap_AMY$species)
sample_info$replicate <- as.factor(REDCap_AMY$serial)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")
sample_info$sample_path = file.path(here::here("processed-data", "01_spaceranger", "first_donor"), sample_info$sample_id, "outs")

head(sample_info)

## Build basic SPE
Sys.time()
#spe <- read10xVisiumWrapper(
#    sample_info$sample_path,
#    sample_info$sample_id,
#    type = "sparse",
#    data = "raw",
#    images = c("lowres", "hires", "detected", "aligned"),
#    load = TRUE,
#    reference_gtf = file.path("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/","genes", "genes.gtf")
#)
Sys.time()
#save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))
load(here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))
spe

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))
spe

## Read in cell counts and segmentation results
#segmentations_list <-
  lapply(sample_info$sample_id, function(sampleid) {
    file <-
      here(
        "processed-data",
        "01_spaceranger",
        "first_donor",
        sampleid,
        "outs",
        "spatial",
        "tissue_spot_counts.csv"
      )
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
  })
#segmentations_list
#
### Merge them (once the these files are done, this could be replaced by an rbind)
#segmentations <-
#  Reduce(function(...) {
#    merge(..., all = TRUE)
#  }, segmentations_list[lengths(segmentations_list) > 0])
#  segmentations
#
### Add the information
#segmentation_match <- match(spe$key, segmentations$key)
#segmentation_info <-
#  segmentations[segmentation_match, -which(
#    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#  )]
#segmentation_match
#
#colData(spe) <- cbind(colData(spe), segmentation_info)
#colData(spe)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 6345
length(no_expr) / nrow(spe) * 100
# [1] 17.33559
spe <- spe[-no_expr, ]


## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
  factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
#spe_raw <- spe

save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))
spe

## Size in Gb
lobstr::obj_size(spe)
# 2.07GB


## Now drop the spots outside the tissue
spe_raw <- spe
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 36601 37298
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
  message("removing spots without counts for spe")
  spe <- spe[, -which(colSums(counts(spe)) == 0)]
  dim(spe)
}

# removing spots without counts for spe
# [1] 30256 37290

lobstr::obj_size(spe)
# 2.04 GB

save(spe, file = here::here("processed-data", "02_build_spe", "spe.Rdata"))
spe

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()