library(here)
library(Matrix)
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(Voyager)
library(VisiumIO)

# ======== Define directories for all 5 samples ========
data1.dir <- here("processed-data", "VisiumHD", "01_spaceranger", 
                   "H1-W369TJK_A1", "outs")
data2.dir <- here("processed-data", "VisiumHD", "01_spaceranger", 
                   "H1-937FVHX_A1", "outs")
data3.dir <- here("processed-data", "VisiumHD", "01_spaceranger", 
                   "H1-937FVHX_D1", "outs")
data4.dir <- here("processed-data", "VisiumHD", "01_spaceranger", 
                   "H1-HW9VGBW_A1", "outs")
data5.dir <- here("processed-data", "VisiumHD", "01_spaceranger", 
                   "H1-HW9VGBW_D1", "outs")

# ======== Helper: read 16µm bins and convert to SPE ========
spe.1 <- TENxVisiumHD(
    spacerangerOut=data1.dir,
    processing="raw", 
    images="lowres", 
    bin_size="016") |>
    import()

# use gene symbols as feature names
gs <- rowData(spe.1)$Symbol
rownames(spe.1) <- make.unique(gs)

spe.2 <- TENxVisiumHD(
    spacerangerOut=data2.dir,
    processing="raw",
    images="lowres", 
    bin_size="016") |>
    import()

gs <- rowData(spe.2)$Symbol
rownames(spe.2) <- make.unique(gs)

spe.3 <- TENxVisiumHD(
    spacerangerOut=data3.dir,
    processing="raw",
    images="lowres", 
    bin_size="016") |>
    import()

gs <- rowData(spe.3)$Symbol
rownames(spe.3) <- make.unique(gs)

spe.4 <- TENxVisiumHD(
    spacerangerOut=data4.dir,
    processing="raw",
    images="lowres", 
    bin_size="016") |>
    import()        

gs <- rowData(spe.4)$Symbol
rownames(spe.4) <- make.unique(gs)

spe.5 <- TENxVisiumHD(
    spacerangerOut=data5.dir,
    processing="raw",
    images="lowres", 
    bin_size="016") |>
    import()

gs <- rowData(spe.5)$Symbol
rownames(spe.5) <- make.unique(gs)

spe.1$sample_id <- "Br9280_CeA"
spe.2$sample_id <- "Br8325_MeA"
spe.3$sample_id <- "Br8325_CeA"
spe.4$sample_id <- "Br9280_ITC"
spe.5$sample_id <- "Br9280_MeA"

# Check for mismatched rowData columns before combining
# The 'Type' column from Space Ranger can differ across runs
for (obj_name in c("spe.1", "spe.2", "spe.3", "spe.4", "spe.5")) {
    obj <- get(obj_name)
    message(obj_name, " rowData cols: ", paste(colnames(rowData(obj)), collapse = ", "))
    if ("Type" %in% colnames(rowData(obj))) {
        message("  Type values: ", paste(unique(rowData(obj)$Type), collapse = ", "))
    }
}

# Drop the 'Type' column from all objects to avoid cbind mismatch
for (obj_name in c("spe.1", "spe.2", "spe.3", "spe.4", "spe.5")) {
    obj <- get(obj_name)
    if ("Type" %in% colnames(rowData(obj))) {
        rowData(obj)$Type <- NULL
        assign(obj_name, obj)
    }
}

# combine
spe.016 <- cbind(spe.1, spe.2, spe.3, spe.4, spe.5)


# ======== Save each separately ========
saveRDS(spe.016, file = here("processed-data", "VisiumHD", "02_build_spe", "spe_016_combined.rds"))