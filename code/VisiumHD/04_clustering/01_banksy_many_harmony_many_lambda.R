library("SpatialFeatureExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("Banksy")
library("harmony")
library("data.table")

# ======= Parse resolution from command line =======
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript 01_banksy_many_harmony_many_lambda.R <resolution>\n",
         "  e.g. Rscript 01_banksy_many_harmony_many_lambda.R 0.4")
}
res <- as.numeric(args[1])
cat("Resolution:", res, "\n")

# set directories
plots_dir <- here("plots", "VisiumHD", "04_clustering")
processed_dir <- here("processed-data", "VisiumHD", "04_clustering", "016_wAI_markers")

# load
spe <- readRDS(here("processed-data", "VisiumHD", "03_quality_control", "spe_qc_016.rds"))

# drop out of tissue for all except Br8325_MeA
spe <- spe[, spe$in_tissue | spe$sample_id == "Br8325_MeA"]


# load Visium domain marker genes
markers <- read.csv(here("processed-data", "Visium", "08_marker_genes", "top100_markers_bs_final_ITC_smoothed_wBLVM.csv"))

# get AI markers. Columns "cluster" and "gene"
ai_markers <- markers %>%
    filter(cluster == "AI") %>%
    pull(gene) %>%
    unique()

# ======= Feature Selection =======
spe <- spe[!grepl("^MT-", rownames(spe)), ]
spe <- logNormCounts(spe)
top <- modelGeneVar(spe)
hvg <- getTopHVGs(top, n = 2000)

# add AI markers to HVG list
hvg <- unique(c(hvg, ai_markers))
spe.sub <- spe[hvg, ]

# ======== Banksy clustering ========
locs <- spatialCoords(spe.sub)
locs <- cbind(locs, sample_id = factor(spe.sub$sample_id))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx - min(sdimx), by = group]
global_max <- max(locs_dt$sdimx) * 1.5
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe.sub)
spatialCoords(spe.sub) <- locs

lambda <- 0.8
k_geom <- c(25, 50)
npcs <- 50
aname <- "logcounts"


set.seed(1000)
spe.sub <- Banksy::computeBanksy(spe.sub, assay_name = aname, 
                                  compute_agf = TRUE, k_geom = k_geom)
spe.sub <- Banksy::runBanksyPCA(spe.sub, use_agf = TRUE, lambda = lambda, npcs = npcs)

reducedDim(spe.sub, "PCA") <- reducedDim(spe.sub, "PCA_M1_lam0.8")
reducedDim(spe.sub, "PCA_M1_lam0.8") <- NULL
spe.sub <- RunHarmony(spe.sub, "sample_id", reduction.save = "HARMONY_M1_lam0.8")

SpatialExperiment::imgData(spe.sub) <- NULL
spe.sub <- Banksy::clusterBanksy(spe.sub, use_agf = TRUE,
                                  lambda = lambda, npcs = npcs, 
                                  resolution = res,
                                  dimred = "HARMONY_M1_lam0.8")

out_file <- file.path(processed_dir,
                       sprintf("Banksy_integrated_lambda_0.8_res%s.rds", res))
saveRDS(spe.sub, out_file)
cat("Saved:", out_file, "\n")