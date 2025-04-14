suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("scran"))

load(here("processed-data", "03_qc_metrics", "spe_local_outliers.Rdata"), verbose = TRUE)
spe

# discard out of tissue spots
spe <- spe[,colData(spe)$in_tissue]

# Discard bad spots
spe <- spe[,!colData(spe)$local_outliers]

# drop mito genes
spe <- spe[!grepl("^MT-", rownames(spe)),]

# drop bad brnum Br9469
spe <- spe[,colData(spe)$brnum != "Br9469"]

# calculate library size factors
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00051  0.35785  0.76570  1.00000  1.31600 11.07102

# I'm having issues of negative size factos here. Not sure the best way to handle.
# Voyager vingette is doing the below. I'll try it. https://pachterlab.github.io/voyager/articles/vig2_visium.html

dim(spe)
#[1]  26758 134288

spe<- spe[, sizeFactors(spe) > 0]
dim(spe)
#[1]  26758 134287


spe <- logNormCounts(spe)

# This actually worked. It looks like only one spot was causing the issue

save(spe, file = here::here("processed-data", "04_normalization", "spe_norm.Rdata"))
