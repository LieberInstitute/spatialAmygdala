suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("scran"))

load(here("processed-data", "03_qc_metrics", "spe_discarded.Rdata"), verbose = TRUE)
spe

# Discard bad spots
spe <- spe[,!colData(spe)$discard]

# calculate library size factors
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00051  0.35785  0.76570  1.00000  1.31600 11.07102

pdf(here("plots", "04_normalization", "hist_size_factors.pdf"), width = 21, height = 10)
hist(sizeFactors(spe), breaks = 20)
dev.off()

spe <- logNormCounts(spe)

save(spe, file = here::here("processed-data", "04_normalization", "spe_norm.Rdata"))