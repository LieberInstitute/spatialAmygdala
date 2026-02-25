
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)
library(scran)

# want to use logcounts for NMF
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy

sce <- as.SingleCellExperiment(rat.amy)

# subset to HVGs
# dec <- modelGeneVar(sce)
# chosen <- getTopHVGs(dec, prop=0.1)

# sce <- sce[chosen, ]

cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(10, 20, 30, 40, 50, 75, 100, 125),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

png(here("plots","snRNAseq","Zhou_cocaine","nmf_CV_results_lower_k.png"),
    height=4,width=8, unit="in",res=300)
plot(cvnmf)
dev.off()

saveRDS(cvnmf, file = here("processed-data", "snRNAseq", "Zhou_cocaine", "nmf_cv_results_lower_k.RDS"))

