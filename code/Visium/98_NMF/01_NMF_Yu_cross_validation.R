# Load libraries
library(RcppML)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(singlet)

# Load data
load(file = here("processed-data", "snRNAseq", "yu_sce_gtf.rda"))
sce <- sce.amy
# want to use logcounts for NMF
sce <- logNormCounts(sce)

cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=c(5,10,20, 30, 40, 50, 75, 100,125),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

saveRDS(cvnmf, file = here("processed-data", "Visium", "98_NMF", "nmf_cv_results_lower_k.RDS"))

png(here("plots","Visium","98_NMF","NMF_Yu","nmf_cv_results_lower_k.png"),
    height=4,width=8, unit="in",res=300)
plot(cvnmf)
dev.off()