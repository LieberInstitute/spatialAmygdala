library("here")
library("Seurat")
library("GeneNMF")
library("RColorBrewer")

# set random seed for reproducibility
set.seed(382)


# load data
rat.amy <- readRDS(here("processed-data","snRNAseq", "GSE212415_seurat.rds"))
rat.amy



# ======= run GeneNMF on samples separately =======
seu.list <- SplitObject(rat.amy, split.by = "sample")

geneNMF.programs <- multiNMF(seu.list, assay="SCT", slot="data", k=10:14, nfeatures = 1000)


#  extract metaprograms
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=15,
                                        weight.explained = 0.7,
                                        max.genes=100)



pdf(here("plots","snRNAseq","Zhou_cocaine","GeneNMF_Zhou_metaprograms_k10-14_k15.pdf"), width=8, height=6)
plotMetaPrograms(geneNMF.metaprograms)
dev.off()