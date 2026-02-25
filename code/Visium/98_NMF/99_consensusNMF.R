library(SpatialExperiment)
library(RcppML)
library(here)
library(Matrix)
library(Seurat)
library(GeneNMF)

plot_dir <- here("plots", "Visium", "98_NMF", "consensusNMF")
processed_dir <- here("processed-data", "snRNAseq")


load(here(processed_dir, "sce_FINAL_human.rda"))
rda.human
sce.libd <- rda.human
sce.libd
# class: SingleCellExperiment 
# dim: 36601 15132 
# metadata(0):
# assays(2): counts logcounts
# rownames(36601): A1BG A1BG-AS1 ... ZYX ZZEF1
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(15132): cell_2 cell_4 ... cell_16967 cell_16969
# colData names(33): orig.ident nCount_originalexp ... fine_celltype
#   ident
# reducedDimNames(1): umap
# mainExpName: NULL
# altExpNames(0):

sce.biccn <- readRDS(here(processed_dir, "biccnAMY_sce.rds"))
sce.biccn
# class: SingleCellExperiment 
# dim: 59357 187167 
# metadata(0):
# assays(1): logcounts
# rownames(59357): ENSG00000259882 ENSG00000273106 ... ENSG00000161405
#   ENSG00000185905
# rowData names(0):
# colnames(187167): 10X355_2:TGTCCACAGGTACCTT 10X355_2:CGTTAGACACTTCTCG
#   ... 10X356_7:AAAGGGCCAGCTACCG 10X205_5:AGGGCCTGTTGGAGGT
# colData names(34): roi organism_ontology_term_id ... ident dataset
# reducedDimNames(2): UMAP TSNE
# mainExpName: RNA
# altExpNames(0):

load(here(processed_dir, "yu_sce_gtf.rda"))
sce.yu <- sce.amy
sce.yu
# class: SingleCellExperiment 
# dim: 31908 91699 
# metadata(0):
# assays(2): counts logcounts
# rownames(31908): MIR1302-2HG FAM138A ... AC240274.1 AC213203.1
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(91699): AAACCCAAGTCGCCAC-1 AAACCCAGTCGAGATG-1 ...
#   TTTGTTGTCACTTCTA-14 TTTGTTGTCCAAATGC-14
# colData names(14): barcode library_id ... orig_anno ident
# reducedDimNames(2): PCA UMAP
# mainExpName: RNA
# altExpNames(0):

# ========== Get common genes ==========

# convert to Ensembl IDs
rownames(sce.libd) <- rowData(sce.libd)$gene_id
rownames(sce.yu) <- rowData(sce.yu)$gene_id

genes.libd <- rownames(sce.libd)
genes.biccn <- rownames(sce.biccn)
genes.yu <- rownames(sce.yu)

# get common
common_genes <- Reduce(intersect, list(genes.libd, genes.biccn, genes.yu))
length(common_genes)
# [1] 31388

# subset to common genes
sce.libd.common <- sce.libd[common_genes, ]
sce.biccn.common <- sce.biccn[common_genes, ]
sce.yu.common <- sce.yu[common_genes, ]

# conver to seurat
obj.libd <- as.Seurat(sce.libd.common, counts="logcounts", data = "logcounts")
obj.biccn <- as.Seurat(sce.biccn.common, counts="logcounts", data = "logcounts")
obj.yu <- as.Seurat(sce.yu.common, counts="logcounts", data = "logcounts")

# make new RNA assay from originalexp in libd
obj.libd[["RNA"]] <- CreateAssayObject(data = GetAssayData(obj.libd, slot = "data"))

# remove objects to save memory
rm(sce.libd)
rm(sce.biccn)
rm(sce.yu)
rm(sce.libd.common)
rm(sce.biccn.common)
rm(sce.yu.common)


# ========= Consensus NMF =========

# make seurat list
seurat_list <- list(LIBD = obj.libd, BICCN = obj.biccn, YU = obj.yu)

set.seed(1234)
geneNMF.programs <- multiNMF(seurat_list, assay="RNA", k=30:40, min.exp = 0.05)


pdf(here(plot_dir, "consensusNMF_gene_programs.pdf"), width=10, height=8)
ph <- plotMetaPrograms(geneNMF.metaprograms,
                       similarity.cutoff = c(0.1,1))
pheatmap
print(ph)
dev.off()

# save results
save(geneNMF.programs, file=here(processed_dir, "Visium","98_NMF", "consensusNMF_gene_programs_k30to40.rda"))