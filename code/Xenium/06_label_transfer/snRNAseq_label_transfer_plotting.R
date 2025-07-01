library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")
library(Voyager)

# save directories
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
plot_dir <- here("plots", "Xenium", "06_label_transfer")


# load xenium data
load(here("processed-data","Xenium", "03_quality_control", "spe_normcounts.Rdata"))
spe

colnames(spe) <- make.unique(colnames(spe), sep = "-")
rownames(spatialCoords(spe)) <- colnames(spe)

sfe <- toSpatialFeatureExperiment(spe)
sfe
# class: SpatialFeatureExperiment 
# dim: 541 1018069 
# metadata(0):
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ENSG00000069431 ENSG00000151388 ...
#   DeprecatedCodeword_0344 DeprecatedCodeword_0373
# rowData names(3): ID Symbol Type
# colnames(1018069): aaaadgkh-1 aaaadlfe-1 ... oihobmbk-1 oihoeehh-1
# colData names(19): cell_id transcript_counts ... cell_area.sf
#   nucleus_area.sf
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

# unit:
# Geometries:
# colGeometries: centroids (POINT) 

# Graphs:
# 20240425__170523__0022862: 
# 20240425__170523__0023004: 
# output-XETG00117__0023153__Region_1__20240329__171204: 
# output-XETG00117__0023154__Region_1__20240329__171203: 

# load label predictions
pred.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "pred_fine_celltype.csv"))

head(pred.broad)
#   scores.Excitatory scores.Inhibitory scores.Non.neuronal       labels
# 1        0.19319302         0.1835533         -0.00430950   Excitatory
# 2        0.06472557         0.0336494          0.15005161 Non-neuronal
# 3        0.36894944         0.1903083          0.14287377   Excitatory
# 4        0.21104344         0.1361702          0.24034109 Non-neuronal
# 5        0.12929172         0.1479399         -0.04223317   Excitatory
# 6        0.29683978         0.1528640          0.20460371   Excitatory
#     delta.next pruned.labels
# 1 2.747675e-01    Excitatory
# 2 8.532605e-02  Non-neuronal
# 3 1.786411e-01    Excitatory
# 4 3.273884e-01  Non-neuronal
# 5 1.110223e-16    Excitatory
# 6 9.223607e-02    Excitatory

head(pred.fine)
#   scores.aBA_ESR1 scores.aBA_GRIK3 scores.ADARB2 scores.Astrocyte
# 1      0.12284806       0.13790068    0.14086617       0.07244408
# 2      0.07269898       0.06013192    0.10152535       0.07743790
# 3      0.13694261       0.14609139    0.14283185       0.10395617
# 4      0.07988761       0.06755231    0.06848542       0.07235379
# 5      0.11494250       0.12203284    0.13516474       0.01514703
# 6      0.12505878       0.14719936    0.12080988       0.09214318
#   scores.BA_MOXD6 scores.BA_MYRIP scores.CALCR_PENK scores.CCK_CNR1
# 1      0.15232349      0.13602711        0.08900065      0.07911925
# 2      0.10343338      0.09154458        0.03211850      0.06784961
# 3      0.13896253      0.14456125        0.10947869      0.14684027
# 4      0.07092909      0.06404432        0.04009634      0.09602650
# 5      0.14064730      0.13404386        0.07204170      0.05104494
# 6      0.12425224      0.13175031        0.07298671      0.03391356
#   scores.DLK1_ZFHX3 scores.Endothelial scores.Ependymal scores.LA_RORB
# 1        0.09761535         0.17050637       0.09386726     0.13203136
# 2        0.06362222         0.20651921       0.09899986     0.07462727
# 3        0.13533000         0.16996673       0.18173545     0.12420277
# 4        0.04450273         0.24743241       0.13928305     0.08122008
# 5        0.06428753         0.07814371      -0.01452194     0.12442988
# 6        0.10827407         0.11935111       0.08485824     0.11862718
#   scores.LA_ZBTB20 scores.LAMP5_EGFR scores.LAMP5_NOS1 scores.Microglia
# 1       0.13755674       0.081538000        0.03380823      0.008917900
# 2       0.08846387      -0.005147238        0.02973724      0.075171067
# 3       0.14758969       0.085630746        0.09235998      0.118661158
# 4       0.06838343       0.057151995        0.06274706      0.101652972
# 5       0.12438038       0.016384171        0.01576153      0.009926586
# 6       0.11727491       0.130513587        0.09798839     -0.001005825
#   scores.Oligodendrocyte scores.PENK_DRD2 scores.PVALB_MYO5B scores.PVALB_ST18
# 1             0.01986781       0.13642741         0.04809575        0.05352965
# 2             0.03287174       0.03617164         0.07546828        0.03530745
# 3             0.10748682       0.12112694         0.07507526        0.10146250
# 4             0.05875594       0.03495307         0.02905279       -0.00740078
# 5             0.00574073       0.05822687         0.10539449        0.05200554
# 6            -0.02938176       0.10255604         0.10170997        0.04115665
#   scores.PVALB_UNC5B scores.RELN_IL1RAPL2 scores.SATB2_RXFP1
# 1         0.03670632           0.13777276         0.11202498
# 2         0.05971944           0.07389180         0.04426244
# 3         0.08351388           0.12094707         0.13985595
# 4         0.07308364           0.13698946         0.06334181
# 5         0.05443218           0.05144707         0.10328338
# 6         0.06345443           0.09547837         0.12985419
#   scores.SATB2_SLC17A8 scores.SATB2_SULF1 scores.SLC17A6_CARTPT
# 1           0.13055180         0.14137710            0.16385014
# 2           0.08508883         0.07699026            0.05608435
# 3           0.14995333         0.14577407            0.13654045
# 4           0.06766420         0.07673315            0.04010701
# 5           0.12061380         0.12886480            0.13513920
# 6           0.14447507         0.13834585            0.14927645
#   scores.SLC17A6_SIM1 scores.SST_NOS1 scores.SST_NXPH2 scores.SST_TAC1
# 1           0.1334108      0.07188033       0.06373150      0.07964873
# 2           0.0822694      0.08577688       0.05271645      0.06285746
# 3           0.1397738      0.12067478       0.09084698      0.11241635
# 4           0.0762196      0.03377172       0.02060851      0.03530189
# 5           0.1216372      0.05666707       0.07772450      0.05189744
# 6           0.1192144      0.10230057       0.07345659      0.09669506
#   scores.SST_TMEM132C scores.TAC1_PPP1R1B scores.TSHZ1_CPNE4 scores.TSHZ1_PRKG1
# 1         0.096136054          0.11959581         0.06374123         0.09875292
# 2         0.034982673          0.03830943         0.04261463         0.02903255
# 3         0.084532915          0.11360230         0.12566796         0.11059062
# 4        -0.009548501          0.03789039         0.02961941         0.02078526
# 5         0.089573192          0.04097391         0.03951455         0.02443286
# 6         0.027326386          0.07574951         0.08235846         0.06583161
#   scores.VIP_NRXN1 scores.VIP_SEMA5A      labels   delta.next pruned.labels
# 1       0.08332581        0.07000219 Endothelial 0.1069711284   Endothelial
# 2       0.05547706        0.06325362 Endothelial 0.1030858335   Endothelial
# 3       0.10002658        0.12005465   aBA_GRIK3 0.0001602253     aBA_GRIK3
# 4       0.05089195        0.06486540 Endothelial 0.1081493607   Endothelial
# 5       0.04719210        0.03414939      ADARB2 0.0049250341        ADARB2
# 6       0.07114378        0.04014417 Endothelial 0.0087818227   Endothelial

# ======== Adding labels back in ========
# NOTE only doing this to match the coldata, need to refun label transfer where I keep the genes
# while dropping the correct sizeFactor < 0 cells.

# keep only rowData Type == Gene Expression
sfe <- sfe[rowData(sfe)$Type == "Gene Expression", ]

# re-normalize the snRNAseq data
sfe <- computeLibraryFactors(sfe)
sfe <- sfe[, sizeFactors(sfe) > 0]
sfe <- logNormCounts(sfe)

# add labels
sfe$pred_broad_celltype <- pred.broad$pruned.labels
sfe$pred_fine_celltype <- pred.fine$pruned.labels



# ==== Plotting the label transfer results ====


# UMAP
png(file=here(plot_dir, "UMAP_broad_labels.png"), width=8, height=8, units="in", res=300)
plotSpots(spe, annotate="pred_broad_celltype", in_tissue=NULL, point_size=0.03, sample_id="sample_id")
dev.off()

png(file=here(plot_dir, "UMAP_fine_labels.png"), width=8, height=8, units="in", res=300)
plotSpots(spe, annotate="pred_fine_celltype", in_tissue=NULL, point_size=0.03, sample_id="sample_id")
dev.off()


rownames(spe) <- rowData(spe)$Symbol

png(file=here(plot_dir, "UMAP_TSHZ1.png"), width=8, height=8, units="in", res=300)
plotSpots(spe, annotate="COL25A1", in_tissue=NULL, point_size=0.01, sample_id="sample_id")
dev.off()
