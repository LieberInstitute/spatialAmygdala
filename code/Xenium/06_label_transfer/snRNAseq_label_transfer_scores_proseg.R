library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("here")
library("dplyr")
library("patchwork")
library("SingleR")

# save directories
processed_dir <- here("processed-data", "Xenium", "06_label_transfer")
plot_dir <- here("plots", "Xenium", "06_label_transfer")


# ===== Load predictions ======
pred.broad <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_broad_celltype.csv"))
pred.fine <- read.csv(here("processed-data", "Xenium", "06_label_transfer", "xenium_5um_pred_fine_celltype.csv"))


head(pred.fine)
#   scores.aBA_ESR1 scores.aBA_GRIK3 scores.ADARB2 scores.Astrocyte
# 1      0.09872651       0.10300679    0.11171737       0.06123702
# 2      0.04286984       0.02774025    0.07434693       0.08539900
# 3      0.13694261       0.14609139    0.14283185       0.10395617
# 4      0.08038352       0.06813568    0.06896444       0.07246244
# 5      0.08206626       0.08818053    0.10113972       0.01423641
# 6      0.11387633       0.13164363    0.11184681       0.08730327
#   scores.BA_MOXD6 scores.BA_MYRIP scores.CALCR_PENK scores.CCK_CNR1
# 1      0.12858343      0.10403016      0.0888968600      0.08034328
# 2      0.06673149      0.05414914     -0.0004774821      0.03829230
# 3      0.13896253      0.14456125      0.1094786875      0.14684027
# 4      0.07148551      0.06423364      0.0405339298      0.09668474
# 5      0.11169287      0.10161190      0.0420377706      0.03856562
# 6      0.12377335      0.13674706      0.0844817234      0.04679598
#   scores.DLK1_ZFHX3 scores.Endothelial scores.Ependymal scores.LA_RORB
# 1        0.08607666         0.18268777       0.08311726     0.09920418
# 2        0.03716161         0.22012139       0.11746365     0.04550745
# 3        0.13533000         0.16996673       0.18173545     0.12420277
# 4        0.04506322         0.24636947       0.13898443     0.08198761
# 5        0.04639713         0.07830833      -0.01492895     0.09792076
# 6        0.12409363         0.13920998       0.10779319     0.11089993
#   scores.LA_ZBTB20 scores.LAMP5_EGFR scores.LAMP5_NOS1 scores.Microglia
# 1       0.11027825       0.064901269       0.044747998      0.014541397
# 2       0.05372013      -0.028612130       0.007335239      0.078843455
# 3       0.14758969       0.085630746       0.092359982      0.118661158
# 4       0.06845703       0.057156127       0.062054689      0.102383662
# 5       0.09833446      -0.003455724       0.001798761     -0.003018766
# 6       0.12446540       0.118930898       0.081430269      0.006411339
#   scores.Oligodendrocyte scores.PENK_DRD2 scores.PVALB_MYO5B scores.PVALB_ST18
# 1            0.030552026       0.12054726         0.04981847       0.054979239
# 2            0.033597508       0.03082014         0.04677568      -0.002371666
# 3            0.107486825       0.12112694         0.07507526       0.101462498
# 4            0.058866257       0.03552298         0.02944326      -0.008028810
# 5           -0.004764702       0.03346073         0.10034628       0.033026799
# 6           -0.021740746       0.11909695         0.12123434       0.050204392
#   scores.PVALB_UNC5B scores.RELN_IL1RAPL2 scores.SATB2_RXFP1
# 1         0.04486442           0.12615417         0.09031730
# 2         0.03023356           0.03429345         0.02805369
# 3         0.08351388           0.12094707         0.13985595
# 4         0.07381069           0.13591516         0.06418643
# 5         0.04385765           0.03328632         0.07753421
# 6         0.08271664           0.09021355         0.12668745
#   scores.SATB2_SLC17A8 scores.SATB2_SULF1 scores.SLC17A6_CARTPT
# 1           0.09900167         0.11400777            0.12946498
# 2           0.06006262         0.04049944            0.05076724
# 3           0.14995333         0.14577407            0.13654045
# 4           0.06787973         0.07608455            0.04104752
# 5           0.09040375         0.08907068            0.10401366
# 6           0.13112850         0.13804452            0.12750950
#   scores.SLC17A6_SIM1 scores.SST_NOS1 scores.SST_NXPH2 scores.SST_TAC1
# 1          0.10539593      0.06406642       0.05095737      0.07139937
# 2          0.04034572      0.06844333       0.02620813      0.03619310
# 3          0.13977375      0.12067478       0.09084698      0.11241635
# 4          0.07740243      0.03426130       0.02078487      0.03563589
# 5          0.08229166      0.04552871       0.06537589      0.04239977
# 6          0.10829146      0.11753641       0.07417622      0.10691057
#   scores.SST_TMEM132C scores.TAC1_PPP1R1B scores.TSHZ1_CPNE4 scores.TSHZ1_PRKG1
# 1          0.08648779          0.08219371         0.05033481        0.079940893
# 2          0.01951014          0.02874554         0.02402173        0.002338244
# 3          0.08453292          0.11360230         0.12566796        0.110590617
# 4         -0.00931028          0.03732195         0.02999443        0.021278575
# 5          0.07452219          0.01209126         0.02326943        0.003431270
# 6          0.03812154          0.09029372         0.08259574        0.072409869
#   scores.VIP_NRXN1 scores.VIP_SEMA5A      labels   delta.next pruned.labels
# 1       0.07510135        0.06336311 Endothelial 0.0532227937   Endothelial
# 2       0.03259651        0.04035029 Endothelial 0.1026577374   Endothelial
# 3       0.10002658        0.12005465   aBA_GRIK3 0.0001602253     aBA_GRIK3
# 4       0.05114876        0.06544312 Endothelial 0.1073850398   Endothelial
# 5       0.02485018        0.01589515 PVALB_MYO5B 0.0005998836   PVALB_MYO5B
# 6       0.07744946        0.05195681  LAMP5_EGFR 0.0028301369    LAMP5_EGFR

# ======= Scores ggplot ======

library(tidyverse)

# Convert your prediction dataframe to tibble
df <- as_tibble(pred.fine)

# Extract just the score that matches the predicted label
df_long <- df %>%
  mutate(cell_id = row_number()) %>% 
  pivot_longer(cols = starts_with("scores."), names_to = "celltype", values_to = "score") %>%
  mutate(celltype = str_remove(celltype, "scores\\.")) %>%
  filter(celltype == labels)

# Now make the boxplot
pdf(file = here(plot_dir, "xenium_5um_label_transfer_scores_fine.pdf"), width = 8, height = 6)
ggplot(df_long, aes(x = celltype, y = score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label")
dev.off()




# Convert your prediction dataframe to tibble
df <- as_tibble(pred.broad)

# Extract just the score that matches the predicted label

df_long <- df %>%
  mutate(cell_id = row_number()) %>%
  pivot_longer(cols = starts_with("scores."), names_to = "celltype", values_to = "score") %>%
  mutate(celltype_clean = str_remove(celltype, "scores\\."),           # remove prefix
         celltype_clean = str_replace_all(celltype_clean, "\\.", "-"), # replace dots with dashes
         label_clean = str_replace_all(labels, "\\.", "-")) %>%        # also clean labels
  filter(celltype_clean == label_clean)


# Now make the boxplot
pdf(file = here(plot_dir, "xenium_5um_label_transfer_scores_broad.pdf"), width = 8, height = 6)
ggplot(df_long, aes(x = celltype, y = score)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted Cell Type", y = "Score for Assigned Label")
dev.off()

