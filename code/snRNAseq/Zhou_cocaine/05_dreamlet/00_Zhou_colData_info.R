library(here)
library(Seurat)
library(SingleCellExperiment)


# load
rat.amy <- readRDS(here("processed-data","snRNAseq","zhou_with_yu_labels.rds"))
sce <- as.SingleCellExperiment(rat.amy)
sce
# class: SingleCellExperiment 
# dim: 17297 163003 
# metadata(0):
# assays(2): counts logcounts
# rownames(17297): AABR07000156.1 Lrp11 ... AABR07043200.1 Pomp
# rowData names(0):
# colnames(163003): AAACCCAAGAAACCCG-1_1 AAACCCACAAAGCACG-1_1 ...
#   TTTGTTGTCTTCGTAT-1_19 TTTGTTGTCTTCTGGC-1_19
# colData names(86): orig.ident nCount_RNA ...
#   yu_space_prediction.score.max ident
# reducedDimNames(0):
# mainExpName: RNA
# altExpNames(2): SCT integrated

colnames(colData(sce))
# [1] "orig.ident"                                         
#  [2] "nCount_RNA"                                         
#  [3] "nFeature_RNA"                                       
#  [4] "sample"                                             
#  [5] "treatment"                                          
#  [6] "addiction.index"                                    
#  [7] "label"                                              
#  [8] "percent.mt"                                         
#  [9] "nCount_SCT"                                         
# [10] "nFeature_SCT"                                       
# [11] "rfid"                                               
# [12] "integrated_snn_res.0.8"                             
# [13] "seurat_clusters"                                    
# [14] "cocaine.low"                                        
# [15] "cocaine.high"                                       
# [16] "X933000320047328"                                   
# [17] "X933000120138592"                                   
# [18] "X933000120138586"                                   
# [19] "X933000320046084"                                   
# [20] "X933000320046077"                                   
# [21] "X933000120138609"                                   
# [22] "X933000320186802"                                   
# [23] "X933000320047225"                                   
# [24] "X933000320046609"                                   
# [25] "X933000320047001"                                   
# [26] "X933000320047132"                                   
# [27] "X933000320186801"                                   
# [28] "X933000320046621"                                   
# [29] "A_933000320046625_JB_257"                           
# [30] "X933000320047104"                                   
# [31] "X933000320045674"                                   
# [32] "Rat_Opioid_HS_1"                                    
# [33] "Rat_Opioid_HS_2"                                    
# [34] "Rat_Amygdala_787A_all_seq"                          
# [35] "batch"                                              
# [36] "yu_cell_predicted.id"                               
# [37] "yu_cell_prediction.score.Micro.Ctss"                
# [38] "yu_cell_prediction.score.Adamts5.St18"              
# [39] "yu_cell_prediction.score.Foxp2.Nxph2"               
# [40] "yu_cell_prediction.score.Isl1"                      
# [41] "yu_cell_prediction.score.Tshz1.Lamp5"               
# [42] "yu_cell_prediction.score.Matn2.Slc24a2"             
# [43] "yu_cell_prediction.score.Peri.Art3"                 
# [44] "yu_cell_prediction.score.Htr3a.Vip"                 
# [45] "yu_cell_prediction.score.Htr3a.Ppp1r1c"             
# [46] "yu_cell_prediction.score.Vwa5b1.Sim1"               
# [47] "yu_cell_prediction.score.Tshz1.Drd1"                
# [48] "yu_cell_prediction.score.Adamts5.Kansl1l"           
# [49] "yu_cell_prediction.score.Astro.Aqp4"                
# [50] "yu_cell_prediction.score.Lamp5.Ndnf"                
# [51] "yu_cell_prediction.score.Vwa5b1.Zbtb7c"             
# [52] "yu_cell_prediction.score.Crim1"                     
# [53] "yu_cell_prediction.score.Endo.Flt1"                 
# [54] "yu_cell_prediction.score.Vgll3"                     
# [55] "yu_cell_prediction.score.Hgf.Grp"                   
# [56] "yu_cell_prediction.score.Calcr"                     
# [57] "yu_cell_prediction.score.Abi3bp"                    
# [58] "yu_cell_prediction.score.Sst.Gldn"                  
# [59] "yu_cell_prediction.score.Htr3a.Col14a1"             
# [60] "yu_cell_prediction.score.Npsr1.Rspo2"               
# [61] "yu_cell_prediction.score.Sst.Drd1"                  
# [62] "yu_cell_prediction.score.Oligo.Mog"                 
# [63] "yu_cell_prediction.score.Rorb"                      
# [64] "yu_cell_prediction.score.Satb2.Ebf2"                
# [65] "yu_cell_prediction.score.Lhx8"                      
# [66] "yu_cell_prediction.score.OPC.Pdgfra"                
# [67] "yu_cell_prediction.score.Hgf"                       
# [68] "yu_cell_prediction.score.Prkcd"                     
# [69] "yu_cell_prediction.score.Pvalb"                     
# [70] "yu_cell_prediction.score.Tfap2c"                    
# [71] "yu_cell_prediction.score.VLMC.Lum"                  
# [72] "yu_cell_prediction.score.Drd2"                      
# [73] "yu_cell_prediction.score.Lamp5.Id2"                 
# [74] "yu_cell_prediction.score.Sst.Nts.Tac2.Crh"          
# [75] "yu_cell_prediction.score.Strip2"                    
# [76] "yu_cell_prediction.score.max"                       
# [77] "yu_space_predicted.id"                              
# [78] "yu_space_prediction.score.non.neuron"               
# [79] "yu_space_prediction.score.Cortical.interneuron.like"
# [80] "yu_space_prediction.score.COA.MEA"                  
# [81] "yu_space_prediction.score.CEA"                      
# [82] "yu_space_prediction.score.IA"                       
# [83] "yu_space_prediction.score.BLA"                      
# [84] "yu_space_prediction.score.NLOT"                     
# [85] "yu_space_prediction.score.max"                      
# [86] "ident"        

# specifically: how many unique rats, and what column identifies them?
unique(sce$ident) 
#  [1] InhNeuron        Astrocytes       ExNeuron         Nos1+           
#  [5] Microglia        Oligodendrocytes Chat+            Cck+/Vip+       
#  [9] OPC              Sst+             Reln+            Endothelial     
# [13] Pvalb+          
# 13 Levels: Cck+/Vip+ Astrocytes Oligodendrocytes InhNeuron OPC ... Pvalb+

unique(sce$yu_cell_predicted.id)
#  [1] "Isl1"             "Astro Aqp4"       "Prkcd"            "Vgll3"           
#  [5] "Hgf"              "Vwa5b1 Zbtb7c"    "Micro Ctss"       "Oligo Mog"       
#  [9] "OPC Pdgfra"       "Npsr1 Rspo2"      "Drd2"             "Tshz1 Drd1"      
# [13] "Lhx8"             "Foxp2 Nxph2"      "Htr3a Col14a1"    "Matn2 Slc24a2"   
# [17] "Satb2 Ebf2"       "Adamts5 St18"     "Hgf Grp"          "Sst Nts/Tac2/Crh"
# [21] "Abi3bp"           "Tfap2c"           "Tshz1 Lamp5"      "Adamts5 Kansl1l" 
# [25] "Endo Flt1"        "Sst Gldn"         "Lamp5 Ndnf"       "Peri Art3"       
# [29] "Htr3a Ppp1r1c"    "Htr3a Vip"        "Sst Drd1"         "Calcr"           
# [33] "Rorb"             "Crim1"            "Lamp5 Id2"        "Pvalb"           
# [37] "VLMC Lum"         "Strip2" 

unique(sce$yu_space_predicted.id)
# [1] "CEA"                       "non-neuron"               
# [3] "BLA"                       "COA/MEA"                  
# [5] "IA"                        "Cortical interneuron-like"

unique(sce$label)
# [1] "cocaine_high" "cocaine_low"  "naive" 

unique(sce$treatment)
# [1] "cocaine" "naive"  

unique(sce$sample)
#  [1] "933000320047328"           "933000120138592"          
#  [3] "933000120138586"           "933000320046084"          
#  [5] "933000320046077"           "933000120138609"          
#  [7] "933000320186802"           "933000320047225"          
#  [9] "933000320046609"           "933000320047001"          
# [11] "933000320047132"           "933000320186801"          
# [13] "933000320046621"           "A_933000320046625_JB_257" 
# [15] "933000320047104"           "933000320045674"          
# [17] "Rat_Opioid_HS_1"           "Rat_Opioid_HS_2"          
# [19] "Rat_Amygdala_787A_all_seq"
