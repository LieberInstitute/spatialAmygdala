## SCRIPT FOR READING LARGE hdf5 DATA SET INTO R

## Load necessary libraries
library(rhdf5)      # For reading hdf5; https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
library(HDF5Array)  # Alternatively option for reading the data matrix with less memory: https://rdrr.io/github/Bioconductor/HDF5Array/
library(data.table) # For fast reading of csv files
library(SingleCellExperiment)
library(here)

## Read in the metadata using fread (fast!)
metadata <- fread(here("processed-data/snRNAseq/abc_atlas/metadata/abc_cell_metadata_extended.csv"))
metadata <- as.data.frame(metadata)
colnames(metadata)
#  [1] "index"                      "cell_label"                
#  [3] "cell_barcode"               "barcoded_cell_sample_label"
#  [5] "library_label"              "feature_matrix_label"      
#  [7] "entity"                     "brain_section_label"       
#  [9] "library_method"             "donor_label"               
# [11] "donor_sex"                  "dataset_label"             
# [13] "x"                          "y"                         
# [15] "cluster_alias"              "region_of_interest_label"  
# [17] "anatomical_division_label"  "abc_sample_id"             
# [19] "subcluster"                 "cluster"                   
# [21] "supercluster"               "neurotransmitter"          
# [23] "subcluster_color"           "cluster_color"             
# [25] "supercluster_color"         "neurotransmitter_color"    
# [27] "region_of_interest_color"  

unique(metadata$anatomical_division_label)
#  [1] "Myelencephalon"     "Pons"               "Cerebellum"        
#  [4] "Midbrain"           "Thalamus"           "Hypothalamus"      
#  [7] "Spinal cord"        "Basal forebrain"    "Hippocampus"       
# [10] "Basal nuclei"       "Amygdaloid complex" "Cerebral cortex"   
# [13] "Extended amygdala"  "Claustrum"  

# get amygdaloid complex samples
samples <- metadata$cell_label[metadata$anatomical_division_label == "Amygdaloid complex"]
length(samples)

unique(metadata$supercluster)
# [1] "Upper rhombic lip"                   "Splatter"                           
#  [3] "Lower rhombic lip"                   "Mammillary body"                    
#  [5] "Thalamic excitatory"                 "Amygdala excitatory"                
#  [7] "Medium spiny neuron"                 "Eccentric medium spiny neuron"      
#  [9] "Miscellaneous"                       "Cerebellar inhibitory"              
# [11] "Midbrain-derived inhibitory"         "CGE interneuron"                    
# [13] "LAMP5-LHX6 and Chandelier"           "MGE interneuron"                    
# [15] "Deep-layer near-projecting"          "Deep-layer corticothalamic and 6b"  
# [17] "Hippocampal CA1-3"                   "Upper-layer intratelencephalic"     
# [19] "Deep-layer intratelencephalic"       "Hippocampal dentate gyrus"          
# [21] "Hippocampal CA4"                     "Oligodendrocyte"                    
# [23] "Committed oligodendrocyte precursor" "Astrocyte"                          
# [25] "Bergmann glia"                       "Oligodendrocyte precursor"          
# [27] "Ependymal"                           "Choroid plexus"                     
# [29] "Fibroblast"                          "Vascular"                           
# [31] "Microglia"

# get Eccentric medium spiny neuron samples
samples <- metadata$cell_label[metadata$supercluster == "Eccentric medium spiny neuron"]
length(samples) 

# get Amygdala excitatory samples
samples <- metadata$cell_label[metadata$supercluster == "Amygdala excitatory"]
length(samples)

# # Note that the order of metadata and counts and the number of cells are DIFFERENT 
# #   (there are 107 cells with no metadata). This will be critical later!

# ## Subsample your data
# ## -- If you want to analyze these data in R and especially using Seurat, you will
# ## --   have a better chance of success if you only work on parts of the data at a
# ## --   time.  We suggest, using information in the metadata (e.g., cell type or
# ## --   brain region columns) to subset.   
# ## -- For this example code, I'm selecting 1000 random cells checking overlap with samples
# set.seed(42)
# use_samples  <- intersect(sample(rownames(metadata),1000),samples)
# read_samples <- sort(match(use_samples,samples))


# ## Read in the count matrix in one of three ways

# ## Strategy #1: Read in the entire data set using h5read. The data set is ~145GB total,
# ##   and as far as I can tell there is no way to convert it into a sparse matrix, transpose
# ##   it, or read it into Seurat without subsetting (at least with the 1TB of RAM my
# ##   computer has.  It may be possible with python or if you have more memory.

# system.time({
#   counts <- h5read("expression_matrix.hdf5", "/data/counts")
# }) 
# #    user  system elapsed 
# # 398.307 146.285 545.006 

## Strategy #2: Read in only a relevant subset of data using h5read. This is the 
##   method that probably works best in most situations.

contents <- h5ls(here("processed-data/snRNAseq/abc_atlas/expression_matrices/WHB-10Xv3/20240330/WHB-10Xv3-Neurons-raw.h5ad"))

system.time({
  counts1 <- h5read(here("processed-data/snRNAseq/abc_atlas/expression_matrices/WHB-10Xv3/20240330/WHB-10Xv3-Neurons-raw.h5ad"),
  name="/X/data",
  index = list(samples,NULL))
  counts1 <- t(counts1)
  subcounts1 <- as(counts1, "dgCMatrix")
})  
#    user  system elapsed 
# 270.724  36.517 307.828    

# ## Strategy #3: Work with DelayedArray format before reading in (via HDF5Array library).  
# ##   This has the advantage of being able to manipulate the data on disk, and can read
# ##   the file directly into sparse matrix format and therefore requires less memory.
# ##   However, this method is much slower as the number of non-sequential cells increases
# ##   and will crash if you have too many cells (I'm not sure what the cutoff is, but I 
# ##   couldn't read the whole data set this way.)

# system.time({
#   counts2 <- HDF5Array("expression_matrix.hdf5", "/data/counts")  # Note that COLUMNS are genes
#   counts2 <- t(counts2)
#   subcounts2 <- as(counts2[,read_samples], "dgCMatrix")
# })
# #    user  system elapsed 
# # 310.024  63.176 375.114 


# ## Assess variable sizes (optional)
# sapply(ls(),function(x){object.size(get(x))})[c("counts","counts1","subcounts1","counts2","subcounts2")]
# #        counts      counts1   subcounts1      counts2   subcounts2
# #  145243576056    124212216     49413768         3368     49413768 


## Add gene and sample names to the data matrix
## -- Note: We'll works with the recommended strategy (#2) for the rest of the script
rownames(subcounts1) <- as.character(genes)
colnames(subcounts1) <- as.character(samples) [read_samples]


## Read the subsetted data and metadata into Seurat
## -- Note that this **WILL NOT WORK** for the full 1 Million+ cell data set!
## -- I would strongly encourage other methods for analysis when dealing with more than ~100,000 cells at a time
## -- Also recall that the order of data and meta do not match, so reorder here
seu <- CreateSeuratObject(counts=subcounts1)          # Put in a Seurat object
met <- as.data.frame(metadata[colnames(subcounts1),]) # Format the metadata
seu <- AddMetaData(seu,met)                           # Add the metadata