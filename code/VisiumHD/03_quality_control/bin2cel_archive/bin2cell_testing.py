import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os

import bin2cell as b2c

#create directory for stardist input/output files
os.makedirs("stardist", exist_ok=True)

root = os.getcwd()

path = "processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/outs/binned_outputs/square_002um/"
source_image_path = "raw-data/images/vis-hd/40x_HE_AMY_s003.tif"

# read in the data and make names unique
adata = b2c.read_visium(path, source_image_path = source_image_path)
adata.var_names_make_unique()
adata

# filter genes and cells (bins)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)
adata

'''
Over the course of the demo, bin2cell will perform a segmentation of both the
H&E image and a gene expression representation of the data. When performing 
segmentation, the resolution of the input images is controlled via the mpp parameter. 
This stands for microns per pixel and translates to how many micrometers are captured in 
each pixel of the input. For example, if using the array coordinates 
(present as .obs["array_row"] and .obs["array_col"]) as an image, each of the 
pixels would have 2 micrometers in it, so the mpp of that particular representation is 2.

In local testing, using an mpp of 0.5 has worked well with both GEX and H&E segmentation.

Since we're already generating a custom resolution H&E image, b2c.scaled_he_image() 
stores it within the object so it can be used for visualisation. The function crops 
the image to an area around the actual spatial grid present in the object, and the new 
coordinates are captured in .obsm["spatial_cropped"]. The new image can be used for 
plotting by providing basis="spatial_cropped" and img_key="0.5_mpp" to sc.pl.spatial(). 
For segmentation purposes, the image needs to be saved to the drive, and the function 
does so to a user-specified save_path.
'''

mpp = 0.5
b2c.scaled_he_image(adata, mpp=mpp, save_path="plots/VisiumHD/03_bin2cell/stardist/AmyHD_HE.tiff")

'''
Visium HD suffers from variable bin sizing. When printing the chips, the 2um bins can have
about 10% variability in their width/height. Inspecting the total counts per spot reveals a 
characteristic striped appearance, with some rows/columns capturing visibly fewer transcripts 
than others.

To overcome this, b2c.destripe() identifies a user-specified quantile (by default 0.99) of 
total counts for each row, then divides the counts of the spots in that row by that value. 
This procedure is then repeated for the columns. .obs[adjusted_counts_key] 
(by default "n_counts_adjusted") is obtained by multiplying the resulting per-spot factor 
by the global quantile of count totals, and the count matrix is by default rescaled to match 
it.
'''

b2c.destripe(adata)

#define a mask to easily pull out this region of the object in the future
mask = ((adata.obs['array_row'] >= 2000) & 
        (adata.obs['array_row'] <= 4000) & 
        (adata.obs['array_col'] >= 1000) & 
        (adata.obs['array_col'] <= 3000)
       )

bdata = adata[mask]

plt.figure(figsize=(10,7.5))
sc.pl.spatial(bdata, color=[None, "n_counts", "n_counts_adjusted"], img_key="0.5_mpp", basis="spatial_cropped")
plt.savefig(root+'/plots/VisiumHD/03_bin2cell/spatial_cropped_sharp.png', dpi=300)

# plt.figure(figsize=(10,7.5))
# sc.pl.spatial(bdata, color=[None, "n_counts", "n_counts_adjusted"])
# plt.savefig(fname = root+'/plots/spatial_cropped.png', dpi=300)

# segmentation of H&E image
b2c.stardist(image_path="plots/VisiumHD/03_bin2cell/stardist/AmyHD_HE.tiff", 
             labels_npz_path="processed-data/VisiumHD/03_bin2cell/stardist/he.npz", 
             stardist_model="2D_versatile_he", 
             prob_thresh=0.01
            )
# took 13 minutes to segement amygdala H&E image: 
# 12/12 [13:10<00:00, 65.84s/it]

# load resulting cell calls
b2c.insert_labels(adata, 
                  labels_npz_path="processed-data/VisiumHD/03_bin2cell/stardist/he.npz", 
                  basis="spatial", 
                  spatial_key="spatial_cropped",
                  mpp=mpp, 
                  labels_key="labels_he"
                 )

# visualize 
bdata = adata[mask]

#0 means unassigned
bdata = bdata[bdata.obs['labels_he']>0]
bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)

plt.figure(figsize=(10,7.5))
sc.pl.spatial(bdata, color=[None, "labels_he"], img_key="0.5_mpp", basis="spatial_cropped")
plt.savefig(root+'/plots/VisiumHD/03_bin2cell/spatial_segemented.png', dpi=300)


# expand labels to rest of cell
b2c.expand_labels(adata, 
                  labels_key='labels_he', 
                  expanded_labels_key="labels_he_expanded"
                 )

bdata = adata[mask]

#the labels obs are integers, 0 means unassigned
bdata = bdata[bdata.obs['labels_he_expanded']>0]
bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)


b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path="plots/VisiumHD/03_bin2cell/stardist/gex.tiff")

b2c.stardist(image_path="plots/VisiumHD/03_bin2cell/stardist/gex.tiff", 
             labels_npz_path="processed-data/VisiumHD/03_bin2cell/stardist/gex.npz", 
             stardist_model="2D_versatile_fluo", 
             prob_thresh=0.05, 
             nms_thresh=0.5
            )

b2c.insert_labels(adata, 
                  labels_npz_path="processed-data/VisiumHD/03_bin2cell/stardist/gex.npz", 
                  basis="array", 
                  mpp=mpp, 
                  labels_key="labels_gex"
                 )

bdata = adata[mask]

#the labels obs are integers, 0 means unassigned
bdata = bdata[bdata.obs['labels_gex']>0]
bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)