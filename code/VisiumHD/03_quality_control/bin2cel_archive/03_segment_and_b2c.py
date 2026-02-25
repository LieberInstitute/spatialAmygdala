# This code was modified from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/code/04_bin2cell
import matplotlib.pyplot as plt
import scanpy as sc
import os
import pandas as pd
from pyhere import here
import session_info
import bin2cell as b2c

sample_info_path = here('raw-data', 'sample_info', 'VisiumHD_sample_info.csv')
stardist_dir = here('processed-data', 'VisiumHD', '03_bin2cell', 'stardist')
plot_dir = here('plots', 'VisiumHD', '03_bin2cell')

mpp = 0.3

os.makedirs(plot_dir, exist_ok=True)

#   Read in sample info and subset to this ID
sample_info = pd.read_csv(sample_info_path)
sample_info = sample_info.iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]

in_path = here(
    'processed-data', 'VisiumHD', '03_bin2cell', f'basic_{sample_info["sample_id"]}.h5ad'
)
out_path = here(
    'processed-data', 'VisiumHD', '03_bin2cell', f'sc_{sample_info["sample_id"]}.h5ad'
)

adata = sc.read(in_path)

#   Region for plots
mask = (
    (adata.obs['array_row'] >= 1000) & 
    (adata.obs['array_row'] <= 1050) & 
    (adata.obs['array_col'] >= 1000) & 
    (adata.obs['array_col'] <= 1050)
)

#   Normalize counts to account for "striping" effect
b2c.destripe(adata)

#   Plot before and after destriping
sc.pl.spatial(
    adata[mask], color=[None, "n_counts", "n_counts_adjusted"],
    img_key=f"{mpp}_mpp", basis="spatial_cropped"
)
plt.savefig(os.path.join(plot_dir, f'{sample_info["sample_id"]}_destripe.png'))
plt.close('all')

#   Segment nuclei on H&E image
b2c.stardist(
    image_path=os.path.join(
        stardist_dir, f'he_{sample_info["sample_id"]}.tiff'
    ),
    labels_npz_path=os.path.join(
        stardist_dir, f'he_{sample_info["sample_id"]}.npz'
    ),
    stardist_model="2D_versatile_he", 
    prob_thresh=0.01
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path=os.path.join(
        stardist_dir, f'he_{sample_info["sample_id"]}.npz'
    ), 
    basis="spatial", 
    spatial_key="spatial_cropped",
    mpp=mpp, 
    labels_key="labels_he"
)

#   Plot some example nuclei
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_he'] > 0]
bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)
sc.pl.spatial(
    bdata, color=[None, "labels_he"], img_key=f"{mpp}_mpp",
    basis="spatial_cropped"
)
plt.savefig(os.path.join(plot_dir, f'{sample_info["sample_id"]}_nuclei.png'))
plt.close('all')

#   Expand labels to attempt to capture cells and not nuclei
b2c.expand_labels(
    adata, 
    labels_key='labels_he', 
    expanded_labels_key="labels_he_expanded"
)

#   Plot example cells
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_he_expanded'] > 0]
bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)
sc.pl.spatial(
    bdata, color=[None, "labels_he_expanded"], img_key=f"{mpp}_mpp",
    basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_info["sample_id"]}_cells_primary.png')
)
plt.close('all')

#   Create an image from gene counts
b2c.grid_image(
    adata,
    "n_counts_adjusted",
    mpp=mpp,
    sigma=5,
    save_path=os.path.join(
        stardist_dir, f'gex_{sample_info["sample_id"]}.tiff'
    )
)

#   Segment cells on the gene-count image
b2c.stardist(
    image_path=os.path.join(
        stardist_dir, f'gex_{sample_info["sample_id"]}.tiff'
    ), 
    labels_npz_path = os.path.join(
        stardist_dir, f'gex_{sample_info["sample_id"]}.npz'
    ), 
    stardist_model="2D_versatile_fluo", 
    prob_thresh=0.05, 
    nms_thresh=0.5
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path = os.path.join(
        stardist_dir, f'gex_{sample_info["sample_id"]}.npz'
    ), 
    basis="array", 
    mpp=mpp, 
    labels_key="labels_gex"
)

#   Plot example cells (segmentation does a very poor job)
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_gex'] > 0]
bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)
sc.pl.spatial(
    bdata, color=[None, "labels_gex"], img_key=f"{mpp}_mpp",
    basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_info["sample_id"]}_cells_secondary.png')
)
plt.close('all')

#   Take the union of cell labels from both segmentation methods
b2c.salvage_secondary_labels(
    adata, 
    primary_label="labels_he_expanded", 
    secondary_label="labels_gex", 
    labels_key="labels_joint"
)

#   Plot union of cell labels
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_joint'] > 0]
bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)
sc.pl.spatial(
    bdata, color=[None, "labels_joint_source", "labels_joint"],
    img_key=f"{mpp}_mpp", basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_info["sample_id"]}_cells_both.png')
)
plt.close('all')

#   Choose to aggregate bins by H&E-segmented cells
adata = b2c.bin_to_cell(
    adata, labels_key="labels_he_expanded", spatial_keys=["spatial", "spatial_cropped"]
)

cell_mask = (
    (adata.obs['array_row'] >= 1450) & 
    (adata.obs['array_row'] <= 1550) & 
    (adata.obs['array_col'] >= 250) & 
    (adata.obs['array_col'] <= 450)
)

#   Plot counts within cells after aggregation of bins
bdata = adata[cell_mask]
sc.pl.spatial(
    bdata, color="bin_count", img_key=f"{mpp}_mpp", basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_info["sample_id"]}_cells_aggregated.png')
)
plt.close('all')

#   Save cell-by-gene AnnData
adata.var_names_make_unique()
sc.write(out_path, adata)

session_info.show()
