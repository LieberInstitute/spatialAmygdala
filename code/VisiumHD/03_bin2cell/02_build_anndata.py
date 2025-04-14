# This code was modified from /dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/code/04_bin2cell/
import scanpy as sc
import os
import pandas as pd
from pyhere import here
import session_info
import bin2cell as b2c

sample_info_path = here('raw-data', 'sample_info', 'VisiumHD_sample_info.csv')
stardist_dir = here('processed-data', 'VisiumHD', '03_bin2cell', 'stardist')
mpp = 0.3

os.makedirs(stardist_dir, exist_ok=True)

#   Read in sample info and subset to this ID
sample_info = pd.read_csv(sample_info_path)
sample_info = sample_info.iloc[int(os.environ['SLURM_ARRAY_TASK_ID']) - 1]

out_path = here(
    'processed-data', 'VisiumHD', '03_bin2cell', f'basic_{sample_info["sample_id"]}.h5ad'
)

#   Read in spaceranger outputs into an AnnData
adata = b2c.read_visium(
    os.path.join(
        sample_info['spaceranger_out'], 'outs', 'binned_outputs', 'square_002um'
    ),
    source_image_path = sample_info['raw_image'],
    spaceranger_image_path = os.path.join(
        sample_info['spaceranger_out'], 'outs', 'spatial'
    )
)

#   Require bins with nonzero counts and genes present in at least 3 bins
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

#   Create a scaled H&E image attached to the object (and for segmentation with
#   stardist)
b2c.scaled_he_image(
    adata,
    mpp = mpp,
    save_path = os.path.join(
        stardist_dir, f'he_{sample_info["sample_id"]}.tiff'
    )
)

sc.write(out_path, adata)

session_info.show()
