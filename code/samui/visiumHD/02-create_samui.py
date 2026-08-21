import os
os.chdir('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')

from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import sys
from loopy.sample import Sample
import tifffile
from PIL import Image
import re
import matplotlib.pyplot as plt

import scanpy as sc
from rasterio import Affine
from loopy.utils.utils import remove_dupes, Url
import re
from scipy.spatial import KDTree
import glob


this_sample = sys.argv[1:]
#with open("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/samui/visiumHD/sample_ids-visium.txt", "r", encoding="utf-8") as f:
#    lines = [line.rstrip("\r\n") for line in f.readlines()]
#file = open("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/samui/visiumHD/sample_ids-visium.txt", "r")
#this_sample = file.readlines()
#this_sample = lines[1]
this_sample = ''.join(this_sample)
#this_sample = 'Br9280_CeA'


#spot_diameter_m = 55e-6 # 5-micrometer diameter for Visium spot
spot_diameter_m = 16e-6 # 16-micrometer diameter for HDbin
#m_per_px = 4.20390369911423e-07
#0.00016÷63.36917768880937 = 0.00000252
#0.00008÷31.684588844404686 = 0.00000252
################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################
spg_path = here("processed-data", "samui", "visiumHD", "spe_visiumHD.h5ad")
spg = sc.read(spg_path)

unique_sample_ids = spg.obs['sample_id'].unique
#unique_capture_ids = spg.obs['capture_id'].unique
#spgP = spg
#path_groups = spg.obs['path_groups'].cat.categories
#spgP = spg[spg.obs['brnum'] == this_sample, :]
#capture_id=spgP.obs['capture_id'].unique()[0]
spgP = spg[spg.obs['sample_id'] == this_sample, :]
sample_id=spgP.obs['sample_id'].unique()[0]

samui_dir = Path(here('processed-data', 'samui', "visiumHD", f"{sample_id}"))
samui_dir.mkdir(parents = True, exist_ok = True)
json_path = Path(here("processed-data", "VisiumHD", "01_spaceranger", "submission_links", sample_id, "outs", "binned_outputs", "square_016um", "spatial", "scalefactors_json.json"))
#json_path = Path(here("processed-data", "VisiumHD", "01_spaceranger", "submission_links", sample_id, "outs", "segmented_outputs", "spatial", "scalefactors_json.json"))
#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']
#m_per_px = spaceranger_json['microns_per_pixel']
#m_per_px = 2.5264872653926174e-07

# spgP.obs.index = spgP.obs.index.str.replace('_'+unique_sample_id , '')
#spgP.obs.index.name = "barcode"


#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spgP.X.toarray(),
    index = spgP.obs.index,
    columns = spgP.var['Symbol']
)
#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()
#gene_df.index.name = None

#gene_df.index = gene_df.index.str.split('_').str[0]

#precast_columns = spgP.obs.filter(like="PRECAST")
#precast_columns = spgP.obs[["spd_label"]].join(precast_columns)
#precast_df = pd.DataFrame(precast_columns)

#sample_df = spgP.obs[['capture_id']].copy()
#spotCalling_df = spgP.obs[['pnn_pos', 'neuropil_pos', 'neun_pos', 'vasc_pos']]
#spotCalling_metrics = spgP.obs[['spg_NDAPI', 'spg_PDAPI', 'spg_IDAPI', 'spg_CNDAPI',
#       'spg_NNeuN', 'spg_PNeuN', 'spg_INeuN', 'spg_CNNeuN', 'spg_NWFA', 'spg_PWFA', 'spg_IWFA', 'spg_CNWFA',
#       'spg_NClaudin5', 'spg_PClaudin5', 'spg_IClaudin5']]

tissue_positions_cols=spgP.obs.filter(like="pxl_")
tissue_positions_df=pd.DataFrame(tissue_positions_cols)
################################################################################
#   Use the Samui API to create the importable directory for this combined "sample"
################################################################################
img_channels = 'rgb'
#default_channels = {'blue': 'DAPI', 'green': 'NeuN', 'yellow': 'Claudin5', 'red': 'WFA', 'white':'segDAPI', 'white':'segNeuN', 'white':'segWFA', 'white':'segClaudin5'}
#img_path = here('processed-data', 'Images', 'VistoSeg', 'Capture_areas', '{}.tif')
img_name = sample_id +'.tif'
#img_path = here('processed-data', 'Images', 'VistoSeg', img_name)
#img_path = here('raw-data', 'Images', img_name)
#img_path = here('raw-data', 'images', 'VisiumHD', img_name)
img_path = here('processed-data', 'samui', 'visiumHD', 'images', img_name)

#tissue_positions_path = Path(here("processed-data", "01_spaceranger", capture_id, "outs", "spatial", "tissue_positions.csv"))
#tissue_positions = pd.read_csv(tissue_positions_path ,index_col = 0).rename({'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions = tissue_positions_df.rename({'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions.index.name = None
tissue_positions = tissue_positions[['x', 'y']].astype(int)
#tp_sub = tissue_positions.reindex(gene_df.index).dropna(how="all")
 
default_gene = 'SNAP25'
assert default_gene in gene_df.columns, "Default gene not in AnnData"

#notes_md_url = Url('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VSPG_image_stitching/feature_notes.md')
#this_sample = Sample(name = samui_dir.name, path = samui_dir, notesMd = notes_md_url)
this_sample = Sample(name = samui_dir.name, path = samui_dir)
this_sample.add_coords(tissue_positions, name = "coords", mPerPx = m_per_px, size = spot_diameter_m)
this_sample.add_image(tiff = img_path, channels = img_channels, scale = m_per_px)

#this_sample.add_csv_feature(precast_df, name = "Domains", coordName = "coords", dataType = "categorical")
#this_sample.add_csv_feature(spotCalling_df, name = "Spot_Calling", coordName = "coords", dataType = "categorical")
#this_sample.add_csv_feature(spotCalling_metrics, name = "Spot_Calling_metrics", coordName = "coords", dataType = "quantitative")
this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")
this_sample.set_default_feature(group = "Genes", feature = default_gene)

#   Add additional requested observational columns (colData columns)
#this_sample.add_csv_feature(
#    spgP.obs["Domain"], name = "Domain", coordName = "coords",
#    dataType = "categorical"
#)
#this_sample.add_csv_feature(
#    spgP.obs["LC_neuromelanin_or_other"], name = "LC neuromelanin or other", coordName = "coords",
#    dataType = "categorical"
#)

#this_sample.add_csv_feature(
#    spgP.obs["N_NMbodies"], name = "N NMbodies", coordName = "coords",
#    dataType = "quantitative"
#)




this_sample.write()


session_info.show()
