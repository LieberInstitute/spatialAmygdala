import os
os.chdir('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/')

from pyhere import here
from pathlib import Path
import session_info
import re
import json

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

out_path = here('processed-data', 'visium_stitching','transformationsHE.csv')
sample_info_path = Path( here('raw-data', 'sample_info_visium', 'Visium_DATA_2023-09-18_1325.csv'))

################################################################################
#   Functions
################################################################################

def theta_from_mat(mat):
    #   Determine theta by using arccos. Here the real angle might by theta_cos
    #   or -1 * theta_cos
    theta_cos = np.arccos(mat[:, 0, 0])
    #   Use 'theta_cos', but negate it when sin(theta) is negative, to account
    #   for the convention that arccos returns positive angles
    return theta_cos * (2 * (mat[:, 1, 0] > 0) - 1)


################################################################################
#   Clean sample info
################################################################################

#   Read in the sample sheet

                       
sample_info = pd.read_csv(sample_info_path).loc[:, ['slide', 'sample_a1', 'region_a1', 'project_a1', 'sample_b1', 'region_b1', 'project_b1', 'sample_c1', 'region_c1', 'project_c1', 'sample_d1', 'region_d1', 'project_d1']]
suffixes = ['_a1', '_b1', '_c1', '_d1']  # Add more suffixes as needed
split_dfs = {suffix: pd.concat([sample_info['slide'], sample_info.filter(like=suffix)], axis=1) for suffix in suffixes}
new_columns = ['Slide', 'brnum', 'region', 'project', 'Array']
for key, df in split_dfs.items():
    # Append the entire 'slide' column with the respective key
    df['array'] = key.replace('_', '').capitalize()
    column_mapping = {col: new_columns[i] for i, col in enumerate(df.columns) if col != 'Slide'}
        # Rename the columns using the mapping
    df.rename(columns=column_mapping, inplace=True)
    
combined_df = pd.concat(split_dfs.values(), axis=0)
sample_info = combined_df[combined_df['project'] == 'spatialAMY_LIBD4125'].reset_index(drop=True)
sample_info.index = sample_info['Slide'] + '_' + sample_info['Array']
#
#import shutil
#directory = here('processed-data/visium_stitching/imageJ/transformationsHE')
#shutil.copytree(here('processed-data/10-image_stitching/imageJ/transformations), directory)
#files = os.listdir(directory)

# Rename each file by prepending 'Br'
#for filename in files:
#    os.rename(os.path.join(directory, filename), os.path.join(directory, f'Br{filename}'))
sample_info['xml_path'] = '/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/visium_stitching/' + sample_info['brnum'] +'/' + sample_info['brnum'] +'_AMY.xml'
################################################################################
#   Determine the path to the full-resolution/raw images and spaceranger output directories
################################################################################

sample_info['spaceranger_dir'] = [
    Path(x)
    for x in here(
        'processed-data', '01_spaceranger', 'first_donor', sample_info.index, 'outs', 'spatial'
    )
]

sample_info['raw_image_path'] = [
    Path(x)
    for x in here(
        'processed-data', 'Images', 'VistoSeg', sample_info.index+'.tif'
    )
]

################################################################################
#   Add the initial transformation estimates for image stitching from ImageJ
################################################################################

transform_df_list = []
finalsample = sample_info['xml_path'].dropna().unique()

for imagej_xml_path in finalsample:
   with open(here(imagej_xml_path)) as f:
       imagej_xml = f.read()
       #   Clean the file; make sure new lines only separate XML elements
   imagej_xml = re.sub('\n', '', imagej_xml)
   imagej_xml = re.sub('\>', '>\n', imagej_xml)
   tree = ET.parse(imagej_xml_path)
   root = tree.getroot()
   titles = [patch.get('title') for patch in root.findall('.//t2_patch')]
   samples = [title.split('_hires.png')[0] for title in titles]
   slide_nums = [title.split("_")[0] for title in samples]
   array_nums = [title.split("_")[1] for title in samples]
   matrices = re.findall(r'matrix\(.*\)".*file_path=', imagej_xml)
   trans_mat = [re.sub(r'matrix\((.*)\).*', '\\1', x).split(',') for x in matrices]
   trans_mat = np.transpose(
       np.array(
           [[float(y) for y in x] for x in trans_mat], dtype = np.float64
       )
           .reshape((-1, 3, 2)),
       axes = [0, 2, 1]
   )
   assert trans_mat.shape[1:3] == (2, 3), "Improperly read ImageJ XML outputs"
   #   Grab sample info for this slide, ordered how the array numbers
   #   appear in the ImageJ output
   this_sample_info = sample_info.loc[[slide_nums[i] + '_' + array_nums[i] for i in range(len(array_nums))]]
   #   Adjust translations to represent pixels in full resolution
   json_path = os.path.join(this_sample_info['spaceranger_dir'].iloc[0],'scalefactors_json.json')
   with open(json_path, 'r') as f:
       spaceranger_json = json.load(f)
   sf = spaceranger_json['tissue_hires_scalef']
   trans_mat[:, :, 2] /= sf
   #   Compile transformation information for this slide in a DataFrame
   transform_df_list.append(
       pd.DataFrame(
           {
               'initial_transform_x': trans_mat[:, 0, 2],
               'initial_transform_y': trans_mat[:, 1, 2],
               'initial_transform_theta': theta_from_mat(trans_mat),
               'Array': array_nums,
               'Slide': slide_nums
           }
       )
   )
#   Combine across donors, then merge into sample_info
transform_df = pd.concat(transform_df_list)
index = sample_info.index
sample_info = pd.merge(sample_info, transform_df, how = 'left', on = ['Slide', 'Array'])
sample_info.index = index

sample_info.to_csv(out_path)

session_info.show()
