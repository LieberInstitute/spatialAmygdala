# following the AbcProjectCache vignette here: https://alleninstitute.github.io/abc_atlas_access/notebooks/getting_started.html

from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

download_base = Path('/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/snRNAseq/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

# list
abc_cache.list_metadata_files('WHB-10Xv3')
# ['anatomical_division_structure_map', 'cell_metadata', 'donor', 'example_genes_all_cells_expression', 'gene', 'region_of_interest_structure_map']

abc_cache.list_metadata_files('WHB-taxonomy')
# ['cluster', 'cluster_annotation_term', 'cluster_annotation_term_set', 'cluster_to_cluster_annotation_membership']

# ===== Merging metadata to create extended cell metadata =====
cell = abc_cache.get_metadata_dataframe(
    directory='WHB-10Xv3',
    file_name='cell_metadata',
    dtype={'cell_label': str}
)
cell.set_index('cell_label', inplace=True)

membership = abc_cache.get_metadata_dataframe(
    directory='WHB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership'
)
membership_groupby = membership.groupby(['cluster_alias', 'cluster_annotation_term_set_name'])
membership.head(5)

term_sets = abc_cache.get_metadata_dataframe(directory='WHB-taxonomy', file_name='cluster_annotation_term_set').set_index('label')
cluster_details = membership_groupby['cluster_annotation_term_name'].first().unstack()
cluster_details = cluster_details[term_sets['name']] # order columns
cluster_details.fillna('Other', inplace=True)
cluster_details.sort_values(['supercluster', 'cluster', 'subcluster'], inplace=True)
cluster_details.head(5)


cluster_colors = membership_groupby['color_hex_triplet'].first().unstack()
cluster_colors = cluster_colors[term_sets['name']]
cluster_colors.sort_values(['supercluster', 'cluster', 'subcluster'], inplace=True)
cluster_colors.head(5)


roi = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='region_of_interest_structure_map')
roi.set_index('region_of_interest_label', inplace=True)
roi.rename(columns={'color_hex_triplet': 'region_of_interest_color'},
           inplace=True)
roi.head(5)


cell_extended = cell.join(cluster_details, on='cluster_alias')
cell_extended = cell_extended.join(cluster_colors, on='cluster_alias', rsuffix='_color')
cell_extended = cell_extended.join(roi[['region_of_interest_color']], on='region_of_interest_label')
cell_extended.head(5)
# [5 rows x 26 columns]

# save extended metadata
cell_extended.reset_index().to_csv(download_base / 'metadata' / 'abc_cell_metadata_extended.csv', index=False)

# download these datasets
abc_cache.list_expression_matrix_files('WHB-10Xv3')
# ['WHB-10Xv3-Neurons/log2', 'WHB-10Xv3-Neurons/raw', 'WHB-10Xv3-Nonneurons/log2', 'WHB-10Xv3-Nonneurons/raw']

abc_cache.get_metadata_path(directory='WHB-10Xv3', file_name='WHB-10Xv3-Neurons/raw')