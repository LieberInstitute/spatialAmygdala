
import anndata as ad
import pandas as pd
import numpy as np
import os

# ------------------------------
# 1. File paths
# ------------------------------
base_dir = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
h5_path = os.path.join(base_dir, "processed-data/snRNAseq/abc_atlas/expression_matrices/WHB-10Xv3/20240330/WHB-10Xv3-Neurons-raw.h5ad")
meta_path = os.path.join(base_dir, "processed-data/snRNAseq/abc_atlas/abc_cell_metadata_extended.csv")
output_path = os.path.join(base_dir, "processed-data/snRNAseq/01_build_sce/WHB_amygdala_ITC_subset.h5ad")

# ------------------------------
# 2. Read metadata
# ------------------------------
print("Loading metadata...")
meta = pd.read_csv(meta_path)

# ------------------------------
# 3. Define target supercluster
# ------------------------------
target_supercluster = "Eccentric medium spiny neuron"
subset_cells = meta.loc[
    meta["supercluster"] == target_supercluster,
    "cell_barcode"
].values

# ------------------------------
# 4. Load AnnData in backed mode
# ------------------------------
print("Loading AnnData in backed mode...")
adata = ad.read_h5ad(h5_path, backed="r")
obs_cells = adata.obs["cell_barcode"].values

# ------------------------------
# 5. Create boolean mask to subset
# ------------------------------
print("Creating subset mask...")
subset_cell_set = set(subset_cells)
mask = [cell in subset_cell_set for cell in obs_cells]
n_keep = np.sum(mask)

# ------------------------------
# 6. Subset and save
# ------------------------------
print("Subsetting and saving to disk...")
adata_subset = adata[mask, :]
adata_subset.write_h5ad(output_path)

# ------------------------------
# 7. Read back in and verify
# ------------------------------
print("Reloading and verifying...")
adata_final = ad.read_h5ad(output_path)

superclusters_present = adata_final.obs["supercluster"].unique()
print("Superclusters in subset:", superclusters_present)
