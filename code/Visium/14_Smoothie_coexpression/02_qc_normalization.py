#!/usr/bin/env python3
"""
02_qc_and_normalize.py

Load per-donor .h5ad, QC filter, CPT + log1p normalize, save.

Input:  processed-data/Smoothie/h5ad/<DonorID>.h5ad
Output: processed-data/Smoothie/normalized/<DonorID>_normalized.h5ad
"""

import os
import anndata as ad
import scanpy as sc

PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
H5AD_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "h5ad")
OUT_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "normalized")
os.makedirs(OUT_DIR, exist_ok=True)

DONORS = ["Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660"]

for donor in DONORS:
    adata = ad.read_h5ad(os.path.join(H5AD_DIR, f"{donor}.h5ad"))
    adata.var_names_make_unique()

    # Ensure spatial coords mapped correctly from zellkonverter
    if "spatial" not in adata.obsm:
        spatial_keys = [k for k in adata.obsm.keys() if "spatial" in k.lower()]
        if spatial_keys:
            adata.obsm["spatial"] = adata.obsm[spatial_keys[0]]
        else:
            raise KeyError(f"{donor}: no spatial coords in obsm. Keys: {list(adata.obsm.keys())}")
    assert adata.obsm["spatial"].shape[1] == 2

    n0, g0 = adata.shape

    # QC (binned/cell-sized tier for standard Visium)
    sc.pp.filter_cells(adata, min_counts=50)
    sc.pp.filter_genes(adata, min_counts=100)
    sc.pp.filter_genes(adata, min_cells=10)

    # CPT + log1p
    sc.pp.normalize_total(adata, target_sum=1e3)
    sc.pp.log1p(adata)

    out_path = os.path.join(OUT_DIR, f"{donor}_normalized.h5ad")
    adata.write_h5ad(out_path)
    print(f"{donor}: {n0}→{adata.shape[0]} spots, {g0}→{adata.shape[1]} genes")