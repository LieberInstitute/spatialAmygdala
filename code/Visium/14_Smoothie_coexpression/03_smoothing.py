#!/usr/bin/env python3
"""
03_smoothing.py

Gaussian smoothing on each normalized donor. Most expensive step.

Input:  processed-data/Smoothie/normalized/<DonorID>_normalized.h5ad
Output: processed-data/Smoothie/smoothed/<DonorID>_smoothed.h5ad
"""

import os
import anndata as ad
import smoothie
smoothie.suppress_warnings()

PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
NORM_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "normalized")
OUT_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "smoothed")
os.makedirs(OUT_DIR, exist_ok=True)

DONORS = ["Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660"]

# Smoothing parameters for standard Visium (pixel coordinates)
# Measured median nearest-neighbor distance: ~277 px (~100 µm spacing)
# So micron_to_unit ≈ 277/100 = 2.77 px/µm
MICRON_TO_UNIT = 2.77       # px/µm (derived from actual spot spacing)
TARGET_MICRONS = 150.0      # ~1.5× spot spacing for standard Visium
GRID_BASED = False          # in-place for cell-sized/binned data
MIN_SPOTS = 3               # low value for Visium's sparse grid (per Smoothie docs)

gaussian_sd = TARGET_MICRONS * MICRON_TO_UNIT
print(f"gaussian_sd = {gaussian_sd:.2f} px ({TARGET_MICRONS} µm × {MICRON_TO_UNIT})")

for i, donor in enumerate(DONORS):
    adata = ad.read_h5ad(os.path.join(NORM_DIR, f"{donor}_normalized.h5ad"))
    sm_adata = smoothie.run_parallelized_smoothing(
        adata,
        grid_based_or_not=GRID_BASED,
        gaussian_sd=gaussian_sd,
        min_spots_under_gaussian=MIN_SPOTS
    )
    sm_adata.write_h5ad(os.path.join(OUT_DIR, f"{donor}_smoothed.h5ad"))
    print(f"{donor} ({i+1}/{len(DONORS)}): {sm_adata.shape[0]} spots × {sm_adata.shape[1]} genes")