#!/usr/bin/env python3
"""
04_correlation_and_network.py

Concatenate smoothed data, compute correlation matrix, select hyperparameters,
build spatial co-expression network.

Input:  processed-data/Smoothie/smoothed/<DonorID>_smoothed.h5ad
Output: processed-data/Smoothie/results/
          pearsonR_mat_concat.npy, gene_names_concat.csv
          hyperparameter_selection/
          network/edge_list.csv, network/node_labels.csv
          modules_df.csv
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import smoothie
smoothie.suppress_warnings()

PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
SM_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "smoothed")
OUT_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "results")
os.makedirs(OUT_DIR, exist_ok=True)

DONORS = ["Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660"]

# Adjust after reviewing hyperparameter plots
PCC_CUTOFF = 0.4
CLUSTERING_POWER = 4

# --- Load smoothed data ---
sm_adata_list = [ad.read_h5ad(os.path.join(SM_DIR, f"{d}_smoothed.h5ad")) for d in DONORS]

# --- Concatenate & correlate ---
sm_X_concat, gene_names = smoothie.concatenate_smoothed_matrices(sm_adata_list)
print(f"Concatenated: {sm_X_concat.shape}, {len(gene_names)} shared genes")

pearsonR_mat, _ = smoothie.compute_correlation_matrix(sm_X_concat)
np.save(os.path.join(OUT_DIR, "pearsonR_mat_concat.npy"), pearsonR_mat)
pd.DataFrame({"gene_name": gene_names}).to_csv(
    os.path.join(OUT_DIR, "gene_names_concat.csv"), index=False)
print(f"Correlation matrix: {pearsonR_mat.shape}")

# --- Hyperparameter selection ---
hp_dir = os.path.join(OUT_DIR, "hyperparameter_selection")
os.makedirs(hp_dir, exist_ok=True)
smoothie.select_clustering_params(
    gene_names=gene_names,
    pearsonR_mat=pearsonR_mat,
    output_folder=hp_dir,
    pcc_cutoffs=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    clustering_powers=[1, 3, 5, 7, 9],
    min_genes_for_module=5
)
print(f"Hyperparameter plots → {hp_dir}/")
print(f"*** Review plots, then adjust PCC_CUTOFF/CLUSTERING_POWER if needed ***")

# --- Build network ---
net_dir = os.path.join(OUT_DIR, "network")
os.makedirs(net_dir, exist_ok=True)
edge_list, node_label_df = smoothie.make_spatial_network(
    pearsonR_mat=pearsonR_mat,
    gene_names=gene_names,
    pcc_cutoff=PCC_CUTOFF,
    clustering_power=CLUSTERING_POWER,
    output_folder=net_dir
)

modules_df = node_label_df.groupby("module_label").filter(lambda x: len(x) >= 2)
modules_df.to_csv(os.path.join(OUT_DIR, "modules_df.csv"), index=False)
print(f"Network: {modules_df['module_label'].nunique()} modules, "
      f"{len(modules_df)} genes → {net_dir}/")