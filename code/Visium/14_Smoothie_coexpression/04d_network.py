#!/usr/bin/env python3
"""
04d_network.py
 
Build spatial co-expression network with chosen parameters.
Run after reviewing 04c hyperparameter plots.
 
Input:  results/pearsonR_mat_spqn.npy (preferred) or pearsonR_mat_concat.npy
Output: results/network/, results/modules_df.csv
"""
 
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import smoothie
smoothie.suppress_warnings()
 
PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
RESULTS_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "results")
 
# Set after reviewing 04c hyperparameter plots
PCC_CUTOFF = 0.9
CLUSTERING_POWER = 9
 
spqn_path = os.path.join(RESULTS_DIR, "pearsonR_mat_spqn.npy")
raw_path = os.path.join(RESULTS_DIR, "pearsonR_mat_concat.npy")
 
if os.path.exists(spqn_path):
    print("Using SpQN-corrected correlation matrix")
    pearsonR_mat = np.load(spqn_path)
else:
    print("Using uncorrected correlation matrix")
    pearsonR_mat = np.load(raw_path)
 
gene_names = pd.read_csv(os.path.join(RESULTS_DIR, "gene_names_concat.csv"))["gene_name"].values
print(f"  {len(gene_names)} genes, matrix: {pearsonR_mat.shape}")
 
net_dir = os.path.join(RESULTS_DIR, "network")
os.makedirs(net_dir, exist_ok=True)
edge_list, node_label_df = smoothie.make_spatial_network(
    pearsonR_mat=pearsonR_mat,
    gene_names=gene_names,
    pcc_cutoff=PCC_CUTOFF,
    clustering_power=CLUSTERING_POWER,
    output_folder=net_dir
)
 
modules_df = node_label_df.groupby("module_label").filter(lambda x: len(x) >= 2)
modules_df.to_csv(os.path.join(RESULTS_DIR, f"modules_df_{PCC_CUTOFF}_{CLUSTERING_POWER}.csv"), index=False)