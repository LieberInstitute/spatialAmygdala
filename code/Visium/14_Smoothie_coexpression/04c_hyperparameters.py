#!/usr/bin/env python3
"""
04c_hyperparameters.py
 
Hyperparameter selection for network clustering.
Review output plots before running 04d_network.py.
 
Input:  results/pearsonR_mat_spqn.npy (preferred) or pearsonR_mat_concat.npy
Output: results/hyperparameter_selection/
"""
 
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import smoothie
smoothie.suppress_warnings()
 
import smoothie.choosing_hyperparameters as _chp
if not hasattr(_chp, "os"):
    import os as _os
    _chp.os = _os
 
PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
RESULTS_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "results")
 
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
 
hp_dir = os.path.join(RESULTS_DIR, "hyperparameter_selection")
os.makedirs(hp_dir, exist_ok=True)
smoothie.select_clustering_params(
    gene_names=gene_names,
    pearsonR_mat=pearsonR_mat,
    output_folder=hp_dir,
    pcc_cutoffs=[0.6, 0.7, 0.8, 0.9],
    clustering_powers=[7, 9, 11, 13],
    min_genes_for_module=10
)
print(f"Hyperparameter plots -> {hp_dir}/")
print("*** Review plots, then run 04d_network.py ***")