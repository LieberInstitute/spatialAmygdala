#!/usr/bin/env python3
"""
05_plot_modules.py

Spatial module and gene expression plots across all donors (PDF output).

Input:  processed-data/Smoothie/smoothed/, results/modules_df.csv, results/network/
Output: processed-data/Smoothie/results/plots/
          module_atlas.pdf, example_genes.pdf, + Smoothie built-in module PDFs
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import smoothie
smoothie.suppress_warnings()

PROJECT_ROOT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
SM_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "smoothed")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "processed-data", "Visium", "14_smoothie_coexpression", "results")
PLOT_DIR = os.path.join(PROJECT_ROOT, "plots", "Visium", "14_smoothie_coexpression")
os.makedirs(PLOT_DIR, exist_ok=True)

DONORS = ["Br9192", "Br9280", "Br2743", "Br6471", "Br8325", "Br6423", "Br6660"]
SPOT_SIZE = 8
FIG_W = 1.5  # width per sample panel


# --- Helpers ---

def spatial_scatter(ax, coords, values, title=None):
    ax.scatter(coords[:, 0], coords[:, 1], c=values, s=SPOT_SIZE,
               cmap="magma", edgecolors="none", rasterized=True)
    ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=6, pad=2)

def gene_vals(sm_adata, gene):
    if gene not in sm_adata.var_names:
        return np.zeros(sm_adata.shape[0])
    idx = list(sm_adata.var_names).index(gene)
    X = sm_adata.X
    return X[:, idx].toarray().flatten() if hasattr(X, "toarray") else np.asarray(X[:, idx]).flatten()

def module_score(sm_adata, genes):
    present = [g for g in genes if g in sm_adata.var_names]
    if not present:
        return np.zeros(sm_adata.shape[0])
    idx = [list(sm_adata.var_names).index(g) for g in present]
    X = sm_adata.X
    vals = X[:, idx].toarray().mean(axis=1) if hasattr(X, "toarray") else X[:, idx].mean(axis=1)
    return np.asarray(vals).flatten()

def donor_row(fig_title, value_fn, donors, sm_list, pdf):
    """Plot one row of donor panels and save to pdf."""
    n = len(donors)
    fig, axes = plt.subplots(1, n, figsize=(FIG_W * n, FIG_W + 0.3))
    if n == 1: axes = [axes]
    for j, (d, sm) in enumerate(zip(donors, sm_list)):
        spatial_scatter(axes[j], sm.obsm["spatial"], value_fn(sm), title=d)
    fig.suptitle(fig_title, fontsize=8, y=1.02)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches="tight", dpi=200)
    plt.close(fig)


# --- Load ---

sm_list = [ad.read_h5ad(os.path.join(SM_DIR, f"{d}_smoothed.h5ad")) for d in DONORS]
modules_df = pd.read_csv(os.path.join(RESULTS_DIR, "modules_df.csv"))
node_label_df = pd.read_csv(os.path.join(RESULTS_DIR, "network", "node_labels.csv"))

mod_labels = sorted([
    m for m in modules_df["module_label"].unique()
    if len(modules_df[modules_df["module_label"] == m]) >= 3
])
print(f"{len(mod_labels)} modules with ≥3 genes")


# --- 1. Smoothie built-in module plots ---

smoothie.plot_modules_multisample(
    sm_list, DONORS, node_label_df,
    output_folder=PLOT_DIR, shared_scaling=False,
    min_genes=3, spot_size=SPOT_SIZE, dpi=300, file_format="pdf"
)


# --- 2. Module atlas ---

with PdfPages(os.path.join(PLOT_DIR, "module_atlas.pdf")) as pdf:
    for m in mod_labels:
        genes = modules_df[modules_df["module_label"] == m]["name"].values
        donor_row(f"Module {m} ({len(genes)} genes)",
                  lambda sm, g=genes: module_score(sm, g), DONORS, sm_list, pdf)
print(f"module_atlas.pdf: {len(mod_labels)} pages")


# --- 3. Example genes (top gene from first 10 modules) ---

with PdfPages(os.path.join(PLOT_DIR, "example_genes.pdf")) as pdf:
    for m in mod_labels[:10]:
        sub = modules_df[modules_df["module_label"] == m]
        top = (sub.sort_values("degree", ascending=False).iloc[0]["name"]
               if "degree" in sub.columns else sub.iloc[0]["name"])
        donor_row(f"{top} (Module {m})",
                  lambda sm, g=top: gene_vals(sm, g), DONORS, sm_list, pdf)
print("example_genes.pdf: done")