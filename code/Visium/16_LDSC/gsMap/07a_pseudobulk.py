#!/usr/bin/env python
"""07a: pseudobulk each donor x compartment from the gsMap h5ad inputs.

neuronal     = AI BL BLD BM CeA CHAT CLA CoA HPC LA MeA PL
non-neuronal = WM.1 WM.2 Endothelial      (user decision: vascular grouped with WM)

Output: /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond/pseudobulk_counts.tsv  (genes x 14)  + spot-count audit
"""
import os, glob
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp

ST, OUT = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap/ST", "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond"
os.makedirs(OUT, exist_ok=True)
NEURONAL = ["AI","BL","BLD","BM","CeA","CHAT","CLA","CoA","HPC","LA","MeA","PL"]
NONNEUR  = ["WM.1","WM.2","Endothelial"]
ANN = "BS_k16_Semisupervised_wAI"

cols, audit, genes = {}, [], None
for f in sorted(glob.glob(f"{ST}/*.h5ad")):
    donor = os.path.basename(f)[:-5]
    a = ad.read_h5ad(f)
    if genes is None: genes = a.var_names.to_numpy()
    assert (a.var_names.to_numpy() == genes).all(), f"gene order differs in {donor}"
    X = a.layers["count"] if "count" in a.layers else a.X
    lab = a.obs[ANN].astype(str).values
    for comp, members in [("neuronal", NEURONAL), ("nonneuronal", NONNEUR)]:
        m = np.isin(lab, members)
        n = int(m.sum())
        audit.append(dict(donor=donor, compartment=comp, n_spots=n))
        if n == 0:
            print(f"WARN {donor} {comp}: 0 spots"); continue
        v = np.asarray(X[m].sum(axis=0)).ravel() if sp.issparse(X) else X[m].sum(axis=0)
        cols[f"{donor}__{comp}"] = np.rint(v).astype(np.int64)
    del a

pb = pd.DataFrame(cols, index=genes); pb.index.name = "gene"
pb.to_csv(f"{OUT}/pseudobulk_counts.tsv", sep="\t")
aud = pd.DataFrame(audit); aud.to_csv(f"{OUT}/pseudobulk_spotcounts.tsv", sep="\t", index=False)

print("pseudobulk shape:", pb.shape)
print(aud.pivot(index="donor", columns="compartment", values="n_spots").to_string())
print("\nlibrary sizes (millions):")
print((pb.sum(0)/1e6).round(2).to_string())
print("\nDONE")
