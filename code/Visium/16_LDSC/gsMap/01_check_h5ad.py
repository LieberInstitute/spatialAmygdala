# 01_check_h5ad.py -- verify the exported h5ad files against gsMap's input contract.
# Full detail for the PILOT section; a one-line contract summary for the rest.
import numpy as np, pandas as pd, scanpy as sc, anndata, glob, os

ST_DIR = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap/ST"
PILOT  = "Br6471"
ANN    = "BS_k16_Semisupervised_wAI"

print("anndata", anndata.__version__, "| scanpy", sc.__version__)

def contract(adata):
    L = adata.layers.get("count")
    d = L.data if hasattr(L, "data") else np.asarray(L).ravel()
    S = adata.obsm.get("spatial")
    S_arr = S.to_numpy() if isinstance(S, pd.DataFrame) else np.asarray(S)
    return {
      "counts in layers['count']"       : "count" in adata.layers,
      "counts integer-valued"           : bool(np.all(d == np.floor(d))) and d.min() >= 0,
      "coords in obsm['spatial']"       : "spatial" in adata.obsm,
      "obsm['spatial'] plain ndarray"   : (not isinstance(S, pd.DataFrame)) and np.issubdtype(S_arr.dtype, np.number),
      "obsm['spatial'] is n_obs x 2"    : S_arr.shape == (adata.n_obs, 2),
      "obsm['spatial'] no NaN"          : not bool(np.isnan(S_arr).any()),
      "annotation present"              : ANN in adata.obs.columns,
      "annotation has no NA"            : int(adata.obs[ANN].isna().sum()) == 0,
      "var_names unique (symbols)"      : bool(adata.var_names.is_unique),
      "obs_names unique"                : bool(adata.obs_names.is_unique),
    }, S_arr, d

# ---------------- pilot, in full ----------------
p = os.path.join(ST_DIR, PILOT + ".h5ad")
print("\n" + "="*62)
print("PILOT:", p)
print("="*62)
adata = sc.read_h5ad(p)
print("\n### adata"); print(adata)

L = adata.layers["count"]
print("\n### layers['count']")
print("  keys        :", list(adata.layers.keys()))
print("  shape       :", L.shape)
print("  dtype       :", L.dtype)
print("  type        :", type(L).__name__)

checks, S_arr, d = contract(adata)
print("  nonzero min :", d.min(), " max:", d.max())
print("  integer     :", bool(np.all(d == np.floor(d))))
print("  total counts:", float(d.sum()))
print("  n nonzero   :", int(d.size))

print("\n### obsm['spatial']")
print("  keys        :", list(adata.obsm.keys()))
print("  python type :", type(adata.obsm["spatial"]).__name__)
print("  shape       :", S_arr.shape)
print("  dtype       :", S_arr.dtype)
print("  x range     :", float(S_arr[:,0].min()), "-", float(S_arr[:,0].max()))
print("  y range     :", float(S_arr[:,1].min()), "-", float(S_arr[:,1].max()))
print("  n NaN       :", int(np.isnan(S_arr).sum()))

print("\n### obs['%s'].value_counts()" % ANN)
print(adata.obs[ANN].value_counts())
print("  dtype       :", adata.obs[ANN].dtype)
print("  n NA        :", int(adata.obs[ANN].isna().sum()))
print("  n categories:", len(adata.obs[ANN].cat.categories))

print("\n### obs (other)")
print("  columns     :", list(adata.obs.columns))
print("  sample_id   :", dict(adata.obs["sample_id"].value_counts()))
if "exclude_overlapping" in adata.obs:
    print("  exclude_overlapping TRUE:",
          int(pd.Series(adata.obs["exclude_overlapping"]).fillna(False).astype(bool).sum()),
          "of", adata.n_obs)

print("\n### var")
print("  n_vars      :", adata.n_vars)
print("  var_names[:8]:", list(adata.var_names[:8]))
print("  unique      :", bool(adata.var_names.is_unique))

print("\n### gsMap input contract -- PILOT")
for k, v in checks.items():
    print(("  PASS  " if v else "  FAIL  ") + k)
print("  ALL PASS:", all(checks.values()))

# ---------------- every section, summarised ----------------
print("\n" + "="*62)
print("ALL SECTIONS")
print("="*62)
print(f"{'sample':10s} {'n_spots':>8s} {'n_genes':>8s} {'n_dom':>6s} {'contract':>9s}")
allok = True
for f in sorted(glob.glob(os.path.join(ST_DIR, "*.h5ad"))):
    a = sc.read_h5ad(f)
    ck, Sa, dd = contract(a)
    ok = all(ck.values()); allok &= ok
    nd = int((a.obs[ANN].value_counts() > 0).sum())
    print(f"{os.path.basename(f)[:-5]:10s} {a.n_obs:8d} {a.n_vars:8d} {nd:6d} {'PASS' if ok else 'FAIL':>9s}")
    if not ok:
        for k, v in ck.items():
            if not v: print("      FAILED:", k)
    del a
print("\nALL FILES PASS:", allok)
