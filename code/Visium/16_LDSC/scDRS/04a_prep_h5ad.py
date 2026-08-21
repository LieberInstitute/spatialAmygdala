#!/usr/bin/env python
"""
04a_prep_h5ad.py -- build the single scDRS input h5ad from the seven per-section
gsMap h5ad files, plus the matching cell-level covariate table.

    python 04a_prep_h5ad.py <st_dir> <samples...> --out <h5ad> --cov <tsv> --annot <col>

WHY REUSE THE gsMap h5ad FILES
    They are already the exact object the other two pipelines analyse: the same
    SPE (spe_harmony_markers_BS_k16_Semisupervised_wAI.rds), the same
    protein-coding gene universe (19,389 symbols), the same BayesSpace domain
    labels. Reading them directly guarantees scDRS sees identical spots and
    identical genes -- no chance of a divergent filtering decision making the
    three methods incomparable.

WHAT scDRS NEEDS THAT gsMap DID NOT
    1. ONE object, not seven. scDRS's control gene sets are matched on
       mean/variance computed across the whole matrix, so all sections must be
       scored in a common feature space. Concatenating also makes the
       donor covariate meaningful.
    2. RAW COUNTS IN .X. The gsMap files carry counts in both X and
       layers['count']; we set X = layers['count'] explicitly and drop the
       layer so `--flag-raw-count True` is unambiguous.
    3. A COVARIATE FILE. scDRS regresses covariates out of the expression
       matrix before scoring. We supply:
         const     intercept (scDRS convention)
         n_genes   detected genes per spot -- the dominant technical axis in
                   Visium, and the covariate the scDRS authors use
         log_umi   log10 total UMI, since Visium spots vary in tissue content
                   far more than dissociated cells do
         donor_*   one-hot dummies for 6 of the 7 donors (the 7th is absorbed
                   by the intercept), which removes section/batch offsets
                   without removing the anatomy we are testing

    Domain labels are NOT covariates -- they enter only at the group-analysis
    step, so the per-spot scores stay annotation-free, exactly as in gsMap.

Spot barcodes are prefixed with the donor to guarantee uniqueness after
concatenation.
"""
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("st_dir")
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cov", required=True)
    ap.add_argument("--annot", required=True)
    a = ap.parse_args()

    adatas = {}
    for s in a.samples:
        path = f"{a.st_dir}/{s}.h5ad"
        d = ad.read_h5ad(path)
        if "count" in d.layers:
            d.X = d.layers["count"]
            del d.layers["count"]
        d.X = sp.csr_matrix(d.X)
        d.obs["donor"] = s
        d.obs_names = [f"{s}_{b}" for b in d.obs_names]
        adatas[s] = d
        print(f"{s}: {d.n_obs} spots x {d.n_vars} genes", flush=True)

    genes = [set(d.var_names) for d in adatas.values()]
    shared = set.intersection(*genes)
    print(f"\ngene universe: {[len(g) for g in genes]} -> shared {len(shared)}")
    if any(len(g) != len(shared) for g in genes):
        raise SystemExit("sections disagree on the gene universe -- inspect before proceeding")

    adata = ad.concat(list(adatas.values()), join="outer", label=None, index_unique=None)
    adata.var = list(adatas.values())[0].var.loc[adata.var_names].copy()

    # --- integrity checks on the counts matrix -------------------------------
    X = adata.X
    assert sp.issparse(X), "X is not sparse"
    dat = X.data
    assert np.all(dat > 0), "negative or zero stored values in X"
    assert np.allclose(dat, np.round(dat)), "X is not integer counts"
    assert adata.obs_names.is_unique, "duplicate spot barcodes after concat"

    n_umi = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()

    adata.obs["n_umi"] = n_umi
    adata.obs["n_genes"] = n_genes

    # --- covariates ----------------------------------------------------------
    cov = pd.DataFrame(index=adata.obs_names)
    cov["const"] = 1
    cov["n_genes"] = n_genes
    cov["log_umi"] = np.log10(n_umi + 1)
    donors = sorted(adata.obs["donor"].unique())
    for d in donors[1:]:                       # drop_first: donors[0] is the reference
        cov[f"donor_{d}"] = (adata.obs["donor"] == d).astype(int).to_numpy()
    cov.index.name = "index"
    cov.to_csv(a.cov, sep="\t")

    adata.write_h5ad(a.out, compression="gzip")

    print(f"\ncombined: {adata.n_obs} spots x {adata.n_vars} genes")
    print(f"total UMI  : {n_umi.sum():,.0f}")
    print(f"n_genes    : min {n_genes.min()}  median {np.median(n_genes):.0f}  max {n_genes.max()}")
    print(f"n_umi      : min {n_umi.min():.0f}  median {np.median(n_umi):.0f}  max {n_umi.max():.0f}")
    print(f"\nspots per donor:\n{adata.obs['donor'].value_counts().to_string()}")
    print(f"\nspots per {a.annot}:\n{adata.obs[a.annot].value_counts().to_string()}")
    print(f"\ncovariates : {list(cov.columns)}")
    print(f"-> {a.out}")
    print(f"-> {a.cov}")


if __name__ == "__main__":
    main()
