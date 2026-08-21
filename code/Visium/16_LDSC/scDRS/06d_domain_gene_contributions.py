#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import scdrs

h5ad_file = sys.argv[1]
cov_file = sys.argv[2]
gs_file = sys.argv[3]
out_dir = sys.argv[4]

annot = "BS_k16_Semisupervised_wAI"

pairs = [
    ("SCZ", "LA"),
    ("MDD", "BM"),
    ("BIP_2024", "CLA"),
]

gm_domains = [
    "CLA", "HPC", "BLD", "LA", "PL", "BL",
    "BM", "CoA", "IA", "CeA", "MeA", "CHAT"
]

os.makedirs(out_dir, exist_ok=True)

# ===== Load and preprocess data =====

adata = scdrs.util.load_h5ad(
    h5ad_file,
    flag_filter_data=True,
    flag_raw_count=True
)

cov = pd.read_csv(cov_file, sep="\t", index_col=0)
cov.index = cov.index.astype(str)

scdrs.preprocess(
    adata,
    cov=cov,
    adj_prop=None,
    n_mean_bin=20,
    n_var_bin=20,
    copy=False
)

print(f"adata: {adata.n_obs} spots x {adata.n_vars} genes")

# clean domain labels
domain = adata.obs[annot].astype(str)

domain = domain.str.replace(
    r"\s*\(\d+/\d+\)\s*$",
    "",
    regex=True
)

domain = domain.replace({"AI": "IA"})

# load gene sets
dict_gs = scdrs.util.load_gs(
    gs_file,
    src_species="human",
    dst_species="human",
    to_intersect=adata.var_names
)

gene_stats = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
param = adata.uns["SCDRS_PARAM"]

cov_list = list(param["COV_MAT"].columns)

# ===== Run trait x domain pairs =====

for trait, target_domain in pairs:

    print()
    print(f"===== {trait} - {target_domain} =====")

    gene_list, magma_weight = dict_gs[trait]

    # match ordering used internally by scDRS
    weight_dict = dict(zip(gene_list, magma_weight))

    gene_list = sorted(weight_dict)
    magma_weight = np.array(
        [weight_dict[g] for g in gene_list],
        dtype=float
    )

    print(f"{trait}: {len(gene_list)} genes")

    # scDRS variance-adjusted gene weights
    var_tech = gene_stats.loc[
        gene_list,
        "var_tech"
    ].values

    scdrs_weight = (
        magma_weight /
        np.sqrt(var_tech + 1e-2)
    )

    scdrs_weight = scdrs_weight / scdrs_weight.sum()

    # target domain
    target_mask = (domain == target_domain).values

    # other gray-matter domains retained for QC / optional contrasts
    other_gm_mask = (
        domain.isin(gm_domains) &
        (domain != target_domain)
    ).values

    print(f"{target_domain} spots : {target_mask.sum()}")
    print(f"other GM spots : {other_gm_mask.sum()}")

    gene_idx = adata.var_names.get_indexer(gene_list)

    beta = param["COV_BETA"].loc[
        gene_list,
        cov_list
    ].values

    gene_mean = param["COV_GENE_MEAN"].loc[
        gene_list
    ].values

    def mean_scdrs_expression(mask):

        x_mean = np.asarray(
            adata.X[mask][:, gene_idx].mean(axis=0)
        ).ravel()

        cov_mean = param["COV_MAT"].loc[
            adata.obs_names[mask],
            cov_list
        ].mean(axis=0).values

        x_mean = (
            x_mean +
            cov_mean @ beta.T +
            gene_mean
        )

        return x_mean

    mean_target = mean_scdrs_expression(target_mask)
    mean_other = mean_scdrs_expression(other_gm_mask)

    # per-gene contribution to raw disease score
    contribution_target = scdrs_weight * mean_target
    contribution_other = scdrs_weight * mean_other

    delta_expression = mean_target - mean_other

    delta_contribution = (
        contribution_target -
        contribution_other
    )

    res = pd.DataFrame({
        "gene": gene_list,
        "trait": trait,
        "domain": target_domain,
        "magma_weight": magma_weight,
        "var_tech": var_tech,
        "scdrs_weight": scdrs_weight,
        "mean_expr_domain": mean_target,
        "mean_expr_otherGM": mean_other,
        "delta_expression": delta_expression,
        "contribution_domain": contribution_target,
        "contribution_otherGM": contribution_other,
        "delta_raw_contribution": delta_contribution
    })

    res = res.sort_values(
        "contribution_domain",
        ascending=False
    )

    res["rank"] = np.arange(1, len(res) + 1)

    # check that the gene-level values exactly reconstruct the raw score
    raw_target = contribution_target.sum()

    print(f"Mean raw {trait} score, {target_domain}: {raw_target:.6f}")
    print(f"Sum of gene contributions          : {contribution_target.sum():.6f}")

    assert np.allclose(
        raw_target,
        contribution_target.sum()
    )

    print()
    print("Top contributors:")

    print(
        res[
            [
                "gene",
                "magma_weight",
                "scdrs_weight",
                "mean_expr_domain",
                "contribution_domain"
            ]
        ].head(20).to_string(index=False)
    )

    out_file = os.path.join(
        out_dir,
        f"{trait}_{target_domain}_gene_contributions.tsv"
    )

    res.to_csv(
        out_file,
        sep="\t",
        index=False
    )

    print(f"\n-> {out_file}")