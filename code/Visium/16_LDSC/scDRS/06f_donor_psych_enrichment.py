#!/usr/bin/env python

import os
import sys
import pandas as pd
import anndata as ad

h5ad_file = sys.argv[1]
score_dir = sys.argv[2]
out_dir = sys.argv[3]

annot = "BS_k16_Semisupervised_wAI"

psych_traits = [
    "SCZ",
    "BIP_2024",
    "MDD",
    "PTSD_F3",
    "ADHD",
    "Autism",
    "Anorexia",
    "OUD_META",
]

gm_domains = [
    "CLA", "HPC", "BLD", "LA", "PL", "BL",
    "BM", "CoA", "IA", "CeA", "MeA", "CHAT"
]

os.makedirs(out_dir, exist_ok=True)

# ===== Load spot metadata =====

adata = ad.read_h5ad(
    h5ad_file,
    backed="r"
)

meta = adata.obs[["donor", annot]].copy()
meta["spot"] = meta.index.astype(str)

meta["domain"] = (
    meta[annot]
    .astype(str)
    .str.replace(r"\s*\(\d+/\d+\)\s*$", "", regex=True)
    .replace({"AI": "IA"})
)

meta = meta[
    meta["domain"].isin(gm_domains)
].copy()

print(meta.groupby(["donor", "domain"]).size())

# ===== Donor x trait x domain mean scDRS score =====

res = []

for trait in psych_traits:

    score_file = os.path.join(
        score_dir,
        f"{trait}.score.gz"
    )

    if not os.path.exists(score_file):
        print(f"missing: {trait}")
        continue

    print(f"reading: {trait}")

    header = pd.read_csv(
        score_file,
        sep="\t",
        nrows=0
    ).columns.tolist()

    index_col = header[0]

    score = pd.read_csv(
        score_file,
        sep="\t",
        usecols=[index_col, "norm_score"]
    )

    score = score.rename(
        columns={index_col: "spot"}
    )

    score["spot"] = score["spot"].astype(str)

    d = meta[
        ["spot", "donor", "domain"]
    ].merge(
        score,
        on="spot",
        how="inner"
    )

    d = (
        d.groupby(
            ["donor", "domain"],
            observed=True
        )
        .agg(
            mean_norm_score=("norm_score", "mean"),
            n_spots=("norm_score", "size")
        )
        .reset_index()
    )

    d["trait"] = trait

    res.append(d)

res = pd.concat(
    res,
    ignore_index=True
)

# ===== Rank domains within each donor x trait =====

res["enrichment_rank"] = (
    res.groupby(
        ["donor", "trait"]
    )["mean_norm_score"]
    .rank(
        method="average",
        pct=True
    )
)

# ===== Average psychiatric enrichment within donor =====

donor_domain = (
    res.groupby(
        ["donor", "domain"],
        observed=True
    )
    .agg(
        mean_psych_rank=("enrichment_rank", "mean"),
        mean_psych_score=("mean_norm_score", "mean"),
        n_traits=("trait", "nunique"),
        n_spots=("n_spots", "first")
    )
    .reset_index()
)

# ===== Summary across donors =====

domain_summary = (
    donor_domain.groupby(
        "domain",
        observed=True
    )
    .agg(
        mean_psych_rank=("mean_psych_rank", "mean"),
        median_psych_rank=("mean_psych_rank", "median"),
        sd_psych_rank=("mean_psych_rank", "std"),
        n_donors=("donor", "nunique")
    )
    .reset_index()
    .sort_values(
        "mean_psych_rank",
        ascending=False
    )
)

print()
print("===== Domain ranking =====")
print(domain_summary.to_string(index=False))

res.to_csv(
    os.path.join(
        out_dir,
        "psych_enrichment_donor_trait_domain.tsv"
    ),
    sep="\t",
    index=False
)

donor_domain.to_csv(
    os.path.join(
        out_dir,
        "psych_enrichment_donor_domain.tsv"
    ),
    sep="\t",
    index=False
)

domain_summary.to_csv(
    os.path.join(
        out_dir,
        "psych_enrichment_domain_summary.tsv"
    ),
    sep="\t",
    index=False
)

print(f"\n-> {out_dir}")