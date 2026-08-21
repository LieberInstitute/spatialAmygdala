#!/usr/bin/env python
"""
07_figures.py -- summary figures and the scDRS vs classic s-LDSC comparison.

    python 07_figures.py --workdir <WORKDIR> --annot <col> \
        --ldsc <../LDSC/ldsc_results.csv> --outdir <FIG_DIR> [--traits SCZ MDD Height]

Produces, for whichever traits have finished:
  <trait>_spatial.png        per-spot normalised scDRS score in tissue space,
                             faceted by donor -- the map that makes the point
  <trait>_domain_box.png     score distribution per BayesSpace domain
  domain_heatmap.png         domain x trait matrix of group-analysis
                             association z-scores, starred at UNCORRECTED
                             MC p < 0.05 (no multiple-testing correction --
                             apply your own across the 15 x 40 grid before
                             calling anything significant)
  method_comparison.png      scDRS domain z vs classic s-LDSC coefficient z,
                             for traits present in both
  domain_summary.tsv         the tidy table behind the heatmap

The comparison panel is scDRS vs classic s-LDSC only -- gsMap results are not
read by this script. The two methods answer the same question with different
units of analysis (domain pseudobulk for s-LDSC, spot-level expression score
for scDRS), so agreement is meaningful evidence. To make it three-way, join
gsMap's cauchy_across_samples/<trait>_cauchy.csv.gz on domain x trait; that is
left as a separate step because gsMap's sweep is still running.
"""
import argparse
import os
import glob
import inspect
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad


def load_scores(score_dir):
    """Per-spot scores for every finished trait.

    Only `norm_score` is ever used downstream (fig_spatial, fig_domain_box), but
    a .full_score.gz carries 1,007 columns -- the score plus 1,000 control-set
    columns and the summary stats -- at ~1 GB compressed each. Reading all of
    them for 42 traits needs ~70 GB of RAM to use one column, so restrict the
    read with usecols. The small .score.gz carries the same `norm_score` over
    the same spots (the control columns are what make the other file large), so
    prefer it when present and fall back to .full_score.gz otherwise.
    """
    # The barcode column is literally named "index" in these files (not blank),
    # so it must be named in usecols or pandas promotes norm_score to the index.
    keep = ("index", "", "norm_score")
    out = {}
    seen = set()
    for pattern, suffix in (("*.score.gz", ".score.gz"),
                            ("*.full_score.gz", ".full_score.gz")):
        for f in sorted(glob.glob(os.path.join(score_dir, pattern))):
            base = os.path.basename(f)
            if suffix == ".score.gz" and base.endswith(".full_score.gz"):
                continue
            trait = base[: -len(suffix)]
            if trait in seen:
                continue
            # index_col=0 keeps the unnamed barcode column as the index.
            out[trait] = pd.read_csv(
                f, sep="\t", index_col=0, usecols=lambda c: c in keep,
            )
            seen.add(trait)
    return out


def load_group(downstream_dir, annot):
    rows = []
    for f in sorted(glob.glob(os.path.join(downstream_dir, f"*.scdrs_group.{annot}"))):
        trait = os.path.basename(f).split(".scdrs_group.")[0]
        g = pd.read_csv(f, sep="\t", index_col=0)
        g["trait"] = trait
        g["domain"] = g.index
        rows.append(g)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def fig_spatial(adata, scores, trait, annot, outpath):
    obs = adata.obs
    common = obs.index.intersection(scores.index)
    v = scores.loc[common, "norm_score"]
    xy = adata[common].obsm["spatial"]
    donors = obs.loc[common, "donor"].to_numpy()
    uniq = sorted(pd.unique(donors))

    ncol = 4
    nrow = int(np.ceil(len(uniq) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 4.0 * nrow))
    axes = np.atleast_1d(axes).ravel()
    lim = float(np.nanpercentile(np.abs(v), 99))

    for ax, d in zip(axes, uniq):
        m = donors == d
        s = ax.scatter(xy[m, 0], -xy[m, 1], c=v[m], s=1.2,
                       cmap="RdBu_r", vmin=-lim, vmax=lim, linewidths=0)
        ax.set_title(d, fontsize=11)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
    for ax in axes[len(uniq):]:
        ax.set_visible(False)

    cb = fig.colorbar(s, ax=axes.tolist(), shrink=0.6, pad=0.02)
    cb.set_label("scDRS normalised score", fontsize=10)
    fig.suptitle(f"{trait} — per-spot disease relevance", fontsize=14)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_domain_box(adata, scores, trait, annot, outpath):
    common = adata.obs.index.intersection(scores.index)
    df = pd.DataFrame({
        "score": scores.loc[common, "norm_score"].to_numpy(),
        "domain": adata.obs.loc[common, annot].astype(str).to_numpy(),
    })
    order = df.groupby("domain")["score"].median().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(9, 4.5))
    # matplotlib renamed boxplot(labels=) -> tick_labels= in 3.9 and removed the
    # old spelling in 3.11 (this env has 3.11.1). Pass whichever the installed
    # version accepts so the script runs on both.
    _bp_label_kw = (
        "tick_labels"
        if "tick_labels" in inspect.signature(plt.Axes.boxplot).parameters
        else "labels"
    )
    ax.boxplot([df.loc[df.domain == d, "score"] for d in order],
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#cfd8e3", edgecolor="#41556b"),
               medianprops=dict(color="#b2182b", linewidth=1.6),
               **{_bp_label_kw: list(order)})
    ax.axhline(0, color="0.5", lw=0.8, ls="--")
    ax.set_ylabel("scDRS normalised score")
    ax.set_title(f"{trait} — score by BayesSpace domain")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def fig_heatmap(grp, outpath):
    zcol = "assoc_mcz" if "assoc_mcz" in grp.columns else "assoc_mcp"
    mat = grp.pivot(index="domain", columns="trait", values=zcol)
    fdr = grp.pivot(index="domain", columns="trait", values="assoc_mcp") \
        if "assoc_mcp" in grp.columns else None

    fig, ax = plt.subplots(figsize=(0.45 * mat.shape[1] + 3, 0.42 * mat.shape[0] + 2))
    lim = float(np.nanmax(np.abs(mat.to_numpy())))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=9)

    if fdr is not None:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                p = fdr.iloc[i, j]
                if pd.notna(p) and p < 0.05:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=9)

    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label(zcol)
    ax.set_title("scDRS domain-level association (* MC p < 0.05)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def fig_method_comparison(grp, ldsc_csv, outpath):
    if not os.path.exists(ldsc_csv):
        print(f"skip comparison: {ldsc_csv} not found")
        return
    ld = pd.read_csv(ldsc_csv)
    ld = ld.rename(columns={"cell": "domain", "Coefficient_z.score": "ldsc_z"})
    # classic pipeline names traits in prose ("Schizophrenia_PGC3"); map the
    # overlap we can resolve unambiguously
    name_map = {
        "Schizophrenia_PGC3": "SCZ", "Schizophrenia": "SCZ_PGC2_CLOZUK",
        "Depression": "MDD", "Depression_ex23andMe": "MDD_ex23andMe",
        "mdd2019edinburgh": "MDD_2019_Edinburgh",
        "Bipolar Disorder": "BIP", "Bipolar Disorder2": "BIP_PGC3",
        "Autism": "Autism", "ADHD": "ADHD", "PTSD": "PTSD",
        "Anorexia": "Anorexia", "Neuroticism": "Neuroticism",
        "Education Years": "EduYears", "Intelligence": "Intelligence",
        "Insomnia": "Insomnia", "Height": "Height", "BMI": "BMI",
        "Type_2_Diabetes": "T2D", "Parkinson Disease": "PD",
        "Alzheimer Disease": "Alzheimer",
    }
    ld["trait"] = ld["trait"].map(name_map)
    ld = ld.dropna(subset=["trait"])

    zcol = "assoc_mcz" if "assoc_mcz" in grp.columns else "assoc_mcp"
    m = grp.merge(ld[["domain", "trait", "ldsc_z"]], on=["domain", "trait"], how="inner")
    if m.empty:
        print("skip comparison: no overlapping domain x trait pairs")
        return

    r = np.corrcoef(m[zcol], m["ldsc_z"])[0, 1]
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.axhline(0, color="0.7", lw=0.8); ax.axvline(0, color="0.7", lw=0.8)
    for t, sub in m.groupby("trait"):
        ax.scatter(sub["ldsc_z"], sub[zcol], s=26, label=t, alpha=0.85)
    ax.set_xlabel("classic s-LDSC coefficient z")
    ax.set_ylabel(f"scDRS group {zcol}")
    ax.set_title(f"scDRS vs stratified LDSC\nn={len(m)} domain×trait pairs, r={r:.2f}")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"comparison: n={len(m)}, pearson r={r:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--ldsc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--score-subdir", default="score")
    ap.add_argument("--downstream-subdir", default="downstream")
    ap.add_argument("--traits", nargs="*", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    score_dir = os.path.join(a.workdir, a.score_subdir)
    down_dir = os.path.join(a.workdir, a.downstream_subdir)

    scores = load_scores(score_dir)
    print(f"scored traits: {sorted(scores)}")
    if not scores:
        raise SystemExit(f"no *.full_score.gz in {score_dir}")

    adata = ad.read_h5ad(a.h5ad, backed="r")

    want = a.traits or sorted(scores)
    for t in want:
        if t not in scores:
            print(f"skip {t}: no score file"); continue
        fig_spatial(adata, scores[t], t, a.annot, os.path.join(a.outdir, f"{t}_spatial.png"))
        fig_domain_box(adata, scores[t], t, a.annot, os.path.join(a.outdir, f"{t}_domain_box.png"))
        print(f"figures for {t}")

    grp = load_group(down_dir, a.annot)
    if not grp.empty:
        grp.to_csv(os.path.join(a.outdir, "domain_summary.tsv"), sep="\t", index=False)
        fig_heatmap(grp, os.path.join(a.outdir, "domain_heatmap.png"))
        fig_method_comparison(grp, a.ldsc, os.path.join(a.outdir, "method_comparison.png"))
    else:
        print(f"no group-analysis files in {down_dir}")

    print(f"-> {a.outdir}")


if __name__ == "__main__":
    main()
