#!/usr/bin/env python
"""
09_controls_compute.py -- Control analyses for the scDRS results paragraph.

Tests four claims, writing one tidy TSV per claim into $WORKDIR/controls/.
No figures here: see 10_controls_figures.py. Nothing is ever deleted or
overwritten outside $WORKDIR/controls/.

  C1  scDRS shows greater domain heterogeneity than gsMap        -> ctrl_heterogeneity.tsv
  C2  BM and LA rank among most-enriched, consistently by donor  -> ctrl_donor_ranks.tsv
  C3  scores are not driven by domain size or library size       -> ctrl_domain_covariates.tsv
                                                                    ctrl_spot_covariates.tsv
  C4  scores are not driven by GWAS power                        -> ctrl_gwas_power.tsv

Run via 09_run_controls.sh (batch; never on a login node).
"""
import os, sys, glob, gzip
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, norm
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------- paths
PROJ    = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
WORKDIR = f"{PROJ}/processed-data/Visium/16_LDSC/scDRS"
GSMAP   = f"{PROJ}/processed-data/Visium/16_LDSC/gsMap_cond_functional/cauchy_across_samples"
H5AD    = f"{WORKDIR}/h5ad/visium_amygdala_scdrs.h5ad"
SCORE   = f"{WORKDIR}/score"
DOWN    = f"{WORKDIR}/downstream"
OUT     = f"{WORKDIR}/controls"
ANNOT   = "BS_k16_Semisupervised_wAI"
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------- EDIT HERE: trait sets
# Domains excluded from display (WM.2 dropped by prior decision).
DROP_DOMAINS = ["WM.2"]

# GWAS with zero genome-wide significant SNPs -- no power to detect anything.
ZERO_POWER = ["Epilepsy_Focal", "OUD_MVP1", "OUD_MVP2", "PTSD"]

# Non-redundant trait set (one freeze per phenotype), as used in the figures.
SELECTED = [
    "ADHD", "Alzheimer_v3", "Anorexia", "Autism", "BIP_2024", "BMI",
    "EduYears", "Epilepsy_GGE", "GSCAN_DrnkWk", "GSCAN_SmkInit", "Height",
    "Insomnia", "Intelligence", "MDD", "Neuroticism", "OUD_META", "PD",
    "PTSD_F3", "SCZ", "Stroke_2022_Any", "T2D",
]

# Psychiatric subset, for the "psychiatric disease scores" claim.
PSYCHIATRIC = [
    "SCZ", "BIP_2024", "MDD", "ADHD", "Autism", "Anorexia",
    "Neuroticism", "PTSD_F3", "Insomnia",
]

# Non-brain controls, for contrast.
CONTROLS = ["Height", "BMI", "T2D"]

FDR_ALPHA = 0.05


def log(m):
    print(f"[controls] {m}", flush=True)


# =====================================================================
# Load scDRS group-level statistics (assoc_mcz / assoc_mcp per domain)
# =====================================================================
def load_scdrs_group():
    rows = []
    for d in sorted(glob.glob(f"{DOWN}/*.scdrs_group.{ANNOT}")):
        trait = os.path.basename(d).split(".scdrs_group.")[0]
        # index column is literally named "group" in scdrs output
        df = pd.read_csv(d, sep="\t")
        df = df.rename(columns={"group": "domain"})
        df["trait"] = trait
        rows.append(df)
    G = pd.concat(rows, ignore_index=True)
    log(f"scDRS group stats: {G.shape[0]} rows, "
        f"{G.trait.nunique()} traits, {G.domain.nunique()} domains")
    return G


# =====================================================================
# Load gsMap Cauchy p-values (conditional-functional arm)
# =====================================================================
def load_gsmap():
    rows = []
    for f in sorted(glob.glob(f"{GSMAP}/*_cauchy.csv.gz")):
        trait = os.path.basename(f).replace("_cauchy.csv.gz", "")
        df = pd.read_csv(f)
        df["trait"] = trait
        rows.append(df)
    M = pd.concat(rows, ignore_index=True)
    M = M.rename(columns={"annotation": "domain"})
    log(f"gsMap cauchy: {M.shape[0]} rows, {M.trait.nunique()} traits")
    return M


# =====================================================================
# C1 -- Heterogeneity across domains: scDRS vs gsMap
# =====================================================================
# Both methods are put on a common footing two ways:
#   (a) concentration -- share of total absolute signal held by the top 3
#       domains. Scale-free, so a z-statistic and a -log10 p are comparable.
#   (b) n_sig -- domains passing BH FDR within the SAME shared grid, so
#       neither method is advantaged by a different multiple-testing family.
def controls_heterogeneity(G, M, keep_dom, shared_traits):
    Sz = G.pivot(index="domain", columns="trait", values="assoc_mcz") \
          .loc[keep_dom, shared_traits]
    Sp = G.pivot(index="domain", columns="trait", values="assoc_mcp") \
          .loc[keep_dom, shared_traits]
    Gp = M.pivot(index="domain", columns="trait", values="p_cauchy") \
          .loc[keep_dom, shared_traits]
    # gsMap p -> two-sided z magnitude, so both methods carry an effect size
    Gz = pd.DataFrame(norm.isf(Gp.clip(lower=1e-320, upper=1 - 1e-16) / 2.0),
                      index=Gp.index, columns=Gp.columns)

    def conc_top3(col):
        v = np.sort(np.abs(col.values))[::-1]
        tot = v.sum()
        return np.nan if tot == 0 else v[:3].sum() / tot

    # BH within each method, across the whole shared grid
    def bh_sig(pmat):
        flat = pmat.values.flatten()
        ok = ~np.isnan(flat)
        out = np.full(flat.shape, np.nan)
        out[ok] = multipletests(flat[ok], alpha=FDR_ALPHA, method="fdr_bh")[1]
        return pd.DataFrame(out.reshape(pmat.shape),
                            index=pmat.index, columns=pmat.columns)

    Sq, Gq = bh_sig(Sp), bh_sig(Gp)

    rows = []
    for t in shared_traits:
        rows.append(dict(
            trait=t,
            scdrs_conc_top3=conc_top3(Sz[t]),
            gsmap_conc_top3=conc_top3(Gz[t]),
            scdrs_spread=Sz[t].max() - Sz[t].min(),
            gsmap_spread=Gz[t].max() - Gz[t].min(),
            scdrs_n_sig=int((Sq[t] < FDR_ALPHA).sum()),
            gsmap_n_sig=int((Gq[t] < FDR_ALPHA).sum()),
            n_domains=len(keep_dom),
        ))
    H = pd.DataFrame(rows)
    H["is_selected"] = H.trait.isin(SELECTED)
    H["is_psychiatric"] = H.trait.isin(PSYCHIATRIC)
    return H, Sz, Gz, Sq, Gq


# =====================================================================
# C2 -- Per-donor domain ranks (BM / LA consistency)
# =====================================================================
# Recomputed from per-spot scores, since the group tables pool donors and
# therefore cannot speak to across-donor consistency.
def controls_donor_ranks(obs, traits):
    rows = []
    for i, t in enumerate(traits, 1):
        f = f"{SCORE}/{t}.score.gz"
        if not os.path.exists(f):
            log(f"  MISSING score file, skipping: {t}")
            continue
        # read only the barcode index + norm_score (full_score has 1007 cols)
        sc = pd.read_csv(f, sep="\t", compression="gzip",
                         usecols=["index", "norm_score"], index_col="index")
        j = obs.join(sc, how="inner")
        if j.shape[0] == 0:
            log(f"  WARNING zero-row join for {t}")
            continue
        g = (j.groupby(["donor", "domain"], observed=True)["norm_score"]
               .agg(["mean", "size"]).reset_index())
        g["trait"] = t
        rows.append(g)
        if i % 10 == 0:
            log(f"  donor means: {i}/{len(traits)} traits")
    R = pd.concat(rows, ignore_index=True).rename(
        columns={"mean": "mean_norm_score", "size": "n_spots"})
    # rank domains within each donor x trait (1 = most enriched)
    R["rank_in_donor"] = (R.groupby(["donor", "trait"])["mean_norm_score"]
                           .rank(ascending=False, method="min").astype(int))
    R["n_domains_in_donor"] = R.groupby(["donor", "trait"])["domain"].transform("size")
    R["is_psychiatric"] = R.trait.isin(PSYCHIATRIC)
    return R


# =====================================================================
# C3 -- Domain size and library size
# =====================================================================
def controls_domain_covariates(obs, Sz, keep_dom):
    cov = (obs[obs.domain.isin(keep_dom)]
           .groupby("domain", observed=True)
           .agg(n_spots=("domain", "size"),
                n_donors=("donor", "nunique"),
                median_n_genes=("n_genes", "median"),
                median_n_umi=("n_umi", "median"),
                mean_n_genes=("n_genes", "mean"),
                mean_n_umi=("n_umi", "mean")))
    psy = [t for t in Sz.columns if t in PSYCHIATRIC]
    cov["mean_mcz_psychiatric"] = Sz[psy].mean(axis=1)
    cov["mean_mcz_all"] = Sz.mean(axis=1)
    cov["max_mcz_psychiatric"] = Sz[psy].max(axis=1)
    return cov.reset_index()


def controls_spot_covariates(obs, traits):
    """Per-trait correlation of spot-level score with library size.

    Covariates (log UMI, n_genes, donor) were regressed out during scoring,
    so these should sit near zero. This is the direct check of that.
    """
    rows = []
    for t in traits:
        f = f"{SCORE}/{t}.score.gz"
        if not os.path.exists(f):
            continue
        sc = pd.read_csv(f, sep="\t", compression="gzip",
                         usecols=["index", "norm_score"], index_col="index")
        j = obs.join(sc, how="inner")
        if j.shape[0] < 100:
            continue
        r_umi, p_umi = spearmanr(j.norm_score, j.n_umi)
        r_gen, p_gen = spearmanr(j.norm_score, j.n_genes)
        rows.append(dict(trait=t, n_spots=j.shape[0],
                         rho_norm_score_vs_n_umi=r_umi, p_n_umi=p_umi,
                         rho_norm_score_vs_n_genes=r_gen, p_n_genes=p_gen))
    S = pd.DataFrame(rows)
    S["is_selected"] = S.trait.isin(SELECTED)
    return S


# =====================================================================
# C4 -- GWAS power
# =====================================================================
def controls_gwas_power(G, sweep, keep_dom):
    Sz = G.pivot(index="domain", columns="trait", values="assoc_mcz").loc[keep_dom]
    Sp = G.pivot(index="domain", columns="trait", values="assoc_mcp").loc[keep_dom]
    per = pd.DataFrame({
        "max_abs_mcz": Sz.abs().max(),
        "max_mcz": Sz.max(),
        "spread_mcz": Sz.max() - Sz.min(),
        "min_mcp": Sp.min(),
    })
    per.index.name = "trait"
    per = per.reset_index().merge(
        sweep[["trait", "n_snp", "gw_sig", "max_N"]], on="trait", how="left")
    per["is_selected"] = per.trait.isin(SELECTED)
    per["is_zero_power"] = per.trait.isin(ZERO_POWER)
    # The staleness sweep predates the two freeze upgrades, so those traits
    # carry no power data and drop out of the C4 correlations. Say so loudly
    # rather than letting the n quietly shrink.
    missing = per.loc[per.gw_sig.isna(), "trait"].tolist()
    if missing:
        log(f"  NOTE no GWAS power data for {len(missing)} trait(s): "
            f"{', '.join(missing)} -> excluded from C4 correlations")
    return per


def corr_block(df, xcol, ycol, label, subsets):
    """Spearman for each subset, reported as its own row."""
    out = []
    for name, mask in subsets:
        d = df.loc[mask, [xcol, ycol]].dropna()
        if d.shape[0] < 4:
            out.append(dict(comparison=label, subset=name, n=d.shape[0],
                            rho=np.nan, p=np.nan))
            continue
        r, p = spearmanr(d[xcol], d[ycol])
        out.append(dict(comparison=label, subset=name, n=d.shape[0],
                        rho=r, p=p))
    return out


# =====================================================================
def main():
    import anndata as ad

    log("reading h5ad obs (backed, obs only)")
    A = ad.read_h5ad(H5AD, backed="r")
    obs = A.obs[[ANNOT, "donor", "n_genes", "n_umi"]].copy()
    obs = obs.rename(columns={ANNOT: "domain"})
    log(f"  obs: {obs.shape[0]} spots, {obs.domain.nunique()} domains, "
        f"{obs.donor.nunique()} donors")

    G = load_scdrs_group()
    M = load_gsmap()
    if not os.path.exists(STALENESS_TSV):
        sys.exit(f"FATAL: missing GWAS power table {STALENESS_TSV}\n"
                 f"  (stage gwas_staleness_sweep.tsv into the code dir first)")
    sweep = pd.read_csv(STALENESS_TSV, sep="\t")
    log(f"GWAS power table: {sweep.shape[0]} traits")

    keep_dom = [d for d in G.domain.unique() if d not in DROP_DOMAINS]
    informative = [t for t in G.trait.unique() if t not in ZERO_POWER]
    shared = sorted(set(informative) & set(M.trait.unique()))
    log(f"keep_dom={len(keep_dom)} informative_traits={len(informative)} "
        f"shared_with_gsmap={len(shared)}")

    # ---- C1
    H, Sz, Gz, Sq, Gq = controls_heterogeneity(G, M, keep_dom, shared)
    H.to_csv(f"{OUT}/ctrl_heterogeneity.tsv", sep="\t", index=False)
    log(f"wrote ctrl_heterogeneity.tsv ({H.shape[0]} traits)")

    # ---- C2
    R = controls_donor_ranks(obs, informative)
    R.to_csv(f"{OUT}/ctrl_donor_ranks.tsv", sep="\t", index=False)
    log(f"wrote ctrl_donor_ranks.tsv ({R.shape[0]} donor x domain x trait rows)")

    # ---- C3
    C = controls_domain_covariates(obs, Sz, keep_dom)
    C.to_csv(f"{OUT}/ctrl_domain_covariates.tsv", sep="\t", index=False)
    log(f"wrote ctrl_domain_covariates.tsv ({C.shape[0]} domains)")

    Sc = controls_spot_covariates(obs, informative)
    Sc.to_csv(f"{OUT}/ctrl_spot_covariates.tsv", sep="\t", index=False)
    log(f"wrote ctrl_spot_covariates.tsv ({Sc.shape[0]} traits)")

    # ---- C4
    P = controls_gwas_power(G, sweep, keep_dom)
    P.to_csv(f"{OUT}/ctrl_gwas_power.tsv", sep="\t", index=False)
    log(f"wrote ctrl_gwas_power.tsv ({P.shape[0]} traits)")

    # ---- summary of every statistical test, one row each
    stats = []
    # C1: paired difference in heterogeneity
    from scipy.stats import wilcoxon
    hh = H.dropna(subset=["scdrs_conc_top3", "gsmap_conc_top3"])
    w, pw = wilcoxon(hh.scdrs_conc_top3, hh.gsmap_conc_top3)
    stats.append(dict(comparison="C1 top3 concentration scDRS vs gsMap",
                      subset="all shared traits", n=hh.shape[0],
                      rho=np.nan, p=pw,
                      extra=f"median scDRS {hh.scdrs_conc_top3.median():.3f} "
                            f"vs gsMap {hh.gsmap_conc_top3.median():.3f} (Wilcoxon)"))
    w2, pw2 = wilcoxon(hh.scdrs_n_sig, hh.gsmap_n_sig)
    stats.append(dict(comparison="C1 n significant domains scDRS vs gsMap",
                      subset="all shared traits", n=hh.shape[0],
                      rho=np.nan, p=pw2,
                      extra=f"median scDRS {hh.scdrs_n_sig.median():.1f} "
                            f"vs gsMap {hh.gsmap_n_sig.median():.1f} (Wilcoxon)"))
    # C3 domain-level
    sub_all = [("all domains", C.index == C.index)]
    for y in ["mean_mcz_psychiatric", "mean_mcz_all"]:
        for x in ["n_spots", "median_n_genes", "median_n_umi", "n_donors"]:
            for r in corr_block(C, x, y, f"C3 {y} vs {x}", sub_all):
                stats.append({**r, "extra": ""})
    # C4 power
    subsets = [("all informative traits", ~P.is_zero_power),
               ("non-redundant selected", P.is_selected)]
    for y in ["max_abs_mcz", "spread_mcz"]:
        for x in ["gw_sig", "max_N", "n_snp"]:
            for r in corr_block(P, x, y, f"C4 {y} vs {x}", subsets):
                stats.append({**r, "extra": ""})
    ST = pd.DataFrame(stats)
    ST.to_csv(f"{OUT}/ctrl_stats_summary.tsv", sep="\t", index=False)
    log(f"wrote ctrl_stats_summary.tsv ({ST.shape[0]} tests)")

    print("\n=== ctrl_stats_summary.tsv ===")
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(ST.to_string(index=False))


STALENESS_TSV = f"{PROJ}/code/Visium/16_LDSC/scDRS/gwas_staleness_sweep.tsv"

if __name__ == "__main__":
    main()
