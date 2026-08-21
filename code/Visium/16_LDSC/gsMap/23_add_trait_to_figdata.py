#!/usr/bin/env python
"""23_add_trait_to_figdata.py -- add one trait to the staged figdata tables.

WHY THIS EXISTS: test2_cauchy_matrix.csv and domain_vs_rest_OR.csv were originally
built in-kernel with no saved code. This script reconstructs BOTH computations so a
newly-run trait (e.g. PTSD_F3) can be appended reproducibly.

  python 23_add_trait_to_figdata.py --trait PTSD_F3 --arm test2
  python 23_add_trait_to_figdata.py --trait PTSD --arm test2 --validate

--validate recomputes an EXISTING trait and diffs against the staged rows. Run it
before trusting any appended row: it is the only check that this reimplementation
matches whatever the original in-kernel code did.

OR CONVENTION -- recovered by grid search against the staged PTSD/test2 rows, and
reproduces all 15 domains to max |dlog2OR| = 0.00000 (and again on a second trait):
  1. significance called per donor at BH-FDR < 0.05;
  2. donors that do not contain the domain at all (a+b == 0) are DROPPED, not
     zero-filled -- CLA exists in one donor, AI in a subset;
  3. +0.5 added to ALL FOUR cells of every retained table (not only zero cells);
  4. Mantel-Haenszel pooling, Robins-Breslow-Greenland variance;
  5. se is reported on the LOG2 scale, i.e. se_ln / ln(2), matching log2OR.
Step 3 is what makes sparse domains (CHAT, Endothelial) finite and step 5 is why an
earlier ln-scale attempt looked "1.44x off" -- that factor was 1/ln(2), not a real
discrepancy. Do not change these without re-running --validate.

METHOD (mirrors what the figures assume):
  * Cauchy row: read straight from gsMap's own cauchy_across_samples output.
  * OR: per donor, spots are called significant at BH-FDR < 0.05 on the per-spot
    p-values, giving a 2x2 (in-domain x significant) table; tables are pooled across
    the 7 donors by Mantel-Haenszel, SE by Robins-Breslow-Greenland.
  * saturated=True when a trait calls <2% or >98% of spots significant in the pooled
    data -- the OR is then undefined everywhere for that trait, not just where a cell
    hit zero, so the figures blank the whole trait.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm, false_discovery_control

CODE = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap"
PD   = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC"
ARM_WORKDIR = {"test1": f"{PD}/gsMap",
               "test2": f"{PD}/gsMap_cond_functional",
               "test3": f"{PD}/gsMap_cond_neuronal"}
FIGDATA = f"{CODE}/figdata"
ANNOT_COL = "BS_k16_Semisupervised_wAI"


def _domains(sample):
    """Spot -> domain label. Only Br6471 has a staged *_domains.csv; for the rest
    the labels come straight from the h5ad obs column the whole pipeline used."""
    csv = f"{PD}/gsMap/{sample}_domains.csv"
    if os.path.exists(csv):
        return pd.read_csv(csv, index_col=0)[ANNOT_COL]
    import anndata as ad
    h5 = f"{PD}/gsMap/ST/{sample}.h5ad"
    if not os.path.exists(h5):
        sys.exit(f"FATAL: no domain source for {sample} (tried {csv} and {h5})")
    obs = ad.read_h5ad(h5, backed="r").obs
    if ANNOT_COL not in obs.columns:
        sys.exit(f"FATAL: {h5} obs lacks {ANNOT_COL}")
    return obs[ANNOT_COL].astype(str)


def mh_or(tabs):
    """Mantel-Haenszel pooled OR + RBG standard error (natural-log scale).
    tabs: list of (a,b,c,d). Callers apply the +0.5 correction before this."""
    num = den = 0.0
    P_ = Q_ = R_ = S_ = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
        P_ += (a + d) / n
        Q_ += (b + c) / n
        R_ += a * d / n
        S_ += b * c / n
    if den == 0 or num == 0:
        return np.nan, np.nan
    orr = num / den
    # Robins-Breslow-Greenland variance of log(OR)
    v = 0.0
    if R_ > 0:
        v += sum((a + d) / n * (a * d / n) for a, b, c, d in tabs
                 for n in [a + b + c + d] if n) / (2 * R_ ** 2)
    if R_ > 0 and S_ > 0:
        v += sum(((a + d) / n * (b * c / n) + (b + c) / n * (a * d / n))
                 for a, b, c, d in tabs for n in [a + b + c + d] if n) / (2 * R_ * S_)
    if S_ > 0:
        v += sum((b + c) / n * (b * c / n) for a, b, c, d in tabs
                 for n in [a + b + c + d] if n) / (2 * S_ ** 2)
    return orr, float(np.sqrt(v)) if v > 0 else np.nan


def load_spots(arm, trait, samples):
    """Per-spot p-values + domain label + donor, concatenated across donors."""
    W = ARM_WORKDIR[arm]
    frames = []
    for s in samples:
        f = f"{W}/{s}/spatial_ldsc/{s}_{trait}.csv.gz"
        if not os.path.exists(f):
            sys.exit(f"FATAL: missing per-spot file {f}")
        d = pd.read_csv(f, usecols=["spot", "p"])
        dom = _domains(s)
        d["domain"] = d["spot"].map(dom)
        d = d.dropna(subset=["domain"])
        # significance is called WITHIN donor, so an unusually strong section cannot
        # drag the calls in the others
        d["sig"] = false_discovery_control(d["p"].to_numpy(), method="bh") < 0.05
        d["sample"] = s
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def compute_or(spots, arm, trait):
    rows = []
    frac_sig = float(spots["sig"].mean())
    saturated = bool(frac_sig < 0.02 or frac_sig > 0.98)
    for dom in sorted(spots["domain"].unique()):
        tabs, nzero = [], 0
        for _, g in spots.groupby("sample", sort=True):
            ind = g["domain"].eq(dom)
            a = int((ind & g["sig"]).sum());  b = int((ind & ~g["sig"]).sum())
            cc = int((~ind & g["sig"]).sum()); dd = int((~ind & ~g["sig"]).sum())
            if a + b == 0:
                continue          # donor lacks this domain entirely -> contributes nothing
            if 0 in (a, b, cc, dd):
                nzero += 1
            tabs.append((a + 0.5, b + 0.5, cc + 0.5, dd + 0.5))
        orr, se_ln = mh_or(tabs) if tabs else (np.nan, np.nan)
        log2or = np.log2(orr) if orr and orr > 0 else np.nan
        se = se_ln / np.log(2) if se_ln and not np.isnan(se_ln) else np.nan
        if se_ln and not np.isnan(se_ln):
            p = 2 * norm.sf(abs(np.log(orr) / se_ln))
        else:
            p = np.nan
        rows.append(dict(arm=arm, trait=trait, domain=dom, OR=orr, log2OR=log2or,
                         se=se, p=p, fdr=np.nan,
                         spots_in=float(spots["domain"].eq(dom).sum()),
                         n_zero_cell=float(nzero), saturated=saturated))
    out = pd.DataFrame(rows)
    ok = out["p"].notna()
    if ok.any():
        out.loc[ok, "fdr"] = false_discovery_control(out.loc[ok, "p"].to_numpy(), method="bh")
    return out, frac_sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--arm", default="test2", choices=list(ARM_WORKDIR))
    ap.add_argument("--cauchy-only", action="store_true",
                    help="append only the Cauchy row, leaving domain_vs_rest_OR.csv alone")
    ap.add_argument("--validate", action="store_true",
                    help="recompute an existing trait and diff, writing nothing")
    a = ap.parse_args()

    if not a.cauchy_only:
        samples = [l.strip() for l in open(f"{CODE}/samples.txt") if l.strip()]
        spots = load_spots(a.arm, a.trait, samples)
        new_or, frac_sig = compute_or(spots, a.arm, a.trait)
        print(f"{a.trait} [{a.arm}]: {len(spots)} spots, {frac_sig:.3%} significant, "
              f"saturated={bool(new_or['saturated'].iloc[0])}")

    or_path = f"{FIGDATA}/domain_vs_rest_OR.csv"

    if a.validate:
        cur = pd.read_csv(or_path)
        old = cur[(cur.arm == a.arm) & (cur.trait == a.trait)].set_index("domain")
        if old.empty:
            sys.exit(f"FATAL: no staged rows for {a.trait}/{a.arm} to validate against")
        cmp = new_or.set_index("domain").join(old, rsuffix="_staged")
        for col in ["log2OR", "se", "spots_in"]:
            d = (cmp[col] - cmp[f"{col}_staged"]).abs()
            print(f"  {col:9s} max abs diff = {d.max():.6g}")
        return

    # ---- append OR rows (idempotent: drop any existing rows for this arm/trait)
    if not a.cauchy_only:
        cur = pd.read_csv(or_path)
        cur = cur[~((cur.arm == a.arm) & (cur.trait == a.trait))]
        pd.concat([cur, new_or], ignore_index=True).to_csv(or_path, index=False)
        print(f"wrote {or_path}  (+{len(new_or)} rows)")

    # ---- append Cauchy row from gsMap's own across-sample output
    cfile = f"{ARM_WORKDIR[a.arm]}/cauchy_across_samples/{a.trait}_cauchy.csv.gz"
    if not os.path.exists(cfile):
        sys.exit(f"FATAL: missing {cfile}")
    cau = pd.read_csv(cfile).set_index("annotation")["p_cauchy"]
    mpath = f"{FIGDATA}/{a.arm}_cauchy_matrix.csv"
    M = pd.read_csv(mpath, index_col=0)
    missing = [c for c in M.columns if c not in cau.index]
    if missing:
        sys.exit(f"FATAL: cauchy output lacks domains {missing}")
    M.loc[a.trait] = [float(cau[c]) for c in M.columns]
    M.to_csv(mpath)
    print(f"wrote {mpath}  (now {len(M)} traits)")


if __name__ == "__main__":
    main()
