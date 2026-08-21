#!/usr/bin/env python
"""27_build_spot_grid_data.py -- stage per-spot z for a trait x donor spot-plot grid.

  python 27_build_spot_grid_data.py --arm test2 --group Psychiatric

Reads the same two sources as 20_spatial_pdfs.py --
  {PROC}/gsMap/ST/{donor}.h5ad                          obsm['spatial'] coordinates
  {arm_workdir}/{donor}/spatial_ldsc/{donor}_{trait}.csv.gz   per-spot z
-- and writes one tidy table so the figure itself can be plain ggplot2.

COORDINATES are normalised per donor: centred, then divided by the larger of the two
half-ranges. Every donor therefore lands in a common [-1, 1] box with its true aspect
ratio intact, which lets the R side use fixed facet scales + coord_fixed and get panels
of equal physical size. Raw pixel ranges differ ~2x between sections, so plotting them
on shared axes without this would shrink some sections to a corner of their panel.

Y IS NOT PRE-FLIPPED here; the R script applies scale_y_reverse, mirroring the
ax.invert_yaxis() in 20_spatial_pdfs.py. Doing it in one place only avoids a
double flip.

OUTPUT: figdata/spot_grid_{arm}_{group}.csv.gz
  trait, donor, x, y, nlp        (x, y in [-1, 1]; nlp = -log10 P, 3 dp)
plus figdata/spot_grid_{arm}_{group}_meta.csv with the shared colour limit.

The plotted quantity is -log10(P), matching the domain-level heatmaps and the
existing spatial figures, so a spot's colour means the same thing in both. P is
floored at 1e-300 before the log, as elsewhere.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import anndata as ad

CODE = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap"
PROC = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC"
ARM_WORKDIR = {"test1": f"{PROC}/gsMap",
               "test2": f"{PROC}/gsMap_cond_functional",
               "test3": f"{PROC}/gsMap_cond_neuronal"}
FIGDATA = f"{CODE}/figdata"

# one GWAS freeze per phenotype, matching the main figures
PREFERRED = {"Psychiatric": ["SCZ", "BIP_PGC3", "MDD", "ADHD", "PTSD_F3", "Autism", "Anorexia"]}


def coords(donor):
    f = f"{PROC}/gsMap/ST/{donor}.h5ad"
    if not os.path.exists(f):
        sys.exit(f"FATAL: no h5ad at {f}")
    a = ad.read_h5ad(f, backed="r")
    xy = pd.DataFrame(np.asarray(a.obsm["spatial"])[:, :2], index=a.obs_names,
                      columns=["x", "y"])
    # centre, then scale by the LARGER half-range so aspect ratio survives
    for c in ("x", "y"):
        xy[c] = xy[c] - (xy[c].max() + xy[c].min()) / 2.0
    half = max(xy.x.abs().max(), xy.y.abs().max())
    return xy / half


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="test2", choices=list(ARM_WORKDIR))
    ap.add_argument("--group", default="Psychiatric")
    ap.add_argument("--traits", nargs="*", default=None,
                    help="explicit trait list, overriding --group")
    a = ap.parse_args()

    W = ARM_WORKDIR[a.arm]
    donors = [l.strip() for l in open(f"{CODE}/samples.txt") if l.strip()]

    if a.traits:
        traits = a.traits
    elif a.group in PREFERRED:
        traits = PREFERRED[a.group]
    else:
        g = pd.read_csv(f"{FIGDATA}/trait_groups.csv", index_col=0)["group"]
        traits = sorted(g.index[g == a.group])
    if not traits:
        sys.exit(f"FATAL: no traits for group {a.group!r}")

    xy_by_donor = {d: coords(d) for d in donors}
    rows = []
    for t in traits:
        for d in donors:
            f = f"{W}/{d}/spatial_ldsc/{d}_{t}.csv.gz"
            if not os.path.exists(f):
                sys.exit(f"FATAL: missing per-spot file {f}")
            sp = pd.read_csv(f, index_col=0)          # not `d` -- that is the donor
            z = -np.log10(sp["p"].clip(lower=1e-300))
            g = xy_by_donor[d]
            common = g.index.intersection(z.index)
            if len(common) < 0.95 * len(z):
                sys.exit(f"FATAL: only {len(common)} of {len(z)} spot ids matched "
                         f"for {d}/{t}; index convention differs")
            rows.append(pd.DataFrame({"trait": t, "donor": d,
                                      "x": g.loc[common, "x"].round(4),
                                      "y": g.loc[common, "y"].round(4),
                                      "nlp": z.loc[common].round(3).values}))
    long = pd.concat(rows, ignore_index=True)

    # Upper limit shared across every panel so one trait's colour is comparable to
    # another's. 99.5th percentile rather than the max: a handful of runaway spots
    # would otherwise wash out every panel. The LOWER end is fixed at 0 -- this is a
    # one-sided significance scale, not a diverging one.
    vlim = float(np.nanpercentile(long.nlp.values, 99.5))

    tag = a.group if not a.traits else "custom"
    out = f"{FIGDATA}/spot_grid_{a.arm}_{tag}.csv.gz"
    long.to_csv(out, index=False, compression="gzip")
    meta = pd.DataFrame([{"arm": a.arm, "group": tag, "vlim": round(vlim, 4),
                          "n_traits": long.trait.nunique(),
                          "n_donors": long.donor.nunique(),
                          "n_spots": len(long),
                          "trait_order": "|".join(traits),
                          "donor_order": "|".join(donors),
                          "pct_clamped": round(100 * float((long.nlp > vlim).mean()), 3),
                          "pct_sig": round(100 * float((long.nlp > -np.log10(0.05)).mean()), 3)}])
    meta.to_csv(f"{FIGDATA}/spot_grid_{a.arm}_{tag}_meta.csv", index=False)
    print(f"wrote {out}")
    print(f"  {long.trait.nunique()} traits x {long.donor.nunique()} donors, "
          f"{len(long)} spot-rows, vlim={vlim:.3f} "
          f"({100 * float((long.nlp > vlim).mean()):.2f}% clamped, "
          f"{100 * float((long.nlp > -np.log10(0.05)).mean()):.1f}% at P<0.05)")


if __name__ == "__main__":
    main()
