#!/usr/bin/env python
"""24_add_trait_power_stats.py -- add one trait to figdata/chi2_single.json and
figdata/trait_groups.csv, the two inputs the power-vs-enrichment figure needs
beyond the Cauchy matrix.

  python 24_add_trait_power_stats.py --trait PTSD --validate          # check first
  python 24_add_trait_power_stats.py --trait PTSD_F3 --group Psychiatric

Stats are computed straight from the LDSC-format sumstats the gsMap run consumed,
so they describe the exact file that produced the enrichment values:
  nsnp      rows with finite Z
  medN      median of the N column
  mean_chi2 mean of Z^2            <- the x-axis of the figure
  lambda_gc median(Z^2) / 0.4549 (rounded chi2_1 median; see CHI2_MEDIAN_NULL)

--validate recomputes an EXISTING entry and diffs, which is the only check that
this matches however the original numbers were produced.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd

CODE = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap"
GWAS = ("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/"
        "processed-data/Visium/16_LDSC/gsMap/GWAS")
CONFIG = f"{CODE}/gwas_config_full.yaml"   # trait label -> sumstats path
FIGDATA = f"{CODE}/figdata"
# The staged values were produced with the 4-dp rounding 0.4549, not the exact
# chi2.ppf(0.5, 1) = 0.4549364231. Kept as-is so recomputed rows match the existing
# ones exactly; the difference is ~8e-05 in lambda_gc and affects nothing plotted.
CHI2_MEDIAN_NULL = 0.4549


def sumstats_path(trait):
    """Trait label -> sumstats file. Filenames in GWAS/ are the ORIGINAL upstream
    names (adhd.gz, height_ldscore.gz), so the label mapping must come from the
    same config gsMap itself was run with -- guessing f"{trait}.gz" is wrong for
    all but a handful of traits."""
    for line in open(CONFIG):
        line = line.split("#")[0].strip()
        if not line or ":" not in line:
            continue
        k, v = line.split(":", 1)
        if k.strip() == trait and v.strip():
            return v.strip()
    direct = f"{GWAS}/{trait}.gz"          # newly added traits may not be in the config yet
    if os.path.exists(direct):
        return direct
    sys.exit(f"FATAL: {trait} not in {CONFIG} and no {direct}")


def stats_for(trait):
    f = sumstats_path(trait)
    if not os.path.exists(f):
        sys.exit(f"FATAL: no sumstats at {f}")
    # header is whitespace-padded in several of these files, not strictly tab-separated
    d = pd.read_csv(f, sep=r"\s+", usecols=["N", "Z"])
    z = pd.to_numeric(d["Z"], errors="coerce")
    n = pd.to_numeric(d["N"], errors="coerce")
    ok = np.isfinite(z)
    chi2 = z[ok] ** 2
    return dict(trait=trait,
                nsnp=int(ok.sum()),
                medN=float(np.median(n[ok])),
                mean_chi2=float(chi2.mean()),
                lambda_gc=float(np.median(chi2) / CHI2_MEDIAN_NULL))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--group", default=None,
                    help="trait_groups.csv label, e.g. Psychiatric (required unless --validate)")
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()

    new = stats_for(a.trait)
    print(f"{a.trait}: nsnp={new['nsnp']} medN={new['medN']:.1f} "
          f"mean_chi2={new['mean_chi2']:.4f} lambda_gc={new['lambda_gc']:.4f}")

    jpath = f"{FIGDATA}/chi2_single.json"
    rows = json.load(open(jpath))

    if a.validate:
        old = next((r for r in rows if r["trait"] == a.trait), None)
        if old is None:
            sys.exit(f"FATAL: {a.trait} not in {jpath}; nothing to validate against")
        for k in ["nsnp", "medN", "mean_chi2", "lambda_gc"]:
            o, nv = old[k], new[k]
            rel = abs(nv - o) / abs(o) if o else abs(nv - o)
            print(f"  {k:10s} staged={o!r:>22} recomputed={nv!r:>22} rel_diff={rel:.3g}")
        return

    if a.group is None:
        sys.exit("FATAL: --group is required when appending")

    rows = [r for r in rows if r["trait"] != a.trait] + [new]
    json.dump(rows, open(jpath, "w"), indent=1)
    print(f"wrote {jpath}  ({len(rows)} entries)")

    gpath = f"{FIGDATA}/trait_groups.csv"
    g = pd.read_csv(gpath, index_col=0)
    g.loc[a.trait] = a.group
    g.to_csv(gpath)
    print(f"wrote {gpath}  ({len(g)} traits)  {a.trait} -> {a.group}")


if __name__ == "__main__":
    main()
