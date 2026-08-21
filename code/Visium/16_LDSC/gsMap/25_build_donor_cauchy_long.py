#!/usr/bin/env python
"""25_build_donor_cauchy_long.py -- stage PER-DONOR domain p-values as a long table
for the supplementary donor-replication figure.

  python 25_build_donor_cauchy_long.py --arm test2 --group Psychiatric

Everything else in figdata/ is pooled across the 7 donors (cauchy_across_samples).
This reads the UNPOOLED per-donor output instead --
  {workdir}/{donor}/cauchy_combination/{donor}_{trait}.Cauchy.csv.gz
-- so a figure can show whether a domain's enrichment replicates donor to donor
or rests on one section.

OUTPUT: figdata/donor_cauchy_long_{arm}.csv  with columns
  trait, donor, domain, p_cauchy, p_median, n_donors_with_domain

Domain coverage is NOT uniform: CLA appears only in Br6660, LA is absent from
Br6660, AI is absent from Br9192. Those (trait, donor, domain) rows are simply
absent from the output rather than zero-filled, and n_donors_with_domain lets the
plot mark them as structurally missing instead of non-significant.
"""
import argparse, os, sys
import pandas as pd

CODE = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap"
PD   = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC"
ARM_WORKDIR = {"test1": f"{PD}/gsMap",
               "test2": f"{PD}/gsMap_cond_functional",
               "test3": f"{PD}/gsMap_cond_neuronal"}
FIGDATA = f"{CODE}/figdata"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="test2", choices=list(ARM_WORKDIR))
    ap.add_argument("--group", default="Psychiatric",
                    help="trait_groups.csv group to include, or 'all'")
    a = ap.parse_args()

    W = ARM_WORKDIR[a.arm]
    donors = [l.strip() for l in open(f"{CODE}/samples.txt") if l.strip()]
    groups = pd.read_csv(f"{FIGDATA}/trait_groups.csv", index_col=0)["group"]
    traits = sorted(groups.index) if a.group == "all" else sorted(groups.index[groups == a.group])
    if not traits:
        sys.exit(f"FATAL: no traits in group {a.group!r}")

    rows, missing = [], []
    for t in traits:
        for d in donors:
            f = f"{W}/{d}/cauchy_combination/{d}_{t}.Cauchy.csv.gz"
            if not os.path.exists(f):
                missing.append(f"{d}/{t}")
                continue
            x = pd.read_csv(f)
            x.insert(0, "donor", d)
            x.insert(0, "trait", t)
            rows.append(x.rename(columns={"annotation": "domain"}))
    if missing:
        # a missing per-donor file means that donor-trait unit never ran; unlike a
        # missing DOMAIN it is a gap in the analysis, so refuse rather than paper over
        sys.exit(f"FATAL: {len(missing)} per-donor Cauchy files absent, e.g. {missing[:5]}")

    long = pd.concat(rows, ignore_index=True)
    n_dom = (long.groupby("domain")["donor"].nunique().rename("n_donors_with_domain"))
    long = long.merge(n_dom, on="domain", how="left")

    out = f"{FIGDATA}/donor_cauchy_long_{a.arm}.csv"
    long.to_csv(out, index=False)
    print(f"wrote {out}")
    print(f"  {long.trait.nunique()} traits x {long.donor.nunique()} donors "
          f"x {long.domain.nunique()} domains = {len(long)} rows")
    part = n_dom[n_dom < len(donors)]
    if len(part):
        print("  domains not present in every donor: "
              + ", ".join(f"{k} ({v}/{len(donors)})" for k, v in part.items()))


if __name__ == "__main__":
    main()
