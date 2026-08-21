#!/usr/bin/env python
"""
02a_sumstats_to_pval.py -- convert an LDSC-format sumstats file to the
SNP / P / N table MAGMA's --pval expects.

    python 02a_sumstats_to_pval.py <in.gz> <out.tsv>

WHY THIS EXISTS
    The shared sumstats in /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/
    are LDSC .sumstats files: they carry a signed Z but no P column, because
    LDSC never needs one. MAGMA's gene analysis is driven by SNP-level P
    values, so we derive P = 2 * Phi(-|Z|) exactly as LDSC's own munge step
    would have done in reverse.

    Two header shapes occur in that directory:
        SNP  N  Z  A1  A2                  (most files)
        SNP  A1  A2  N  CHISQ  Z           (the PASS_* / UKB_460K files)
    Both are handled; the column names are matched case-insensitively.

    A THIRD format is auto-detected: the PGC sumstats-VCF used by newer
    freezes (e.g. ptsd2024 / PGC-PTSD Freeze 3, downloaded into
    $WORKDIR/gwas_new/). Those files begin with '##' meta lines and carry

        #CHROM  ID  POS  A1  A2  FREQ  NEFF  Z  P  DIRE

    They already have a P column, so no Z -> P conversion happens; ID maps to
    SNP and NEFF (effective sample size, the correct quantity for a
    case/control meta-analysis) maps to N. Detection is on the leading '##',
    so no flag is needed -- pass the path and the right branch runs.

NUMERICAL NOTE
    Two-sided P from Z is computed in log space (scipy.stats.norm.logsf) and
    exponentiated, so |Z| up to ~38 stays representable. Anything that still
    underflows is floored at 1e-300 -- MAGMA rejects P == 0. The number of
    floored SNPs is reported; for these sumstats (max |Z| ~ 9.6 for SCZ) it
    should be zero.

    Rows with a missing/non-finite Z or N are dropped and counted.

Writes a headed, tab-delimited file:  SNP  P  N
"""
import gzip
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

P_FLOOR = 1e-300


def _is_pgc_vcf(src):
    """True if src is a PGC sumstats-VCF (## meta lines, then a #CHROM header).

    The PGC distributes newer freezes (e.g. ptsd2024 / Freeze 3) in this format
    rather than LDSC .sumstats. It already carries a P column, so no Z -> P
    conversion is needed -- but the columns are named differently:

        #CHROM  ID  POS  A1  A2  FREQ  NEFF  Z  P  DIRE

    ID -> SNP, P -> P, NEFF -> N.  NEFF (effective sample size) is the right
    quantity for MAGMA's ncol= on a case/control meta-analysis.
    """
    op = gzip.open if str(src).endswith(".gz") else open
    with op(src, "rt") as fh:
        first = fh.readline()
    return first.startswith("##")


def _read_pgc_vcf(src):
    """Return a DataFrame with SNP / P / N columns from a PGC sumstats-VCF."""
    op = gzip.open if str(src).endswith(".gz") else open
    with op(src, "rt") as fh:
        n_meta = 0
        for line in fh:
            if line.startswith("##"):
                n_meta += 1
            else:
                break
    df = pd.read_csv(src, sep="\t", skiprows=n_meta, engine="c",
                     na_values=["NA", "."], low_memory=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    cols = {c.upper(): c for c in df.columns}
    for need in ("ID", "P"):
        if need not in cols:
            raise SystemExit(f"{src}: PGC-VCF missing '{need}'; found {list(df.columns)}")
    ncol = "NEFF" if "NEFF" in cols else ("N" if "N" in cols else None)
    if ncol is None:
        raise SystemExit(f"{src}: PGC-VCF has neither NEFF nor N; found {list(df.columns)}")
    return pd.DataFrame({
        "SNP": df[cols["ID"]].astype(str),
        "P": pd.to_numeric(df[cols["P"]], errors="coerce"),
        "N": pd.to_numeric(df[cols[ncol]], errors="coerce"),
    }), ncol


def main_vcf(src, dst):
    raw, ncol = _read_pgc_vcf(src)
    n_in = len(raw)
    keep = np.isfinite(raw["P"].to_numpy(dtype=float)) & np.isfinite(raw["N"].to_numpy(dtype=float))
    n_drop = int((~keep).sum())
    out = raw[keep].copy()

    n_floor = int((out["P"] < P_FLOOR).sum())
    out["P"] = out["P"].clip(lower=P_FLOOR, upper=1.0)

    n_before = len(out)
    out = out.drop_duplicates(subset="SNP", keep="first")
    out.to_csv(dst, sep="\t", index=False, float_format="%.6g")

    n_gw = int((out["P"] < 5e-8).sum())
    print(
        f"{src}  [PGC sumstats-VCF]\n"
        f"  input rows      : {n_in}\n"
        f"  dropped (NA P/N): {n_drop}\n"
        f"  duplicate SNPs  : {n_before - len(out)}\n"
        f"  written         : {len(out)}\n"
        f"  min P           : {out['P'].min():.3g}\n"
        f"  GW-sig (P<5e-8) : {n_gw}\n"
        f"  floored at 1e-300: {n_floor}\n"
        f"  N column used   : {ncol}\n"
        f"  median N        : {np.median(out['N']):.0f}\n"
        f"  -> {dst}",
        flush=True,
    )


def main(src, dst):
    if _is_pgc_vcf(src):
        return main_vcf(src, dst)

    df = pd.read_csv(src, sep=r"\s+", engine="c")
    cols = {c.upper(): c for c in df.columns}

    if "SNP" not in cols:
        raise SystemExit(f"{src}: no 'SNP' column; found {list(df.columns)}")

    # --- sample size -----------------------------------------------------
    # Preference order matters. Ricopili "daner" files (PGC bip2024 and other
    # newer freezes) carry no N column; they report Neff_half, which is HALF
    # the effective sample size summed across cohorts. Using it raw would
    # understate N two-fold, so it is doubled. Nca/Nco are a last resort:
    # pooling case/control totals across cohorts overstates Neff whenever the
    # case:control ratio varies between them.
    if "N" in cols:
        n = pd.to_numeric(df[cols["N"]], errors="coerce").to_numpy(dtype=float)
        n_src = "N"
    elif "NEFF" in cols:
        n = pd.to_numeric(df[cols["NEFF"]], errors="coerce").to_numpy(dtype=float)
        n_src = "NEFF"
    elif "NEFF_HALF" in cols:
        n = 2.0 * pd.to_numeric(df[cols["NEFF_HALF"]], errors="coerce").to_numpy(dtype=float)
        n_src = "2 x NEFF_HALF"
    elif "NCA" in cols and "NCO" in cols:
        nca = pd.to_numeric(df[cols["NCA"]], errors="coerce").to_numpy(dtype=float)
        nco = pd.to_numeric(df[cols["NCO"]], errors="coerce").to_numpy(dtype=float)
        n = 4.0 / (1.0 / nca + 1.0 / nco)
        n_src = "4/(1/Nca+1/Nco)"
    else:
        raise SystemExit(f"{src}: no N / Neff / Neff_half / Nca+Nco; found {list(df.columns)}")

    # --- test statistic --------------------------------------------------
    # A file that already carries P (daner format) is used as-is; deriving P
    # from OR/SE would only reintroduce rounding.
    if "Z" not in cols and "CHISQ" not in cols and "P" in cols:
        p_direct = pd.to_numeric(df[cols["P"]], errors="coerce").to_numpy(dtype=float)
        snp = df[cols["SNP"]].astype(str).to_numpy()
        keep = np.isfinite(p_direct) & np.isfinite(n) & (p_direct > 0)
        n_drop = int((~keep).sum())
        p = np.clip(p_direct[keep], P_FLOOR, 1.0)
        out = pd.DataFrame({"SNP": snp[keep], "P": p, "N": n[keep]})
        n_before = len(out)
        out = out.drop_duplicates(subset="SNP", keep="first")
        out.to_csv(dst, sep="\t", index=False, float_format="%.6g")
        print(
            f"{src}  [P column used directly]\n"
            f"  input rows      : {len(df)}\n"
            f"  dropped (NA/<=0): {n_drop}\n"
            f"  duplicate SNPs  : {n_before - len(out)}\n"
            f"  written         : {len(out)}\n"
            f"  min P           : {out['P'].min():.3g}\n"
            f"  GW-sig (P<5e-8) : {int((out['P'] < 5e-8).sum())}\n"
            f"  N column used   : {n_src}\n"
            f"  median N        : {np.median(out['N']):.0f}\n"
            f"  -> {dst}",
            flush=True,
        )
        return

    if "Z" in cols:
        z = pd.to_numeric(df[cols["Z"]], errors="coerce").to_numpy(dtype=float)
    elif "CHISQ" in cols:
        z = np.sqrt(pd.to_numeric(df[cols["CHISQ"]], errors="coerce").to_numpy(dtype=float))
    else:
        raise SystemExit(f"{src}: neither Z nor CHISQ; found {list(df.columns)}")

    snp = df[cols["SNP"]].astype(str).to_numpy()

    keep = np.isfinite(z) & np.isfinite(n)
    n_drop = int((~keep).sum())

    # two-sided p in log space, then exponentiate
    logp = np.log(2.0) + norm.logsf(np.abs(z[keep]))
    p = np.exp(logp)
    n_floor = int((p < P_FLOOR).sum())
    p = np.clip(p, P_FLOOR, 1.0)

    out = pd.DataFrame({"SNP": snp[keep], "P": p, "N": n[keep]})
    out = out.drop_duplicates(subset="SNP", keep="first")
    out.to_csv(dst, sep="\t", index=False, float_format="%.6g")

    print(
        f"{src}\n"
        f"  input rows      : {len(df)}\n"
        f"  dropped (NA Z/N): {n_drop}\n"
        f"  duplicate SNPs  : {len(snp[keep]) - len(out)}\n"
        f"  written         : {len(out)}\n"
        f"  max |Z|         : {np.nanmax(np.abs(z[keep])):.3f}\n"
        f"  min P           : {out['P'].min():.3g}\n"
        f"  floored at 1e-300: {n_floor}\n"
        f"  median N        : {np.median(out['N']):.0f}\n"
        f"  -> {dst}",
        flush=True,
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
