#!/usr/bin/env python
"""
10_controls_figures.py -- Figures for the control analyses.

Reads only the TSVs written by 09_controls_compute.py, so this is cheap to
re-run while adjusting figure aesthetics. No cluster resources needed beyond
a couple of minutes on a compute node.

Outputs into $WORKDIR/controls/:
    fig_ctrl_heterogeneity.png     scDRS vs gsMap domain heterogeneity
    fig_ctrl_donor_consistency.png BM / LA ranks across donors
    fig_ctrl_confounds.png         domain size, library size, GWAS power

ALL STYLE KNOBS ARE IN THE "STYLE" BLOCK BELOW -- edit and re-run.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

PROJ    = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala"
WORKDIR = f"{PROJ}/processed-data/Visium/16_LDSC/scDRS"
OUT     = f"{WORKDIR}/controls"

# ============================================================ STYLE
DPI        = 300
FS_TICK    = 13
FS_AXIS    = 15
FS_TITLE   = 17
FS_SUP     = 20
FS_ANNOT   = 12
PT_SIZE    = 90          # scatter marker area
PT_BIG     = 170         # highlighted markers
GREY       = "#D5D5D5"
GREY_EDGE  = "#AEAEAE"
C_SCDRS    = "#4C72B0"
C_GSMAP    = "#DD8452"
C_BM       = "#55A868"
C_LA       = "#4C72B0"
C_CLA      = "#C44E52"
C_NULL     = "#BAB0AC"
HIGHLIGHT  = {"BM": C_BM, "LA": C_LA, "CLA": C_CLA}
PSY_COLOR  = "#8172B3"
CTRL_COLOR = "#BAB0AC"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "axes.labelsize": FS_AXIS,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})
# ===================================================================


def sub(ax, text):
    """Italic grey subtitle above an axes."""
    ax.text(0, 1.030, text, transform=ax.transAxes, fontsize=FS_ANNOT + 1,
            color="#666", style="italic", va="bottom")


def panel_title(ax, letter, text):
    ax.set_title(f"{letter}   {text}", loc="left", pad=28, fontsize=FS_TITLE)


# =====================================================================
def fig_heterogeneity(H):
    """scDRS resolves more domain-to-domain structure than gsMap."""
    from scipy.stats import wilcoxon
    h = H.dropna(subset=["scdrs_conc_top3", "gsmap_conc_top3"]).copy()
    fig = plt.figure(figsize=(15.5, 6.4))
    gs = fig.add_gridspec(1, 3, wspace=0.34, left=0.055, right=0.985,
                          top=0.780, bottom=0.155)

    # a: paired concentration
    ax = fig.add_subplot(gs[0, 0])
    for _, r in h.iterrows():
        ax.plot([0, 1], [r.gsmap_conc_top3, r.scdrs_conc_top3],
                color=GREY, lw=0.9, zorder=2)
    ax.scatter(np.zeros(len(h)), h.gsmap_conc_top3, s=PT_SIZE, c=C_GSMAP,
               edgecolor="white", lw=0.8, zorder=4, label="gsMap")
    ax.scatter(np.ones(len(h)), h.scdrs_conc_top3, s=PT_SIZE, c=C_SCDRS,
               edgecolor="white", lw=0.8, zorder=4, label="scDRS")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["gsMap", "scDRS"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("signal concentration\n(top-3 domains / all domains)")
    _, p_conc = wilcoxon(h.scdrs_conc_top3, h.gsmap_conc_top3)
    panel_title(ax, "a", "Signal concentration")
    sub(ax, f"n = {len(h)} traits · Wilcoxon p = {p_conc:.1e}")
    ax.set_box_aspect(1)

    # b: number of significant domains
    ax = fig.add_subplot(gs[0, 1])
    mx = int(max(h.scdrs_n_sig.max(), h.gsmap_n_sig.max())) + 1
    ax.plot([0, mx], [0, mx], color="#BBB", lw=1.0, ls="--", zorder=1)
    ax.scatter(h.gsmap_n_sig, h.scdrs_n_sig, s=PT_SIZE, c=GREY,
               edgecolor=GREY_EDGE, lw=0.7, zorder=3)
    psy = h[h.is_psychiatric]
    ax.scatter(psy.gsmap_n_sig, psy.scdrs_n_sig, s=PT_BIG, c=PSY_COLOR,
               edgecolor="white", lw=1.2, zorder=5, label="psychiatric")
    ax.set_xlabel("gsMap  domains at FDR < 0.05")
    ax.set_ylabel("scDRS  domains at FDR < 0.05")
    ax.legend(frameon=False, fontsize=FS_ANNOT, loc="upper left")
    _, p_nsig = wilcoxon(h.scdrs_n_sig, h.gsmap_n_sig)
    panel_title(ax, "b", "Domains called significant")
    sub(ax, f"median {h.gsmap_n_sig.median():.0f} (gsMap) vs "
            f"{h.scdrs_n_sig.median():.0f} (scDRS) · p = {p_nsig:.1e}")
    ax.set_box_aspect(1)

    # c: spread comparison
    ax = fig.add_subplot(gs[0, 2])
    ax.scatter(h.gsmap_spread, h.scdrs_spread, s=PT_SIZE, c=GREY,
               edgecolor=GREY_EDGE, lw=0.7, zorder=3)
    ax.scatter(psy.gsmap_spread, psy.scdrs_spread, s=PT_BIG, c=PSY_COLOR,
               edgecolor="white", lw=1.2, zorder=5)
    rho, p = spearmanr(h.gsmap_spread, h.scdrs_spread)
    ax.set_xlabel("gsMap  domain spread (z)")
    ax.set_ylabel("scDRS  domain spread (mcz)")
    panel_title(ax, "c", "Domain spread, method vs method")
    sub(ax, f"Spearman ρ = {rho:+.3f}, p = {p:.1e}")
    ax.set_box_aspect(1)

    # Headline is DERIVED from the tests, not asserted. Direction and
    # significance both come from the paired comparisons computed above.
    alpha = 0.05
    d_conc = h.scdrs_conc_top3.median() - h.gsmap_conc_top3.median()
    d_nsig = h.scdrs_n_sig.median() - h.gsmap_n_sig.median()
    # Word the headline for the metric actually tested (top-3 concentration),
    # and only generalize to "more heterogeneity" when the independent
    # domain-count test agrees in direction. Otherwise report both.
    conc_sig = p_conc < alpha
    nsig_sig = p_nsig < alpha
    if conc_sig and d_conc > 0 and nsig_sig and d_nsig < 0:
        head = ("scDRS concentrates disease signal in fewer domains "
                "than gsMap")
    elif conc_sig and d_conc > 0:
        head = ("scDRS signal is more concentrated in the top domains "
                "than gsMap")
    elif conc_sig and d_conc < 0:
        head = ("gsMap signal is more concentrated in the top domains "
                "than scDRS")
    else:
        head = ("Top-domain signal concentration does not differ detectably "
                "between scDRS and gsMap")
    if conc_sig and nsig_sig and (d_conc > 0) == (d_nsig > 0):
        # both tests significant but pointing the same way on both axes:
        # more concentrated AND more domains called -> flag, don't smooth over
        head += " (note: domain counts differ in the same direction)"
    fig.suptitle(head, fontsize=FS_SUP, x=0.055, ha="left", y=0.950)
    print(f"  [C1] concentration: median diff {d_conc:+.3f}, p = {p_conc:.2e}"
          f" | n_sig median diff {d_nsig:+.1f}, p = {p_nsig:.2e}")
    f = f"{OUT}/fig_ctrl_heterogeneity.png"
    fig.savefig(f, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {f}")


# =====================================================================
def fig_donor_consistency(R):
    """BM and LA rank near the top in every donor, for psychiatric traits."""
    psy = R[R.is_psychiatric].copy()
    fig = plt.figure(figsize=(15.5, 6.4))
    gs = fig.add_gridspec(1, 3, wspace=0.36, left=0.058, right=0.985,
                          top=0.780, bottom=0.230)

    # a: mean rank per domain (lower = more enriched)
    ax = fig.add_subplot(gs[0, 0])
    mr = (psy.groupby("domain")["rank_in_donor"]
             .agg(["mean", "std", "size"]).sort_values("mean"))
    ypos = np.arange(len(mr))
    cols = [HIGHLIGHT.get(d, GREY) for d in mr.index]
    ax.barh(ypos, mr["mean"], color=cols, edgecolor="none", height=0.72)
    ax.set_yticks(ypos); ax.set_yticklabels(mr.index)
    ax.invert_yaxis()
    ax.set_xlabel("mean rank across donor × trait\n(1 = most enriched)")
    panel_title(ax, "a", "Mean domain rank")
    sub(ax, f"{psy.trait.nunique()} psychiatric traits × "
            f"{psy.donor.nunique()} donors")
    ax.set_box_aspect(1)

    # b: per-donor rank of BM / LA / CLA
    ax = fig.add_subplot(gs[0, 1])
    donors = sorted(psy.donor.unique())
    xpos = np.arange(len(donors))
    for dom, col in HIGHLIGHT.items():
        sel = psy[psy.domain == dom]
        if sel.empty:
            continue
        med = [sel[sel.donor == d]["rank_in_donor"].median() for d in donors]
        ax.plot(xpos, med, "-o", color=col, lw=2.0, ms=9,
                label=dom, zorder=4)
    ax.set_xticks(xpos)
    ax.set_xticklabels(donors, rotation=90)
    ax.set_ylabel("median rank (1 = most enriched)")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=FS_ANNOT, loc="lower right")
    panel_title(ax, "b", "Consistency across donors")
    sub(ax, "median over psychiatric traits, within donor")
    ax.set_box_aspect(1)

    # c: how often each domain lands in the top 3
    ax = fig.add_subplot(gs[0, 2])
    top3 = (psy.assign(in_top3=psy.rank_in_donor <= 3)
               .groupby("domain")["in_top3"].mean().sort_values(ascending=False))
    ypos = np.arange(len(top3))
    cols = [HIGHLIGHT.get(d, GREY) for d in top3.index]
    ax.barh(ypos, 100 * top3.values, color=cols, edgecolor="none", height=0.72)
    ax.set_yticks(ypos); ax.set_yticklabels(top3.index)
    ax.invert_yaxis()
    ax.set_xlabel("% of donor × trait combinations\nin the top 3 domains")
    panel_title(ax, "c", "Top-3 frequency")
    sub(ax, "CLA present in one donor only")
    ax.set_box_aspect(1)

    # Headline derived: state the actual per-donor worst case for BM and LA
    # rather than asserting "every donor".
    bits = []
    for dom in ["BM", "LA"]:
        s = psy[psy.domain == dom]
        if s.empty:
            continue
        per_donor_med = s.groupby("donor")["rank_in_donor"].median()
        worst_med = per_donor_med.max()          # no truncation
        worst_any = int(s.rank_in_donor.max())   # worst single trait x donor
        ndon = s.donor.nunique()
        # State the median bound as a median, and carry the true worst case.
        bits.append(f"{dom} median rank ≤ {worst_med:.1f} in all {ndon} donors")
        print(f"  [C2] {dom}: worst per-donor median = {worst_med:.1f}; "
              f"worst single trait x donor rank = {worst_any}; "
              f"n_donors = {ndon}; top-3 rate "
              f"{100 * (s.rank_in_donor <= 3).mean():.0f}%")
    head = "; ".join(bits) if bits else "Per-donor domain ranks"
    fig.suptitle(head, fontsize=FS_SUP, x=0.058, ha="left", y=0.950)
    f = f"{OUT}/fig_ctrl_donor_consistency.png"
    fig.savefig(f, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {f}")


# =====================================================================
def fig_confounds(C, S, P):
    """Scores track neither domain size, library size, nor GWAS power."""
    fig = plt.figure(figsize=(15.5, 10.6))
    gs = fig.add_gridspec(2, 3, wspace=0.36, hspace=0.52, left=0.058,
                          right=0.985, top=0.885, bottom=0.075)

    def scat_dom(ax, xcol, xlabel, letter, title, logx=False):
        cols = [HIGHLIGHT.get(d, GREY) for d in C.domain]
        ax.scatter(C[xcol], C.mean_mcz_psychiatric, s=PT_BIG, c=cols,
                   edgecolor="white", lw=1.1, zorder=4)
        for _, r in C.iterrows():
            if r.domain in HIGHLIGHT:
                ax.annotate(r.domain, (r[xcol], r.mean_mcz_psychiatric),
                            textcoords="offset points", xytext=(10, 4),
                            fontsize=FS_ANNOT, color=HIGHLIGHT[r.domain],
                            fontweight="bold")
        if logx:
            ax.set_xscale("log")
            # default log minor ticks collide badly at this panel width
            from matplotlib.ticker import LogLocator, NullFormatter
            ax.xaxis.set_major_locator(LogLocator(base=10, numticks=4))
            ax.xaxis.set_minor_formatter(NullFormatter())
        rho, p = spearmanr(C[xcol], C.mean_mcz_psychiatric)
        ax.axhline(0, color="#CCC", lw=0.8, zorder=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("mean assoc_mcz\n(psychiatric traits)")
        panel_title(ax, letter, title)
        sub(ax, f"ρ = {rho:+.3f}, p = {p:.3f}  (n = {len(C)} domains)")
        ax.set_box_aspect(1)

    scat_dom(fig.add_subplot(gs[0, 0]), "n_spots", "spots per domain",
             "a", "Domain size", logx=True)
    scat_dom(fig.add_subplot(gs[0, 1]), "median_n_umi", "median UMI per spot",
             "b", "Library size")
    scat_dom(fig.add_subplot(gs[0, 2]), "median_n_genes",
             "median genes per spot", "c", "Detected genes")

    # d: per-spot correlation with library size, per trait
    ax = fig.add_subplot(gs[1, 0])
    s = S.sort_values("rho_norm_score_vs_n_umi")
    ypos = np.arange(len(s))
    ax.barh(ypos, s.rho_norm_score_vs_n_umi, color=GREY,
            edgecolor=GREY_EDGE, lw=0.5, height=0.8)
    ax.axvline(0, color="#333", lw=1.0)
    ax.set_yticks([]); ax.set_ylabel(f"{len(s)} traits")
    ax.set_xlabel("Spearman ρ:\nspot score vs spot UMI")
    # Scale to the data, with a reference band at |rho| = 0.1 so the reader
    # sees how far below any meaningful effect these sit. A fixed +/-0.5 axis
    # renders every bar invisible and looks like missing data.
    lim = max(0.02, 1.35 * s.rho_norm_score_vs_n_umi.abs().max())
    ax.set_xlim(-lim, lim)
    panel_title(ax, "d", "Per-spot library size")
    sub(ax, f"max |ρ| = {s.rho_norm_score_vs_n_umi.abs().max():.4f} "
            f"(axis spans ±{lim:.3f})")
    ax.set_box_aspect(1)

    # e/f: GWAS power
    def scat_pow(ax, xcol, xlabel, letter, title):
        p_ = P[~P.is_zero_power].dropna(subset=[xcol, "max_abs_mcz"])
        selm = p_.is_selected
        ax.scatter(p_.loc[~selm, xcol], p_.loc[~selm, "max_abs_mcz"],
                   s=PT_SIZE, c=GREY, edgecolor=GREY_EDGE, lw=0.6,
                   zorder=3, label="redundant freeze")
        ax.scatter(p_.loc[selm, xcol], p_.loc[selm, "max_abs_mcz"],
                   s=PT_BIG, c=C_SCDRS, edgecolor="white", lw=1.1,
                   zorder=5, label="non-redundant")
        ax.set_xscale("symlog")
        r_all, pa = spearmanr(p_[xcol], p_.max_abs_mcz)
        d = p_[selm]
        r_sel, ps = spearmanr(d[xcol], d.max_abs_mcz)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("max |assoc_mcz|")
        ax.legend(frameon=False, fontsize=FS_ANNOT - 1, loc="upper left")
        panel_title(ax, letter, title)
        # two lines: one subtitle per trait set, so neither runs off the panel
        sub(ax, f"all n={len(p_)}: ρ={r_all:+.3f} p={pa:.3f}\n"
                f"selected n={len(d)}: ρ={r_sel:+.3f} p={ps:.3f}")
        # Traits with no power data cannot appear here -- name them on the
        # panel so their absence is visible, not silent.
        nomiss = P.loc[P[xcol].isna() & ~P.is_zero_power, "trait"].tolist()
        if nomiss:
            ax.text(0.985, 0.03, "no power data: " + ", ".join(nomiss),
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=FS_ANNOT - 2, color="#999", style="italic")
        ax.set_box_aspect(1)

    scat_pow(fig.add_subplot(gs[1, 1]), "gw_sig",
             "genome-wide significant SNPs", "e", "GWAS power (loci)")
    scat_pow(fig.add_subplot(gs[1, 2]), "max_N",
             "GWAS sample size (max N)", "f", "GWAS power (N)")

    # Headline derived from the tests actually run. A non-significant rho on
    # ~15 domains / ~35 traits is ABSENCE OF EVIDENCE, not evidence of no
    # association, so the wording says "no detectable" and the panel
    # subtitles carry every rho and n for the reader to judge.
    tests = []
    for xc in ["n_spots", "median_n_umi", "median_n_genes"]:
        r_, p_ = spearmanr(C[xc], C.mean_mcz_psychiatric)
        tests.append((f"domain {xc}", r_, p_, len(C)))
    pw = P[~P.is_zero_power].dropna(subset=["gw_sig", "max_abs_mcz"])
    for xc in ["gw_sig", "max_N"]:
        r_, p_ = spearmanr(pw[xc], pw.max_abs_mcz)
        tests.append((f"GWAS {xc}", r_, p_, len(pw)))
    sel = pw[pw.is_selected]
    for xc in ["gw_sig", "max_N"]:
        r_, p_ = spearmanr(sel[xc], sel.max_abs_mcz)
        tests.append((f"GWAS {xc} (selected)", r_, p_, len(sel)))
    umax = S.rho_norm_score_vs_n_umi.abs().max()

    sig = [t for t in tests if t[2] < 0.05]
    worst = max(tests, key=lambda t: abs(t[1]))
    if not sig:
        head = (f"No detectable association with domain size, library size, "
                f"or GWAS power (all p > 0.05; max |ρ| = {abs(worst[1]):.2f})")
    else:
        names = ", ".join(t[0] for t in sig)
        head = (f"Association detected for: {names} "
                f"(max |ρ| = {abs(worst[1]):.2f})")
    fig.suptitle(head, fontsize=FS_SUP - 2, x=0.058, ha="left", y=0.960)
    print("  [C3/C4] confound tests:")
    for nm, r_, p_, n_ in tests:
        flag = "  <-- p < 0.05" if p_ < 0.05 else ""
        print(f"    {nm:26s} rho={r_:+.3f} p={p_:.4f} n={n_}{flag}")
    print(f"    per-spot |rho| vs UMI, max over traits = {umax:.4f}")
    miss = P.loc[P.gw_sig.isna(), "trait"].tolist()
    if miss:
        print(f"    NOTE excluded from GWAS-power panels (no power data): "
              f"{', '.join(miss)}")
    f = f"{OUT}/fig_ctrl_confounds.png"
    fig.savefig(f, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {f}")


# =====================================================================
def main():
    H = pd.read_csv(f"{OUT}/ctrl_heterogeneity.tsv", sep="\t")
    R = pd.read_csv(f"{OUT}/ctrl_donor_ranks.tsv", sep="\t")
    C = pd.read_csv(f"{OUT}/ctrl_domain_covariates.tsv", sep="\t")
    S = pd.read_csv(f"{OUT}/ctrl_spot_covariates.tsv", sep="\t")
    P = pd.read_csv(f"{OUT}/ctrl_gwas_power.tsv", sep="\t")
    fig_heterogeneity(H)
    fig_donor_consistency(R)
    fig_confounds(C, S, P)


if __name__ == "__main__":
    main()
