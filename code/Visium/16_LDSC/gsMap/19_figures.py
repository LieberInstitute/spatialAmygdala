#!/usr/bin/env python3
"""
19_figures.py -- regenerate the gsMap domain-enrichment figure set.

Reads the staged summary tables in figdata/ and writes PNGs to figures/.
Self-contained: no cluster access, no artifact store, no kernel state.

Usage
-----
    python 19_figures.py                 # all figures
    python 19_figures.py cond_P          # one figure by key
    python 19_figures.py --list          # show keys

Inputs (figdata/)
-----------------
    test1_cauchy_matrix.csv    40 traits x 15 domains, Cauchy P, default baseline (arm 0)
    test2_cauchy_matrix.csv    same, functional-conditioned arm (arm A / "test 2")
    domain_vs_rest_OR.csv      Mantel-Haenszel domain-vs-rest OR, both arms, long format
    trait_groups.csv           trait -> phenotype group
    neuro_pct_all.csv          per-donor, per-domain neuronal-marker % of counts
    gwas_selection_table.csv   25 candidate GWAS covering the 10 phenotypes that had more
                               than one, with lambda_GC and mean chi2; selected==True marks
                               the 10 winners. NOT the analysis subset on its own -- the
                               15 single-GWAS phenotypes are absent and load() adds them.

Provenance of the two derived statistics
---------------------------------------
P     Cauchy combination across spots within a donor, then across donors (gsMap's own
      run_cauchy_combination). Absolute enrichment: "is this domain enriched at all?"
OR    Per donor, spots are called trait-associated at FDR<0.05, then a 2x2 table
      (in-domain x called) is built and pooled across donors by Mantel-Haenszel with a
      Haldane-Anscombe 0.5 correction; variance is Robins-Breslow-Greenland. Relative
      preference: "are called spots concentrated HERE versus the rest of the amygdala?"
      Computed upstream; this script only plots the saved estimates.

Caveats encoded in the figures
------------------------------
  * Traits whose FDR call saturates (<2% or >98% of spots) have undefined OR and are
    drawn as white cells. Across the 25-trait subset that is six traits in the conditional
    arm (Alzheimer_v3, Anorexia, Autism, Height, Stroke_2022_Any, T2D) and two in the
    default (Stroke_2022_Any, T2D).
  * CLA is present in one donor only; AI and LA in six of seven. Marked "n/7" in teal.
  * The P colour scale is anchored at P=0.05 (white), so the grey/orange boundary is the
    significance threshold rather than a midpoint of the data range.
"""

import argparse
import json
import os
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from scipy.stats import linregress, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "figdata")
OUT = os.path.join(HERE, "figures")

# ---------------------------------------------------------------- style

FS_TITLE, FS_LABEL, FS_TICK, FS_NOTE = 8, 6.5, 6, 5.5
THR = -np.log10(0.05)          # 1.301 -- the white point of the P scale

GCOL = {"Psychiatric": "#3B6FB6", "Substance use": "#4E9F50",
        "Cognitive/behavioral": "#8E6FBF", "Neurological": "#D97C2B",
        "Non-brain control": "#6E6E6E"}
GLAB = {"Psychiatric": "Psychiatric", "Substance use": "Substance use",
        "Cognitive/behavioral": "Cognitive/ behavioral", "Neurological": "Neurological",
        "Non-brain control": "Non-brain control"}

# light grey below threshold -> white at P=0.05 -> orange-red above
CMAP_P = LinearSegmentedColormap.from_list(
    "lightgrey_white_orangered",
    ["#D2D2D2", "#E4E4E4", "#FFFFFF", "#FDBE85", "#F16913", "#B30000"])

GRIDC = "#BFBFBF"       # cell border / spine colour
TEALC = "#00707A"       # reduced donor-coverage marks
SIGC = "#B2182B"        # forest, significant
NSC = "#AAAAAA"         # forest, not significant


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300,
        "font.size": FS_TICK, "axes.titlesize": FS_TITLE, "axes.labelsize": FS_LABEL,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": None, "figure.autolayout": False,
    })


# ---------------------------------------------------------------- data

def load():
    """Load every staged table and derive the shared row/column orderings."""
    d = {}
    d["P1"] = pd.read_csv(f"{DATA}/test1_cauchy_matrix.csv", index_col=0)
    d["P2"] = pd.read_csv(f"{DATA}/test2_cauchy_matrix.csv", index_col=0)
    d["OR"] = pd.read_csv(f"{DATA}/domain_vs_rest_OR.csv")
    d["grp"] = pd.read_csv(f"{DATA}/trait_groups.csv", index_col=0)["group"]
    d["gsel"] = pd.read_csv(f"{DATA}/gwas_selection_table.csv")
    npc = pd.read_csv(f"{DATA}/neuro_pct_all.csv")

    # The 25-trait best-powered subset = 10 winners from the replicate sets (chosen on
    # genome-wide mean chi2; gwas_selection_table.csv documents that choice) + the 15
    # phenotypes that had only one GWAS and so needed no selection.
    sel = d["gsel"]
    picked = sel.loc[sel["selected"].astype(bool), "trait"].tolist()
    singles = ["ADHD", "Autism", "Anorexia", "PTSD", "Insomnia", "Neuroticism",
               "EduYears", "Intelligence", "PD", "GSCAN_AgeSmk", "GSCAN_CigDay",
               "GSCAN_SmkCes", "Height", "BMI", "T2D"]
    traits = picked + [t for t in singles if t not in picked]
    have = [t for t in traits if t in d["P1"].index and t in d["P2"].index]
    if len(have) != len(traits):
        print(f"    WARNING: {sorted(set(traits) - set(have))} absent from a cauchy matrix")
    assert len(have) >= 20, f"expected ~25 subset traits, got {len(have)}"
    d["traits"] = have

    DOM = list(d["P2"].columns)
    d["DOM"] = DOM
    d["ndon"] = npc.groupby("domain")["sample"].nunique().reindex(DOM).fillna(0).astype(int).to_dict()
    d["npct"] = npc.groupby("domain")["pct_neuronal"].mean().reindex(DOM)

    # -log10 P, conditional arm, on the subset
    d["N1"] = -np.log10(d["P1"].loc[traits, DOM].clip(lower=1e-300))
    d["N2"] = -np.log10(d["P2"].loc[traits, DOM].clip(lower=1e-300))

    # orderings come from the CONDITIONAL arm: strongest mean signal first.
    # (Recomputing per arm matters -- the default-arm ordering differs.)
    d["torder"] = d["N2"].mean(axis=1).sort_values(ascending=False).index.tolist()
    d["dorder"] = d["N2"].mean(axis=0).sort_values(ascending=False).index.tolist()

    # OR matrices, saturated traits masked to NaN (undefined, not zero)
    for arm, key in [("test1", "O1"), ("test2", "O2")]:
        sub = d["OR"][d["OR"].arm == arm]
        M = sub.pivot(index="trait", columns="domain", values="log2OR")
        sat = sub.groupby("trait")["saturated"].any()
        M.loc[[t for t in sat[sat].index if t in M.index]] = np.nan
        d[key] = M.reindex(index=traits, columns=DOM)

    # GWAS power. Two files are needed and it is easy to get this wrong:
    # chi2_full.json covers the 25 CANDIDATE GWAS of the 10 replicate phenotypes, while
    # chi2_single.json covers the 15 single-GWAS phenotypes. Both are 25-and-15, and
    # gwas_selection_table.csv is also 25 rows -- using it alone silently yields the wrong
    # 25 traits (it lacks Height/PTSD/ADHD/... entirely).
    pw_full = pd.DataFrame(json.load(open(f"{DATA}/chi2_full.json"))).set_index("trait")
    pw_sing = pd.DataFrame(json.load(open(f"{DATA}/chi2_single.json"))).set_index("trait")
    POW = pd.concat([pw_full.loc[[t for t in picked if t in pw_full.index]], pw_sing])
    missing = [t for t in have if t not in POW.index]
    assert not missing, f"no power metrics for {missing}"
    d["POW"] = POW.loc[have]
    # mean genome-wide chi2 is the power axis, NOT median N: MDD has 685k samples at
    # chi2 1.85 while SCZ has 126k at 2.03, and N correlates with enrichment far weaker
    # (rho +0.37, P 0.07) than chi2 does (+0.77, P 7e-06).
    d["power"] = d["POW"]["mean_chi2"]
    return d


# ---------------------------------------------------------------- helpers

def white_cells(ax, M):
    """Undefined (saturated) cells: white with a thin border, never grey fill --
    grey would read as a low value on the P scale."""
    ys, xs = np.where(~np.isfinite(M.values))
    for yy, xx in zip(ys, xs):
        ax.add_patch(plt.Rectangle((xx - .5, yy - .5), 1, 1, facecolor="white",
                                   edgecolor=GRIDC, lw=0.25, zorder=3))


def heat_axes(ax, M, dorder, torder, grp, ndon, ylabels=True):
    ax.set_xticks(range(len(torder)))
    ax.set_xticklabels(torder, rotation=90, fontsize=FS_TICK)
    for tl in ax.get_xticklabels():
        tl.set_color(GCOL.get(grp.get(tl.get_text()), "#333"))
    ax.set_yticks(range(len(dorder)))
    if ylabels:
        ax.set_yticklabels([f"{d}  {ndon[d]}/7" if ndon[d] < 7 else d for d in dorder],
                           fontsize=FS_LABEL)
        for tl in ax.get_yticklabels():
            if "/7" in tl.get_text():
                tl.set_color(TEALC)
    else:
        ax.tick_params(labelleft=False)
    ax.tick_params(length=1.6, pad=1.6)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_edgecolor(GRIDC); sp.set_linewidth(0.5)


def group_legend(fig, x, yfig):
    hs = [plt.Line2D([], [], color=GCOL[g], lw=2.6, label=GLAB[g]) for g in GCOL]
    fig.legend(handles=hs, loc="lower left", bbox_to_anchor=(x, yfig), ncol=5,
               frameon=False, fontsize=FS_NOTE, handlelength=1.1, columnspacing=1.1)


def check_overlaps(fig, name):
    """Fail loudly on colliding text -- adjacent tick labels are exempt."""
    FigureCanvasAgg(fig); fig.canvas.draw(); r = fig.canvas.get_renderer()
    tt = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    ticks = set()
    for ax in fig.axes:
        ticks |= set(ax.get_xticklabels()) | set(ax.get_yticklabels())
    ov = [(a.get_text()[:16], b.get_text()[:16])
          for i, (a, ba) in enumerate(tt) for b, bb in tt[i + 1:]
          if ba.overlaps(bb) and not (a in ticks and b in ticks)]
    print(f"    overlaps: {len(ov)}" + (f"  {ov[:3]}" if ov else ""))
    return ov


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, f"{name}.png")
    fig.savefig(p, dpi=300)
    print(f"    wrote {p}  ({os.path.getsize(p):,} B)")
    plt.close(fig)


# ---------------------------------------------------------------- figures

def fig_cond_P(d):
    """Standalone conditional-arm P heatmap. Square cells; margins sized to the
    rendered label extents rather than guessed."""
    apply_style()
    N2, dorder, torder = d["N2"], d["dorder"], d["torder"]
    M = N2.reindex(index=torder, columns=dorder).T          # domains x traits
    nr, nc = M.shape
    norm = TwoSlopeNorm(vmin=0, vcenter=THR, vmax=float(np.ceil(np.nanmax(M.values))))

    # probe pass: measure how much room the labels actually need
    probe = plt.figure(figsize=(6, 5)); pax = probe.add_subplot(111)
    pax.set_xticks(range(nc)); pax.set_xticklabels(torder, rotation=90, fontsize=FS_TICK)
    pax.set_yticks(range(nr)); pax.set_yticklabels(
        [f"{x}  {d['ndon'][x]}/7" if d["ndon"][x] < 7 else x for x in dorder], fontsize=FS_LABEL)
    FigureCanvasAgg(probe); probe.canvas.draw(); pr = probe.canvas.get_renderer()
    need_b = max(t.get_window_extent(pr).height for t in pax.get_xticklabels()) / probe.dpi
    need_l = max(t.get_window_extent(pr).width for t in pax.get_yticklabels()) / probe.dpi
    plt.close(probe)

    cell = 0.20
    w_ax, h_ax = nc * cell, nr * cell
    L, R, T, B = need_l + 0.10, 0.60, 0.34, need_b + 0.52
    # the title, legend and caption are wider than a narrow grid; widen so nothing clips
    w_ax = max(w_ax, 4.35 - L)
    fig = plt.figure(figsize=(L + w_ax + R, T + h_ax + B))
    W, H = fig.get_size_inches()
    ax = fig.add_axes([L / W, B / H, w_ax / W, h_ax / H])
    im = ax.imshow(M.values, cmap=CMAP_P, norm=norm, aspect="equal")
    white_cells(ax, M)
    heat_axes(ax, M, dorder, torder, d["grp"], d["ndon"])

    cax = fig.add_axes([(L + w_ax + 0.10) / W, B / H, 0.10 / W, h_ax / H])
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([0, THR, 5, 10, 15]); cb.set_ticklabels(["0", "1.3", "5", "10", "15"])
    cb.ax.tick_params(labelsize=FS_NOTE, length=1.6, pad=1.4)
    cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(GRIDC)
    fig.text((L + w_ax + 0.15) / W, (B + h_ax + 0.05) / H, r"$-\log_{10}P$",
             fontsize=FS_LABEL, ha="center", va="bottom")
    fig.text(L / W, (B + h_ax + 0.07) / H,
             "Domain enrichment under functional conditioning",
             fontsize=FS_TITLE, ha="left", va="bottom")
    group_legend(fig, L / W, 0.30 / H)
    fig.text(L / W, 0.09 / H,
             "Cauchy-combined P per domain, 7 donors. White = P 0.05; grey below, "
             "orange-red above. Teal n/7: domain in fewer than 7 donors.",
             fontsize=FS_NOTE, color="#555", ha="left")
    check_overlaps(fig, "cond_P")
    return fig


def fig_cond_3panel(d):
    """P heatmap | SCZ forest | OR heatmap.

    The forest is the reading key: it is literally one column of panel c drawn as
    estimates with 95% CIs. Rows and columns are IDENTICAL across a and c -- asserted
    below, not assumed, because the two matrices are built from different sources.
    """
    apply_style()
    dorder, torder, grp, ndon = d["dorder"], d["torder"], d["grp"], d["ndon"]
    PD = d["N2"].reindex(index=torder, columns=dorder).T
    OD = d["O2"].reindex(index=torder, columns=dorder).T
    forest = (d["OR"][(d["OR"].arm == "test2") & (d["OR"].trait == "SCZ")]
              .set_index("domain").reindex(dorder))

    assert list(PD.index) == list(OD.index) == list(forest.index) == dorder, "row axes differ"
    assert list(PD.columns) == list(OD.columns) == torder, "column axes differ"

    normP = TwoSlopeNorm(vmin=0, vcenter=THR, vmax=float(np.ceil(np.nanmax(PD.values))))
    olim = float(np.ceil(np.nanpercentile(np.abs(OD.values[np.isfinite(OD.values)]), 98)))

    fig = plt.figure(figsize=(7.4, 5.6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[3.0, 0.95, 3.0], wspace=0.085,
                           left=0.105, right=0.925, top=0.855, bottom=0.245)
    aP, aF, aO = (fig.add_subplot(gs[0]),
                  fig.add_subplot(gs[1]),
                  fig.add_subplot(gs[2]))
    aF.sharey(aP); aO.sharey(aP)

    iP = aP.imshow(PD.values, cmap=CMAP_P, norm=normP, aspect="auto")
    white_cells(aP, PD)
    aP.set_title(r"Conditional enrichment  $-\log_{10}P$", loc="left", pad=5, fontsize=7)

    y = np.arange(len(dorder))
    sig = forest["fdr"].values < 0.05
    aF.axvline(0, color="#888", lw=0.8, ls="--", zorder=1)
    for m, col in [(sig, SIGC), (~sig, NSC)]:
        if m.any():
            aF.errorbar(forest["log2OR"].values[m], y[m], xerr=1.96 * forest["se"].values[m],
                        fmt="o", ms=3.0, lw=0, elinewidth=0.85, capsize=1.5, color=col, zorder=3)
    aF.set_title("SCZ", loc="left", pad=5, fontsize=7)
    aF.set_xlabel(r"$\log_2$OR", fontsize=FS_LABEL)
    aF.tick_params(labelsize=FS_NOTE)
    aF.tick_params(labelleft=False)   # sharey would otherwise echo panel a's domain labels

    iO = aO.imshow(OD.values, cmap="PRGn", vmin=-olim, vmax=olim, aspect="auto")
    white_cells(aO, OD)
    aO.set_title(r"Spatial specificity  $\log_2$OR", loc="left", pad=5, fontsize=7)

    heat_axes(aP, PD, dorder, torder, grp, ndon, ylabels=True)
    heat_axes(aO, OD, dorder, torder, grp, ndon, ylabels=False)
    aP.set_ylim(len(dorder) - .5, -.5)

    c1 = fig.colorbar(iP, ax=aP, fraction=0.026, pad=0.014)
    c1.set_ticks([0, THR, 5, 10, 15]); c1.set_ticklabels(["0", "1.3", "5", "10", "15"])
    c2 = fig.colorbar(iO, ax=aO, fraction=0.026, pad=0.014)
    for cb in (c1, c2):
        cb.ax.tick_params(labelsize=FS_NOTE, length=1.6, pad=1.4)
        cb.outline.set_linewidth(0.5); cb.outline.set_edgecolor(GRIDC)

    ytop = max(a.get_position().y1 for a in (aP, aF, aO))
    for lab, a in zip("abc", (aP, aF, aO)):
        fig.text(a.get_position().x0 - 0.020, ytop + 0.085, lab,
                 fontsize=9, fontweight="bold", ha="left", va="top")
    group_legend(fig, 0.105, 0.048)
    fig.text(0.105, 0.012,
             "Functional-conditioned arm, 7 donors. Rows (domains) and columns (traits) are "
             "identical in a and c; b plots the SCZ column of c with 95% CI.\nWhite cells: OR "
             "undefined (trait saturated at FDR<0.05). Teal n/7: domain in fewer than 7 donors.",
             fontsize=FS_NOTE, color="#555", ha="left", linespacing=1.5)

    # colourbar labels must be placed AFTER a draw -- positions are stale before it
    fig.canvas.draw()
    for cb, lab in [(c1, r"$-\log_{10}P$"), (c2, r"$\log_2$OR")]:
        b = cb.ax.get_position()
        fig.text(b.x0 + b.width / 2, b.y1 + 0.012, lab,
                 ha="center", va="bottom", fontsize=FS_NOTE)
    check_overlaps(fig, "cond_3panel")
    return fig


def fig_OR_scatter(d):
    """Default vs conditional log2 OR.

    The headline is the INTERCEPT, not the slope: conditioning applies a near-constant
    downward shift to every domain rather than redistributing signal among them. The
    right panel therefore plots the residual after removing that uniform offset -- the
    raw per-domain medians are all negative and would just show the offset 15 times.
    """
    apply_style()
    O1, O2, dorder = d["O1"], d["O2"], d["dorder"]
    J = (O1.stack().rename("d1").to_frame()
         .join(O2.stack().rename("d2"), how="inner").reset_index())
    J.columns = ["trait", "domain", "d1", "d2"]
    J = J.dropna()
    J["group"] = J.trait.map(d["grp"])
    lr = linregress(J.d1, J.d2)

    fig, axs = plt.subplots(1, 2, figsize=(7.9, 3.8),
                            gridspec_kw=dict(wspace=0.32, left=0.080, right=0.985,
                                             top=0.855, bottom=0.205))
    a, b = axs
    lim = (min(J.d1.min(), J.d2.min()) - .6, max(J.d1.max(), J.d2.max()) + .6)
    a.plot(lim, lim, color="#AAA", lw=.85, ls="--", zorder=1)
    for g in GCOL:
        gr = J[J.group == g]
        if len(gr):
            a.scatter(gr.d1, gr.d2, s=8, color=GCOL[g], alpha=.75, lw=0, zorder=3,
                      label=GLAB[g])
    xs = np.linspace(*lim, 50)
    a.plot(xs, lr.intercept + lr.slope * xs, color=SIGC, lw=1.05, zorder=4)
    a.set_xlim(lim); a.set_ylim(lim); a.set_aspect("equal")
    a.set_xlabel(r"$\log_2$OR, default baseline"); a.set_ylabel(r"$\log_2$OR, conditional")
    a.set_title("Conditioning preserves spatial preference", loc="left", pad=6, fontsize=7)
    a.annotate(f"slope {lr.slope:.2f}\nintercept {lr.intercept:+.2f}\n"
               f"r {lr.rvalue:.3f}\nn {len(J)} cells",
               xy=(.04, .96), xycoords="axes fraction", va="top",
               fontsize=FS_NOTE, linespacing=1.3)
    a.legend(loc="lower right", frameon=False, fontsize=5, handlelength=.8,
             labelspacing=.22, borderpad=.15)

    med = (J.assign(delta=J.d2 - J.d1).groupby("domain")["delta"].median().reindex(dorder))
    resid = med - med.median()
    b.barh(range(len(dorder)), resid.values, height=.7, lw=0,
           color=["#1B7837" if v > 0 else "#762A83" for v in resid.values])
    b.axvline(0, color="#666", lw=.9)
    b.set_yticks(range(len(dorder))); b.set_yticklabels(dorder, fontsize=FS_LABEL)
    b.set_ylim(len(dorder) - .5, -.5)
    b.set_xlabel(r"$\Delta\log_2$OR relative to the uniform shift")
    b.set_title("Departure from uniform attenuation", loc="left", pad=6, fontsize=7)
    b.tick_params(labelsize=FS_NOTE)
    fig.text(0.080, 0.045,
             f"{len(J)} domain x trait cells, {J.trait.nunique()} traits with defined OR in both "
             f"arms. Conditioning shifts every domain down by a near-constant\n"
             f"{abs(float(med.median())):.2f} $\\log_2$ units (intercept {lr.intercept:+.2f}); "
             f"bars show the residual after removing that offset.",
             fontsize=FS_NOTE, color="#555", ha="left", linespacing=1.5)
    check_overlaps(fig, "OR_scatter")
    return fig


def fig_two_arm_P(d):
    """Default and conditional P side by side, shared colour scale.

    The scale must span BOTH arms: the default arm reaches -log10 P 26 while the
    conditional tops out near 16, so a vmax taken from the conditional alone
    saturates the default panel.
    """
    apply_style()
    dorder, torder, grp, ndon = d["dorder"], d["torder"], d["grp"], d["ndon"]
    A = d["N1"].reindex(index=torder, columns=dorder).T
    B = d["N2"].reindex(index=torder, columns=dorder).T
    both = np.r_[A.values.ravel(), B.values.ravel()]
    vmax = float(np.ceil(np.nanpercentile(both, 99.5)))
    norm = TwoSlopeNorm(vmin=0, vcenter=THR, vmax=vmax)

    fig, axs = plt.subplots(1, 2, figsize=(7.8, 6.4), sharey=True,
                            gridspec_kw=dict(wspace=0.07, left=0.20, right=0.885,
                                             top=0.885, bottom=0.265))
    for ax, M, ttl in [(axs[0], A, "Default baseline"),
                       (axs[1], B, "+ functional conditioning")]:
        im = ax.imshow(M.values, cmap=CMAP_P, norm=norm, aspect="auto")
        white_cells(ax, M)
        ax.set_title(ttl, loc="left", pad=6)
    heat_axes(axs[0], A, dorder, torder, grp, ndon, ylabels=True)
    heat_axes(axs[1], B, dorder, torder, grp, ndon, ylabels=False)
    axs[0].set_ylim(len(dorder) - .5, -.5)

    cb = fig.colorbar(im, ax=axs.tolist(), fraction=0.030, pad=0.02)
    cb.set_ticks([0, THR, 6, 12, 18, vmax])
    cb.set_ticklabels(["0", "1.3\n(P=0.05)", "6", "12", "18", rf"$\geq${int(vmax)}"])
    cb.ax.tick_params(labelsize=FS_NOTE)
    cb.set_label(r"$-\log_{10}$($P$ value)")
    fig.suptitle("Domain enrichment, traits ordered by mean significance",
                 x=0.045, ha="left", y=0.965)
    group_legend(fig, 0.20, 0.062)
    nd = ", ".join(f"{k} {v}/7" for k, v in sorted(ndon.items()) if v < 7)
    fig.text(0.20, 0.020,
             f"Teal n/7 under a domain: present in fewer than 7 donors ({nd}); those "
             f"columns average over fewer sections.", fontsize=FS_NOTE, color="#555", ha="left")
    check_overlaps(fig, "two_arm_P")
    return fig


def fig_two_arm_OR(d):
    """Default and conditional OR side by side, shared symmetric scale."""
    apply_style()
    dorder, torder, grp, ndon = d["dorder"], d["torder"], d["grp"], d["ndon"]
    A = d["O1"].reindex(index=torder, columns=dorder).T
    B = d["O2"].reindex(index=torder, columns=dorder).T
    both = np.r_[A.values.ravel(), B.values.ravel()]
    olim = float(np.ceil(np.nanpercentile(np.abs(both[np.isfinite(both)]), 98)))

    fig, axs = plt.subplots(1, 2, figsize=(7.8, 6.4), sharey=True,
                            gridspec_kw=dict(wspace=0.07, left=0.20, right=0.885,
                                             top=0.885, bottom=0.265))
    for ax, M, ttl in [(axs[0], A, "Default baseline"),
                       (axs[1], B, "+ functional conditioning")]:
        im = ax.imshow(M.values, cmap="PRGn", vmin=-olim, vmax=olim, aspect="auto")
        white_cells(ax, M)
        ax.set_title(ttl, loc="left", pad=6)
    heat_axes(axs[0], A, dorder, torder, grp, ndon, ylabels=True)
    heat_axes(axs[1], B, dorder, torder, grp, ndon, ylabels=False)
    axs[0].set_ylim(len(dorder) - .5, -.5)

    cb = fig.colorbar(im, ax=axs.tolist(), fraction=0.030, pad=0.02)
    cb.set_label(r"$\log_2$ OR (domain vs rest of amygdala)")
    cb.ax.tick_params(labelsize=FS_NOTE)
    fig.suptitle("Spatial specificity: trait-associated spots in domain vs rest",
                 x=0.045, ha="left", y=0.965)
    group_legend(fig, 0.20, 0.062)
    nd = ", ".join(f"{k} {v}/7" for k, v in sorted(ndon.items()) if v < 7)
    fig.text(0.20, 0.020,
             "Mantel-Haenszel OR pooled across donors (Haldane 0.5); spots called at FDR<0.05. "
             f"White cells: OR undefined. Teal n/7: fewer than 7 donors ({nd}).",
             fontsize=FS_NOTE, color="#555", ha="left")
    check_overlaps(fig, "two_arm_OR")
    return fig


def fig_forest_OR(d):
    """SCZ forest as a reading key, then BOTH arms' OR heatmaps.

    Same idea as cond_3panel but showing default and conditional side by side, so the
    forest belongs to the DEFAULT arm here (it is the boxed column of the middle panel).
    """
    apply_style()
    dorder, torder, grp, ndon = d["dorder"], d["torder"], d["grp"], d["ndon"]
    A = d["O1"].reindex(index=torder, columns=dorder).T
    B = d["O2"].reindex(index=torder, columns=dorder).T
    forest = (d["OR"][(d["OR"].arm == "test1") & (d["OR"].trait == "SCZ")]
              .set_index("domain").reindex(dorder))
    assert list(A.index) == list(B.index) == list(forest.index) == dorder

    both = np.r_[A.values.ravel(), B.values.ravel()]
    olim = float(np.ceil(np.nanpercentile(np.abs(both[np.isfinite(both)]), 98)))

    fig = plt.figure(figsize=(9.4, 5.8))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.05, 3.0, 3.0], wspace=0.10,
                           left=0.135, right=0.90, top=0.855, bottom=0.275)
    aF, a1, a2 = (fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2]))
    a1.sharey(aF); a2.sharey(aF)

    y = np.arange(len(dorder))
    sig = forest["fdr"].values < 0.05
    aF.axvline(0, color="#888", lw=.8, ls="--", zorder=1)
    for m, col in [(sig, SIGC), (~sig, NSC)]:
        if m.any():
            aF.errorbar(forest["log2OR"].values[m], y[m], xerr=1.96 * forest["se"].values[m],
                        fmt="o", ms=3.2, lw=0, elinewidth=.9, capsize=1.6, color=col, zorder=3)
    aF.set_title("SCZ, default arm", loc="left", pad=5, fontsize=7)
    aF.set_xlabel(r"$\log_2$OR", fontsize=FS_LABEL)
    aF.set_yticks(y)
    aF.set_yticklabels([f"{x}  {ndon[x]}/7" if ndon[x] < 7 else x for x in dorder],
                       fontsize=FS_LABEL)
    for tl in aF.get_yticklabels():
        if "/7" in tl.get_text():
            tl.set_color(TEALC)
    aF.tick_params(labelsize=FS_NOTE)
    for sp in ("top", "right"):
        aF.spines[sp].set_visible(False)

    for ax, M, ttl in [(a1, A, "Default baseline"), (a2, B, "+ functional conditioning")]:
        im = ax.imshow(M.values, cmap="PRGn", vmin=-olim, vmax=olim, aspect="auto")
        white_cells(ax, M)
        ax.set_title(ttl, loc="left", pad=5, fontsize=7)
        heat_axes(ax, M, dorder, torder, grp, ndon, ylabels=False)
    aF.set_ylim(len(dorder) - .5, -.5)

    # box the SCZ column the forest expands
    if "SCZ" in torder:
        j = torder.index("SCZ")
        a1.add_patch(plt.Rectangle((j - .5, -.5), 1, len(dorder), fill=False,
                                   edgecolor=SIGC, lw=1.1, zorder=5))

    cb = fig.colorbar(im, ax=[a1, a2], fraction=0.024, pad=0.015)
    cb.set_label(r"$\log_2$ OR (domain vs rest)", fontsize=FS_NOTE)
    cb.ax.tick_params(labelsize=FS_NOTE)
    hs = [plt.Line2D([], [], marker="o", color=SIGC, lw=0, ms=3.2, label="FDR < 0.05"),
          plt.Line2D([], [], marker="o", color=NSC, lw=0, ms=3.2, label="n.s.")]
    aF.legend(handles=hs, loc="upper left", bbox_to_anchor=(0, -0.115), frameon=False,
              fontsize=5, handlelength=.8, labelspacing=.25)
    ytop = max(a.get_position().y1 for a in (aF, a1, a2))
    for lab, a in zip("abc", (aF, a1, a2)):
        fig.text(a.get_position().x0 - 0.022, ytop + 0.085, lab,
                 fontsize=9, fontweight="bold", ha="left", va="top")
    group_legend(fig, 0.135, 0.048)
    nd = ", ".join(f"{k} {v}/7" for k, v in sorted(ndon.items()) if v < 7)
    fig.text(0.135, 0.012,
             "Left: SCZ column of the default-arm heatmap (boxed), with 95% CI -- every heatmap "
             f"cell is one such estimate. White cells: OR undefined.\nTeal n/7: domain in fewer "
             f"than 7 donors ({nd}); intervals there rest on fewer strata.",
             fontsize=FS_NOTE, color="#555", ha="left", linespacing=1.5)
    check_overlaps(fig, "forest_OR")
    return fig


def fig_power(d):
    """GWAS power vs domain enrichment, functional-conditioned arm only.

    One panel. The question it answers: is the enrichment we report just a readout of how
    well-powered each GWAS is? Power predicts the MAGNITUDE of -log10P strongly, so the
    answer for absolute significance is partly yes -- which is why this plot exists and
    why the OR panels, being scale-free, carry the spatial claim instead.
    """
    apply_style()
    tr = d["traits"]
    x = d["power"].loc[tr]
    y = d["N2"].loc[tr].mean(axis=1)          # conditional arm only
    grp, DOM = d["grp"], d["DOM"]
    rho, p = spearmanr(x, y)

    HL = ["SCZ", "MDD", "PTSD", "Height"]
    missing = [t for t in HL if t not in tr]
    assert not missing, f"highlight traits absent from subset: {missing}"

    fig, ax = plt.subplots(figsize=(4.7, 4.3))
    # non-highlighted traits recede; the four requested ones are focal (figure-style 4.2)
    rest = [t for t in tr if t not in HL]
    ax.scatter(x[rest], y[rest], c=[GCOL.get(grp.get(t), "#999") for t in rest],
               s=20, lw=0, alpha=0.45, zorder=2)
    ax.scatter(x[HL], y[HL], c=[GCOL.get(grp.get(t), "#999") for t in HL],
               s=58, lw=1.0, edgecolor="white", zorder=4)

    ax.set_xscale("log")
    ax.set_xticks([1, 2, 4, 8]); ax.set_xticklabels(["1", "2", "4", "8"])
    ax.set_xlabel(r"GWAS power   (mean genome-wide $\chi^2$)")
    ax.set_ylabel(r"mean $-\log_{10}P$ across the 15 domains")
    ax.set_title(f"Enrichment magnitude tracks GWAS power\n"
                 f"Spearman $\\rho$ = {rho:+.2f}, P = {p:.0e}, n = {len(tr)}", loc="left")
    ax.margins(0.11)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # leader lines placed away from the local point cloud, then verified below
    OFF = {"SCZ": (16, 6), "MDD": (16, -4), "PTSD": (16, 9), "Height": (-14, 14)}
    for t in HL:
        dx, dy = OFF[t]
        ax.annotate(t, (x[t], y[t]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=FS_LABEL, fontweight="bold",
                    color=GCOL.get(grp.get(t), "#333"),
                    ha="right" if dx < 0 else "left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.7, color="#888",
                                    shrinkA=0, shrinkB=3), zorder=5)

    seen = [g for g in GCOL if g in set(grp.get(t) for t in tr)]
    hs = [plt.Line2D([], [], marker="o", lw=0, ms=4, color=GCOL[g], label=GLAB[g])
          for g in seen]
    ax.legend(handles=hs, loc="upper left", frameon=False, fontsize=FS_NOTE,
              handlelength=.8, labelspacing=.28, borderpad=.2,
              bbox_to_anchor=(0.005, 0.985))
    fig.text(0.012, 0.035,
             f"Functional-conditioned arm (test 2), {len(tr)}-trait best-powered\n"
             f"subset. Colour = trait group. Power explains magnitude but NOT\n"
             f"which domains rank highest: 0 of {len(DOM)} domains significant\n"
             f"after Bonferroni (see power_vs_domain_shape.csv).",
             fontsize=FS_NOTE, color="#555", ha="left", linespacing=1.5)
    fig.subplots_adjust(left=0.145, right=0.985, top=0.855, bottom=0.285)

    # verify each leader ends nearer its own point than any other (figure-style 6.9)
    fig.canvas.draw()
    for t in HL:
        px, py = ax.transData.transform((x[t], y[t]))
        near = min(tr, key=lambda o: (ax.transData.transform((x[o], y[o]))[0] - px) ** 2
                   + (ax.transData.transform((x[o], y[o]))[1] - py) ** 2)
        assert near == t, f"leader for {t} lands nearest {near}"
    check_overlaps(fig, "power")
    return fig


FIGS = {
    "cond_P": (fig_cond_P, "fig_conditional_P"),
    "forest_OR": (fig_forest_OR, "fig_forest_OR"),
    "power": (fig_power, "fig_power_vs_enrichment"),
    "cond_3panel": (fig_cond_3panel, "fig_conditional_3panel"),
    "OR_scatter": (fig_OR_scatter, "fig_OR_scatter"),
    "two_arm_P": (fig_two_arm_P, "fig_tall_P"),
    "two_arm_OR": (fig_two_arm_OR, "fig_heatmaps_OR"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("keys", nargs="*", help="figure keys (default: all)")
    ap.add_argument("--list", action="store_true", help="list keys and exit")
    args = ap.parse_args()

    if args.list:
        for k, (_, fn) in FIGS.items():
            print(f"  {k:14s} -> figures/{fn}.png")
        return

    missing = [f for f in os.listdir(DATA)] if os.path.isdir(DATA) else []
    if not missing:
        sys.exit(f"no input tables in {DATA}/ -- see the docstring for the expected files")

    d = load()
    print(f"loaded: {len(d['traits'])} traits x {len(d['DOM'])} domains")
    print(f"trait order (conditional arm, strongest first): {', '.join(d['torder'][:4])} ...")
    print(f"domain order: {', '.join(d['dorder'][:4])} ...")

    keys = args.keys or list(FIGS)
    bad = [k for k in keys if k not in FIGS]
    if bad:
        sys.exit(f"unknown key(s) {bad}; use --list")
    for k in keys:
        fn, out = FIGS[k]
        print(f"\n[{k}]")
        save(fn(d), out)


if __name__ == "__main__":
    main()
