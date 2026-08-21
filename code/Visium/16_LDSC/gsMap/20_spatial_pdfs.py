#!/usr/bin/env python3
"""
20_spatial_pdfs.py -- one PDF per trait: per-spot gsMap z across all donors and arms.

Each PDF is a 4 x 7 grid.
  row 1  anatomical domains (categorical) -- the reference for reading the rows below
  row 2  arm 0, default gsMap baseline
  row 3  arm A, + full functional conditioning
  row 4  arm B, + functional + neuronal axis
  columns  the 7 donors, in samples.txt order

Colour: per-spot z, diverging, symmetric, on a scale SHARED across all 21 z panels of
that trait, so panels are comparable within a PDF but not between PDFs (traits differ in
power by orders of magnitude -- a shared global scale would flatten the weak ones).
Scatter is rasterized, text stays vector, so the PDFs are small enough to email.

Usage
-----
    python 20_spatial_pdfs.py                  # all 25 subset traits
    python 20_spatial_pdfs.py SCZ MDD          # named traits
    python 20_spatial_pdfs.py --all-traits      # all 40, not just the subset
    python 20_spatial_pdfs.py --outdir DIR

Inputs (all on /dcs04, read-only)
---------------------------------
    {WORK}/ST/{donor}.h5ad                              obsm['spatial'] coordinates and
                                                        obs['BS_k16_Semisupervised_wAI']
    {arm_workdir}/{donor}/spatial_ldsc/{donor}_{trait}.csv.gz    per-spot beta/se/z/p

Notes
-----
  * Spot IDs are identical across the three sources (verified: 36,876/36,876 for Br6471).
    The script asserts the intersection is non-empty and reports any shortfall per panel.
  * Visium obsm['spatial'] is pixel space with y increasing downward, so the y axis is
    inverted for anatomical orientation.
  * Domain coverage is computed per PDF from the annotation itself, not hardcoded. As of
    this writing: CLA 1/7 donors, LA 6/7, everything else 7/7. Note AI is ANNOTATED in all
    seven but has only 2 spots in Br9192, so tables that require a minimum spot count
    (e.g. neuro_pct_all.csv, and hence the heatmap captions in 19_figures.py) report AI as
    6/7. Both are correct: annotated in 7, usefully quantified in 6.
  * Absent domains simply do not appear in that donor's reference panel -- sampling, not
    missing data.
"""

import argparse
import os
import sys
from importlib.util import module_from_spec, spec_from_file_location

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
PROC = ("/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/"
        "processed-data/Visium/16_LDSC")
WORK = f"{PROC}/gsMap"
ANNOT = "BS_k16_Semisupervised_wAI"   # obs column gsMap was run with
ARMS = [("arm 0  default baseline", f"{PROC}/gsMap"),
        ("arm A  + functional", f"{PROC}/gsMap_cond_functional"),
        ("arm B  + functional + neuronal axis", f"{PROC}/gsMap_cond_neuronal")]

# 15 domains; colours chosen so the grey-matter block reads warm and the
# non-neuronal block (CHAT/Endothelial/WM) reads cool-grey
DOM_COL = {
    "CLA": "#8C510A", "BM": "#BF812D", "PL": "#DFC27D", "CoA": "#F6E8C3",
    "LA": "#C7522B", "AI": "#E17C51", "HPC": "#F2A578", "BLD": "#D98E5A",
    "MeA": "#9E4A6E", "CeA": "#C4739A", "BL": "#E0A7C0",
    "CHAT": "#80CDC1", "Endothelial": "#35978F", "WM.1": "#01665E", "WM.2": "#003C30",
}


def subset_traits():
    """Reuse 19_figures.load()'s subset definition rather than restating it."""
    spec = spec_from_file_location("figs", os.path.join(HERE, "19_figures.py"))
    m = module_from_spec(spec)
    sys.modules["figs"] = m
    spec.loader.exec_module(m)
    return m.load()["traits"]


def donors():
    p = os.path.join(HERE, "samples.txt")
    return [l.strip() for l in open(p) if l.strip()]


def load_geometry(donor, cache={}):
    """Coordinates + domain labels for one donor. Read once, reused across traits."""
    if donor in cache:
        return cache[donor]
    import anndata as ad
    a = ad.read_h5ad(f"{WORK}/ST/{donor}.h5ad", backed="r")
    xy = pd.DataFrame(np.asarray(a.obsm["spatial"])[:, :2],
                      index=list(a.obs_names), columns=["x", "y"])
    # Domain labels come from the ST h5ad itself -- the same obs column gsMap was run
    # with. (A {donor}_domains.csv exists for Br6471 only, a pilot leftover; reading
    # that would have worked for one donor and failed for the other six.)
    assert ANNOT in a.obs.columns, f"{donor}: {ANNOT} not in obs ({list(a.obs.columns)[:8]}...)"
    dom = a.obs[ANNOT].astype(str).replace({"nan": np.nan})
    g = xy.join(dom.rename("domain"), how="left")
    cache[donor] = g
    return g


def load_z(armdir, donor, trait):
    f = f"{armdir}/{donor}/spatial_ldsc/{donor}_{trait}.csv.gz"
    if not os.path.exists(f):
        return None
    return pd.read_csv(f, index_col=0)["z"]


def panel(ax, g, vals, vlim, cmap="RdBu_r", s=0.9):
    ax.scatter(g.x, g.y, c=vals, cmap=cmap, vmin=-vlim, vmax=vlim,
               s=s, lw=0, marker=".", rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def make_pdf(trait, dons, geo, outdir):
    # gather every panel's z first, so the colour scale can span the whole trait
    Z = {}
    for lab, armdir in ARMS:
        for d in dons:
            z = load_z(armdir, d, trait)
            if z is not None:
                Z[(lab, d)] = z
    if not Z:
        print(f"  {trait}: no per-spot output in any arm -- skipped")
        return None
    allz = np.concatenate([v.values for v in Z.values()])
    vlim = float(np.nanpercentile(np.abs(allz), 99.5))
    vlim = max(vlim, 1e-6)

    nrow, ncol = 1 + len(ARMS), len(dons)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.05 * ncol, 2.25 * nrow))
    axs = np.atleast_2d(axs)

    # row 1: anatomical domains
    present = set()
    for j, d in enumerate(dons):
        g = geo[d]
        cols = g.domain.map(DOM_COL).fillna("#DDDDDD")
        present |= set(g.domain.dropna().unique())
        ax = axs[0, j]
        ax.scatter(g.x, g.y, c=cols, s=0.9, lw=0, marker=".", rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{d}\n{len(g):,} spots", fontsize=6.5, pad=3, y=1.005)
    axs[0, 0].set_ylabel("domains", fontsize=7)

    # rows 2..4: one arm each
    short = {}
    for i, (lab, armdir) in enumerate(ARMS, start=1):
        for j, d in enumerate(dons):
            ax = axs[i, j]
            z = Z.get((lab, d))
            g = geo[d]
            if z is None:
                ax.text(.5, .5, "not run", ha="center", va="center",
                        fontsize=6, color="#999", transform=ax.transAxes)
                ax.set_xticks([]); ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                continue
            common = g.index.intersection(z.index)
            # A zero intersection means the spot-id conventions diverged between the
            # h5ad and the ldsc output -- that renders as an empty panel, which looks
            # like a null result rather than a bug. Fail loudly instead.
            assert len(common) > 0, (
                f"{trait} / {d} / {lab}: no shared spot ids between "
                f"{WORK}/ST/{d}.h5ad ({len(g)} spots) and the ldsc output "
                f"({len(z)} spots); first ids {list(g.index[:2])} vs {list(z.index[:2])}")
            gg = g.loc[common]
            panel(ax, gg, z.loc[common].values, vlim)
            if len(common) < len(g):
                short.setdefault(lab, []).append(f"{d} {len(g)-len(common)}")
        axs[i, 0].set_ylabel(lab.replace("  ", "\n"), fontsize=6.5, linespacing=1.4)

    ndon_dom = {k: sum(k in set(geo[d].domain.dropna()) for d in dons) for k in DOM_COL}
    partial = {k: v for k, v in ndon_dom.items() if 0 < v < len(dons)}

    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                               norm=plt.Normalize(vmin=-vlim, vmax=vlim))
    cb = fig.colorbar(sm, cax=fig.add_axes([0.945, 0.20, 0.008, 0.42]))
    cb.set_label("per-spot z", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)

    hs = [Line2D([], [], marker="o", lw=0, ms=3.2, color=DOM_COL[k],
                 label=k if ndon_dom[k] == len(dons) else f"{k} ({ndon_dom[k]}/{len(dons)})")
          for k in DOM_COL if k in present]
    fig.legend(handles=hs, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=min(len(hs), 8), frameon=False, fontsize=5.5,
               handlelength=.8, columnspacing=1.0)
    fig.suptitle(f"{trait} -- per-spot enrichment across donors and conditioning arms",
                 x=0.012, ha="left", fontsize=9)
    fig.text(0.012, 0.045,
             f"Colour: per-spot z, symmetric scale shared across all {len(Z)} z panels of "
             f"this trait (|z| <= {vlim:.2f} at the 99.5th percentile). Scale is per-trait, "
             f"NOT comparable between PDFs.\nRow 1 shows anatomical domains for reference. "
             f"y axis inverted for anatomical orientation. Domains not in every donor: "
             + (", ".join(f"{k} {v}/{len(dons)}" for k, v in sorted(partial.items()))
                if partial else "none") + "."
             + ("" if not short else "\nSpots missing from the ldsc output: "
                + "; ".join(f"{k.split('  ')[0]}: {', '.join(v)}" for k, v in short.items())),
             fontsize=5.5, color="#555", ha="left", linespacing=1.5)
    fig.subplots_adjust(left=0.055, right=0.935, top=0.885, bottom=0.125,
                        wspace=0.04, hspace=0.18)

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"spatial_{trait}.pdf")
    with PdfPages(out) as pp:
        pp.savefig(fig, dpi=200)
    plt.close(fig)
    print(f"  {trait}: {len(Z)} panels, vlim {vlim:.2f} -> {out} "
          f"({os.path.getsize(out)/1e6:.1f} MB)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("traits", nargs="*")
    ap.add_argument("--all-traits", action="store_true",
                    help="all 40 traits rather than the 25-trait subset")
    ap.add_argument("--outdir", default=os.path.join(HERE, "figures", "spatial"))
    args = ap.parse_args()

    dons = donors()
    if args.traits:
        traits = args.traits
    elif args.all_traits:
        traits = [l.strip() for l in open(os.path.join(HERE, "traits.txt")) if l.strip()]
    else:
        traits = subset_traits()
    print(f"{len(traits)} traits x {len(dons)} donors x {len(ARMS)} arms -> {args.outdir}")

    print("loading geometry")
    geo = {d: load_geometry(d) for d in dons}
    for d in dons:
        n_dom = geo[d].domain.nunique()
        print(f"  {d}: {len(geo[d]):,} spots, {n_dom} domains")

    ok = 0
    for t in traits:
        if make_pdf(t, dons, geo, args.outdir):
            ok += 1
    print(f"\n{ok}/{len(traits)} PDFs written to {args.outdir}")


if __name__ == "__main__":
    main()
