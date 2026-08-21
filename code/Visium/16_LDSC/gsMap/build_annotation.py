
import os, sys, gzip, json, glob
import numpy as np, pandas as pd
import pyranges as pr
from gsMap.generate_ldscore import load_gtf, overlaps_gtf_bim

# cache the 1.1GB GTF read so load_gtf can be called twice cheaply
import pyranges as _pr
_gtf_cache = {}
_orig_read_gtf = _pr.read_gtf
def _cached_read_gtf(f, *a, **k):
    if f not in _gtf_cache:
        _gtf_cache[f] = _orig_read_gtf(f, *a, **k)
    return _gtf_cache[f].copy()
_pr.read_gtf = _cached_read_gtf
import gsMap.generate_ldscore as _gl
_gl.pr.read_gtf = _cached_read_gtf


COND = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond"
WORK = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap"
RES  = "/users/mtotty/claude_scratch/gsmap/gsMap_resource"
ANNOT= "/users/mtotty/claude_scratch/gsmap/additional_annotation/gsMap_additional_annotation"
GTF  = RES + "/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf"
OUT  = COND + "/annotation_neuronal"
QC   = COND + "/annotation_qc"
WINDOW = 50000
os.makedirs(OUT, exist_ok=True); os.makedirs(QC, exist_ok=True)

log = lambda *a: (print(*a, flush=True))

# ---------- STEP 1: gene universe + affine rescale ----------
feathers = sorted(glob.glob(f"{WORK}/*/latent_to_gene/*_gene_marker_score.feather"))
usets = {f.split("/")[-1].split("_gene")[0]: set(pd.read_feather(f, columns=["HUMAN_GENE_SYM"])["HUMAN_GENE_SYM"])
         for f in feathers}
U_int = sorted(set.intersection(*usets.values()))
U_uni = sorted(set.union(*usets.values()))
log(f"[universe] per-sample sizes: { {k:len(v) for k,v in usets.items()} }")
log(f"[universe] intersection={len(U_int)} union={len(U_uni)}")

ax = pd.read_csv(f"{COND}/gene_neuronal_axis.tsv", sep="\t")
log(f"[edgeR] genes with delta: {len(ax)}  dup={ax.gene.duplicated().sum()}")
delta_map = ax.set_index("gene")["delta"]

def build_axis(universe, delta_series, tag):
    s = pd.Series(0.0, index=pd.Index(universe, name="gene_name"))
    common = s.index.intersection(delta_series.index)
    s.loc[common] = delta_series.loc[common].values
    lo, hi = float(s.min()), float(s.max())
    rescaled = (s - lo) / (hi - lo)
    info = dict(tag=tag, n_universe=len(s), n_with_delta=int(len(common)),
                n_zero_pre=int(len(s) - len(common)), delta_min=lo, delta_max=hi,
                offset_b=float((0.0 - lo) / (hi - lo)))
    log(f"[axis:{tag}] universe={len(s)} with_delta={len(common)} zeroed={len(s)-len(common)} "
        f"delta range [{lo:.4f},{hi:.4f}] offset_b={info['offset_b']:.4f}")
    return rescaled.to_frame("neuronal_axis"), info

axis_main, info_main = build_axis(U_int, delta_map, "main_intersection")
axis_union, info_union = build_axis(U_uni, delta_map, "union")

# abundance-sensitivity: drop bottom logCPM quartile of edgeR genes
q25 = ax["logCPM"].quantile(0.25)
delta_hi = delta_map[ax.set_index("gene")["logCPM"] > q25]
axis_abund, info_abund = build_axis(U_int, delta_hi, "abund_top75")
log(f"[abund] logCPM q25={q25:.4f}  genes retained={len(delta_hi)} of {len(delta_map)}")

# ---------- STEP 2: gsMap gene windows (imported, not reimplemented) ----------
log("[gtf] loading via gsMap.load_gtf ...")
gtf_pr, mk_used = load_gtf(GTF, axis_main, WINDOW)
gtf_pr_u, _     = load_gtf(GTF, axis_union, WINDOW)
log(f"[gtf] genes in windows: main={len(mk_used)} union={len(gtf_pr_u.df)}")
gtf_df = gtf_pr.df
gtf_df[["Start","End","TSS","TED"]].to_csv(f"{QC}/gene_windows_check.csv.gz", index=False)

# ---------- per-chromosome mapping + write ----------
per_chrom, hist_bins = [], np.linspace(0, 1, 101)
hist_tot = np.zeros(100, dtype=np.int64)
spear_x, spear_y = [], []
hdr_ref = None
for chrom in range(1, 23):
    src = f"{ANNOT}/baseline.{chrom}.annot.gz"
    dst = f"{OUT}/baseline.{chrom}.annot.gz"
    hdr = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
    assert hdr[:4] == ["CHR","BP","SNP","CM"], f"unexpected header {hdr[:4]}"
    assert "neuronal_axis" not in hdr, "column already present"
    if hdr_ref is None: hdr_ref = hdr
    else: assert hdr == hdr_ref, f"header mismatch chr{chrom}"

    bim = pd.read_csv(src, sep="\t", usecols=["CHR","BP","SNP"], dtype={"CHR":str})
    # replicate PlinkBEDFile.convert_bim_to_pyrange exactly (BIM 1-based -> PyRanges 0-based)
    bp = pd.DataFrame({"Chromosome": "chr" + bim["CHR"].astype(str),
                       "SNP": bim["SNP"], "Start": bim["BP"].astype(int) - 1,
                       "End": bim["BP"].astype(int)})
    bim_pr = pr.PyRanges(bp)

    def map_scores(gpr, axis_df):
        nb = overlaps_gtf_bim(gpr, bim_pr)              # nearest-TSS, one row per SNP
        m = nb.set_index("SNP")["gene_name"]
        sc = m.map(axis_df["neuronal_axis"])
        return bim["SNP"].map(sc).fillna(0.0).to_numpy(), m

    vals, mapped = map_scores(gtf_pr, axis_main)
    vals_u, mapped_u = map_scores(gtf_pr_u, axis_union)
    vals_a, _ = map_scores(gtf_pr, axis_abund)

    n = len(bim); n_map = int(mapped.index.isin(bim["SNP"]).sum())
    per_chrom.append(dict(chrom=chrom, n_snp=n, n_mapped=n_map,
                          frac_mapped=n_map/n, n_mapped_union=len(mapped_u),
                          mean=float(vals.mean()), median=float(np.median(vals)),
                          nonzero=int((vals>0).sum())))
    hist_tot += np.histogram(vals, bins=hist_bins)[0]
    idx = np.random.default_rng(chrom).choice(n, size=min(60000, n), replace=False)
    spear_x.append(vals[idx]); spear_y.append(vals_a[idx])

    # ---------- STEP 3: byte-identical append ----------
    with gzip.open(src, "rt") as fi, gzip.open(dst, "wt", compresslevel=6) as fo:
        first = fi.readline().rstrip("\n")
        fo.write(first + "\tneuronal_axis\n")
        for i, line in enumerate(fi):
            fo.write(line.rstrip("\n") + "\t" + f"{vals[i]:.6f}" + "\n")
        assert i + 1 == n, f"row count mismatch chr{chrom}: {i+1} vs {n}"
    log(f"[chr{chrom}] snps={n} mapped={n_map} ({n_map/n:.3f}) nonzero={(vals>0).sum()} -> {dst}")

pc = pd.DataFrame(per_chrom); pc.to_csv(f"{QC}/per_chrom_stats.csv", index=False)

# ---------- STEP 4: QC ----------
from scipy.stats import spearmanr
sx, sy = np.concatenate(spear_x), np.concatenate(spear_y)
rho, pval = spearmanr(sx, sy)
allvals = []
for chrom in range(1, 23):
    allvals.append(pd.read_csv(f"{OUT}/baseline.{chrom}.annot.gz", sep="\t",
                               usecols=["neuronal_axis"])["neuronal_axis"].to_numpy())
allvals = np.concatenate(allvals)
qs = [0,1,5,10,25,50,75,90,95,99,100]
summary = dict(
    universe=dict(per_sample={k: len(v) for k, v in usets.items()},
                  intersection=len(U_int), union=len(U_uni)),
    axis_main=info_main, axis_union=info_union, axis_abund=info_abund,
    logCPM_q25=float(q25), n_genes_edgeR=int(len(ax)),
    n_snp_total=int(len(allvals)),
    frac_snp_mapped=float(pc.n_mapped.sum() / pc.n_snp.sum()),
    frac_snp_mapped_union=float(pc.n_mapped_union.sum() / pc.n_snp.sum()),
    quantiles={f"q{q}": float(np.percentile(allvals, q)) for q in qs},
    mean=float(allvals.mean()), sd=float(allvals.std()),
    frac_zero=float((allvals == 0).mean()),
    spearman_vs_abund_filtered=float(rho), spearman_n=int(len(sx)),
)
json.dump(summary, open(f"{QC}/annotation_qc_summary.json", "w"), indent=2)
pd.DataFrame({"bin_left": hist_bins[:-1], "bin_right": hist_bins[1:], "count": hist_tot}
             ).to_csv(f"{QC}/axis_histogram.csv", index=False)
np.savez_compressed(f"{QC}/spearman_sample.npz", main=sx.astype(np.float32), abund=sy.astype(np.float32))
axis_main.join(pd.Series(delta_map, name="delta")).to_csv(f"{QC}/gene_axis_rescaled.csv.gz")
log("[done] " + json.dumps({k: summary[k] for k in
    ["frac_snp_mapped","frac_zero","mean","spearman_vs_abund_filtered"]}))
