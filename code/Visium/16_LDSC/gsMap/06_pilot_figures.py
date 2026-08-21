
import os, glob
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import spearmanr
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

W="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap"; LC="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/LDSC"; OUT="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/pilot_figs"
os.makedirs(OUT, exist_ok=True)
BASE, MID, SMALL = 9, 8, 7
mpl.rcParams.update({"figure.dpi":110,"savefig.dpi":300,"savefig.bbox":"tight",
    "font.size":BASE,"axes.titlesize":BASE,"axes.labelsize":BASE,"legend.fontsize":MID,
    "xtick.labelsize":SMALL,"ytick.labelsize":SMALL,"axes.spines.top":False,
    "axes.spines.right":False,"xtick.direction":"out","ytick.direction":"out",
    "legend.frameon":False,"axes.titlelocation":"left","axes.axisbelow":True})

SAMPLE="Br6471"
dom = pd.read_csv(W+"/Br6471_domains.csv", index_col=0); dom.columns=["domain"]
a = ad.read_h5ad(f"{W}/ST/{SAMPLE}.h5ad", backed="r")
xy = pd.DataFrame(np.asarray(a.obsm["spatial"]), index=a.obs_names, columns=["x","y"])

# per-spot results
traits={}
for f in sorted(glob.glob(f"{W}/{SAMPLE}/spatial_ldsc/{SAMPLE}_*.csv.gz")):
    t=os.path.basename(f).replace(f"{SAMPLE}_","").replace(".csv.gz","")
    traits[t]=pd.read_csv(f,index_col=0)

# gsMap's OFFICIAL domain-level Cauchy output (authoritative)
gs_rows=[]
for f in sorted(glob.glob(f"{W}/{SAMPLE}/cauchy_combination/{SAMPLE}_*.Cauchy.csv.gz")):
    t=os.path.basename(f).replace(f"{SAMPLE}_","").replace(".Cauchy.csv.gz","")
    d=pd.read_csv(f)
    d["trait"]=t; d=d.rename(columns={"annotation":"domain"})
    gs_rows.append(d)
gs=pd.concat(gs_rows)
gs["gsmap_neglog10p"]=-np.log10(gs.p_cauchy.clip(lower=1e-300))
nsp=dom["domain"].value_counts().rename("n_spots")
gs=gs.merge(nsp,left_on="domain",right_index=True,how="left")
gs.to_csv(OUT+"/gsmap_domain_cauchy.csv",index=False)
print("traits:", sorted(gs.trait.unique()))

# classic LDSC
cl=pd.read_csv(LC+"/ldsc_results.csv",index_col=0)
cl["classic_neglog10p"]=-np.log10(cl.p_zcore.clip(lower=1e-300))
MAP={"SCZ":"Schizophrenia","MDD":"Depression","Height":"Height"}
NOTE={"MDD":"different GWAS (classic=mdd2019edinburgh)"}
cr=[]
for gt,ct in MAP.items():
    if gt not in gs.trait.values: continue
    g=gs[gs.trait==gt].set_index("domain")
    c=cl[cl.trait==ct].set_index("cell")
    j=g.join(c[["Coefficient_z.score","p_zcore","FDR","Enrichment","classic_neglog10p"]],how="inner")
    j["classic_trait"]=ct; j["caveat"]=NOTE.get(gt,"exact sumstats match")
    cr.append(j.reset_index())
cmp=pd.concat(cr)
if "domain" not in cmp.columns and "index" in cmp.columns:
    cmp=cmp.rename(columns={"index":"domain"})
print("\n=== SPEARMAN (gsMap official Cauchy vs classic coefficient z) ===")
for t in cmp.trait.unique():
    d=cmp[cmp.trait==t]
    r=spearmanr(d["Coefficient_z.score"],d["gsmap_neglog10p"])
    cmp.loc[cmp.trait==t,"spearman_rho_vs_z"]=r.statistic
    cmp.loc[cmp.trait==t,"spearman_p"]=r.pvalue
    print(f"{t:8s} n={len(d):2d}  rho={r.statistic:+.4f}  p={r.pvalue:.3g}")
cmp.to_csv(OUT+"/gsmap_vs_classic_ldsc.csv",index=False)

ORDER=[t for t in ["SCZ","MDD","Height"] if t in traits]

# ============ FIG 1 : spot maps ============
n=len(ORDER)+1
fig,axes=plt.subplots(1,n,figsize=(3.15*n,3.7))
axd=axes[0]
doms=[d for d in sorted(dom.domain.unique()) if (dom.domain==d).sum()>0]
cmap20=plt.get_cmap("tab20"); colmap={d:cmap20(i%20) for i,d in enumerate(doms)}
cent={}
for d in doms:
    idx=dom.index[dom.domain==d]
    axd.scatter(xy.loc[idx,"x"],-xy.loc[idx,"y"],s=0.7,lw=0,color=colmap[d],rasterized=True)
    cent[d]=(xy.loc[idx,"x"].median(),-xy.loc[idx,"y"].median())
span=abs(max(v[1] for v in cent.values())-min(v[1] for v in cent.values()))
minsep=span*0.062; placed=[]
for d in sorted(cent,key=lambda k:cent[k][1]):
    x,y=cent[d]
    while any(abs(y-py_)<minsep and abs(x-px)<span*0.30 for px,py_ in placed): y+=minsep*0.65
    placed.append((x,y))
    axd.annotate(d,(x,y),fontsize=SMALL,ha="center",va="center",zorder=5,
        path_effects=[pe.withStroke(linewidth=2.2,foreground="white")])
axd.set_title("Domains (BayesSpace k16)"); axd.set_aspect("equal"); axd.axis("off")
vmax=float(np.nanmax([np.nanpercentile(-np.log10(traits[t].p.clip(lower=1e-300)),99) for t in ORDER]))
for ax,t in zip(axes[1:],ORDER):
    df=traits[t].join(dom,how="inner").join(xy,how="inner")
    v=-np.log10(df.p.clip(lower=1e-300)); o=np.argsort(v.values)
    sc=ax.scatter(df.x.values[o],-df.y.values[o],c=v.values[o],s=0.7,lw=0,
                  cmap="magma_r",vmin=0,vmax=vmax,rasterized=True)
    ax.set_title(t+("   (non-brain comparator)" if t=="Height" else ""))
    ax.set_aspect("equal"); ax.axis("off")
    cb=fig.colorbar(sc,ax=ax,fraction=0.042,pad=0.02)
    cb.set_label(r"$-\log_{10}P$ per spot",fontsize=SMALL); cb.ax.tick_params(labelsize=SMALL)
fig.savefig(OUT+"/fig1_gsmap_spotplot.png"); plt.close(fig)

# ============ FIG 2 : gsMap vs classic ============
WM={"WM.1","WM.2"}
OFF={"BM":(7,-4),"LA":(-4,10),"BLD":(-17,7),"PL":(7,4),"CoA":(-21,-9),"BL":(7,-9),
     "CeA":(7,-4),"MeA":(8,3),"AI":(-14,-10),"HPC":(8,-2),"Endothelial":(-30,9),
     "CHAT":(8,-2),"WM.1":(8,3),"WM.2":(8,3)}
tls=[t for t in ORDER if t in cmp.trait.values]
fig,axes=plt.subplots(1,len(tls),figsize=(4.5*len(tls),4.2),squeeze=False)
for ax,t in zip(axes[0],tls):
    d=cmp[cmp.trait==t]; iswm=d.domain.isin(WM)
    ax.scatter(d.loc[~iswm,"Coefficient_z.score"],d.loc[~iswm,"gsmap_neglog10p"],
               s=44,color="#3b4cc0",zorder=3)
    ax.scatter(d.loc[iswm,"Coefficient_z.score"],d.loc[iswm,"gsmap_neglog10p"],
               s=44,facecolor="white",edgecolor="#3b4cc0",lw=1.5,zorder=3)
    for _,r in d.iterrows():
        ax.annotate(r.domain,(r["Coefficient_z.score"],r["gsmap_neglog10p"]),
            textcoords="offset points",xytext=OFF.get(r.domain,(6,3)),fontsize=SMALL,
            zorder=6,path_effects=[pe.withStroke(linewidth=2.2,foreground="white")])
    ax.axvline(0,color="0.78",lw=0.8,zorder=1)
    rho=d["spearman_rho_vs_z"].iloc[0]; pv=d["spearman_p"].iloc[0]
    ax.set_xlabel("classic stratified LDSC coefficient z   (negative = depletion)")
    ax.set_ylabel(r"gsMap  $-\log_{10}P$  (Cauchy, domain)")
    ttl=f"{t}   Spearman $\\rho$={rho:+.2f} (P={pv:.2g})"
    if t=="MDD": ttl+="\ndifferent GWAS than classic"
    if t=="Height": ttl+="\nnon-brain comparator"
    ax.set_title(ttl); ax.margins(0.20)
    ax.set_ylim(bottom=min(0,ax.get_ylim()[0]))
h=[plt.Line2D([],[],marker="o",ls="",mfc="#3b4cc0",mec="#3b4cc0",ms=7,label="neuronal / other domain"),
   plt.Line2D([],[],marker="o",ls="",mfc="white",mec="#3b4cc0",mew=1.5,ms=7,label="white matter")]
axes[0][0].legend(handles=h,loc="upper left",borderaxespad=0.6,handletextpad=0.5,labelspacing=0.35)
fig.savefig(OUT+"/fig2_gsmap_vs_classic.png"); plt.close(fig)

# ============ FIG 3 : trait x domain heatmap ============
piv=gs.pivot_table(index="domain",columns="trait",values="gsmap_neglog10p")
piv=piv[[t for t in ORDER if t in piv.columns]]
piv=piv.loc[piv.mean(axis=1).sort_values(ascending=False).index]
fig,ax=plt.subplots(figsize=(1.5+0.95*piv.shape[1],0.34*piv.shape[0]+1.5))
im=ax.imshow(piv.values,cmap="magma_r",aspect="auto")
ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns)
ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v=piv.values[i,j]
        ax.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=SMALL,
                color="white" if v>piv.values.max()*0.55 else "0.15")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03)
cb.set_label(r"$-\log_{10}P$ (Cauchy)",fontsize=SMALL); cb.ax.tick_params(labelsize=SMALL)
ax.set_title("Domain enrichment by trait, "+SAMPLE)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
fig.savefig(OUT+"/fig3_trait_domain_heatmap.png"); plt.close(fig)
print("\nFILES:",sorted(os.listdir(OUT)))
