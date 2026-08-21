# scDRS — per-spot disease relevance scores for the Visium amygdala data

Third pipeline in `16_LDSC/`, alongside classic stratified LDSC (`../LDSC/`) and
gsMap (`../gsMap/`). Same Visium data, same GWAS, same BayesSpace domains — a
different statistical route from GWAS to spatial map.

**Status (2026-08-09): ANALYSIS COMPLETE. Steps 00–06 and 08 done — 42/42
traits scored and 42/42 group tables written. Step 07 (the in-repo figure
script) has not been run; publication figures were built separately.**
See `RUNLOG.md` for the authoritative execution record.

Step 06 is now run as **`06b_score_array.sh` + `06c_downstream_array.sh`**
(one array task per trait, 20-wide), NOT the original single-job
`06_score_full.sh`, which is retained for reference but is ~14x slower and
effectively single-threaded. Request **32 G** per array task, not the 24 G
currently in the headers — measured peaks reached 24.15 GiB against a 24 GiB
allocation. See RUNLOG "step 06 redesigned as an array".

- `00_setup_env.sh` — **done**. Conda env at
  `/users/mtotty/claude_scratch/scdrs/envs/scdrs` (scdrs 1.0.2, python 3.11.15).
  **This env needs a one-line patch** — see "NumPy 2.x" below.
- `01_prep_magma_ref.sh` — **done** (SLURM `34696538`). 9,997,231 variants
  merged; 18,547 / 18,575 genes annotated at window 10,1.5.
- `02_magma_array.sh` — **done, 42 traits.** `traits.tsv` grew from 40 to 42:
  `PTSD_F3` (row 41) and `BIP_2024` (row 42) are current PGC freezes downloaded
  to `$WORKDIR/gwas_new/`, replacing shared files that had 0 and 50
  genome-wide significant SNPs respectively.
- `03_munge_gs.sh` — **done**, 42 traits in `magma_top1000.gs` (~26 s).
- `04_prep_h5ad.sh` — **done** (3 m 29 s). 224,021 spots × 19,389 genes.
- `05_score_pilot.sh` — **done, gate PASSED.** SCZ / MDD / Height.
  Height gives **0** significant domains vs 5 each for SCZ and MDD; SCZ
  correlates with the existing s-LDSC ranking at **rho +0.842, p = 0.0002**.
  Peak RSS **24.6 GB**.
- `06_score_full.sh` — **submitted** (SLURM `34721738`), 42 traits, ~35–40 h.
  `--mem` reduced 64 G → 32 G on the pilot's measured peak.
- `07_run_figures.sh` — not started; runs after `06`.
- `08_capture_provenance.sh` — **done**. Nine files in `provenance/`.
  Worth re-running after `06` to capture the patched env.

**NumPy 2.x incompatibility (must reapply on any env rebuild).** scdrs 1.0.2
calls `np.float_`, removed in NumPy 2.0, at `scdrs/method.py` lines 1055 and
1062. The pilot completed all scoring and then crashed in `perform-downstream`.
Fix applied to the scratch env — the original is preserved as
`method.py.orig-numpy1`:

    sed -i 's/np\.float_/np.float64/g' \
      /users/mtotty/claude_scratch/scdrs/envs/scdrs/lib/python3.11/site-packages/scdrs/method.py

Alternative: pin `numpy<2` in `00_setup_env.sh`. Note `06` also calls
`perform-downstream`, so the patch must be in place before the sweep runs.

Files created outside this directory: the conda env, the MAGMA reference panel
(both under `/users/mtotty/claude_scratch/scdrs/`), and `magma/SCZ.genes.out`
under `processed-data/Visium/16_LDSC/scDRS/`. Nothing existing anywhere has been
modified or deleted.

---

## How the three methods differ

| | `../LDSC/` | `../gsMap/` | this directory |
|---|---|---|---|
| GWAS input | SNP sumstats + LD scores | SNP sumstats + LD scores | **MAGMA gene-level Z** |
| gene→trait link | binary top-10% specific genes, ±100 kb | continuous GSS mapped to SNPs | top-1000 genes weighted by MAGMA Z |
| unit of analysis | domain pseudobulk | individual spot | **individual spot** |
| null model | LD-score regression | LD-score regression per spot | **1000 matched control gene sets, MC empirical p** |
| domain result | coefficient z / p | Cauchy combination of spot p | group analysis over spot scores |
| runtime | hours | ~310 h of step-4 compute | **hours, single job** |

The key practical difference: scDRS never touches SNPs after MAGMA. It asks
whether a spot's *expression* is unusually high across the trait's top genes,
relative to control gene sets matched on mean expression and expression
variance. That makes it fast and makes the control model explicit, but it also
means all LD/annotation modelling is compressed into MAGMA's gene test upstream.
Where all three methods agree, the result is genuinely triangulated.

---

## Pipeline

```
00_setup_env.sh          conda env in scratch, scdrs 1.0.2                 (once)
01_prep_magma_ref.sh     merge 1000G EUR panel, build SNP/gene loc, annotate (once)
02_magma_array.sh        40-way array: Z→P, then MAGMA gene analysis per trait
  └ 02a_sumstats_to_pval.py     LDSC sumstats → SNP/P/N
03_munge_gs.sh           gene×trait Z matrix → scdrs .gs gene sets
  └ 03a_build_zscore_matrix.py
04_prep_h5ad.sh          7 gsMap h5ad → 1 concatenated h5ad + covariates
  └ 04a_prep_h5ad.py
05_score_pilot.sh        PILOT: SCZ, MDD, Height + domain group analysis
06_score_full.sh         FULL: all 40 traits + domain group analysis
07_run_figures.sh        spatial maps, domain boxplots, heatmap, method comparison
  └ 07_figures.py
08_capture_provenance.sh env versions + input checksums → provenance/  (any time)
config.sh                every path, in one place
traits.tsv               40 traits → sumstats filenames
RUNLOG.md                what was actually run, when, with which job id
```

Run in order; each `sbatch` must finish before the next. `02` is an array
(`sbatch 02_magma_array.sh`, 20 concurrent); everything else is a single job.

---

## Why MAGMA is needed, and what feeds it

scDRS consumes **gene sets with weights**, not sumstats. The published scDRS
gene sets come from MAGMA gene analysis, so we reproduce that step on the same
40 GWAS the other two pipelines use.

Everything MAGMA needs already exists on JHPCE — nothing is downloaded:

| input | path | note |
|---|---|---|
| MAGMA 1.10 | `module load magma/1.10` | LIBD module; compute nodes only |
| LD panel | `…/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/` | 22 per-chr plink files, merged by `01` |
| gene locations | `…/statsgen/PRS/ref_EUR/NCBI37.3.gene.loc` | hg19/NCBI37, 19,427 genes |
| GWAS sumstats | `…/LDSC/base/gwas_brain/` | the same 40 files `../gsMap/03b_stage_gwas_symlinks.sh` links |

**Build consistency.** `NCBI37.3.gene.loc` is hg19, matching `gene_meta_hg19.txt`
used by `../LDSC/03_bed.R` and the `gencode.v46lift37` GTF gsMap uses. All three
pipelines are on the same coordinate system.

**Gene IDs.** `01_prep_magma_ref.sh` rewrites the gene-loc file with the **symbol**
in column 1, so MAGMA emits symbols in `.genes.out` and they join directly to
`var_names` in the h5ad. No Entrez translation step, no silent ID loss.

**The one real conversion.** These sumstats are LDSC-format: `SNP N Z A1 A2`,
with no P column (two files instead carry `SNP A1 A2 N CHISQ Z`). MAGMA needs
P. `02a_sumstats_to_pval.py` computes `P = 2·Φ(−|Z|)` in log space and reports
max |Z|, min P, and how many values hit the 1e-300 floor. For these files max
|Z| ≈ 9.6 (SCZ), so nothing should floor — if the log says otherwise, stop.

**Window.** `10,1.5` kb (MAGMA default, and what the scDRS authors used). Note
this is *much* tighter than the ±100 kb window `../LDSC/03_bed.R` uses. That is a
real methodological difference between the pipelines, not an oversight — MAGMA's
gene test and s-LDSC's annotation are doing different things. Worth a
sensitivity run at a wider window if a trait disagrees between methods.

---

## The spatial input

`04a_prep_h5ad.py` reads **the same seven h5ad files gsMap uses**
(`processed-data/Visium/16_LDSC/gsMap/ST/`) — 224,021 spots × 19,389
protein-coding genes, `BS_k16_Semisupervised_wAI` domain labels. Reusing them
means scDRS cannot silently diverge from gsMap on gene universe or spot
filtering.

Three changes are required for scDRS:

1. **One object, not seven.** scDRS matches control gene sets on mean and
   variance computed across the whole matrix, so all sections must share a
   feature space. Barcodes are prefixed with donor to stay unique.
2. **Raw counts in `.X`.** gsMap needed `layers['count']`; scDRS wants counts in
   `.X` with `--flag-raw-count True` and does its own normalisation + log1p.
3. **A covariate file** — `const`, `n_genes`, `log_umi`, and one-hot donor
   dummies (6 of 7; the intercept absorbs the reference donor). Depth is the
   dominant technical axis in Visium and donor removes section offsets.

**Domain labels are not covariates.** They enter only at the group-analysis
step, so per-spot scores stay annotation-free — the same design choice gsMap
makes, and what lets "the anatomy came back out" count as a result.

---

## Decisions to confirm before running

1. **Output locations are outside this folder.** Scripts write results to
   `processed-data/Visium/16_LDSC/scDRS/` (mirroring where gsMap writes) and
   heavy intermediates to `/users/mtotty/claude_scratch/scdrs/`. Both are new,
   previously non-existent paths; nothing existing is touched. If you want
   everything inside this code directory instead, change `WORKDIR`/`SCRATCH` at
   the top of `config.sh` — nothing else needs editing.
2. **`/dcs04/lieber` is 99% full (18 T free).** Footprint here is modest: the
   merged plink panel ~1.5 GB and per-trait P files ~30 MB each go to *scratch*;
   the project tree gets `.genes.out` (~1 MB × 40), the concatenated h5ad
   (~1.5 GB), and score files (~50 MB × 40 with control scores retained).
3. **Gene window `10,1.5` kb** vs the ±100 kb the classic pipeline uses — see above.
4. **`--n-ctrl 1000`** gives MC p-values resolved to ~1e-3 with an extended
   analytic tail. Raising it costs linearly.
5. **Pilot first.** `05_score_pilot.sh` runs SCZ + MDD (positive controls, both
   significant in the classic results) and Height (negative control). If Height
   shows domain structure, the covariate model or the gene sets are wrong and
   the full sweep should not launch.

---

## Expected runtime

**Measured** where a step has run (see `RUNLOG.md`); the rest are estimates from
the data dimensions and the gsMap runs on the same sections.

| step | estimate | **measured** | basis / note |
|---|---|---|---|
| `00` env build | ~10 min | **~9 min** | |
| `01` merge + annotate | 1–3 h | **~2 min** | far faster than estimated — the plink merge streams; `magma --annotate` alone took 8 s |
| `02` per trait | 5–15 min | **~5 min** (SCZ) | ~1.1 M SNPs × 19 k genes |
| `03` munge-gs | < 5 min | — | text only |
| `04` h5ad concat | 20–40 min | — | 224 k × 19 k sparse, gzip write |
| `05` pilot, 3 traits | 1–3 h | — | 1000 control sets × 224 k spots |
| `06` full, 40 traits | 8–24 h | — | matrix loaded once, per-trait control sampling |
| `07` figures | 10–20 min | — | |
| `08` provenance | ~2 min | **~1 min** | md5 over ~5 GB of inputs |

Queue wait dominates wall-clock, not compute: step `01` sat ~50 min behind the
concurrent gsMap array before running for 2 minutes. Keeping `--mem` tight is
what shortens that wait — dropping `01` from 48 G to 24 G moved it from
`Reason=Priority` to schedulable immediately.

Memory is the number to watch: the concatenated matrix is ~600 M non-zeros but
stays sparse throughout (see *Resource requests*), so `05` requests 48 G and
`06` 64 G. The pilot's `/usr/bin/time -v` line gives the real peak — right-size
`06` from it rather than guessing.

---

## Resource requests

Sized to schedule fast, not to be safe by inflation — on `shared`, an oversized
`--mem` is the main reason a job waits. Cores are cheap here and RAM is the
scarce dimension, so cores are spent where they help and memory is kept tight.

| step | CPU | mem | walltime | why |
|---|---|---|---|---|
| `00` env | 4 | 12 G | 2 h | pip/conda solve; disk-bound |
| `01` magma ref | 4 | 24 G | 12 h | plink merge of ~9 M SNPs is the peak |
| `02` magma | 2 | 8 G | 3 d | **per array task**; MAGMA is single-threaded |
| `03` munge-gs | 1 | 8 G | 1 h | text manipulation on a 19 k × 40 table |
| `04` h5ad | 4 | 48 G | 3 d | concat holds inputs + result transiently |
| `05` pilot | 8 | 48 G | 1 d | sparse throughout — see below |
| `06` full | 8 | 64 G | 3 d | pilot + headroom; **re-size from pilot** |
| `07` figures | 2 | 16 G | 2 h | h5ad read in `backed="r"` mode |
| `08` provenance | 1 | 4 G | 1 h | md5 streaming, I/O-bound |

**The scoring jobs were over-requested by ~4× in the first draft.** Checking the
scdrs 1.0.2 source (`scdrs/pp.py`) settled it: when `.X` is sparse *and* a
`--cov-file` is supplied, scDRS uses **implicit covariate correction** — it
keeps the covariate betas and applies them chunk-wise rather than materialising
a corrected dense matrix. The 224 k × 19 k matrix stays sparse (~600 M
non-zeros, ~7 GB as float32 CSR) instead of becoming a ~17 GB dense array. The
original 192 G was sized for a densification that never happens.

**MAGMA is single-threaded**, so `02` gets 2 cores (one for MAGMA, one for the
Z→P conversion) rather than 4, and runs 20-wide instead of 10 — small jobs
schedule quickly, so 40 traits should clear in about two waves.

**Right-size `06` from the pilot.** Both `05` and `06` wrap `compute-score` in
`/usr/bin/time -v`; read `Maximum resident set size` from the pilot log. Memory
does not scale with trait count (traits are scored in sequence), so the pilot's
peak predicts the sweep's. If the pilot peaks under 32 G, lower `06` to match.

---

## Outputs

```
processed-data/Visium/16_LDSC/scDRS/
├── magma/<trait>.genes.out            MAGMA gene-level Z and P
├── gs/magma_zscore_matrix.tsv         gene × trait Z
├── gs/magma_top1000.gs                scDRS gene sets
├── h5ad/visium_amygdala_scdrs.h5ad    224 k spots, raw counts
├── h5ad/cov.tsv                       covariates
├── score/<trait>.score.gz             per-spot raw/norm score, MC p, FDR
├── score/<trait>.full_score.gz        + control scores
├── downstream/<trait>.scdrs_group.BS_k16_Semisupervised_wAI
└── figures/                           maps, boxplots, heatmap, comparison
```

The deliverable that answers the original question is
`downstream/<trait>.scdrs_group.*`: per domain, the association statistic and
its MC p-value — directly comparable to `../LDSC/ldsc_results.csv` and to gsMap's
`cauchy_across_samples/`. `07_figures.py` plots that comparison.

---

## Reproducibility

Everything needed to re-run this later lives in **this directory**. Software and
bulk intermediates live in scratch, which is not backed up — so the things that
would be expensive to reconstruct are captured here as text:

```
provenance/        <- sbatch 08_capture_provenance.sh
├── env_spec.yml           full conda environment export
├── env_pip_freeze.txt     pip freeze (exact pinned versions)
├── env_conda_list.txt     conda list
├── env_versions.txt       scdrs / scanpy / anndata / numpy / ... versions
├── tool_versions.txt      magma and plink versions, module names
├── input_manifest.txt     md5 of every GWAS file, plink .bim, h5ad, gene.loc
├── code_checksums.txt     md5 of every script in this directory
├── resolved_config.txt    config.sh with all variables expanded
└── history.txt            append-only log of when snapshots were taken
RUNLOG.md          <- what was actually run, when, with which job id
logs/              <- stdout/stderr of every step
```

Run `sbatch 08_capture_provenance.sh` **now** (to pin the inputs) and **again
after the sweep finishes** (to pin the environment that produced the results).
It overwrites the snapshot files and appends a dated line to `history.txt`, so
you always have a current snapshot plus a record of every time one was taken.

The md5 pass over the GWAS files and h5ad inputs takes a few minutes — that cost
buys the ability to prove, later, that the inputs were the same bytes.

**Logging caveat.** Each script declares its own
`#SBATCH --output=.../logs/<step>.%j.log`, which works when you run
`sbatch <script>` directly. Agent-dispatched runs go through a wrapper that
calls `bash <script>`, bypassing that directive — those logs are copied into
`logs/` after the fact. `RUNLOG.md` records which is which.

---

## Guardrails

Every script only creates files under this directory, `WORKDIR`, or `SCRATCH`.
Nothing deletes, moves, or overwrites anything, anywhere — `02`, `04`, and `06`
skip work whose output already exists rather than clobbering it. Shared inputs
(`gwas_brain/`, `referencefiles/`, `PRS/ref_EUR/`, `gsMap/ST/`) are opened
read-only. All nine job scripts run on `shared`. Walltimes are set per step with
headroom over the expected runtime rather than uniformly — see the *Resource
requests* table above. `02`, `04`, and `06` keep `--time=3-00:00:00` per house
convention because their runtime is least predictable; the rest are capped
lower. No script may be run on a login node.
