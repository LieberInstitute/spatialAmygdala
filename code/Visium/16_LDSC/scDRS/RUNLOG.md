# Run log

Append-only record of what was actually executed, when, and what came out.
Every entry should let a future reader reconstruct the run without guessing.
Add a row when you submit; fill in the outcome when it lands.

Log files for each step live in `logs/`. Provenance snapshots (package
versions, input checksums) live in `provenance/` — regenerate with
`sbatch 08_capture_provenance.sh`.

---

## 2026-08-07

| step | how | job id | outcome |
|---|---|---|---|
| scaffold created | agent, create-only copy | — | 15 files + `logs/` written to `code/Visium/16_LDSC/scDRS/`. Two stray wrapper files (`cmd.sh`, `job.sh`) also landed; inert, left in place. |
| `00_setup_env.sh` | agent-dispatched, `shared` 4 CPU / 32 G | `d841b831` | **success**, ~9 min. Env at `/users/mtotty/claude_scratch/scdrs/envs/scdrs`, python 3.11, scdrs 1.0.2. CLI verified (`compute_score`, `munge_gs`, `perform_downstream` all present). Log: `logs/00_setup_env.agentjob_d841b831.log`. |
| doc corrections | agent, copy | `ef04dd71` | `07_figures.py` + `README.md` updated in place. |
| `01_prep_magma_ref.sh` | agent-dispatched, `shared` 4 CPU | `38946cec` (SLURM `34696538`) | submitted at 48 G / 3 d, then **right-sized in place to 24 G / 12 h** via `scontrol update JobId=34696538 MinMemoryNode=24576 TimeLimit=12:00:00` (note: `scontrol` wants MB, not `24G`). **success**, ~2 min once it started (waited ~50 min in queue behind the gsMap array). Merged **9,997,231 variants / 489 EUR individuals** across 22 autosomes, no merge conflicts. Gene-loc keyed on symbol: **18,575 autosomal genes** after dedup. `magma --annotate` (window 10,1.5) mapped **45.01% of SNPs** to ≥1 gene and annotated **18,547 / 18,575 genes** (28 empty — genes in SNP-poor regions, expected). The `plink --version` guard confirmed PLINK v1.90b7 loaded, not PuTTY. Log: `logs/01_prep_magma_ref.agentjob_38946cec.log`. |
| `08_capture_provenance.sh` | agent-dispatched, `shared` 2 CPU / 8 G | `6308847b` (SLURM `34696707`) | **success**, ~1 min. Wrote 9 files to `provenance/`. Captured env (scdrs 1.0.2, scanpy 1.11.5, anndata 0.12.19, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1), tools (MAGMA v1.10 linux/s, PLINK v1.90b7 64-bit 16 Jan 2023), and md5 for all 40 GWAS files, 7 h5ad inputs, 22 plink `.bim`, gene.loc, source SPE, `ldsc_results.csv`. |
| `02_magma_array.sh` | user, manual `sbatch` | SLURM `34696995` | **FAILED — all 40 tasks, exit 1 in 0-1 s.** Cause: `$0` does not resolve to the script's real path under `sbatch`; see *FAILURE + FIX* below. All 9 scripts patched to use `$CODE_DIR`. Fix verified on task 34 (SCZ), which completed and wrote `magma/SCZ.genes.out`. **Awaiting resubmission** (`sbatch 02_magma_array.sh`) for the remaining 39 traits; step 02 is idempotent so SCZ will skip. |

### FAILURE + FIX: `$0` does not resolve under `sbatch` (job 34696995)

The user submitted `sbatch 02_magma_array.sh`. **All 40 tasks failed in 0-1 s**,
exit code 1, each with:

```
/var/spool/slurm/d/job34696995/slurm_script: line 39:
/var/spool/slurm/d/job34696995/config.sh: No such file or directory
```

**Cause.** Every script sourced its config with
`source "$(dirname "$(readlink -f "$0")")/config.sh"`. Under `sbatch`, SLURM
copies the batch script to `/var/spool/slurm/d/jobNNNN/slurm_script` and runs
that copy, so `$0` is the spool path — `dirname` yields the spool directory,
where `config.sh` does not exist.

This was latent: steps `00`, `01`, and `08` all succeeded because they were
agent-dispatched as `bash /full/path/to/script.sh`, where `$0` *is* the real
path. The bug only appears under direct `sbatch`, which is the documented way
to run these scripts.

**Fix.** All 9 scripts now use an explicit, overridable code directory:

```bash
CODE_DIR=${CODE_DIR:-/dcs04/.../code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"
```

**Verified** by running the array script with `SLURM_ARRAY_TASK_ID=34` (SCZ):

| check | value |
|---|---|
| input SNPs | 1,150,814 (11 duplicates dropped, 0 missing Z/N) |
| max \|Z\| | 9.582 — matches the value measured directly from the sumstats |
| min P | 9.53e-22, **0 floored** at 1e-300 |
| median N | 126,281 |
| genes tested | 17,703 of 18,547 annotated |
| top genes | DPYD, C10orf32, FES, CACNB2, CHRM4, DGKZ, NRGN — canonical PGC3 loci |
| known risk genes | DRD2 Z=5.87, CACNA1C Z=6.73, TCF4 Z=6.22, CACNB2 Z=7.50, ZNF804A Z=4.86, GRIN2A Z=4.64 |

`SCZ.genes.out` is therefore already complete; the array will skip task 34 on
resubmission because step `02` is idempotent.

Note: `MAGMA` output is **space-delimited**, not tab — `grep -P "^GENE\t"` finds
nothing. Match on field 1 with awk instead.

### Right-sizing pass (same day, after the above)

All nine scripts had their `--mem` cut, most substantially. Two source-level
checks drove it rather than guesswork:

- **scDRS does not densify.** In `scdrs/pp.py` (v1.0.2), a sparse `.X` plus a
  `--cov-file` triggers *implicit* covariate correction — betas are stored and
  applied chunk-wise, never materialising a corrected dense matrix. The
  224k × 19k matrix stays sparse (~600 M non-zeros, ~7 GB CSR) instead of
  becoming a ~17 GB dense array. `05`/`06` were sized 192 G for a densification
  that never happens; now 48 G / 64 G.
- **MAGMA is single-threaded** (no thread/batch flag in `--help`), so `02`
  dropped from 4 cores to 2 and from 32 G to 8 G, and went 20-wide instead of
  10 — small jobs schedule faster.

| step | before | after |
|---|---|---|
| `00` | 4 CPU / 32 G / 3 h | 4 / 12 G / 2 h |
| `01` | 4 / 48 G / 3 d | 4 / 24 G / 12 h |
| `02` | 4 / 32 G / `%10` | 2 / 8 G / `%20` |
| `03` | 2 / 16 G / 2 h | 1 / 8 G / 1 h |
| `04` | 4 / 128 G / 3 d | 4 / 48 G / 3 d |
| `05` | 8 / 192 G / 3 d | 8 / 48 G / 1 d |
| `06` | 8 / 192 G / 3 d | 8 / 64 G / 3 d |
| `07` | 2 / 64 G / 4 h | 2 / 16 G / 2 h |
| `08` | 2 / 8 G / 4 h | 1 / 4 G / 1 h |

Also added to `01_prep_magma_ref.sh`: a guard asserting `plink --version`
starts with `PLINK`. On login nodes `/usr/bin/plink` is PuTTY's SSH client, and
the genetics plink module refuses to load off-compute — without the guard a
login-node run would silently call the wrong binary.

**Concurrent load.** The gsMap 40-trait array (`34696402`, plus `34696332`) was
running throughout, occupying most of the user's `shared` slots, so scDRS jobs
queued rather than starting immediately. This is expected and harmless — the two
pipelines share no files.

---

## 2026-08-08 — GWAS freeze upgrades, pilot gate, full sweep

### Two traits upgraded to newer published freezes

Audit of the 40 shared sumstats found **4 files with zero genome-wide
significant SNPs** (Epilepsy_Focal, OUD_MVP1, OUD_MVP2, PTSD) and OUD_META with
exactly 1. See `gwas_staleness_sweep.tsv` for the full per-trait table.

PTSD and BIP were replaced with current PGC freezes. **The shared
`gwas_brain/` directory was NOT modified** — it is owned by another user and
the existing s-LDSC/gsMap results depend on it. New files were downloaded to
`$WORKDIR/gwas_new/` and the old traits were kept alongside the new ones so the
power difference is measurable.

| trait | source | figshare | GW-sig SNPs | Bonferroni genes |
|---|---|---|---|---|
| `PTSD` (old) | shared `PTSD.gz` | — | 0 | 0 |
| `PTSD_F3` | ptsd2024, PMID 38637617 | 26349322 | 7,108 | **194** |
| `BIP` (old) | shared `bp_ldscore.gz` | — | 50 | 49 |
| `BIP_PGC3` | shared `bp3_ldscore.gz` | — | 479 | 117 |
| `BIP_2024` | bip2024, PMID 39843750 | 27216117 | — | **215** |

Both new files are **hg19** (verified in headers), matching the MAGMA
reference — no liftover needed. `traits.tsv` now has **42 data rows**;
`PTSD_F3` is row 41, `BIP_2024` row 42.

### Converter extended for two new formats (backward compatible)

`02a_sumstats_to_pval.py` auto-detects format; no flags were added and the
40 original LDSC-format traits process exactly as before.

- **PGC sumstats-VCF** (ptsd2024): detected on a leading `##`. Columns
  `#CHROM ID POS A1 A2 FREQ NEFF Z P DIRE`; maps `ID`→SNP, `P`→P, `NEFF`→N.
  Already has P, so no Z→P conversion.
- **Ricopili daner** (bip2024): has `P` but **no `N` column**. Sample size comes
  from `Neff_half`, which is HALF the summed effective N — it is **doubled**.
  Using it raw would understate N two-fold. Verified: 2 × 81,683.64 = 163,367,
  correctly below the pooled `4/(1/Nca+1/Nco)` estimate of 220,416.
  Fallback order is `N` → `NEFF` → `2×NEFF_HALF` → `4/(1/Nca+1/Nco)`.

`02_magma_array.sh` now accepts an **absolute path** in `traits.tsv` column 2
(`case "$SS" in /*) SRC=$SS ;;`), which is how the self-downloaded sumstats are
reached without touching `$GWAS_SRC`.

### scdrs 1.0.2 is incompatible with NumPy 2.x — patched

The pilot (`34715015`) ran **2 h 43 m, completed all scoring, then crashed** in
`perform-downstream`:

    AttributeError: `np.float_` was removed in the NumPy 2.0 release.

`scdrs/method.py` calls the removed `np.float_` alias at lines 1055 and 1062
inside `gearys_c`. The env has NumPy 2.4.6.

**Fix applied** in the scratch env (NOT in this code directory):

    sed -i 's/np\.float_/np.float64/g' $ENVDIR/lib/python3.11/site-packages/scdrs/method.py

The unpatched original is preserved as `method.py.orig-numpy1` (44,886 bytes
vs 44,888 patched). Nothing was deleted. **A rebuilt environment must reapply
this patch**, or pin `numpy<2`. Verified by re-running the group analysis on the
already-computed pilot scores — all 3 traits produced output, ~45 min.

Note `06_score_full.sh` also calls `perform-downstream`, so submitting the
sweep before patching would have wasted ~35 h and failed at the same line.

### Pilot gate: PASSED

| criterion | result |
|---|---|
| Height spatially uniform | **PASS** — mean-score spread 0.216 vs SCZ 0.817; Kruskal H 565 vs 5,574; **0 domains at mcp<0.05** vs 5 for SCZ and 5 for MDD |
| SCZ agrees with s-LDSC | **PASS** — Spearman rho **+0.842**, p = 0.0002 over 14 domains. Height control: rho −0.244 (ns) |
| Peak RSS | **24.6 GB** measured against 48 G requested |

Height is the informative control precisely because its gene-level signal is
the *strongest* in the set (mean gene Z 3.12, 3× SCZ's) — the test is spatial
uniformity, not weakness.

Top SCZ domains by scdrs `assoc_mcz`: CLA 3.83, HPC 3.76, LA 2.99, BM 2.79,
WM.2 2.44. Most depleted: CHAT −2.04, WM.1 −1.96, MeA −1.78, Endothelial −1.74.
MDD is *not* a weaker SCZ — it peaks at LA 3.36 / BM 3.02 and reaches
significance at CoA 2.08 and AI 2.03 where SCZ is flat or negative.

### Step 06 right-sized and submitted

`--mem` dropped **64 G → 32 G** on the measured 24.6 GB peak; the previous
header is saved as `06_score_full.sh.bak-64G`. Submitted as SLURM
**34721738**, 42 traits. Estimated 35–40 h from the pilot's 3-traits-in-2h43m.

### Bug found and fixed in this session's own QC

The first staleness sweep reported `max|Z| = 0.000` for `Alzheimer_v2`
*alongside 424 genome-wide significant SNPs* — self-contradictory. Cause: the
awk guard `z+0 == z+0` is a tautology, so two non-numeric Z entries in
`Alzheimer_ldscore2.gz` entered the comparison and `a > mx` then compared
lexically. True value is **max|Z| = 39.626, 422 hits**. Rescanned all 40 with a
strict regex numeric guard; only that one file was affected, and it is the only
file in the directory containing non-numeric Z entries. The corrected table
carries an `n_bad_Z` column so this failure mode is visible rather than silent.

---

## 2026-08-09 — step 06 redesigned as an array; pipeline completed 42/42

### The single-job design was ~14x slower than it needed to be

`06_score_full.sh` loops all traits in one job, justified in its header by
amortizing the h5ad load. Measured on the running job (`34721738`):

| quantity | measured |
|---|---|
| setup (load h5ad, normalize, gene-match) | ~7 min, once |
| per-trait control sampling (1,000 sets) | ~32-39 min |
| AveCPU vs elapsed | **02:08:48 / 02:10:15 = 1.0 core of 8** |

So setup is only ~15-20% of ONE trait's cost, and `scdrs compute-score` is
effectively **single-threaded** (no BLAS thread vars are set in `config.sh`).
The amortization argument does not hold. Projected 42-trait cost: ~48 h serial
(the cancelled job had saved 2 traits in 2 h 19 m) versus ~2.5 h for a 20-wide
array. Replaced by `06b_score_array.sh` (scoring) + `06c_downstream_array.sh`
(group analysis), both `--array=1-42%20`, 2 CPU each.

**The array is per TRAIT, not per sample** — all 7 donors are already merged
into one h5ad (224,021 spots), unlike gsMap's sample x trait fan-out.

### Memory: 24 G was too tight

Extrapolated from two points (pilot 3 traits = 14.1 GB, sweep 42 traits =
22.0 GB) to ~13.5 GB baseline + ~200 MB/trait, predicting a one-trait peak of
13.7 GB, and requested 24 G. **The real per-task range was 13-24 GiB**, with
task 12 at 25,311,944K = 24.15 GiB against a 24 GiB allocation (101%, no OOM
kill). Peak memory also depends on how much of the trait's gene set matches the
19,389 measured genes, which trait count alone does not capture.
**A rerun should request 32 G.**

Note SLURM reports MaxRSS in KiB — dividing by 1e6 (not 1024^2) overstates it.

### Idempotency must test gzip integrity, not file size

`06b`'s original guard was `[ -s "$SCORE_DIR/$TRAIT.score.gz" ]`. When the
serial sweep was cancelled mid-write, `AUD_EA_MVP.full_score.gz` was left
854 MB but truncated; the size check accepted it, the array skipped rescoring,
and `06c` then died with:

    EOFError: Compressed file ended before the end-of-stream marker was reached

Fixed — both outputs must now pass `gzip -t` before a trait is skipped. An
integrity scan of all 84 score files found exactly one corrupt file (that one);
it was rescored and now passes.

### 06c walltime raised 6 h -> 18 h

Three traits (EduYears, Smoking, Stroke_2022_Any) hit the original 6 h wall and
were killed having written nothing. All three then completed in **26-30 min**
on resubmission, so the wall was not a per-trait cost problem (see below).
Completed-task elapsed across the array spans 2 s to 6 h, median 44 min.

### UNRESOLVED: tasks that run for hours without creating their log file

Six occurrences across both arrays. Signature, every time:

- `sacct` reports the job and its `.batch`/`.extern` steps as RUNNING;
- the `#SBATCH --output` file is **never created**, while sibling tasks write
  their first line within seconds and ~4.7-93 KB total;
- no `AveCPU`/`MaxRSS` is ever recorded;
- **cancel + resubmit completes normally in 20-90 min.**

Affected: GSCAN_AgeSmk (9 h 42 m -> 1 h 10 m), EduYears (4 h 14 m -> 26 min),
AUD_EA_MVP (3 h 28 m -> 23 min), plus the three "timeouts" above. Nodes differed
each time (156, 099, 109, others), so excluding individual nodes is not a fix.
`sstat` cannot be used to diagnose it — it fails to resolve array task IDs on
this cluster (`sstat -j 34723020_18` returns "No steps running for job
34723020"; the `.batch` suffix resolves to the wrong StepId). Compute nodes are
not reachable by ssh from the login node. **Cause not identified** — worth
raising with the cluster admins; it is not specific to this pipeline.

### Merged spatial coordinates — verified, not assumed

The h5ad `obsm['spatial']` was checked against every candidate column in the
source SpatialExperiment (`spe_harmony_markers_BS_k16_Semisupervised_wAI.rds`):

| candidate | max deviation (px) |
|---|---|
| **`pxl_col_in_fullres` / `pxl_row_in_fullres`** | **0.0 / 0.0 on all 224,021 spots** |
| `*_rounded` | 160 / 139 |
| `*_original` | 58,006 / 50,902 |
| `array_col` / `array_row` | 62,288 / 74,413 |

Barcodes differ in form between the two objects (the h5ad prefixes the donor);
stripping `^Br\d+_` joins 224,021/224,021. Using `_original` would have stacked
every section of a donor on top of itself.

Capture areas within a donor **do overlap** in merged coordinates: at a 277 px
spot pitch, 57,707 spots (25.8%) have a cross-area neighbour within 139 px,
ranging from 1.4% (Br9280) to 41.5% (Br6423). This is expected and plots are
built on the merged coordinates regardless.

### CLA enrichment is confounded with a single donor

CLA is present in **1 of 7 donors** (Br6660, 4,834 spots); LA in 6 of 7 (absent
in Br6660); the other 13 domains in 7/7. Three checks:

1. Br6660 is not globally elevated — per-donor global means all within +/-0.005,
   and it ranks 7th of 7 for SCZ.
2. Within Br6660 alone, CLA still ranks 1st of its 14 domains for MDD,
   BIP_PGC3, Insomnia, ADHD and BMI (2nd for SCZ, behind BLD).
3. Br6660 agrees with the other donors — Spearman rho +0.918 (p<1e-4) on the 13
   domains present in all 7; Br8325 is highest at +0.934.

So the effect is real within that donor but **n=1 donor with 4,834
pseudo-replicate spots**. It cannot be separated from a Br6660-specific effect
and should be reported with that caveat, not as a population-level finding.
Partial density confound: domain mean SCZ score correlates with domain median
`n_genes` at rho 0.495 (though WM.2 argues against a pure density story).

### Final state

**42/42 traits scored, 42/42 group tables.** 38 domain x trait pairs reach
FDR<0.05 of 630 tests; 21 of 42 traits have at least one significant domain.
Signal concentrates in LA, BM, CLA and WM.2.

Negative controls behave: **Height ranks last of 38 informative traits on
domain-mean spread (0.216)**, has zero significant domains, and shows the
largest negative residual against GWAS power. Four GWAS with zero genome-wide
significant SNPs (Epilepsy_Focal, OUD_MVP1, OUD_MVP2, and the old PTSD file)
are excluded from interpretation rather than deleted.

**PTSD's 2024 freeze is well powered (7,108 gw-sig SNPs, N=641,533) and still
has no significant amygdala domain** — its flatness is a property of the trait
in this tissue, not a power limitation.

Caveat for interpretation: **22 of the significant pairs sit at the Monte Carlo
p-value floor** of 1/(n_ctrl+1) = 0.000999. Within that group, ranking is by
`assoc_mcz`, not by p-value. Raising `--n-ctrl` would resolve them further.

---

## Notes for reproducing

- **Recommended invocation is plain `sbatch <script>` from this directory.** Each
  script carries its own `#SBATCH --output=.../logs/<step>.%j.log`, so logs land
  in `logs/` automatically.
- **Agent-dispatched runs are different.** The compute channel wraps the script
  (`bash <script>`), which means the inner `#SBATCH --output` directive is *not*
  honoured — the log goes to the job workdir instead and must be copied into
  `logs/` afterwards. If a step's log is missing from `logs/`, that is why.
- Steps are idempotent: `02`, `04`, and `06` skip work whose output already
  exists, so a partial run can be resubmitted as-is.
