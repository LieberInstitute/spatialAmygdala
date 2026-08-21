# gsMap — spatial GWAS enrichment for the Visium amygdala data

Spot-level GWAS heritability enrichment with [gsMap](https://github.com/JianYang-Lab/gsMap)
(Song et al., *Nature* 2025), run alongside the classic stratified-LDSC pipeline in `../LDSC/`.

Both pipelines use the **same BayesSpace domains** (`BS_k16_Semisupervised_wAI`) and the
**same GWAS sumstats** (`/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/`), so results are
directly comparable. The difference is the unit of analysis:

| | `../LDSC/` (classic) | this directory (gsMap) |
|---|---|---|
| unit | domain pseudobulk | individual spot |
| annotation | binary, top-10% specific genes ±100 kb | continuous gene specificity score (GSS) |
| domain result | LDSC coefficient z / p | Cauchy combination of per-spot p-values |
| scale | 15 domains × 40 traits | 36k spots/section × 40 traits, then aggregated |

gsMap learns a latent representation of each spot with a graph VAE, derives a GSS from each
spot's homogeneous neighbours, maps GSS onto SNPs, and runs stratified LDSC per spot. Domain
labels enter **only** at the final Cauchy step — the per-spot statistics are annotation-free,
so recovery of anatomical structure is an emergent result rather than an input.

---

## Layout

    16_LDSC/
    ├── LDSC/                      classic stratified-LDSC pipeline (pre-existing)
    │   └── ldsc_results.csv       15 domains × 40 traits, the comparison target
    └── gsMap/                     <- this directory
        ├── 00_setup_env.sh              build the conda env (run once)
        ├── 00_download_resources.sh     fetch + verify the gsMap resource bundle (run once)
        ├── 01a_inspect_spe.R            read-only SPE audit -> 01a_spe_summary.txt
        ├── 01_spe_to_h5ad.R             SPE -> one h5ad per capture area
        ├── 01_check_h5ad.py             assert the gsMap input contract
        ├── 02_gsmap_quick_mode.sh       PILOT: one section × 3 traits
        ├── 03_create_slice_mean.sh      cross-section gene-rank normalisation
        ├── 03b_stage_gwas_symlinks.sh   symlink the 40 sumstats into the workdir
        ├── 04a_gsmap_cache_pass.sh      SWEEP stage 1: cache steps 1-3 + Height (7-way array)
        ├── 04b_gsmap_trait_array.sh     SWEEP stage 2: step 4 per trait (40-way array)
        ├── 04_gsmap_array.sh            SUPERSEDED single-stage sweep -- see note below
        ├── 05_cauchy_across_samples.sh  domain p-values pooled over sections
        ├── 06_pilot_figures.py          spot map + gsMap-vs-classic comparison
        ├── gwas_config_pilot.yaml       3 traits (SCZ, MDD, Height)
        ├── gwas_config_full.yaml        all 40 traits
        └── samples.txt                  the 7 capture areas, one per line

**Software and reference data live in scratch, not in the project tree** (`/dcs04/lieber` is
~99% full):

    /users/mtotty/claude_scratch/gsmap/envs/gsmap      conda env, gsMap 1.73.7, python 3.11 (CPU-only)
    /users/mtotty/claude_scratch/gsmap/gsMap_resource  reference bundle (842 MB)

Outputs go to `../../../processed-data/Visium/16_LDSC/gsMap/`:

    ST/<sample>.h5ad                     gsMap inputs (7 sections)
    GWAS/                                symlinks to the sumstats (never copies)
    <sample>/find_latent_representations GVAE latent embedding
    <sample>/latent_to_gene              gene specificity scores
    <sample>/spatial_ldsc/<sample>_<trait>.csv.gz   per-spot beta/se/z/p
    <sample>/cauchy_combination          domain-level p-values
    <sample>/report                      HTML report with diagnostics

Activate the environment with:

    module load conda/3-24.3.0
    source activate /users/mtotty/claude_scratch/gsmap/envs/gsmap

---

## Input conventions (why the h5ad looks the way it does)

- **Gene universe**: protein-coding only, deduplicated on `gene_name`, symbols as `var_names`
  (gsMap matches genes by symbol). 19,389 genes; gsMap retains ~17,717 after intersecting with
  the GTF. This mirrors `../LDSC/01_aggregate.R`'s convention, minus `filterByExpr` — gsMap does
  its own expression handling.
- **Counts in `layers['count']`** (singular). gsMap's `--data_layer` default is `counts`
  (plural), so **every command must pass `--data_layer count` explicitly.**
- **Coordinates in `obsm['spatial']`** — stitched full-resolution pixel coordinates.
- **Species is human**, so `--homolog_file` is never passed.
- **Domain filter**: domains with <70 spots genome-wide are dropped, matching the classic
  pipeline. 15 domains survive; `CLA` exists only in Br6660, so most sections carry 14.
- **Genome build**: the bundle GTF is `gencode.v46lift37` = **hg19**, consistent with the
  classic pipeline's `gene_meta_hg19.txt`.

Capture areas and spot counts (`samples.txt`):

    Br6471  36876      Br6660  36756      Br6423  34211      Br2743  32556
    Br8325  29481      Br9192  28306      Br9280  25835

---

## How to run

### Once, to set up

    bash 00_setup_env.sh          # inside a compute session, not the login node
    bash 00_download_resources.sh

### Once, to build inputs

    sbatch 01a_inspect_spe.R      # via a wrapper; login nodes force-quit Rscript
    sbatch 01_spe_to_h5ad.R

### Pilot (validated; see Results below)

    sbatch 02_gsmap_quick_mode.sh

### Full sweep — two stages

Each stage must finish before the next starts.

    bash   03b_stage_gwas_symlinks.sh              # light file op, no scheduler
    sbatch 03_create_slice_mean.sh                 # cross-section gene-rank normalisation
    sbatch 04a_gsmap_cache_pass.sh                 # stage 1: 7-way array, ~100 min each
    sbatch --array=1-21,23-40 04b_gsmap_trait_array.sh   # stage 2: 40-way array, ~7.7 h each

**Why two stages.** Steps 1–2 (GVAE latent + gene specificity scores) cost ~32 min per section
and are identical for every trait, so `04a` pays them once per section and caches the result.
`04b` then fans out over traits, running only step 4 (`run_spatial_ldsc`) against those caches.
One trait × 7 sections is ~7.7 h, so the whole sweep finishes in about that time if your SLURM
concurrency limit lets all 40 tasks run at once — versus ~45 h per job for the single-stage
layout.

**What `04a` actually costs.** It invokes `quick_mode`, which has no "stop after step 3" flag, so
it necessarily runs step 4 for whichever trait drives it. That trait is Height — the sweep's
negative control — so nothing is wasted, but the per-task cost is
**~32 min (steps 1–2) + ~67 min (Height step 4) ≈ 100 min**, not the 32 min of caching alone.
Because Height is therefore already complete, stage 2 skips task 22
(`--array=1-21,23-40`); plain `--array=1-40` just recomputes an identical Height.

Both stages skip sample-trait pairs whose output already exists, so re-running after a partial
failure is safe. To re-run one trait: `sbatch --array=34 04b_gsmap_trait_array.sh` (34 = SCZ,
per `traits.txt`).

`05_cauchy_across_samples.sh` is only needed if you want to re-pool domain p-values separately —
`04b` already runs the across-sample Cauchy combination per trait into
`cauchy_across_samples/<trait>_cauchy.csv.gz`.

**Concurrency safety** (verified against the gsMap 1.73.7 source, not assumed): step 4 writes
only `{workdir}/{sample}/spatial_ldsc/{sample}_{trait}.csv.gz` — unique per trait — reads the
shared per-sample intermediates read-only, and logs to stdout via `StreamHandler` rather than a
shared file. 40 concurrent tasks over the same 7 sections therefore do not collide.

**`04_gsmap_array.sh` is superseded.** It runs the sweep as 7 jobs of 40 traits (~45 h each,
~4.5 days wall-clock with `%3` throttling). It is kept for reference and still works, but
`04a` + `04b` is the recommended path.

---

## Runtime (measured on Br6471, 36,876 spots, 8 CPU, no GPU)

| step | cost | scope |
|---|---|---|
| 1 find_latent_representations (GVAE) | 13 min | per sample |
| 2 latent_to_gene (GSS) | 19 min | per sample |
| 3 generate_ldscore | 0 min — quick_mode uses the pre-built SNP×gene matrix | — |
| 4 spatial_ldsc | **66.6 min** | **per sample × trait** |

Wall-clock for a job is roughly `(samples × traits) × 66 min`; steps 1–2 are a one-time
32 min per sample and are cached for reuse across traits. Step 4 dominates everything.

**Planning the sweep.** 40 traits × 7 sections = 280 sample-trait pairs ≈ 310 h of step-4
compute. How that maps to wall-clock depends entirely on how it is split:

| layout | jobs | per-job wall-clock |
|---|---|---|
| 7 samples × 40 traits (as `04_gsmap_array.sh` is written) | 7 | ~45 h |
| 1 trait × 7 samples | 40 | ~7.7 h |
| 1 trait × 3–4 samples | 80 | ~3.3–4.4 h |
| 1 trait × 1 sample | 280 | ~1.1 h |

Traits batched into one job share only the LD-score disk load, not the regression, so batching
saves far less than it appears to. Match the job count to your SLURM concurrency limit —
submitting more jobs than you have slots just queues them.

Inside step 4, spots are processed in serial chunks of 1,000, parallelised *within* a chunk by
`thread_map(max_workers=--max_processes)`. Because those are Python threads, scaling past ~8
cores is unverified; the serial outer loop and its per-chunk disk reads do not benefit from
more cores at all. Test with `--chunk_range` before assuming a larger allocation helps.

---

## Results — pilot (Br6471, SCZ)

Domain-level Cauchy p-values (gsMap's own `run_cauchy_combination` output) against
`../LDSC/ldsc_results.csv`:

| trait | Spearman ρ | P | n domains | note |
|---|---|---|---|---|
| SCZ | +0.76 | 0.0024 | 13 | exact sumstats match |
| MDD | +0.89 | 0.00011 | 12 | different GWAS (classic = mdd2019edinburgh) |
| Height | +0.58 | 0.039 | 13 | non-brain comparator |

Over the 13 domains shared with the classic table, the domains classic LDSC ranked highest are
also gsMap's strongest, and both methods place white matter last:

| domain | gsMap Cauchy P | gsMap −log₁₀P | classic z | classic FDR |
|---|---|---|---|---|
| BM | 3.00e-24 | 23.5 | +3.20 | 0.029 |
| LA | 2.50e-22 | 21.6 | +2.86 | 0.054 |
| BLD | 1.82e-19 | 18.7 | +2.35 | 0.126 |
| PL | 3.96e-18 | 17.4 | +1.33 | 0.477 |
| AI | 2.16e-16 | 15.7 | +0.89 | 0.651 |
| MeA | 5.55e-16 | 15.3 | −1.36 | 0.477 |
| BL | 1.67e-15 | 14.8 | +0.96 | 0.614 |
| CeA | 1.88e-13 | 12.7 | +0.87 | 0.663 |
| HPC | 3.53e-11 | 10.5 | +2.82 | 0.055 |
| Endothelial | 2.02e-05 | 4.7 | −3.64 | 0.013 |
| CHAT | 9.70e-01 | 0.0 | −0.69 | 0.725 |
| WM.2 | 1.00e+00 | 0.0 | −1.14 | 0.562 |
| WM.1 | 1.00e+00 | 0.0 | −2.62 | 0.077 |

(CoA is gsMap's 3rd strongest domain at 4.28e-21 but has no row in the classic table, so it is
absent from the merge.)


Per-spot signal follows the anatomy — continuous across the basolateral complex and cortical-like
nuclei, dropping sharply at white-matter boundaries — despite domain labels playing no part in
the per-spot regression.

### Interpretation caveats

1. **Per-spot p-values are not independent.** Neighbouring spots share GSS through the spatial
   KNN graph; 77% of Br6471 spots have p < 0.05 for SCZ. Rankings are meaningful; per-spot
   "discoveries" are not.
2. **The top domains span a wide range** (BM 3.0e-24 to HPC 3.5e-11), so they are not tied —
   but the weakest domains floor at exactly p = 1, which is a boundary, not a measurement.
3. **Depletion is signed in classic LDSC but not in gsMap.** Endothelial's classic z = −3.64 is
   *significant depletion*; comparing it against an unsigned −log₁₀P misreads the disagreement.
4. **Trait-specificity: resolved for the pilot.** SCZ enriches nearly every neuronal domain,
   which could have been genuine polygenicity or a bias toward high-expression grey matter.
   Height discriminates them: its basolateral values are p ≈ 0.004–0.016, nine to twenty orders
   of magnitude weaker than SCZ, and its two strongest domains are Endothelial (1.4e-15) and
   WM.1 (3.5e-08) — exactly the domains SCZ and MDD rank last. The grey/white contrast is
   therefore not a generic expression artifact.
5. **Trait-name mismatch with the classic table.** Classic "Depression" is `mdd2019edinburgh`;
   the pilot's MDD is `MDD_ldscore_PGC_UKB_23andme`. Only **SCZ** and **Height** are exact
   sumstats matches. `gwas_config_full.yaml` keeps both MDD variants under distinct names.
6. **All values above come from gsMap's own `run_cauchy_combination`.** An earlier
   hand-rolled Cauchy implementation gave the same domain *ranking* but different
   magnitudes (e.g. SCZ/BM 1.0e-16 vs the official 3.0e-24, ρ = 0.79 vs the official
   0.76); those numbers are superseded and should not be quoted.
7. **Height is a non-brain comparator, not a null control.** It correlates with classic
   LDSC (ρ = +0.58, P = 0.039), so it is not orthogonal to the classic results. What
   distinguishes it is *which* domains it favours: Endothelial (1.4e-15) and WM.1
   (3.5e-08) are its strongest, exactly the domains SCZ and MDD rank last.

---

## Guardrails observed

- Everything written here stays in `code/Visium/16_LDSC/gsMap/`,
  `processed-data/Visium/16_LDSC/gsMap/`, or `/users/mtotty/claude_scratch/gsmap/`.
- GWAS sumstats are **symlinked, never copied**, and never modified.
- No script deletes anything.
- Array jobs are submitted manually, not dispatched programmatically.

---

# Three-arm conditional design (added Aug 8 2026)

The sweep documented above is **arm 0**. It is complete but conditions on only gsMap's
default baseline, which holds **two annotations** (`all_gene`, `base`) - verified by reading
`quick_mode/baseline/baseline.1.l2.ldscore.feather` (3 columns incl. SNP) and confirmed to be
gsMap's design on BOTH the quick_mode and step-by-step code paths (`generate_ldscore.py:497`,
`columns=["all_gene","base"]`). That is one gene-level annotation, not the ~97-annotation
baselineLD model classic stratified LDSC uses.

| arm | conditioning set | annotation cols | LD scores | status |
|---|---|---|---|---|
| **0** | gsMap default (`all_gene`, `base`) | 3 | `quick_mode/baseline` (precomputed) | **complete**, 40 traits x 7 samples |
| **A** | 52 functional annotations (baselineLD family) | 53 | `gsMap_cond_functional/` | LD scores complete; pilot running |
| **B** | those 52 + `neuronal_axis` | 54 | `gsMap_cond_neuronal/` | LD scores complete; pilot running |

## Why arm B exists

Arm 0 has one dominant axis: **PC1 = 62.35%** of trait x domain variance (66.84% on the 24
traits with dynamic range >= 5), contrasting neuronal domains (HPC +0.254, LA +0.237, CoA
+0.214) against vasculature and white matter (Endothelial -0.669, WM.1 -0.454).

Removing that axis post-hoc (residualizing on PC1) is **circular** - it conditions on a
component estimated from the same data. Arm B instead puts the axis **in the model** as a
continuous SNP annotation, so the regression accounts for it rather than the analyst
subtracting it afterwards.

**Read PC1 as a gradient, not enriched-vs-null.** Across the 37 brain traits only WM.1
(median p 0.997) and WM.2 (0.995) are near-null. CHAT is significant in 31 of 37 traits
(median p 3e-04) and Endothelial in 23 of 37 (median 0.011). A negative PC1 loading means
"less enriched than the grey-matter average", NOT "unenriched".

## The neuronal_axis annotation

Paired pseudobulk contrast, `~ donor + compartment`, 7 donors x 2 compartments, TMM
normalisation, shrunken logFC via `predFC(prior.count=5)`.

- **neuronal** (12): AI BL BLD BM CeA CHAT CLA CoA HPC LA MeA PL
- **non-neuronal** (3): WM.1 WM.2 **Endothelial** - Endothelial included by explicit design
  decision (non-neuronal, and its PC1 loading matches WM.1's)

A single **signed** axis, not two one-sided annotations: a two-column {delta, |delta|} form
spans the same space as the pair and over-conditions. Signed values are safe for the LD-score
path but note `generate_ldscore.py:453` sums annotation values for the M file, where negatives
would cancel - hence the affine rescale to [0,1] before writing.

QC (all hard-stop checks passed): 14,992 of 19,389 genes kept by `filterByExpr`; all 10
myelin markers negative (CLDN11 -2.72 ... OLIG1 -1.28); all 8 neuronal markers positive
(SLC17A7 +1.65 ... GAD2 +0.89); all 5 endothelial markers negative; donor consistency mean
off-diagonal r = **0.789**. Abundance sensitivity: rebuilding without the bottom logCPM
quartile gives SNP-level Spearman **0.964** against the shipped version.

Gene -> SNP mapping **imports gsMap's own** `load_gtf` / `overlaps_gtf_bim` rather than
reimplementing them. This matters: `load_gtf` swaps TSS/TED for negative-strand genes
(lines 84-87) *after* building the window, so a hand-rolled "nearest TSS" would mis-assign
every minus-strand SNP. Window is **gene body +/-50 kb** (lines 74-79), not TSS +/-50 kb.
60.69% of 9,997,231 SNPs map to a gene; the rest are annotated 0.

## The silent-unconditioned trap

`config.py:1185-1207` - `process_additional_baseline_annotation` sets
`use_additional_baseline_annotation = False` **silently** when
`{ldscore_save_dir}/additional_baseline/` is absent. No error, no warning. The run then
produces **unconditioned** results and **exits 0**.

`09_cond_spatial_ldsc.sh` therefore guards explicitly:
`exit 4` annotation dir missing, `exit 5` log did not confirm the baseline was used,
`exit 6` expected output absent. Never trust the exit code alone.

## Measured step 3 cost (`run_generate_ldscore`)

quick_mode **skips step 3 entirely** - it symlinks the precomputed baseline. The conditional
arms must build it. Two chromosomes measured on Br6471, arm A:

| | chr22 | chr2 |
|---|---|---|
| elapsed | 1m53s | 11m00s |
| sacct MaxRSS | 15.75 GiB | 36.37 GiB |

Runtime and memory scale with **chromosome SNP count** - chromosomes are not interchangeable
budget units. gsMap's own "Peak memory usage" log line **under-reports sacct MaxRSS by
1.55-1.78x**; size `--mem` from sacct only. Tiered submission used in practice:
`--array=1-6 --mem=60G`, `--array=7-13 --mem=45G`, `--array=14-22 --mem=35G`.

Full build measured: **308/308 chromosome-builds** (22 chrom x 7 samples x 2 arms), ~1 wall
hour per sample-arm, ~0.84 TB total.

## Rerunning the conditional arms

```bash
# step 3 - LD scores (ALREADY DONE for both arms; only if rebuilding)
for S in Br6471 Br6660 Br6423 Br2743 Br8325 Br9192 Br9280; do
  mkdir -p ${WORK}_cond_${ARM}/$S/latent_to_gene
  ln -sfn ${WORK}/$S/latent_to_gene/${S}_gene_marker_score.feather \\
          ${WORK}_cond_${ARM}/$S/latent_to_gene/${S}_gene_marker_score.feather
  sbatch --array=1-6   --mem=60G --time=02:00:00 --cpus-per-task=4 \\
         --export=ALL,ARM=functional,SAMPLE=$S --exclude=compute-058,compute-175 \\
         08_cond_generate_ldscore.sh
  # repeat --array=7-13 --mem=45G and --array=14-22 --mem=35G; ARM=neuronal for arm B
done

# step 4 - spatial LDSC per trait, per arm
sbatch --mem=48G --exclude=compute-058,compute-175 \\
       --export=ALL,ARM=functional,SAMPLE=Br6471,TRAIT=SCZ 09_cond_spatial_ldsc.sh
```

The `latent_to_gene` symlink is a **required prerequisite** -
`run_generate_ldscore` reads `{sample}_gene_marker_score.feather` and fails without it.
`08_cond_generate_ldscore.sh` does not create it.

## Memory, corrected

`--mem=24G` caused OUT_OF_MEMORY kills in arm 0 step 4 (EduYears, T2D, Neuroticism, Autism).
EduYears exited OUT_OF_MEMORY **with complete results** - the clearest illustration of why
success is judged by output files, never exit codes. Use **48G** for step 4.

## Arm 0 final numbers

40 traits x 15 domains, 280 of 280 section files, all 40 cross-sample Cauchy files present.
Slowest trait SCZ_PGC2_CLOZUK at 6h49m (peak 15.4 GB); typical 1-3.5 h.

Residual structure after removing PC1 is **+0.004 above the simulated null** - essentially
nothing. The null for k trait-vectors orthogonal to a common direction is **-1/(k-1)**, not
zero (simulated: -0.099 for k=10). Any residual-correlation claim must use that null.

## Cross-workdir prerequisites (the failure class that bit us three times)

gsMap assumes every step shares ONE workdir. Running conditional variants in separate workdirs
(`gsMap_cond_functional`, `gsMap_cond_neuronal`) breaks that assumption, and it never warns —
the expensive step succeeds and a downstream cheap step dies, or worse, exits 0 having done nothing.

Link these from the arm-0 workdir into EVERY sample directory of EVERY conditional workdir
BEFORE submitting. Symlinks only (`ln -sfn`); never move or delete the arm-0 originals.

| needed by | file | symptom if missing |
|---|---|---|
| step 3 `run_generate_ldscore` | `<S>/latent_to_gene/<S>_gene_marker_score.feather` | errors immediately |
| step 4 `run_cauchy_combination` | `<S>/find_latent_representations/<S>_add_latent.h5ad` | **regression runs 2.5 h, writes per-spot output, THEN FileNotFoundError** |

```bash
WORK=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap
for ARM in functional neuronal; do
  for S in Br6471 Br6660 Br6423 Br2743 Br8325 Br9192 Br9280; do
    mkdir -p ${WORK}_cond_${ARM}/$S/latent_to_gene \
             ${WORK}_cond_${ARM}/$S/find_latent_representations
    ln -sfn ${WORK}/$S/latent_to_gene/${S}_gene_marker_score.feather \
            ${WORK}_cond_${ARM}/$S/latent_to_gene/${S}_gene_marker_score.feather
    ln -sfn ${WORK}/$S/find_latent_representations/${S}_add_latent.h5ad \
            ${WORK}_cond_${ARM}/$S/find_latent_representations/${S}_add_latent.h5ad
  done
done
```

Re-running the Cauchy step alone after linking takes ~3 s and 1.1 GB, so units whose regression
already completed are cheap to repair — link, re-run Cauchy only, never recompute the regression.

A third instance of the same class: gsMap's step 6 report generation raised `KeyError: 'Z'`
(`diagnosis.py:23`) in the original arm-0 pilot with all science complete. Cosmetic, but it is why
**exit codes are never the success criterion here** — judge by output files.

## Measured step 4 cost, conditioned (pilot, Br6471, 8 CPU / 48G)

| arm | trait | node | wall | MaxRSS |
|---|---|---|---|---|
| A | SCZ | compute-053 | 5.76 h | 26.6 GB |
| A | MDD | compute-160 | 2.56 h | 26.2 GB |
| A | Height | compute-160 | 2.41 h | 23.4 GB |
| B | SCZ | compute-153 | 2.66 h | 27.0 GB |
| B | MDD | compute-162 | 2.62 h | 26.6 GB |
| B | Height | compute-153 | 2.56 h | 23.4 GB |

Five of six at **2.41-2.66 h (median 2.56)**; the sixth ran at half chunk-rate on **compute-053**
(6.4 vs ~14 chunks/h) — a slow node, not a heavier workload. **Add compute-053 to `--exclude`**
alongside compute-058 and compute-175.

**The unit is one donor-trait, not one trait.** 40 traits x 7 donors = 280 tasks per arm.
At the 2.56 h median that is ~717 task-h per arm, ~1,434 task-h for both. This is ~7x arm 0's
per-trait cost per arm: arm 0 covered all seven donors in 1-3.5 h with 2 annotation columns, while
the conditional regression carries 52-53 columns and costs 2.5-3 h *per donor-trait*.

MaxRSS 23.4-27.0 GB means the script's original `--mem=24G` would have OOM-killed all six tasks.

## Pilot result (Br6471 only; single-sample Cauchy, NOT the 7-donor values)

Spearman rank correlation of the 14 domain p-values (CLA absent from Br6471):

| trait | test1 vs test2 | test1 vs test3 | test2 vs test3 |
|---|---|---|---|
| SCZ | 0.991 | 0.974 | 0.982 |
| MDD | 0.916 | 0.890 | 0.969 |
| Height | 0.521 | 0.525 | 0.996 |

**Rankings hold for brain traits; magnitudes attenuate heavily.** Median shift across the 10
enriched domains: SCZ +6.7 orders (test 2) / +10.7 (test 3); MDD +4.2 / +6.0. SCZ's BM goes
3.0e-24 -> 1.1e-16 -> 1.6e-10 and 10 of 14 domains stay under p<0.05 in both conditional arms.
Test 3 attenuates consistently MORE than test 2 while preserving the same order — the neuronal
axis is partly collinear with the domain structure, as expected since both derive from the same
expression data.

**Height is the informative case.** Unconditioned it shows p<0.05 in **13 of 14 domains**
(Endothelial 1.4e-15) — not credible for a body-height GWAS in amygdala tissue, and exactly the
baseline inadequacy conditional analysis exists to catch. Conditioned, the neuronal and
white-matter domains collapse to p ~ 0.92-1.0 (CHAT 0.31 test 2 / 0.48 test 3) and only
Endothelial (1.7e-04 / 2.3e-04) and WM.1 (0.030 / 0.103) survive. Endothelial behaves the same way
for SCZ: 2.0e-05 unconditioned -> 0.95 / 0.36 conditioned, i.e. its arm-0 enrichment was largely
functional-annotation artifact while the neuronal domains' was not.

**Implication not yet checked:** the unconditioned 40-trait results likely carry inflated
enrichment for the other non-brain traits (T2D, BMI, Stroke).

**Verdict: outcome (b)** — rankings hold, p-values weaken, which is the collinearity signature and
a real anatomical result: amygdala domain identity is substantially but not wholly captured by
functional genomic annotation. Full 40-trait sweeps for BOTH conditional arms are running per the
user's decision to publish all three tests.


## Figures (19_figures.py)

All domain-enrichment figures are generated by a single script from staged summary
tables. No cluster access, no kernel state, no artifact store -- it runs anywhere
pandas/matplotlib/scipy are available.

    python 19_figures.py            # all figures -> figures/
    python 19_figures.py --list     # show available keys
    python 19_figures.py cond_P     # regenerate one

| key | output | content |
|---|---|---|
| `cond_P` | fig_conditional_P.png | conditional-arm P heatmap, square cells |
| `forest_OR` | fig_forest_OR.png | SCZ forest (default arm) + both arms' OR heatmaps |
| `power` | fig_power_vs_enrichment.png | GWAS power vs enrichment, conditional arm only |
| `cond_3panel` | fig_conditional_3panel.png | P + SCZ forest + OR, shared axes |
| `OR_scatter` | fig_OR_scatter.png | default vs conditional OR concordance |
| `two_arm_P` | fig_tall_P.png | both arms, P, shared colour scale |
| `two_arm_OR` | fig_heatmaps_OR.png | both arms, OR, shared symmetric scale |

### Inputs (figdata/)

| file | content |
|---|---|
| test1_cauchy_matrix.csv | 40 traits x 15 domains, Cauchy P, default baseline |
| test2_cauchy_matrix.csv | same, functional-conditioned arm |
| domain_vs_rest_OR.csv | Mantel-Haenszel domain-vs-rest OR, both arms, long format |
| trait_groups.csv | trait -> phenotype group |
| neuro_pct_all.csv | per-donor, per-domain neuronal-marker % of counts |
| chi2_single.json | power metrics for the 15 single-GWAS phenotypes |
| chi2_full.json | power metrics for the 25 candidate GWAS of the 10 replicate phenotypes |
| gwas_selection_table.csv | 25 candidate GWAS for the 10 phenotypes that had more than one; `selected==True` marks the 10 winners. The 15 single-GWAS phenotypes are absent -- `load()` adds them to form the 25-trait subset |

### The two statistics

**P** -- Cauchy combination across spots within a donor, then across donors (gsMap's own
`run_cauchy_combination`). Absolute enrichment: *is this domain enriched at all?*

**OR** -- per donor, spots are called trait-associated at FDR<0.05; a 2x2 table
(in-domain x called) is pooled across donors by Mantel-Haenszel with a Haldane-Anscombe
0.5 correction and Robins-Breslow-Greenland variance. Relative preference: *are called
spots concentrated HERE versus the rest of the amygdala?* The OR panel shows structure
the P panel cannot express -- white matter is not merely non-significant but strongly
depleted (median log2 OR -7.7 for WM.1).

### Conventions the script enforces

  * The 25-trait subset = 10 best-powered picks from the replicate sets (per
    `gwas_selection_table.csv`) + the 15 phenotypes with only one GWAS. The selection
    table alone contains only the 10 -- both halves are needed.
  * Row and column orderings are recomputed **per arm** from mean signal. The default and
    conditional arms give different orderings; do not reuse one for the other.
  * Traits whose FDR call saturates (<2% or >98% of spots called) have undefined OR and
    are drawn as **white** cells with a thin border -- never grey, which would read as a
    low value. Six such traits in the conditional arm, two in the default.
  * CLA is present in **one donor**; AI and LA in six of seven. Marked "n/7" in teal and
    named in every caption.
  * The P colour scale is anchored at P=0.05 (white), so the grey/orange boundary is the
    significance threshold rather than a midpoint of the data range.
  * The shared two-arm P scale spans **both** arms (default reaches -log10 P 26, the
    conditional ~16); a vmax from the conditional alone saturates the default panel.
  * Every figure runs a text-collision check and prints the overlap count. The one
    residual overlap in `cond_3panel` is a long trait label meeting the group legend.

### GWAS power: use mean chi2, not N

`fig_power_vs_enrichment` plots mean genome-wide chi2, NOT median sample size. N is a poor
proxy here: MDD has 685k samples at chi2 1.85 while SCZ has 126k at 2.03, and N correlates
with enrichment far more weakly (rho +0.37, P 0.07) than chi2 does (+0.77, P 7e-06,
default arm).

The two chi2 files are a trap worth naming. `chi2_full.json` covers the 25 *candidate* GWAS
of the 10 replicate phenotypes; `chi2_single.json` covers the 15 single-GWAS phenotypes.
Both files, and `gwas_selection_table.csv`, involve the numbers 25 and 15 in different
roles -- using the selection table alone to define the analysis subset silently yields the
wrong 25 traits (it contains no Height, PTSD, ADHD, ...). `load()` unions the correct two
halves and asserts power metrics exist for every subset trait.

What the figure shows: power predicts enrichment MAGNITUDE (rho +0.51 conditional arm,
+0.77 default), but not WHICH domains rank highest -- after removing each trait's scale,
0 of 15 domains correlate with power after Bonferroni. That asymmetry is why the OR panels,
being scale-free, carry the spatial claim while the P panels carry significance.

### Adding the neuronal-axis arm (test 3)

Export its Cauchy matrix as `figdata/test3_cauchy_matrix.csv` and append its OR rows to
`domain_vs_rest_OR.csv` with `arm='test3'`; then extend `load()` with a `P3`/`O3` pair
following the existing `test1`/`test2` blocks.


## Spatial PDFs (20_spatial_pdfs.py)

One PDF per trait, showing per-spot enrichment across all donors and all three
conditioning arms side by side.

    python 20_spatial_pdfs.py              # the 25-trait subset
    python 20_spatial_pdfs.py SCZ MDD      # named traits
    python 20_spatial_pdfs.py --all-traits # all 40
    sbatch 20_run_spatial_pdfs.sh          # full set on the cluster, ~15 s/trait

Output: `figures/spatial/spatial_<TRAIT>.pdf`, ~4 MB each. Each is a 4 x 7 grid --
anatomical domains on the top row for reference, then arm 0 (default), arm A
(functional), arm B (functional + neuronal axis); one column per donor in
`samples.txt` order.

### Inputs (read directly from /dcs04, nothing staged)

| path | content |
|---|---|
| `{WORK}/ST/{donor}.h5ad` | `obsm['spatial']` coordinates + `obs['BS_k16_Semisupervised_wAI']` |
| `{arm}/{donor}/spatial_ldsc/{donor}_{trait}.csv.gz` | per-spot beta/se/z/p |

Domain labels come from the ST h5ad, NOT from `{donor}_domains.csv` -- that file exists
for Br6471 only (a pilot leftover) and reading it fails for the other six donors.

### Colour scale

Per-spot z on a symmetric diverging scale, clipped at the 99.5th percentile of |z|,
**shared across all 21 z panels within a trait but not between traits**. Traits differ in
power by orders of magnitude (SCZ vlim 9.14 vs Height 5.44), so one global scale would
render everything below the top few traits blank. Panels are therefore comparable within a
PDF and NOT between PDFs; every caption states this.

### Domain coverage, and why two numbers circulate

Coverage is derived per PDF from the annotation, never hardcoded. From the h5ad obs
column: **CLA 1/7 donors, LA 6/7, all others 7/7**.

`AI` is a deliberate discrepancy worth knowing about. It is *annotated* in all seven
donors, but Br9192 has only **2 AI spots**, so any table with a minimum-spot-count
threshold drops it -- including `neuro_pct_all.csv`, and hence the "AI 6/7" marks in the
`19_figures.py` heatmap captions. Both numbers are correct about different things:
annotated in 7, usefully quantified in 6. The spatial PDFs plot raw annotation, so they
report 7/7.

### Guards

  * Asserts the annotation column exists in every donor's h5ad.
  * Asserts a non-empty spot-id intersection per panel -- a divergence in id conventions
    would otherwise render as an empty panel that reads like a null result.
  * Panels with a partial intersection are summarized in the caption (donor + count of
    missing spots) rather than marked per panel.
  * Traits with no per-spot output in any arm are skipped with a message, not a crash.
