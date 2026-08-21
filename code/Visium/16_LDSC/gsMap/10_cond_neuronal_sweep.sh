#!/bin/bash
#SBATCH --job-name=cond_neuronal
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_neuro_%A_%a.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_neuro_%A_%a.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --exclude=compute-058,compute-175,compute-053
# -----------------------------------------------------------------------------
# 10_cond_neuronal_sweep.sh -- TEST 3 full sweep: 40 traits x 7 samples = 280 units.
# ARM=neuronal: 52 functional annotations + neuronal_axis = 53 columns.
#
# Array index K in 1..280 maps trait-major:
#   TRAIT = line ceil(K/7) of traits.txt ; SAMPLE = line ((K-1)%7)+1 of samples.txt
# so a contiguous block of 7 completes one whole trait (enables early Cauchy).
#
# THE SILENT-UNCONDITIONED TRAP (gsMap config.py:1185-1207): when the additional
# baseline dir is absent, gsMap sets use_additional_baseline_annotation=False
# with NO warning and produces UNCONDITIONED numbers at exit 0. For this arm a
# false null is the worst possible outcome, so we guard three ways:
#   exit 4 = annotation dir missing, or not 22 chroms, or not 53 annotations
#   exit 5 = log did not confirm baseline used (or explicitly denied it)
#   exit 6 = expected output file absent/empty
# Success is judged by OUTPUT FILE existence + size, NEVER by exit code.
# -----------------------------------------------------------------------------
set -uo pipefail

PROJ=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala
CODE=$PROJ/code/Visium/16_LDSC/gsMap
ARM=neuronal
WORKDIR=$PROJ/processed-data/Visium/16_LDSC/gsMap_cond_${ARM}
SCRATCH=/users/mtotty/claude_scratch/gsmap
ANNOT=BS_k16_Semisupervised_wAI
EXPECT_ANNOT=53
NS=7

K=${SLURM_ARRAY_TASK_ID:?no array task id}
TI=$(( (K - 1) / NS + 1 ))
SI=$(( (K - 1) % NS + 1 ))
TRAIT=$(sed -n "${TI}p" $CODE/traits.txt)
SAMPLE=$(sed -n "${SI}p" $CODE/samples.txt)
[ -n "$TRAIT" ]  || { echo "FATAL: no trait at line $TI";  exit 3; }
[ -n "$SAMPLE" ] || { echo "FATAL: no sample at line $SI"; exit 3; }

GWAS=$(grep -E "^${TRAIT}:[[:space:]]" $CODE/gwas_config_full.yaml | head -1 | sed 's/^[^:]*:[[:space:]]*//')
[ -n "$GWAS" ] && [ -e "$GWAS" ] || { echo "FATAL: sumstats for $TRAIT not found ('$GWAS')"; exit 3; }

CONFDIR=$WORKDIR/cond_confirm
mkdir -p "$CONFDIR" "$CODE/logs"
CONF=$CONFDIR/${SAMPLE}_${TRAIT}.conf
SPOT=$WORKDIR/$SAMPLE/spatial_ldsc/${SAMPLE}_${TRAIT}.csv.gz
CAU=$WORKDIR/$SAMPLE/cauchy_combination/${SAMPLE}_${TRAIT}.Cauchy.csv.gz

echo "**** Job starts ****"; date
echo "K=$K ARM=$ARM TRAIT=$TRAIT SAMPLE=$SAMPLE"
echo "sumstats: $GWAS"
echo "host: $(hostname)  cpus: ${SLURM_CPUS_PER_TASK:-NA}  mem: ${SLURM_MEM_PER_NODE:-NA}"

# --- guard 1: the conditioning inputs must really be there, and be THIS arm ---
ABDIR=$WORKDIR/$SAMPLE/generate_ldscore/additional_baseline
[ -d "$ABDIR" ] || { echo "FATAL(4): $ABDIR absent -> gsMap would SILENTLY run unconditioned"; exit 4; }
NM=0; NBAD=0
for ch in $(seq 1 22); do
    f=$ABDIR/baseline.${ch}.l2.M
    g=$ABDIR/baseline.${ch}.l2.ldscore.feather
    if [ -s "$f" ] && [ -s "$g" ]; then
        NM=$((NM+1))
        nf=$(awk '{print NF; exit}' "$f")
        [ "$nf" = "$EXPECT_ANNOT" ] || { echo "chr$ch has $nf annotations, expected $EXPECT_ANNOT"; NBAD=$((NBAD+1)); }
    else
        echo "chr$ch missing .M or .ldscore.feather"; NBAD=$((NBAD+1))
    fi
done
echo "additional_baseline: $NM/22 chroms present, $NBAD bad"
if [ "$NM" -ne 22 ] || [ "$NBAD" -ne 0 ]; then
    echo "FATAL(4): additional_baseline incomplete or wrong annotation count for ARM=$ARM"; exit 4
fi

# --- guard 1b: Cauchy's prerequisite must exist BEFORE the 2.5 h regression ---
# run_cauchy_combination reads {workdir}/{sample}/find_latent_representations/
# {sample}_add_latent.h5ad, which gsMap does NOT create in a conditional workdir.
# The pilot lost 2 of 6 tasks to this: regression succeeded, Cauchy then died with
# FileNotFoundError. Checking it up front turns a 2.5 h loss into an instant exit.
LATENT=$WORKDIR/$SAMPLE/find_latent_representations/${SAMPLE}_add_latent.h5ad
[ -e "$LATENT" ] || { echo "FATAL(7): $LATENT absent -> regression would succeed and Cauchy would die. Symlink it from the arm-0 workdir first."; exit 7; }
echo "cauchy prerequisite present: $LATENT -> $(stat -Lc%s "$LATENT") bytes"

module load conda/3-24.3.0
source activate $SCRATCH/envs/gsmap
export PATH="$SCRATCH/envs/gsmap/bin:$PATH"
echo "gsmap: $(which gsmap)"

RUNLOG=$WORKDIR/$SAMPLE/cond_${ARM}_${TRAIT}_run.log
mkdir -p "$(dirname "$RUNLOG")"

if [ -s "$SPOT" ]; then
    echo "spatial_ldsc output already present, skipping regression: $SPOT"
    LDSC_RC=0
else
    # No pipe into grep: a short-circuiting reader would SIGPIPE gsmap mid-run.
    gsmap run_spatial_ldsc \
        --workdir "$WORKDIR" \
        --sample_name "$SAMPLE" \
        --trait_name "$TRAIT" \
        --sumstats_file "$GWAS" \
        --num_processes 8 \
        --use_additional_baseline_annotation True 2>&1 | tee "$RUNLOG"
    LDSC_RC=${PIPESTATUS[0]}
    echo "run_spatial_ldsc rc=$LDSC_RC (non-fatal; outputs decide)"

    # --- guard 2: the log must AFFIRM conditioning and never deny it ---
    if grep -Fq "Baseline annotation is not provided" "$RUNLOG"; then
        echo "FATAL(5): log says baseline annotation NOT provided -> UNCONDITIONED"; exit 5
    fi
    if ! grep -Fq "Using additional baseline annotations" "$RUNLOG"; then
        echo "FATAL(5): log never confirmed additional baseline -> UNCONDITIONED"; exit 5
    fi
    echo "CONFIRMED: conditioning active"
fi

# --- guard 3: output must exist and be non-empty ---
[ -s "$SPOT" ] || { echo "FATAL(6): no per-spot output $SPOT (rc=$LDSC_RC)"; exit 6; }
echo "per-spot output present: $SPOT ($(stat -c%s "$SPOT") bytes)"

CAU_RC=0
if [ -s "$CAU" ]; then
    echo "per-sample Cauchy already present: $CAU"
else
    gsmap run_cauchy_combination \
        --workdir "$WORKDIR" --sample_name "$SAMPLE" \
        --trait_name "$TRAIT" --annotation "$ANNOT" || CAU_RC=$?
fi
[ -s "$CAU" ] || { echo "FATAL(6): no Cauchy output $CAU (rc=$CAU_RC)"; exit 6; }

# --- machine-readable per-unit conditioning confirmation ---
printf 'sample=%s\ttrait=%s\tarm=%s\tannot=%s\tchroms=%s\tconfirmed=%s\tspot_bytes=%s\tcauchy_bytes=%s\tldsc_rc=%s\tcauchy_rc=%s\tnode=%s\tjob=%s\n' \
    "$SAMPLE" "$TRAIT" "$ARM" "$EXPECT_ANNOT" "$NM" \
    "$(grep -Fq 'Using additional baseline annotations' "$RUNLOG" && echo yes || echo prior_run)" \
    "$(stat -c%s "$SPOT")" "$(stat -c%s "$CAU")" "$LDSC_RC" "$CAU_RC" \
    "$(hostname)" "${SLURM_ARRAY_JOB_ID:-NA}_${K}" > "$CONF"

echo "OK  K=$K $SAMPLE / $TRAIT"
echo "**** Job ends ****"; date
