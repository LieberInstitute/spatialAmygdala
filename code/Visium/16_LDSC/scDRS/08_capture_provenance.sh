#!/bin/bash
#SBATCH --job-name=scdrs_prov
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS/logs/08_provenance.%j.log
# -----------------------------------------------------------------------------
# 08_capture_provenance.sh -- freeze everything needed to reproduce this
# pipeline: exact package versions, exact input files with checksums, exact
# tool versions, and checksums of the code itself.
#
#   sbatch 08_capture_provenance.sh
#
# Re-runnable at any time; each run OVERWRITES the files in provenance/ with a
# current snapshot and APPENDS a dated stamp to provenance/history.txt, so the
# record of when each snapshot was taken is never lost.
#
# WHY THIS EXISTS
#   The conda env and the MAGMA reference panel live in scratch, which is not
#   backed up and not part of the project tree. If scratch is ever cleared,
#   provenance/env_spec.yml + provenance/env_pip_freeze.txt are what let you
#   rebuild a byte-comparable environment, and input_manifest.txt is what lets
#   you prove the inputs were the same files.
#
# Writes ONLY to $CODE/provenance/. Reads everything else read-only.
# -----------------------------------------------------------------------------
set -uo pipefail
# NOTE: under sbatch, SLURM copies this script to a spool directory, so "$0"
# does NOT resolve to this file's real location -- $(dirname $(readlink -f $0))
# yields /var/spool/slurm/... and config.sh is not found there. Hardcode the
# code directory instead. Override with CODE_DIR=... if you relocate the tree.
CODE_DIR=${CODE_DIR:-/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/scDRS}
source "$CODE_DIR/config.sh"

PROV=$CODE/provenance
mkdir -p "$PROV" "$CODE/logs"

echo "**** Job starts ****"; date; hostname

STAMP=$(date '+%Y-%m-%d %H:%M:%S %Z')

# --- 1. environment ----------------------------------------------------------
echo "[1/5] environment"
if [ -x "$PY" ]; then
    activate_env
    conda env export -p "$ENVDIR"        > "$PROV/env_spec.yml"        2>/dev/null
    conda list      -p "$ENVDIR"         > "$PROV/env_conda_list.txt"  2>/dev/null
    "$ENVDIR/bin/pip" freeze             > "$PROV/env_pip_freeze.txt"  2>/dev/null
    "$PY" - > "$PROV/env_versions.txt" <<'PY'
import platform
print("python  ", platform.python_version())
for m in ("scdrs","scanpy","anndata","numpy","pandas","scipy","statsmodels","sklearn","matplotlib"):
    try:
        mod = __import__(m)
        print(f"{m:12s}", getattr(mod, "__version__", "?"))
    except Exception as e:
        print(f"{m:12s} NOT IMPORTABLE ({e.__class__.__name__})")
PY
    echo "  env captured"
else
    echo "  env not built yet ($PY missing) -- run 00_setup_env.sh first" \
        | tee "$PROV/env_versions.txt"
fi

# --- 2. external tools -------------------------------------------------------
echo "[2/5] external tool versions"
{
    echo "captured: $STAMP"
    echo "host    : $(hostname)"
    echo
    echo "== modules referenced by this pipeline =="
    echo "conda module : $CONDA_MODULE"
    echo "magma        : magma/1.10 -> /jhpce/shared/libd/core/magma/1.10/magma"
    echo "plink        : plink/1.90b"
    echo
    module load magma/1.10 2>/dev/null
    echo "== magma --version =="
    magma --version 2>&1 | head -5
    echo
    module load plink/1.90b 2>/dev/null
    echo "== plink --version =="
    plink --version 2>&1 | head -3
} > "$PROV/tool_versions.txt" 2>&1
cat "$PROV/tool_versions.txt"

# --- 3. input manifest -------------------------------------------------------
# size + mtime + md5 for every file this pipeline reads. md5 over ~3 GB takes
# a couple of minutes; that is the price of being able to prove input identity.
echo "[3/5] input manifest (md5 -- this is the slow step)"
{
    echo "captured: $STAMP"
    echo
    echo "== GWAS sumstats ($GWAS_SRC) =="
    while IFS=$'\t' read -r trait ss; do
        [ "$trait" = "TRAIT" ] && continue
        f=$GWAS_SRC/$ss
        if [ -r "$f" ]; then
            printf '%s\t%s\t%s\n' "$trait" "$ss" "$(md5sum "$f" | cut -d' ' -f1)"
        else
            printf '%s\t%s\tUNREADABLE\n' "$trait" "$ss"
        fi
    done < "$TRAITS_TSV"

    echo
    echo "== MAGMA gene locations =="
    ls -l --time-style=long-iso "$GENE_LOC_SRC"
    md5sum "$GENE_LOC_SRC"

    echo
    echo "== 1000G EUR Phase3 plink (source, per chromosome) =="
    for chr in $(seq 1 22); do
        b=$PLINK_SRC/${PLINK_PREFIX}.${chr}.bim
        [ -r "$b" ] && printf 'chr%-3s %12s bytes  %s\n' \
            "$chr" "$(stat -c%s "$b")" "$(md5sum "$b" | cut -d' ' -f1)"
    done

    echo
    echo "== spatial inputs (gsMap ST h5ad) =="
    for s in Br6471 Br6660 Br6423 Br2743 Br8325 Br9192 Br9280; do
        f=$ST_DIR/$s.h5ad
        if [ -r "$f" ]; then
            printf '%-8s %12s bytes  %s  %s\n' "$s" "$(stat -c%s "$f")" \
                "$(date -r "$f" '+%Y-%m-%d')" "$(md5sum "$f" | cut -d' ' -f1)"
        fi
    done

    echo
    echo "== source SPE that produced those h5ad =="
    # 2.06 GB -- the md5 takes ~1 min but pins the object the h5ad derive from
    SPE=$PROJ/processed-data/Visium/07_clustering/BayesSpace/MarkerGenes/spe_harmony_markers_BS_k16_Semisupervised_wAI.rds
    ls -l --time-style=long-iso "$SPE" 2>&1
    md5sum "$SPE" 2>&1

    echo
    echo "== comparison target =="
    ls -l --time-style=long-iso "$PROJ/code/Visium/16_LDSC/LDSC/ldsc_results.csv" 2>&1
    md5sum "$PROJ/code/Visium/16_LDSC/LDSC/ldsc_results.csv" 2>&1
} > "$PROV/input_manifest.txt" 2>&1
echo "  wrote $PROV/input_manifest.txt"

# --- 4. code checksums -------------------------------------------------------
echo "[4/5] code checksums"
{
    echo "captured: $STAMP"
    echo
    ( cd "$CODE" && md5sum ./*.sh ./*.py ./*.tsv ./*.md 2>/dev/null )
} > "$PROV/code_checksums.txt"
cat "$PROV/code_checksums.txt"

# --- 5. resolved configuration ----------------------------------------------
echo "[5/5] resolved config"
{
    echo "captured: $STAMP"
    echo
    echo "PROJ         = $PROJ"
    echo "CODE         = $CODE"
    echo "WORKDIR      = $WORKDIR"
    echo "SCRATCH      = $SCRATCH"
    echo "ENVDIR       = $ENVDIR"
    echo "GWAS_SRC     = $GWAS_SRC"
    echo "GENE_LOC_SRC = $GENE_LOC_SRC"
    echo "PLINK_SRC    = $PLINK_SRC"
    echo "ST_DIR       = $ST_DIR"
    echo "ANNOT_COL    = $ANNOT_COL"
    echo "MAGMA_WINDOW = $MAGMA_WINDOW"
    echo "N_CTRL       = $N_CTRL"
    echo "GS_NMAX      = $GS_NMAX"
} > "$PROV/resolved_config.txt"
cat "$PROV/resolved_config.txt"

# --- history stamp -----------------------------------------------------------
echo "$STAMP  snapshot taken on $(hostname) (job ${SLURM_JOB_ID:-manual})" \
    >> "$PROV/history.txt"

echo
ls -la "$PROV"
echo "**** Job ends ****"; date
