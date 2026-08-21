#!/bin/bash
#SBATCH --job-name=cond_cauchy_fix
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_cauchy_fix_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/cond_cauchy_fix_%j.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --exclude=compute-058,compute-175
set -uo pipefail
export PATH="/users/mtotty/claude_scratch/gsmap/envs/gsmap/bin:$PATH"
WD=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap_cond_functional
for T in MDD Height; do
  gsmap run_cauchy_combination --workdir "$WD" --sample_name Br6471 \
    --trait_name "$T" --annotation BS_k16_Semisupervised_wAI
  echo "rc=$? for $T"
  ls -l "$WD/Br6471/cauchy_combination/Br6471_$T.Cauchy.csv.gz" || echo "MISSING $T"
done

