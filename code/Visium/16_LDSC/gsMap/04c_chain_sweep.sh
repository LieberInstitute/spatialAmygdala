#!/bin/bash
#SBATCH --job-name=gsmap_chain
#SBATCH --output=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/chain_%j.out
#SBATCH --error=/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/logs/chain_%j.err
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --dependency=afterok:34693837
# ---------------------------------------------------------------------------
# Submits the 40-trait sweep only after the 7-way cache pass succeeds.
# SLURM enforces the ordering, so nothing has to poll.
# ---------------------------------------------------------------------------
set -uo pipefail
cd /dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap
echo "cache pass finished; verifying caches before submitting sweep"; date
MISSING=0
while read -r s; do
  [ -z "$s" ] && continue
  if [ ! -d "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/code/Visium/16_LDSC/gsMap/../../../../processed-data/Visium/16_LDSC/gsMap/$s/latent_to_gene" ]; then
    echo "MISSING latent_to_gene: $s"; MISSING=1
  fi
done < samples.txt
if [ "$MISSING" -ne 0 ]; then
  echo "ERROR: caches incomplete -- NOT submitting sweep"; exit 1
fi
echo "all 7 caches present; submitting 04b"
sbatch 04b_gsmap_trait_array.sh
date
