#!/bin/bash
# 21_convert_ptsd_f3.sh -- PGC-PTSD 2024 freeze -> LDSC/gsMap sumstats.
#
# PROVENANCE: identifiers below are quoted from the user's own written record of the
# download session (relayed 2026-08-15); NOT verified against the file or the source here.
# Confirm before citing in a manuscript: PGC "ptsd2024", figshare 26349322,
# PubMed 38637617, Nature Genetics. File landed in scDRS/gwas_new/ in that earlier session.
#
# WHY THIS EXISTS: the PTSD.gz that the 40-trait sweeps consumed is a stale freeze
# (1,169,745 SNPs, mean chi2 1.079, ZERO genome-wide significant loci). PTSD's weak
# gsMap result is an input artifact, not biology. The 2024 freeze has 7,193,585 SNPs
# and 7,108 genome-wide significant SNPs.
#
# The source file is named .vcf.gz but is NOT VCF -- it is a tab-delimited table:
#   #CHROM  ID  POS  A1  A2  FREQ  NEFF  Z  P  DIRE
# gsMap's regression_read.py selects columns BY NAME (usecols=["SNP","Z","N"]), so we
# emit exactly the header shape the existing 40 files use: SNP N Z A1 A2.
#   SNP <- ID,  N <- NEFF (effective N, already per-SNP),  Z <- Z
#
# DO NOT write into /dcs04/lieber/shared/statsgen/... -- gsMap's PTSD.gz is a symlink
# into that shared tree and other groups' pipelines read it. Output stays in our own
# GWAS dir under a NEW name so both freezes remain side by side.
# NOTE: no `pipefail`. Several steps below are `zcat ... | head -1`, where head exits
# early and SIGPIPEs zcat (rc 141). With pipefail+errexit that aborts the whole script.
set -e

SRC="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/scDRS/gwas_new/eur_ptsd_pcs_v4_aug3_2021.vcf.gz"
OUTDIR="/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/Visium/16_LDSC/gsMap/GWAS"
OUT="$OUTDIR/PTSD_F3.gz"

[ -f "$SRC" ] || { echo "FATAL: source missing: $SRC" >&2; exit 1; }
[ -e "$OUT" ] && { echo "FATAL: $OUT exists; refusing to overwrite" >&2; exit 1; }

# Resolve column positions from the header rather than assuming order.
hdr=$(zcat "$SRC" 2>/dev/null | grep -v '^##' | { IFS= read -r l; printf '%s' "$l"; } )
idx () { echo "$hdr" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
c_id=$(idx ID); c_a1=$(idx A1); c_a2=$(idx A2); c_n=$(idx NEFF); c_z=$(idx Z)
for v in c_id c_a1 c_a2 c_n c_z; do
  [ -n "${!v}" ] || { echo "FATAL: could not locate column for $v" >&2; exit 1; }
done
echo "columns: ID=$c_id A1=$c_a1 A2=$c_a2 NEFF=$c_n Z=$c_z"

zcat "$SRC" | grep -v '^##' | awk -v OFS='\t' \
  -v i="$c_id" -v a="$c_a1" -v b="$c_a2" -v n="$c_n" -v z="$c_z" '
  NR==1 { print "SNP","N","Z","A1","A2"; next }
  $i=="" || $z=="" || $n=="" { skip++; next }        # drop incomplete rows
  $z=="NA" || $n=="NA" || $z+0!=$z { skip++; next }   # drop non-numeric Z
  { print $i,$n,$z,$a,$b; kept++ }
  END { printf("kept=%d skipped=%d\n", kept, skip) > "/dev/stderr" }
' | gzip > "$OUT.tmp"

mv "$OUT.tmp" "$OUT"
echo "wrote $OUT"
# QC: row count, header, and the two numbers that prove this is the new freeze.
zcat "$OUT" 2>/dev/null | { head -2 || true; }
echo "n_snp: $(zcat "$OUT" | tail -n +2 | wc -l)"
zcat "$OUT" | tail -n +2 | awk '{s+=$3*$3; if($3*$3>m)m=$3*$3; n++}
  END{printf("mean_chi2=%.4f max_chi2=%.1f median_N~%s\n", s/n, m, "see NEFF")}'
