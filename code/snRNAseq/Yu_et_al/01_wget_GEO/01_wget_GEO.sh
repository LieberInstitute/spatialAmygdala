#!/bin/bash
#SBATCH --job-name=dl_Yu_fastq
#SBATCH --output=logs/dl_Yu_fastq.txt
#SBATCH --error=logs/dl_Yu_fastq.txt
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=mtotty2@jh.edu

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

# ======== Setup ========
# Yu et al. (2023) Cell Discovery — GSE195445
# "Molecular and cellular evolution of the amygdala across species
#  analyzed by single-nucleus transcriptome profiling"
# Human amygdala snRNA-seq (10x Chromium v3, NovaSeq 6000)
#   AMH1: GSM5836863 / SRX13924613 (8 runs)
#   AMH2: GSM5836864 / SRX13924614 (3 runs)
#   AMH3: GSM5836865 / SRX13924615 (3 runs)

module load sra-toolkit
module list

OUTDIR="../../../../raw-data/snRNAseq/Yu_et_al"
mkdir -p ${OUTDIR}
mkdir -p logs

THREADS=4
TMPDIR="${OUTDIR}/tmp"
mkdir -p ${TMPDIR}

# ======== SRR accessions per sample ========
AMH1_RUNS="SRR17762195 SRR17762196 SRR17762197 SRR17762198 SRR17762199 SRR17762200 SRR17762201 SRR17762202"
AMH2_RUNS="SRR17762203 SRR17762204 SRR17762205"
AMH3_RUNS="SRR17762206 SRR17762207 SRR17762208"

# ======== Download FASTQs ========
for SAMPLE in AMH1 AMH2 AMH3; do
    echo ""
    echo "======== Downloading ${SAMPLE} ========"
    date

    SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
    mkdir -p ${SAMPLE_DIR}

    # Get the run list for this sample
    RUNS_VAR="${SAMPLE}_RUNS"
    RUNS="${!RUNS_VAR}"

    for SRR in ${RUNS}; do
        # Skip if compressed FASTQs already exist for this run
        if ls ${SAMPLE_DIR}/${SRR}*.fastq.gz 1>/dev/null 2>&1; then
            echo "  Skipping ${SRR} — FASTQs already exist"
            continue
        fi

        echo "  Prefetching ${SRR}..."
        prefetch ${SRR} --max-size 50G

        echo "  Extracting FASTQs for ${SRR}..."
        fasterq-dump ${SRR} \
            --outdir ${SAMPLE_DIR} \
            --temp ${TMPDIR} \
            --threads ${THREADS} \
            --mem ${THREADS}G \
            --split-files \
            --include-technical

        if ls ${SAMPLE_DIR}/${SRR}*.fastq 1>/dev/null 2>&1; then
            echo "  Compressing ${SRR} FASTQs..."
            gzip ${SAMPLE_DIR}/${SRR}*.fastq
        else
            echo "  WARNING: No FASTQ files found for ${SRR} — fasterq-dump may have failed"
        fi

        # Clean up prefetch cache to save disk space
        rm -rf ${SRR}

        echo "  Done with ${SRR}"
    done

    echo "  Files for ${SAMPLE}:"
    ls -lh ${SAMPLE_DIR}/
done

rm -rf ${TMPDIR}

echo ""
echo "**** All downloads complete ****"
echo "**** Final file listing ****"
ls -lhR ${OUTDIR}/

echo "**** Job ends ****"
date
