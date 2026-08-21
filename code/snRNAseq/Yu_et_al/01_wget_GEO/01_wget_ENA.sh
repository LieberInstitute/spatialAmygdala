#!/bin/bash
#SBATCH --partition=transfer
#SBATCH --job-name=dl_Yu_fastq_ena
#SBATCH --output=logs/dl_Yu_fastq_ena.txt
#SBATCH --error=logs/dl_Yu_fastq_ena.txt
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
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

OUTDIR="../../../../processed-data/snRNAseq/Yu_et_al"
mkdir -p ${OUTDIR}
mkdir -p logs

# ======== SRR accessions per sample ========
AMH1_RUNS="SRR17762195 SRR17762196 SRR17762197 SRR17762198 SRR17762199 SRR17762200 SRR17762201 SRR17762202"
AMH2_RUNS="SRR17762203 SRR17762204 SRR17762205"
AMH3_RUNS="SRR17762206 SRR17762207 SRR17762208"

# ======== Download FASTQs directly from ENA ========
# ENA hosts pre-built FASTQ.gz files — no SRA conversion needed
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

        echo "  Querying ENA for ${SRR} FASTQ URLs..."
        FASTQ_URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp" \
            | tail -1 | cut -f2)

        if [ -z "${FASTQ_URLS}" ]; then
            echo "  WARNING: No FASTQ URLs found on ENA for ${SRR}"
            continue
        fi

        for URL in $(echo ${FASTQ_URLS} | tr ';' ' '); do
            echo "  Downloading ftp://${URL}..."
            wget -q -P ${SAMPLE_DIR} "ftp://${URL}"
        done

        echo "  Done with ${SRR}"
    done

    echo "  Files for ${SAMPLE}:"
    ls -lh ${SAMPLE_DIR}/
done

echo ""
echo "**** All downloads complete ****"
echo "**** Final file listing ****"
ls -lhR ${OUTDIR}/

echo "**** Job ends ****"
date