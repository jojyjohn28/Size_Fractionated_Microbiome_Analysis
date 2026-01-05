#!/bin/bash
#SBATCH --job-name=cutadapt_day1
#SBATCH --partition=camplab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=cutadapt_%j.out
#SBATCH --error=cutadapt_%j.err

set -euo pipefail

# ---- Load module ----
module load cutadapt/4.9

# ---- Paths ----
IN_DIR="/project/your_project/raw_fastq"
OUT_DIR="/project/your_project/qc/cutadapt"
THREADS="${SLURM_CPUS_PER_TASK}"

# Illumina adapters
ADAPT_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPT_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

mkdir -p "${OUT_DIR}"/{trimmed,logs}

echo "[$(date)] Cutadapt started"

shopt -s nullglob
for R1 in "${IN_DIR}"/*_R1.fastq.gz; do
    sample=$(basename "$R1" _R1.fastq.gz)
    R2="${IN_DIR}/${sample}_R2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "[WARN] Missing R2 for ${sample}, skipping"
        continue
    fi

    echo "[$(date)] Processing ${sample}"

    cutadapt \
      -j "${THREADS}" \
      -a "${ADAPT_R1}" \
      -A "${ADAPT_R2}" \
      -q 15,15 \
      --minimum-length 36 \
      -o "${OUT_DIR}/trimmed/${sample}_R1.trimmed.fastq.gz" \
      -p "${OUT_DIR}/trimmed/${sample}_R2.trimmed.fastq.gz" \
      "$R1" "$R2" \
      > "${OUT_DIR}/logs/${sample}.cutadapt.log"

done
shopt -u nullglob

echo "[$(date)] Cutadapt finished"

