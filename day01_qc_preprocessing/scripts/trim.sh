#!/bin/bash
# Usage:
#   bash trimmomatic.sh sample_R1.fastq.gz sample_R2.fastq.gz

set -euo pipefail

FORWARD=$1
REVERSE=$2

BASE=$(basename "${FORWARD}" _R1.fastq.gz)

trimmomatic PE -threads 8 \
  "${FORWARD}" "${REVERSE}" \
  "paired_${BASE}_R1.fastq.gz" "unpaired_${BASE}_R1.fastq.gz" \
  "paired_${BASE}_R2.fastq.gz" "unpaired_${BASE}_R2.fastq.gz" \
  ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36
#This is a simple wrapper script you can keep for interactive runs
