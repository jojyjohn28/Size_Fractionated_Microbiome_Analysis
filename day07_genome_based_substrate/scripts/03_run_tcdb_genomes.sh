#!/bin/bash
#SBATCH --job-name=tcdb_genomes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --array=1-1800%30
#SBATCH -o /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/tcdb_%A_%a.out
#SBATCH -e /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/tcdb_%A_%a.err

module load diamond

BASE=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes
TCDB_DB=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/databases/tcdb/tcdb
LIST=${BASE}/01_lists/genomes_fasta.list

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST})
GENOME_ID=$(basename ${GENOME})
GENOME_ID=${GENOME_ID%.fna}
GENOME_ID=${GENOME_ID%.fa}
GENOME_ID=${GENOME_ID%.fasta}

PROTEINS=${BASE}/02_proteins/${GENOME_ID}/${GENOME_ID}.faa
OUTDIR=${BASE}/04_tcdb

mkdir -p ${OUTDIR}

echo "Running TCDB search for ${GENOME_ID}"
echo "Proteins: ${PROTEINS}"

if [ ! -s ${PROTEINS} ]; then
    echo "ERROR: protein file missing or empty"
    exit 1
fi

diamond blastp \
  --db ${TCDB_DB} \
  --query ${PROTEINS} \
  --out ${OUTDIR}/${GENOME_ID}_tcdb.tsv \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --evalue 1e-5 \
  --query-cover 50 \
  --subject-cover 50 \
  --threads 16 \
  --sensitive \
  --block-size 4.0 \
  --index-chunks 1

echo "Hits found:"
wc -l ${OUTDIR}/${GENOME_ID}_tcdb.tsv
