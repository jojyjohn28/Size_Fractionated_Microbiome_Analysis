#!/bin/bash
#SBATCH --job-name=prodigal_genomes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --array=1-1800%50
#SBATCH -o /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/prodigal_%A_%a.out
#SBATCH -e /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/prodigal_%A_%a.err

module load prodigal

BASE=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes
LIST=${BASE}/01_lists/genomes_fasta.list

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST})
GENOME_ID=$(basename ${GENOME})
GENOME_ID=${GENOME_ID%.fna}
GENOME_ID=${GENOME_ID%.fa}
GENOME_ID=${GENOME_ID%.fasta}

OUTDIR=${BASE}/02_proteins/${GENOME_ID}
mkdir -p ${OUTDIR}

echo "Processing genome: ${GENOME_ID}"
echo "Input: ${GENOME}"

prodigal \
  -i ${GENOME} \
  -a ${OUTDIR}/${GENOME_ID}.faa \
  -d ${OUTDIR}/${GENOME_ID}.ffn \
  -o ${OUTDIR}/${GENOME_ID}.gff \
  -p meta \
  -f gff

echo "Protein count:"
grep -c "^>" ${OUTDIR}/${GENOME_ID}.faa
