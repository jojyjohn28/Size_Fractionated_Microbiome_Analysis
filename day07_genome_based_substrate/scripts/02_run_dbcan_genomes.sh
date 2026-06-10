#!/bin/bash
#SBATCH --job-name=dbcan_genomes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --array=1-1824%20
#SBATCH -o /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/dbcan_%A_%a.out
#SBATCH -e /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/00_logs/dbcan_%A_%a.err

source ~/.bashrc
conda activate dbcan_env

BASE=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes
DBCAN_DB=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/databases/dbcan
LIST=${BASE}/01_lists/genomes_fasta.list

GENOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST})
GENOME_ID=$(basename ${GENOME})
GENOME_ID=${GENOME_ID%.fna}
GENOME_ID=${GENOME_ID%.fa}
GENOME_ID=${GENOME_ID%.fasta}

PROTEINS=${BASE}/02_proteins/${GENOME_ID}/${GENOME_ID}.faa
OUTDIR=${BASE}/03_dbcan/${GENOME_ID}

mkdir -p ${OUTDIR}

echo "Running dbCAN for ${GENOME_ID}"
echo "Proteins: ${PROTEINS}"

if [ ! -s ${PROTEINS} ]; then
    echo "ERROR: protein file missing or empty"
    exit 1
fi

run_dbcan \
  ${PROTEINS} \
  protein \
  --out_dir ${OUTDIR} \
  --db_dir ${DBCAN_DB} \
  --hmm_cpu 16 \
  --tools hmmer

if [ -f ${OUTDIR}/overview.txt ]; then
    echo "CAZyme hits:"
    grep -v "^#" ${OUTDIR}/overview.txt | wc -l
else
    echo "WARNING: overview.txt not created"
fi
