#!/bin/bash
#SBATCH --job-name=tcdb_search
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --array=1-36
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/tcdb_%A_%a.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/tcdb_%A_%a.err

module load diamond

SAMPLE_LIST="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/sample_list.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

PROTEINS="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/prodigal_MGs/${SAMPLE}/${SAMPLE}.faa"
TCDB_DB="/project/bcampb7/camplab/AL_JJ_oct23/databases/tcdb/tcdb"
OUTDIR="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/06_function_sets/tcdb_results"

mkdir -p ${OUTDIR}

echo "================================================"
echo "Processing sample: ${SAMPLE}"
echo "Input: ${PROTEINS}"
echo "Protein count: $(grep -c '^>' ${PROTEINS})"
echo "================================================"
date

diamond blastp \
  --db ${TCDB_DB} \
  --query ${PROTEINS} \
  --out ${OUTDIR}/${SAMPLE}_tcdb.tsv \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  --evalue 1e-5 \
  --query-cover 50 \
  --subject-cover 50 \
  --threads 16 \
  --sensitive \
  --block-size 4.0 \
  --index-chunks 1

echo "Sample ${SAMPLE} complete"
echo "Hits found: $(wc -l < ${OUTDIR}/${SAMPLE}_tcdb.tsv)"
date
EOF
