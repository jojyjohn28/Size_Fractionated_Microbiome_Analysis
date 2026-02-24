#!/bin/bash
#SBATCH --job-name=dbcan_search
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=48:00:00
#SBATCH --array=1-36
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/dbcan_%A_%a.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/dbcan_%A_%a.err

source ~/.bashrc
conda activate dbcan_env

# Get sample name
SAMPLE_LIST="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/06_function_sets/sample_list.txt"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

PROTEINS="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/prodigal_MGs/${SAMPLE}/${SAMPLE}.faa"
DBCAN_DB="/project/bcampb7/camplab/Jojy/MTM2_RoL/dbcan"
OUTDIR="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/06_function_sets/dbcan_results/${SAMPLE}"

mkdir -p ${OUTDIR}

echo "================================================"
echo "Processing sample: ${SAMPLE}"
echo "Input: ${PROTEINS}"
echo "Protein count: $(grep -c '^>' ${PROTEINS} 2>/dev/null || echo 'N/A')"
echo "================================================"
date

# Run dbCAN
run_dbcan \
  ${PROTEINS} \
  protein \
  --out_dir ${OUTDIR} \
  --db_dir ${DBCAN_DB} \
  --hmm_cpu 16 \
  --tools hmmer

echo "Sample ${SAMPLE} complete"
if [ -f ${OUTDIR}/overview.txt ]; then
    echo "CAZymes found: $(grep -v "^#" ${OUTDIR}/overview.txt | wc -l)"
else
    echo "WARNING: overview.txt not created"
fi
date
