#!/bin/bash
#SBATCH --job-name=eggnog_catalog
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=camplab
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/eggnog_catalog_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/eggnog_catalog_%j.err

export EGGNOG_DATA_DIR=/project/bcampb7/camplab/AL_JJ_oct23/databases/eggnog

# Load eggnog-mapper module or activate conda environment
# module load eggnog-mapper
# conda activate eggnog_env

OUTDIR="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/05_eggnog"
mkdir -p ${OUTDIR}

emapper.py \
  -i /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/MG_gene_catalog90.faa \
  --itype proteins \
  -o ${OUTDIR}/MG_catalog \
  --output_dir ${OUTDIR} \
  --cpu 32 \
  --data_dir ${EGGNOG_DATA_DIR}

echo "eggNOG annotation complete"
