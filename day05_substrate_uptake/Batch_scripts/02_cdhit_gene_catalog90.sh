#!/bin/bash
#SBATCH --job-name=cdhit90
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --partition=camplab
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/cdhit90_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/cdhit90_%j.err

module load cd-hit/4.8.1

cd /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog

cd-hit-est \
  -i all_MG_genes.fna \
  -o MG_gene_catalog90.fna \
  -c 0.90 \
  -n 8 \
  -G 0 \
  -aS 0.90 \
  -T 32 \
  -M 0
