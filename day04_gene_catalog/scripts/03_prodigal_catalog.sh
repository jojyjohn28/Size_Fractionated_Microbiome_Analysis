#!/bin/bash
#SBATCH --job-name=prodigal_catalog
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250
#SBATCH --time=336:00:00
#SBATCH --partition=camplab
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/prodigal_catalog_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/prodigal_catalog_%j.err

module load prodigal/2.6.3

cd /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog

prodigal -i MG_gene_catalog90.fna \
  -a MG_gene_catalog90.faa \
  -d MG_gene_catalog90.genes.fna \
  -f gff \
  -o MG_gene_catalog90.gff \
  -p meta

echo "Prodigal complete on catalog"
