#!/bin/bash
#SBATCH --job-name=prodigal_MG
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=camplab
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/prodigal_MG_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/prodigal_MG_%j.err

module load prodigal/2.6.3

MG_ASM="/project/bcampb7/camplab/0_Raw_Data/Metagenome_estury/Assembly_All/MGs"
OUT="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/prodigal_MGs"
mkdir -p "$OUT"

for fa in ${MG_ASM}/*_contigs.fasta; do
  sample=$(basename "$fa" _contigs.fasta)
  od="${OUT}/${sample}"
  mkdir -p "$od"

  prodigal -i "$fa" -p meta \
    -a "${od}/${sample}.faa" \
    -d "${od}/${sample}.fna" \
    -f gff -o "${od}/${sample}.gff"
done
