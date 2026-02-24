#!/bin/bash
#SBATCH --job-name=cat_PE_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=350G
#SBATCH --time=72:00:00
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/cat_PE_map_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/cat_PE_map_%j.err

set -euo pipefail

module load bowtie2
module load samtools

IDX="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/bowtie_index_unique/MG_gene_catalog90.unique"

READDIR="/project/bcampb7/camplab/0_Raw_Data/Metatrasnscriptome_Estury/metatranscriptome_all"
R1_SUFFIX="_RNA_1.fq"
R2_SUFFIX="_RNA_2.fq"

OUTDIR="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/03_mapping_MT_paired_unique_noGFF"
mkdir -p "$OUTDIR"/{bam,idxstats,logs}

samples=$(ls ${READDIR}/*${R1_SUFFIX} | sed "s/${R1_SUFFIX}$//" | xargs -n1 basename)
echo "Paired samples: $(echo "$samples" | wc -l)"

for s in $samples; do
  R1="${READDIR}/${s}${R1_SUFFIX}"
  R2="${READDIR}/${s}${R2_SUFFIX}"

  echo "==== Mapping PE: $s ===="
  bowtie2 --very-sensitive -x "$IDX" -1 "$R1" -2 "$R2" \
    -p ${SLURM_CPUS_PER_TASK} --no-unal \
    2> "$OUTDIR/logs/${s}.bowtie2.log" \
  | samtools view -@ ${SLURM_CPUS_PER_TASK} -bS - \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "$OUTDIR/bam/${s}.sorted.bam"

  samtools index "$OUTDIR/bam/${s}.sorted.bam"

  # Per-reference mapped read counts (good for CPM/RPKM later)
  samtools idxstats "$OUTDIR/bam/${s}.sorted.bam" > "$OUTDIR/idxstats/${s}.idxstats.tsv"
done

echo "DONE (paired mapping + idxstats)"