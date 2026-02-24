#!/bin/bash
#SBATCH --job-name=MT_PE_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=350G
#SBATCH --time=72:00:00
#SBATCH -o /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/MT_PE_map_%j.out
#SBATCH -e /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/00_logs/MT_PE_map_%j.err

set -euo pipefail

module load bowtie2 samtools

IDX="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog/bowtie_index_unique/MG_gene_catalog90.unique"

READDIR="/project/bcampb7/camplab/0_Raw_Data/Metatrasnscriptome_Estury/MT_trimmed_50/links_all_MT/after_repair_sh"

OUTDIR="/project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/03_mapping_MT_repaired_paired_unique"
mkdir -p "$OUTDIR"/{bam,idxstats,logs,tmp}
export TMPDIR="$OUTDIR/tmp"

# Find unique sample bases from repaired R1 files (RNA1 + RNA2)
# Example base: CP_Spr01G08_RNA1  (from R1.tr50_CP_Spr01G08_RNA1.fq)
bases=$(ls -1 "${READDIR}"/R1.tr50_*_RNA*.fq 2>/dev/null \
  | xargs -n1 basename \
  | sed 's/^R1\.tr50_//' \
  | sed 's/\.fq$//' \
  | sort -u)

echo "Found bases: $(echo "$bases" | wc -l)"

for b in $bases; do
  R1="${READDIR}/R1.tr50_${b}.fq"
  R2="${READDIR}/R2.tr50_${b}.fq"
  [[ -s "$R1" ]] || { echo "Missing R1: $R1"; continue; }
  [[ -s "$R2" ]] || { echo "Missing R2: $R2"; continue; }

  echo "==== Mapping MT PE: $b ===="
  bowtie2 --very-sensitive -x "$IDX" -1 "$R1" -2 "$R2" \
    -p ${SLURM_CPUS_PER_TASK} --no-unal \
    2> "$OUTDIR/logs/${b}.bowtie2.log" \
  | samtools view -@ ${SLURM_CPUS_PER_TASK} -bS - \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -T "$TMPDIR/${b}" \
      -o "$OUTDIR/bam/${b}.sorted.bam"

  samtools index "$OUTDIR/bam/${b}.sorted.bam"
  samtools idxstats "$OUTDIR/bam/${b}.sorted.bam" > "$OUTDIR/idxstats/${b}.idxstats.tsv"

  tail -n 3 "$OUTDIR/logs/${b}.bowtie2.log" || true
done

echo "DONE: MT repaired paired mapping"

