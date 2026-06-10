#!/bin/bash

BASE=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes
GENOMES=/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/struo2_jan20/humann3_genome_analysis/genomes_top2000/genomes_fasta

mkdir -p ${BASE}/{00_logs,01_lists,02_proteins,03_dbcan,04_tcdb,05_summary}

find ${GENOMES} -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) \
  | sort > ${BASE}/01_lists/genomes_fasta.list

wc -l ${BASE}/01_lists/genomes_fasta.list
