#### Day 1 — From Raw Reads to Clean Data

This folder contains the preprocessing workflow for the **Size-Fractionated Microbiome Analysis** series.

**Goal:** convert raw paired-end shotgun metagenomic reads into clean, analysis-ready FASTQ files for downstream profiling (Kaiju, MetaPhlAn, mOTUs) and optional genome-resolved workflows.

#### Steps

**Quality control (FastQC)**
➡️ See: fastqc.md

**Trimming & filtering (Trimmomatic or Cutadapt)**
➡️ See: trimming.md

**Post-trimming QC (FastQC again)**
Confirm adapters and low-quality tails are reduced.

#### Scripts

All SLURM and bash scripts are in: scripts/

fastqc_loop.sh

trimmomatic.slurm

cutadapt.slurm
