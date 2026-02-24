#### Day 05 – Substrate Uptake & Utilization Analysis

This directory contains scripts for profiling substrate degradation and uptake machinery in size-fractionated estuarine microbiomes (FL vs PA). The workflow integrates gene prediction, functional annotation, read mapping, and heatmap visualization.

#### 📁 Directory Structure

1. Batch_scripts/

**8HPC job scripts for large-scale processing:**

1. 01_prodigal_MGs.sh – Gene prediction from metagenomes (Prodigal)

2. 02_cdhit_gene_catalog90.sh – Gene catalog clustering at 90% identity (CD-HIT)

3. 01_map_paired_to_catalog.slurm – Mapping metagenomic reads to gene catalog

4. run_MT_repaired_PE_map.slurm – Mapping metatranscriptomic reads

5. 01_tcdb_search_array.slurm – TCDB search for ABC transporters

6. 02_dbcan_search_array.slurm – dbCAN search for CAZymes

7. Main Scripts

8. make_counts_matrix.py
   Generates gene count matrices from mapping outputs.

9. TonB_heatmap.R
   Produces ordered heatmaps (MG + MT) for TonB-dependent transporters.

#### 🔬 Functional Targets

This day focuses on:

● CAZymes (dbCAN)

● ABC transporters (TCDB)

● TonB-dependent transporters

● Expression vs genomic potential comparison

● FL (<0.8 μm) vs PA (>0.8 μm) structuring

#### 🎯 Goal

To quantify and visualize substrate degradation and uptake strategies across Chesapeake and Delaware estuaries, integrating metagenomic potential with metatranscriptomic expression.

Read today's blog for more deatails : https://jojyjohn28.github.io/blog/day5-substrate-uptake-patterns/
