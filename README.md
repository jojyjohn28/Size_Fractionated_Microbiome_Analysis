🧬 Size_Fractionated_Microbiome_Analysis

An ongoing, step-by-step bioinformatics workflow for analyzing size-fractionated microbial communities (Particle-Attached vs Free-Living) using shotgun metagenomics and metatranscriptomics.

This repository documents the analysis day by day, following real research progress rather than a polished, post-hoc pipeline.
I update this README daily as new steps, scripts, and results are added.

🎯 Project Goals

Compare Particle-Attached (PA) vs Free-Living (FL) microbial communities

Quantify:

Taxonomic composition (total & active communities)

Functional potential (DNA)

Functional activity (RNA)

Integrate:

Read-based profiling

Functional annotation

Ecological statistics

Co-occurrence and redundancy analyses

📅 Analysis Timeline (Living Document)
✅ Day 1 — From Raw Reads to Clean Data

Title: Size Fractionated Microbiome Analysis — Day 1: From Raw Reads to Clean Data

Focus:

Raw FASTQ quality assessment

Adapter trimming and quality filtering

Preparing high-quality reads for downstream analysis

Key tools:

FastQC / NanoPlot

Trimmomatic / fastp

Basic QC summaries

📂 Folder: day01_qc_preprocessing/

✅ Day 2 — Kaiju Classification & Read Extraction

Title: Size Fractionated Microbiome Analysis — Day 2: Kaiju Classification and Extraction of Bacterial & Archaeal Reads

Focus:

Protein-level taxonomic classification with Kaiju

Generating taxonomy summary tables

Extracting bacterial and archaeal reads for focused downstream analyses

Key tools:

Kaiju

Custom parsing scripts

Read extraction workflows

📂 Folder: day02_kaiju_taxonomy/

✅ Day 3 — Species-Level Profiling with mOTUs

Title: Size-Fractionated Microbiome Analysis — Day 3: Species-Level Profiling with mOTUs

Focus:

Marker-gene–based species-level profiling

Total community (DNA) vs active community (RNA)

Batch processing of samples

Manual merging of mOTUs profiles

Visualization-ready abundance tables

Key tools:

mOTUs v3

Bash-based merging

R (heatmaps, stacked barplots)

📂 Folder: day03_motus_profiling/

🔜 Upcoming Steps (Planned)

Day 4: Total vs Active community comparisons (DNA vs RNA)

Day 5: Functional profiling using CAZymes & transporters

Day 6: Non-redundant gene catalog construction

Day 7: DNA:RNA ratios and activity-based niche modeling

Day 8: Co-occurrence networks (season, bay, fraction)

Day 9: Functional redundancy modeling

(Timeline may adapt as analysis evolves — this reflects real research, not a fixed tutorial.)
