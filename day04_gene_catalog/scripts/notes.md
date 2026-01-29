# Day 4: SLURM Scripts Documentation

## 📋 Overview

This document describes each SLURM script used in the gene catalog construction workflow, including inputs, outputs, and usage.

---

## Script 1: `01_prodigal_MGs.sh`

### Purpose

Predict genes from individual metagenome assemblies using Prodigal in metagenome mode.

### Input Files

- **Assembly files:** `/project/bcampb7/camplab/0_Raw_Data/Metagenome_estury/Assembly_All/MGs/<sample>/final.contigs.fa`
- **Sample list:** `mg_assembly_list.txt` (one sample name per line)

### Output Files

- `01_catalog/prodigal_MGs/<sample>/<sample>.fna` — Gene sequences (nucleotide)
- `01_catalog/prodigal_MGs/<sample>/<sample>.faa` — Protein sequences (amino acid)
- `01_catalog/prodigal_MGs/<sample>/<sample>.gff` — Gene annotations

### Usage

```bash
# Create sample list
ls /path/to/assemblies > mg_assembly_list.txt

# Edit array size in script: --array=1-N (N = number of samples)
sbatch scripts/01_prodigal_MGs.sh
```

### Key Parameters

- `--array=1-N` — Array job for multiple samples
- `-p meta` — Metagenome mode
- `--cpus-per-task=4` — 4 CPUs per sample
- `--mem=16G` — 16GB RAM per sample

---

## Script 2: `02_cdhit_catalog90.sh`

### Purpose

Cluster all predicted genes at 90% sequence identity to create a non-redundant gene catalog.

### Input Files

- `01_catalog/all_MG_genes.fna` — Concatenated genes from all samples (67,609,507 sequences)

### Output Files

- `01_catalog/MG_gene_catalog90.fna` — Non-redundant gene catalog (40,619,901 sequences)
- `01_catalog/MG_gene_catalog90.fna.clstr` — Cluster information

### Usage

```bash
# First concatenate all genes
cat 01_catalog/prodigal_MGs/*/*.fna > 01_catalog/all_MG_genes.fna

# Submit clustering job
sbatch scripts/02_cdhit_catalog90.sh
```

### Key Parameters

- `-c 0.90` — 90% sequence identity threshold
- `-aS 0.90` — 90% alignment coverage
- `-G 0` — Global sequence identity
- `--cpus-per-task=32` — 32 CPUs
- `--mem=200G` — 200GB RAM (critical for large datasets)
- `--time=48:00:00` — 48-hour time limit

---

## Script 3: `03_prodigal_catalog.sh`

### Purpose

Predict protein sequences from the non-redundant gene catalog.

### Input Files

- `01_catalog/MG_gene_catalog90.fna` — Non-redundant gene catalog

### Output Files

- `01_catalog/MG_gene_catalog90.faa` — Protein sequences (for annotation)
- `01_catalog/MG_gene_catalog90.genes.fna` — CDS nucleotide sequences
- `01_catalog/MG_gene_catalog90.gff` — Gene annotations

### Usage

````bash
sbatch scripts/03_prodigal_catalog.sh

### Key Parameters
- `-p meta` — Metagenome mode
- `--cpus-per-task=4` — 4 CPUs
- `--mem=32G` — 32GB RAM
- `--time=6:00:00` — 6-hour time limit

---

## Script 4: `04_eggnog_catalog.sh`

### Purpose
Functionally annotate proteins using eggNOG-mapper database.

### Input Files
- `01_catalog/MG_gene_catalog90.faa` — Protein sequences (40,619,901 proteins)
- **Database:** `$EGGNOG_DATA_DIR` (must be downloaded beforehand)

### Output Files
- `05_eggnog/MG_catalog.emapper.annotations` — Main annotation table (KEGG, COG, GO terms)
- `05_eggnog/MG_catalog.emapper.seed_orthologs` — Ortholog assignments
- `05_eggnog/MG_catalog.emapper.hits` — DIAMOND search results

### Usage
```bash
# One-time database setup
export EGGNOG_DATA_DIR=/project/bcampb7/camplab/AL_JJ_oct23/databases/eggnog
download_eggnog_data.py --data_dir ${EGGNOG_DATA_DIR}

# Submit annotation job
sbatch scripts/04_eggnog_catalog.sh
````

### Key Parameters

- `--itype proteins` — Input is protein sequences
- `--cpu 32` — 32 CPUs
- `--mem=100G` — 100GB RAM
- `--time=48:00:00` — 48-hour time limit (may need more for 40M proteins)

---

## Workflow Execution Order

```
1. 01_prodigal_MGs.sh     → Predict genes from assemblies
2. [Manual concatenation]      → Combine all genes
3. 02_cdhit_catalog90.sh   → Cluster at 90% identity
4. 03_prodigal_catalog.sh → Predict proteins on catalog
5. 04_eggnog_catalog.sh   → Annotate proteins
```

---

## File Locations Summary

```
Input:
  /project/bcampb7/camplab/0_Raw_Data/Metagenome_estury/Assembly_All/MGs/

Intermediate:
  01_catalog/prodigal_MGs/          # Gene predictions per sample
  01_catalog/all_MG_genes.fna       # Concatenated genes (67M)

Output:
  01_catalog/MG_gene_catalog90.fna  # Non-redundant catalog (40M)
  01_catalog/MG_gene_catalog90.faa  # Protein sequences
  05_eggnog/MG_catalog.emapper.annotations  # Functional annotations

Logs:
  00_logs/                          # SLURM job logs
```

---

**Last updated:** January 29, 2025
