# Day 4: Gene Catalog Construction & Functional Annotation

## 📋 Overview

This directory contains scripts and documentation for building a non-redundant prokaryotic gene catalog from metagenome assemblies and annotating it with functional databases.

**Objective:** Create a unified, non-redundant reference gene catalog from all metagenome assemblies to enable:

- Functional profiling across samples
- Direct comparison of DNA (potential) vs RNA (expression)
- Substrate uptake transporter analysis (ABC, TBDT, CAZymes)
- Functional redundancy calculations

## 🎯 Key Outputs

1. **Non-redundant gene catalog** (90% identity clustering)
2. **Protein sequences** for functional annotation
3. **eggNOG functional annotations** (KEGG, COG, GO terms)
4. **Substrate-specific gene sets** (ABC transporters, TBDTs)

## 📊 Statistics

| Metric                           | Value            |
| -------------------------------- | ---------------- |
| Original genes (all assemblies)  | 67,609,507       |
| Non-redundant genes (90% CD-HIT) | 40,619,901       |
| Redundancy removed               | ~40% (27M genes) |
| Protein sequences for annotation | 40,619,901       |

## 📁 Directory Structure

```
day04_gene_catalog/
├── README.md                           # This file
├── NOTES.md                            # Analysis notes and troubleshooting
├── scripts/
   ├── 01_prodigal_MGs.sh         # Gene prediction from assemblies
   ├── 02_cdhit_catalog90.sh      # CD-HIT clustering at 90%
   ├── 03_prodigal_catalog.sh     # Protein prediction on catalog
   └── 04_eggnog_catalog.sh      # Functional annotation

```

## 🚀 Quick Start

### Prerequisites

```bash
# Required modules/tools
module load prodigal/2.6.3
module load cd-hit/4.8.1
# eggnog-mapper (conda environment or module)
```

### Step-by-Step Execution

#### 1. Gene Prediction from Assemblies

```bash
# Create sample list
ls /path/to/assemblies > mg_assembly_list.txt

# Edit array size in script (--array=1-N)
N=$(wc -l < mg_assembly_list.txt)

# Submit job
sbatch scripts/01_prodigal_MGs.sh
```

**Output:** `01_catalog/prodigal_MGs/<sample>/<sample>.fna/.faa/.gff`

#### 2. Concatenate All Genes

```bash
cd /project/bcampb7/camplab/AL_JJ_oct23/substrate_FRed/01_catalog
cat prodigal_MGs/*/*.fna > all_MG_genes.fna

# Verify
grep -c "^>" all_MG_genes.fna
```

#### 3. CD-HIT Clustering

```bash
sbatch scripts/02_cdhit_catalog90.sh
```

**Output:** `MG_gene_catalog90.fna` (non-redundant catalog)

**Parameters:**

- `-c 0.90` — 90% sequence identity
- `-aS 0.90` — 90% alignment coverage
- `-G 0` — Global alignment

#### 4. Protein Prediction on Catalog

```bash
sbatch scripts/03_prodigal_catalog.sh
```

**Output:** `MG_gene_catalog90.faa` (protein sequences)

#### 5. Functional Annotation with eggNOG

```bash
# One-time database setup
export EGGNOG_DATA_DIR=/path/to/databases/eggnog
download_eggnog_data.py --data_dir ${EGGNOG_DATA_DIR}

# Run annotation
sbatch scripts/04_eggnog_catalog.sh
```

**Output:** `05_eggnog/MG_catalog.emapper.annotations`

**Expected runtime:** 24-48+ hours for ~40M proteins

## 📈 Next Steps (Day 5)

1. Build Bowtie2 index from gene catalog
2. Map MG reads → catalog (DNA abundance)
3. Map MT reads → catalog (RNA abundance)
4. Calculate gene abundance (TPM/RPKM)
5. Extract substrate-specific transporters (ABC, TBDT)
6. Compute functional redundancy (FRed)

## 🔗 Dependencies

| Tool          | Version | Purpose                 |
| ------------- | ------- | ----------------------- |
| Prodigal      | 2.6.3   | Gene/protein prediction |
| CD-HIT        | 4.8.1   | Sequence clustering     |
| eggNOG-mapper | 2.1+    | Functional annotation   |
| Python        | 3.7+    | Scripting               |
| SLURM         | -       | Job scheduling          |

## 📚 References

- **Prodigal:** Hyatt et al. (2010) BMC Bioinformatics
- **CD-HIT:** Fu et al. (2012) Bioinformatics
- **eggNOG:** Huerta-Cepas et al. (2019) Nucleic Acids Research

## 💡 Tips

- **Memory requirements:** CD-HIT needs substantial RAM (200GB recommended)
- **Time estimates:**
  - Prodigal per sample: 10-30 min
  - CD-HIT: 5+ days depending on your file size
  - eggNOG: 24-48+ hours
- **Checkpoint files:** Keep intermediate files for troubleshooting
- **Verify outputs:** Always check sequence counts after each step

## 📧 Contact

For questions about this workflow:

- **Blog:** https://jojyjohn28.github.io/blog/size-fractionated-microbiome-analysis-day4/

---

**Last updated:** January 29, 2025
