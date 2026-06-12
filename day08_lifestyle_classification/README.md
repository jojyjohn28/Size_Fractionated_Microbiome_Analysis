# FL vs PA Genome Lifestyle Classification

**Blog post:** [fl-pa-genome-lifestyle-classification](https://jojyjohn.github.io/blog/fl-pa-genome-lifestyle-classification/)

---

## Overview

This repository contains the R and Python scripts for classifying ~2,000 representative estuarine genomes as **free-living (FL)** or **particle-associated (PA)** based on their abundance patterns across 36 paired size-fractionated metagenomes from Chesapeake Bay and Delaware Bay.

The lifestyle classification is used as a genome-level metadata layer for downstream comparisons of:

- Genome size
- CAZyme (carbohydrate-active enzyme) content
- ABC transporter content
- TonB-dependent transporter (TBDT) content
- Other functional traits

---

## Biological rationale

Genomes do not come with lifestyle labels. A genome binned from a PA library might be genuinely particle-adapted — or it might be a free-living organism captured in the PA fraction by chance. The correct approach is to look at how a genome's **abundance behaves across many paired samples**: does it consistently appear at higher levels in PA fractions, or FL fractions?

Three complementary metrics are used:

| Metric                    | What it measures                                                      |
| ------------------------- | --------------------------------------------------------------------- |
| Mean abundance comparison | Which fraction has higher average abundance?                          |
| Log₂ fold-change          | How large is the PA vs FL difference?                                 |
| Paired dominance count    | In how many matched PA-FL pairs does the genome favour each fraction? |

---

## Sample naming convention

```
CP_Spr01G08  →  G08 suffix = Particle-Associated (PA)
CP_Spr01L08  →  L08 suffix = Free-Living (FL)
CP_Spr01     →  base name = matched pair identifier
```

---

## Repository structure

```
FL_vs_PA_lifestyle_classification/
├── README.md
├── R/
│   └── classify_genome_lifestyle.R        # Full R workflow
├── Python/
│   ├── merge_cazy_lifestyle.py            # Merge CAZyme counts with lifestyle
│   └── merge_transporters_lifestyle.py    # Merge transporter counts with lifestyle
├── input/                                 # Input files (not tracked — see Data availability)
│   ├── Abd_updated_new.xlsx               # Genome abundance table
│   ├── cazy_genome.xlsx                   # CAZyme counts per genome
│   ├── life_style_reference.xlsx          # Lifestyle reference (output of R script)
│   └── transporters.tsv                   # Transporter counts per genome
└── output/                                # Generated output files
    ├── Abd_updated_new_with_lifestyle.xlsx
    ├── Abd_updated_new_with_lifestyle.csv
    ├── genome_metadata_lifestyle.csv       # ← primary output for downstream analyses
    ├── CAZy_with_lifestyle.xlsx
    └── transporters_with_lifestyle.csv
```

---

## Input files

### 1. `Abd_updated_new.xlsx` — genome abundance table

Tab-delimited Excel file. Required columns:

| Column            | Description                                              |
| ----------------- | -------------------------------------------------------- |
| `Genome`          | Genome identifier (e.g., `GB_GCA_000153445.1`)           |
| `Genome_size_Mbp` | Assembled genome size in megabases                       |
| `Num_contigs`     | Number of contigs in the assembly                        |
| `*G08` columns    | Relative abundance in PA samples (one column per sample) |
| `*L08` columns    | Relative abundance in FL samples (one column per sample) |

Generated with CoverM in genome mode:

```bash
coverm genome \
  --bam-files *.bam \
  --genome-fasta-directory /path/to/genomes/ \
  --methods relative_abundance \
  --min-covered-fraction 0 \
  --output-file genome_relative_abundance.tsv
```

### 2. `cazy_genome.xlsx` — CAZyme counts per genome

Output of dbCAN annotation, summarized per genome. Required columns: `Genome`, one column per CAZyme family class (GH, GT, PL, CE, AA, CBM).

### 3. `transporters.tsv` — transporter counts per genome

Output of DIAMOND TCDB annotation, summarized per genome. Required columns: `Genome`, transporter family count columns (e.g., `ABC_gene_count`, `TBDT_count`, `MFS_count`).

---

## Running the workflow

### Step 1 — Classify genomes (R)

```r
Rscript R/classify_genome_lifestyle.R
```

**Output:** `genome_metadata_lifestyle.csv` — the primary output used by all downstream scripts.

### Step 2 — Merge CAZyme counts with lifestyle (Python)

```bash
python Python/merge_cazy_lifestyle.py
```

**Output:** `CAZy_with_lifestyle.xlsx`

### Step 3 — Merge transporter counts with lifestyle (Python)

```bash
python Python/merge_transporters_lifestyle.py
```

**Output:** `transporters_with_lifestyle.csv`

---

## Output columns

| Column                     | Description                                                          |
| -------------------------- | -------------------------------------------------------------------- |
| `Genome`                   | Genome identifier                                                    |
| `Genome_size_Mbp`          | Total assembled genome size (Mbp)                                    |
| `Num_contigs`              | Number of contigs                                                    |
| `Mean_PA`                  | Mean relative abundance across all PA samples                        |
| `Mean_FL`                  | Mean relative abundance across all FL samples                        |
| `Prev_PA`                  | Number of PA samples where genome was detected (abundance > 0)       |
| `Prev_FL`                  | Number of FL samples where genome was detected                       |
| `Total_Prev`               | Total detection count across both fractions                          |
| `log2_PA_FL`               | log₂(mean_PA / mean_FL) — positive = PA-biased, negative = FL-biased |
| `Lifestyle`                | Final classification: `PA_associated` or `FL_associated`             |
| `PA_Higher_Count`          | Matched pairs where PA abundance > FL abundance                      |
| `FL_Higher_Count`          | Matched pairs where FL abundance > PA abundance                      |
| `Equal_Count`              | Matched pairs where abundances were equal                            |
| `Total_Paired_Comparisons` | Total matched PA-FL pairs (17 in this dataset)                       |
| `PA_Dominance_Percent`     | % of matched pairs where genome was PA-dominant                      |
| `FL_Dominance_Percent`     | % of matched pairs where genome was FL-dominant                      |

---

## Classification results (this dataset)

| Lifestyle     | Genomes   |
| ------------- | --------- |
| PA_associated | 767       |
| FL_associated | 233       |
| **Total**     | **1,000** |

---

## Software

| Tool     | Version | Purpose                          |
| -------- | ------- | -------------------------------- |
| R        | ≥ 4.3   | Lifestyle classification         |
| readxl   | —       | Read Excel input                 |
| dplyr    | —       | Data manipulation                |
| openxlsx | —       | Write Excel output               |
| Python   | ≥ 3.9   | Merging functional annotations   |
| pandas   | —       | DataFrame operations and merging |
| openpyxl | —       | Read/write Excel in Python       |

Install R packages:

```r
install.packages(c("readxl", "dplyr", "openxlsx"))
```

Install Python packages:

```bash
pip install pandas openpyxl
```

---

## Common issues

**ID mismatch after merge**  
Run the diagnosis block in either Python script. Common causes: trailing whitespace (fixed by `.str.strip()`), version suffix mismatch (`GCA_000153445.1` vs `GCA_000153445`), or Excel converting accessions to floats.

**Unequal PA and FL column counts**  
Normal and expected. The mean calculation (`rowMeans`) handles unequal column counts correctly. Only matched pairs (same base name) are used for the dominance count analysis.

**`log2_PA_FL = NaN`**  
Both `Mean_PA` and `Mean_FL` were 0 for that genome. The pseudocount `1e-9` prevents this — check that the pseudocount line was not accidentally removed.

---

## Data availability

Raw abundance tables and genome FASTA files are stored on the Palmetto HPC cluster:

```
/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/
```

Data will be deposited to NCBI SRA upon manuscript submission.

---

## Part of the FL vs PA genome analysis series

| Day       | Analysis                                          | Repository                                                                                 |
| --------- | ------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| Day 6     | CAZyme and transporter annotation (1,824 genomes) | [FL_vs_PA_genome_analysis](https://github.com/JojyJohn/FL_vs_PA_genome_analysis)           |
| Day 7     | Co-occurrence network analysis                    | [cooccurrence_network_analysis](https://github.com/JojyJohn/cooccurrence_network_analysis) |
| **Day 8** | **Genome lifestyle classification**               | **this repo**                                                                              |

---
