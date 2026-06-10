# Day 7\_ FL vs PA Genome-Based CAZyme and Transporter Analysis

## Overview

This repository contains all scripts for annotating **carbohydrate-active enzymes (CAZymes)** and **membrane transporters (TCDB)** across **1,824 representative estuarine genomes**, and linking those functional profiles to observed free-living (FL) vs particle-associated (PA) abundance patterns in Chesapeake Bay and Delaware Bay metagenomes.

The core biological question: do FL and PA genomes carry different functional toolkits — and if so, which enzymes and transporters define each lifestyle?

**Blog post:** [blog](https://jojyjohn.github.io/blog/fl-pa-cazyme-transporter-genome-analysis/)

---

## Biological context

| Lifestyle                    | Ecological niche                                               | Expected functional signature                                        |
| ---------------------------- | -------------------------------------------------------------- | -------------------------------------------------------------------- |
| **Particle-associated (PA)** | Attached to organic particles; access to concentrated polymers | Large genomes · Many GH/PL CAZymes · TonB-dependent transporters     |
| **Free-living (FL)**         | Suspended in water column; dilute dissolved organic matter     | Streamlined genomes · Fewer CAZymes · High-affinity ABC transporters |
| **Generalist**               | Detected in both fractions                                     | Intermediate profiles                                                |

---

## Dataset

- **1,824 genomes** — representative MAGs selected from Chesapeake Bay and Delaware Bay estuarine metagenomes
- **30 samples** — paired PA (G08) and FL (L08) metagenomes across two bays and multiple seasons
- Genome classification based on log₂(mean PA abundance / mean FL abundance):
  - PA-associated: log₂FC > 1 → **337 genomes**
  - FL-associated: log₂FC < -1 → **254 genomes**
  - Generalist: −1 ≤ log₂FC ≤ 1 → **1,233 genomes**

---

## Repository structure

```
FL_vs_PA_genome_analysis/
├── scripts/
│   ├── 00_make_genome_list.sh               # Build genome path list
│   ├── 01_predict_proteins_prodigal.sh      # SLURM array: Prodigal gene prediction
│   ├── 02_run_dbcan_genomes.sh              # SLURM array: dbCAN CAZyme annotation
│   ├── 03_run_tcdb_genomes.sh               # SLURM array: DIAMOND TCDB search
│   ├── run_dbcan_hmm_genome.sh              # Single-genome dbCAN wrapper (for testing)
│   ├── calculate_genome_sizes.py            # Genome size from FASTA
│   ├── concatenate_dbcan_tcdb.py            # Merge per-genome annotation outputs
│   └── extract_1824_abundance_and_classify.py  # Abundance extraction and FL/PA/Gen classification
├── README.md
└── environment/
    └── environment_notes.md                 # Software versions and conda environments
```

> **Note:** Raw data (genome FASTA files, abundance tables) are stored on the Palmetto HPC cluster and are not included in this repository. See **Data availability** below.

---

## Pipeline overview

```
1,824 genome FASTA files
        │
        ▼
[00] Build genome list
        │
        ▼
[01] Prodigal — protein prediction
        │ .faa per genome
        ├────────────────────────────┐
        ▼                            ▼
[02] dbCAN                      [03] DIAMOND
     CAZyme annotation               TCDB transporter search
     (HMM-based)                     (sequence similarity)
        │                            │
        ▼                            ▼
 all_dbcan_overview_with_genome.tsv  all_tcdb_with_genome.tsv
        │                            │
        └──────────────┬─────────────┘
                       ▼
        [Abundance extraction + FL/PA classification]
                       │
                       ▼
        genome_integrated_FL_PA_analysis.tsv
                       │
                       ▼
        [Statistical comparisons in R]
        Kruskal-Wallis · Wilcoxon · Visualization
```

---

## Scripts

### `00_make_genome_list.sh`

Generates a sorted list of all genome FASTA paths. Used as the input manifest for all SLURM array jobs.

```bash
bash scripts/00_make_genome_list.sh
```

Output: `00_lists/genomes_fasta.list`

---

### `01_predict_proteins_prodigal.sh`

SLURM array job — runs Prodigal gene prediction on each genome independently.

```bash
sbatch scripts/01_predict_proteins_prodigal.sh
```

- Mode: `-p meta` (metagenome mode — appropriate for MAGs of any size)
- Array: `1-1824%50` (50 simultaneous jobs)
- Output per genome: `.faa` (proteins), `.ffn` (nucleotide CDS), `.gff` (coordinates)

---

### `02_run_dbcan_genomes.sh`

SLURM array job — runs dbCAN HMM-based CAZyme annotation on each genome's predicted proteins.

```bash
sbatch scripts/02_run_dbcan_genomes.sh
```

- Tool: HMMER only (`--tools hmmer`)
- Database: dbCAN-HMMdb-V14
- Array: `1-1824%20` (20 simultaneous jobs; memory-intensive)
- Key output: `overview.txt` per genome

---

### `03_run_tcdb_genomes.sh`

SLURM array job — runs DIAMOND blastp against the TCDB protein database.

```bash
sbatch scripts/03_run_tcdb_genomes.sh
```

- E-value cutoff: `1e-5`
- Coverage thresholds: `--query-cover 50 --subject-cover 50`
- Mode: `--sensitive`
- Array: `1-1824%30`
- Output: `{genome_id}_tcdb.tsv` per genome

---

### `run_dbcan_hmm_genome.sh`

Single-genome dbCAN wrapper — useful for testing parameters or re-running individual failed jobs without submitting the full array.

```bash
bash scripts/run_dbcan_hmm_genome.sh <genome_id>
```

---

### `calculate_genome_sizes.py`

Calculates genome size (total bp, Mbp) and contig count for all genomes from the FASTA files.

```bash
python scripts/calculate_genome_sizes.py
```

Output: `genome_sizes.tsv`

---

### `concatenate_dbcan_tcdb.py`

Merges per-genome annotation outputs into single combined tables.

```bash
python scripts/concatenate_dbcan_tcdb.py
```

Outputs:

- `all_dbcan_overview_with_genome.tsv` — all CAZyme annotations with genome ID
- `all_tcdb_with_genome.tsv` — all TCDB DIAMOND hits with genome ID

---

### `extract_1824_abundance_and_classify.py`

Extracts the 1,824 annotated genomes from the full 25,838-genome abundance table, calculates mean PA and FL abundances, computes log₂(PA/FL), and assigns lifestyle classifications.

```bash
python scripts/extract_1824_abundance_and_classify.py
```

Output: `genome_lifestyle_classification.tsv`

Classification thresholds:

- **PA-associated:** log₂(PA/FL) > 1
- **FL-associated:** log₂(PA/FL) < −1
- **Generalist:** −1 ≤ log₂(PA/FL) ≤ 1

---

## Software and versions

| Tool              | Version  | Purpose                |
| ----------------- | -------- | ---------------------- |
| Prodigal          | 2.6.3    | Gene prediction        |
| dbCAN / run_dbcan | 4.x      | CAZyme annotation      |
| dbCAN-HMMdb       | V14      | CAZyme HMM database    |
| DIAMOND           | 2.1.x    | TCDB similarity search |
| TCDB database     | May 2026 | Transporter reference  |
| Python            | 3.9+     | Data processing        |
| R                 | 4.3+     | Statistical analysis   |
| pandas            | —        | Python data handling   |
| BioPython         | —        | FASTA parsing          |

Full conda environment files will be added to `environment/` as the analysis is finalized.

---

## Planned downstream analyses

- [ ] Integrated table: genome size + CAZyme counts + transporter counts + lifestyle
- [ ] Kruskal-Wallis tests across PA / FL / Generalist groups
- [ ] Pairwise Wilcoxon tests with BH correction
- [ ] CAZyme family composition comparison (GH, GT, PL, CE, AA, CBM)
- [ ] TBDT vs ABC transporter enrichment analysis
- [ ] Visualization: boxplots, heatmaps, PCA by lifestyle
- [ ] Specific CAZyme family targets (e.g., alginate lyases, cellulose-active GHs)

---

## Data availability

Raw genome FASTA files and metagenome abundance tables are stored on the Palmetto HPC cluster:

```
Genomes:   /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/struo2_jan20/humann3_genome_analysis/genomes_top2000/genomes_fasta
Databases: /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/databases/
Analysis:  /project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes
```

Data will be deposited to a public repository (NCBI/Zenodo) upon manuscript submission.

---

## Citation

If you use these scripts, please cite this repository:

> John J. (2026). FL vs PA Genome-Based CAZyme and Transporter Analysis. GitHub. https://github.com/JojyJohn/Fl_vs_PA_genome_analysis

---

## Contact

**Jojy John**  
Campbell Lab, Clemson University  
GitHub: [@JojyJohn](https://github.com/JojyJohn)  
Blog: [jojyjohn.github.io](https://jojyjohn.github.io)

_Analysis ongoing — results and scripts will be updated as the project progresses._
