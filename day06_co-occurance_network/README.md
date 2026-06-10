# Co-occurrence Network Analysis — FL vs PA Estuarine Metagenomes

**Blog post:** [jojyjohn.github.io/blog/co-occurrence-network-analysis-genome-resolved](https://jojyjohn.github.io/blog/co-occurrence-network-analysis-genome-resolved/)

---

## Overview

This repository contains the complete R workflow for building genome-resolved co-occurrence networks from estuarine metagenome relative abundance data. Networks are constructed separately for four ecological groups — Chesapeake Bay free-living (CP_FL), Chesapeake Bay particle-associated (CP_PA), Delaware Bay free-living (DE_FL), and Delaware Bay particle-associated (DE_PA) — using Spearman correlation with centered log-ratio (CLR) normalization.

**Core outputs:**

- Publication-ready network figures (phylum-coloured, Zi-Pi, module-salinity)
- Network topology statistics per group
- Keystone taxa identified by Zi-Pi role classification
- Module ecology tests against salinity and size fraction

---

## Biological questions

| Question                                                   | Method                                                  |
| ---------------------------------------------------------- | ------------------------------------------------------- |
| Do FL and PA communities show different network structure? | Network statistics (modularity, density, pos:neg ratio) |
| Which genomes are keystone taxa?                           | Zi-Pi role classification                               |
| Which network modules respond to salinity?                 | Spearman ρ per module                                   |
| Do module abundances differ between FL and PA?             | Wilcoxon test per module                                |

---

## Dataset

- **Genomes:** ~1,000 representative MAGs from GTDB r220
- **Samples:** 36 estuarine metagenomes — Chesapeake Bay and Delaware Bay
- **Groups:** CP_FL · CP_PA · DE_FL · DE_PA
- **Correlation method:** Spearman + CLR normalization (via NetCoMi)
- **Edge threshold:** |ρ| ≥ 0.30

---

## Repository structure

```
cooccurrence_network_analysis/
├── cooccurrence_network_analysis.R   # ← main analysis script (this repo)
├── README.md
├── input/                            # input files (not tracked in git — see below)
│   ├── genome_relative_abundance.tsv
│   ├── metadata_updated.csv
│   └── gtdb_taxonomy.xlsx
├── results_figures/                  # generated figures (PDF + TIFF)
│   ├── Figure3A_network_phylum_panel.pdf
│   ├── Figure3B_ZiPi_role_plot.pdf
│   └── Figure3C_module_salinity_correlations.pdf
└── results_tables/                   # generated tables (CSV)
    ├── network_stats_final.csv
    ├── node_table_with_taxonomy.csv
    ├── top_hubs.csv
    ├── zipi_roles.csv
    ├── keystone_taxa.csv
    ├── module_membership.csv
    ├── module_phylum_summary.csv
    ├── module_salinity_stats.csv
    ├── module_fraction_stats.csv
    ├── Figure3B_ZiPi_role_summary.csv
    └── sessionInfo.txt
```

> Raw abundance tables and genome FASTA files are stored on the Palmetto HPC cluster and are not included in this repository. See **Data availability** below.

---

## Input files

### 1. Relative abundance table — `genome_relative_abundance.tsv`

Tab-delimited. Genomes as rows, samples as columns. Values are relative abundances (0–1, summing to 1 per sample).

Generated with CoverM in genome mode:

```bash
coverm genome \
  --bam-files *.bam \
  --genome-fasta-directory /path/to/genomes/ \
  --methods relative_abundance \
  --min-covered-fraction 0 \
  --output-file genome_relative_abundance.tsv
```

### 2. Metadata — `metadata_updated.csv`

Tab-delimited. One row per sample. Required columns:

| Column          | Description                                    | Values                             |
| --------------- | ---------------------------------------------- | ---------------------------------- |
| `sample`        | Sample ID (must match abundance table columns) | e.g. `CP_S1_FL`                    |
| `Bay`           | Estuary                                        | `Chesapeake`, `Delaware`           |
| `size_fraction` | Ecological compartment                         | `Free Living`, `Particle Attached` |
| `season`        | Season of sampling                             | `Spring`, `Summer`, `Fall`         |
| `Salinity`      | Salinity measurement                           | numeric                            |

### 3. Taxonomy — `gtdb_taxonomy.xlsx`

One row per genome. Required columns: `genome` (matching row names of abundance table), `gtdb_taxonomy` (full GTDB string: `d__Bacteria;p__...;g__Genus;s__Species`).

---

## Running the analysis

### Install dependencies

```r
install.packages(c("tidyverse", "igraph", "tidygraph", "ggraph",
                   "ggrepel", "readxl", "viridis", "patchwork"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("microeco")

install.packages("NetCoMi",
  repos = c("https://stefpeschel.r-universe.dev", "https://cloud.r-project.org"))
```

### Run

```r
source("cooccurrence_network_analysis.R")
```

All figures and tables are written to `results_figures/` and `results_tables/` automatically.

---

## Key parameters (adjust as needed)

| Parameter                       | Default              | Description                                            |
| ------------------------------- | -------------------- | ------------------------------------------------------ |
| Global prevalence filter        | ≥ 30% of all samples | Genomes present in fewer samples are excluded globally |
| Global mean abundance filter    | > 0.05%              | Ultra-rare genomes excluded                            |
| Within-group prevalence filter  | ≥ 40%                | Tighter — within a single group of ~9 samples          |
| Within-group abundance filter   | > 0.1%               | Applied after group subsetting                         |
| Correlation method              | Spearman + CLR       | `normMethod = "clr"` in NetCoMi                        |
| Edge threshold                  | \|ρ\| ≥ 0.30         | `thresh = 0.30` in NetCoMi                             |
| Zi threshold (module hub)       | Zi > 2.5             | Standard Olesen et al. 2007 cutoff                     |
| Pi threshold (connector)        | Pi > 0.62            | Standard Olesen et al. 2007 cutoff                     |
| Top edges shown in network plot | 300                  | `subset_top_edges(top_n = 300)`                        |

---

## Figures produced

### Figure 3A — Network visualisation (2×2 panel)

Four co-occurrence networks (CP_FL, CP_PA, DE_FL, DE_PA) plotted with Fruchterman-Reingold force-directed layout. Nodes coloured by phylum, sized by degree. Positive edges green, negative edges red. Top 15 hubs labelled by genus.

### Figure 3B — Zi-Pi role classification

Each genome plotted at its (Pi, Zi) coordinates, coloured by role. Dashed lines mark the Zi = 2.5 and Pi = 0.62 thresholds. Faceted by ecological group.

- **Peripheral** (green): Zi ≤ 2.5, Pi ≤ 0.62 — most taxa
- **Connector** (red): Zi ≤ 2.5, Pi > 0.62 — bridges between modules
- **Module hub** (blue): Zi > 2.5, Pi ≤ 0.62 — dominant within one module
- **Network hub** (purple): Zi > 2.5, Pi > 0.62 — globally influential (rare)

### Figure 3C — Module–salinity correlations

Bar chart of Spearman ρ between module abundance and salinity, per network module per group. Significance stars from BH-corrected p-values. Diverging blue–red colour scale.

---

## Output tables

| File                             | Description                                                                      |
| -------------------------------- | -------------------------------------------------------------------------------- |
| `network_stats_final.csv`        | Nodes, edges, density, modularity, path length, pos:neg ratio per group          |
| `node_table_with_taxonomy.csv`   | Degree, betweenness, closeness, eigenvector centrality + taxonomy for every node |
| `top_hubs.csv`                   | Top 20 nodes by degree per group                                                 |
| `zipi_roles.csv`                 | Full Zi-Pi table with role assignments for every node                            |
| `keystone_taxa.csv`              | Connectors, module hubs, and network hubs only                                   |
| `module_membership.csv`          | Module assignment per genome per group                                           |
| `module_phylum_summary.csv`      | Phylum composition of each module                                                |
| `module_salinity_stats.csv`      | Spearman ρ + BH-corrected p-value for module–salinity correlation                |
| `module_fraction_stats.csv`      | Wilcoxon FL vs PA test per module                                                |
| `Figure3B_ZiPi_role_summary.csv` | Count of nodes per role per group                                                |
| `sessionInfo.txt`                | R and package versions                                                           |

---

## Software versions

| Package                | Purpose                                                |
| ---------------------- | ------------------------------------------------------ |
| `NetCoMi`              | Network construction (Spearman + CLR, thresholding)    |
| `igraph`               | Graph operations, centrality metrics, module detection |
| `tidygraph` / `ggraph` | Tidy graph manipulation and ggplot2-style plotting     |
| `microeco`             | Microtable container for abundance + metadata          |
| `ggrepel`              | Non-overlapping text labels                            |
| `patchwork`            | Multi-panel figure assembly                            |
| `viridis`              | Perceptually uniform colour palettes                   |

Full version details in `results_tables/sessionInfo.txt` after running the script.

---

## Common issues

**`stopifnot()` fails at sample alignment**  
Run `trimws()` on both `rownames(abund)` and `meta$sample`. The most common cause is trailing whitespace from Excel export.

**Network too dense (density > 0.1)**  
Increase the edge threshold: change `thresh = 0.30` to `thresh = 0.40` or `0.50` in Section 8.

**No network hubs in Zi-Pi plot**  
Expected for sparse environmental networks. The Zi > 2.5 / Pi > 0.62 thresholds were proposed for gut microbiome data and may not be reached in environmental communities with fewer samples. Document this as an ecological observation.

**`netConstruct` memory error on large matrices**  
The within-group prevalence filter (≥40%) should reduce the matrix to a manageable size. If the group matrix still has >500 taxa, increase the filter to 0.50.

---
