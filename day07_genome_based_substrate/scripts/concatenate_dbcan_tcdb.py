#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

BASE = Path("/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes")

DBCAN_DIR = BASE / "03_dbcan"
TCDB_DIR  = BASE / "04_tcdb"
OUTDIR    = BASE / "05_summary"

OUTDIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# dbCAN: subfolder/overview.txt
# -----------------------------
dbcan_dfs = []

for file in sorted(DBCAN_DIR.glob("*/overview.txt")):
    genome = file.parent.name

    try:
        df = pd.read_csv(file, sep="\t")
        df.insert(0, "Genome", genome)
        dbcan_dfs.append(df)
    except Exception as e:
        print(f"Skipping dbCAN {file}: {e}")

if dbcan_dfs:
    dbcan_all = pd.concat(dbcan_dfs, ignore_index=True)
    dbcan_out = OUTDIR / "all_dbcan_overview_with_genome.tsv"
    dbcan_all.to_csv(dbcan_out, sep="\t", index=False)
    print(f"dbCAN saved: {dbcan_out}")
    print(f"dbCAN rows: {len(dbcan_all)}")
    print(f"dbCAN genomes: {dbcan_all['Genome'].nunique()}")
else:
    print("No dbCAN overview.txt files found")


# -----------------------------
# TCDB: genome_tcdb.tsv files
# -----------------------------
tcdb_cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
]

tcdb_dfs = []

for file in sorted(TCDB_DIR.glob("*_tcdb.tsv")):
    genome = file.name.replace("_tcdb.tsv", "")

    try:
        df = pd.read_csv(file, sep="\t", header=None, names=tcdb_cols)
        df.insert(0, "Genome", genome)
        tcdb_dfs.append(df)
    except Exception as e:
        print(f"Skipping TCDB {file}: {e}")

if tcdb_dfs:
    tcdb_all = pd.concat(tcdb_dfs, ignore_index=True)
    tcdb_out = OUTDIR / "all_tcdb_with_genome.tsv"
    tcdb_all.to_csv(tcdb_out, sep="\t", index=False)
    print(f"TCDB saved: {tcdb_out}")
    print(f"TCDB rows: {len(tcdb_all)}")
    print(f"TCDB genomes: {tcdb_all['Genome'].nunique()}")
else:
    print("No TCDB files found")
