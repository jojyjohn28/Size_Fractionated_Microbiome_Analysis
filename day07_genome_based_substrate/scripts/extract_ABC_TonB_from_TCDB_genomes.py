#!/usr/bin/env python3

import pandas as pd
import os

BASE = "/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes"

tcdb_file = f"{BASE}/05_summary/all_tcdb_with_genome.tsv"

outdir = f"{BASE}/05_summary/transporter_gene_lists"
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(tcdb_file, sep="\t")

# Make lowercase searchable text from TCDB title and subject ID
df["search_text"] = (
    df["sseqid"].astype(str) + " " +
    df["stitle"].astype(str)
).str.lower()

# ABC transporters
abc = df[
    df["search_text"].str.contains("abc", na=False) |
    df["search_text"].str.contains("atp-binding cassette", na=False)
].copy()

# TonB / TBDT transporters
tonb = df[
    df["search_text"].str.contains("tonb", na=False) |
    df["search_text"].str.contains("ton-b", na=False) |
    df["search_text"].str.contains("tonb-dependent", na=False) |
    df["search_text"].str.contains("tonb dependent", na=False) |
    df["search_text"].str.contains("tbdt", na=False) |
    df["search_text"].str.contains("outer membrane receptor", na=False)
].copy()

# Save full hit tables
abc.to_csv(f"{outdir}/ABC_tcdb_hits.tsv", sep="\t", index=False)
tonb.to_csv(f"{outdir}/TonB_TBDT_tcdb_hits.tsv", sep="\t", index=False)

# Save only protein/gene IDs
abc["qseqid"].drop_duplicates().to_csv(
    f"{outdir}/ABC_gene_ids.txt",
    index=False,
    header=False
)

tonb["qseqid"].drop_duplicates().to_csv(
    f"{outdir}/TonB_gene_ids.txt",
    index=False,
    header=False
)

# Save genome-level count summaries
abc_genome_counts = abc.groupby("Genome")["qseqid"].nunique().reset_index()
abc_genome_counts.columns = ["Genome", "ABC_gene_count"]

tonb_genome_counts = tonb.groupby("Genome")["qseqid"].nunique().reset_index()
tonb_genome_counts.columns = ["Genome", "TonB_TBDT_gene_count"]

abc_genome_counts.to_csv(
    f"{outdir}/ABC_genome_counts.tsv",
    sep="\t",
    index=False
)

tonb_genome_counts.to_csv(
    f"{outdir}/TonB_TBDT_genome_counts.tsv",
    sep="\t",
    index=False
)

print("Done.")
print(f"ABC hits: {len(abc)}")
print(f"ABC genes: {abc['qseqid'].nunique()}")
print(f"ABC genomes: {abc['Genome'].nunique()}")

print(f"TonB/TBDT hits: {len(tonb)}")
print(f"TonB/TBDT genes: {tonb['qseqid'].nunique()}")
print(f"TonB/TBDT genomes: {tonb['Genome'].nunique()}")

print(f"Output directory: {outdir}")
