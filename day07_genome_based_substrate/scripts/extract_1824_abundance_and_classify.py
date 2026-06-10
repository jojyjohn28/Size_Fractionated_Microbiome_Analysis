#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os

BASE = "/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/struo2_jan20/humann3_genome_analysis"
OUTDIR = "/project/bcampb7/camplab/Jojy/Fl_vs_PA_2026_Jojy/CAZyme_TCDB_genomes/05_summary"

os.makedirs(OUTDIR, exist_ok=True)

abund_file = f"{BASE}/genome_abundance_table_approachA.tsv"
genome_list_file = f"{BASE}/annotated_1824_genomes.txt"

abund = pd.read_csv(abund_file, sep="\t")
abund = abund.rename(columns={abund.columns[0]: "Sample"})

genomes = pd.read_csv(genome_list_file, header=None)[0].tolist()

# Match accession before species name, e.g.
# RS_GCF_000016185.1|s__Something
matched_cols = ["Sample"]
matched_genomes = []

for g in genomes:
    hits = [c for c in abund.columns if c.startswith(g + "|")]
    if hits:
        matched_cols.extend(hits)
        matched_genomes.append(g)

abund_1824 = abund[matched_cols]

out_abund = f"{OUTDIR}/genome_abundance_1824.tsv"
abund_1824.to_csv(out_abund, sep="\t", index=False)

print(f"Input annotated genomes: {len(genomes)}")
print(f"Recovered genomes in abundance table: {len(matched_genomes)}")
print(f"Saved: {out_abund}")

# Long format classification
records = []

for col in matched_cols[1:]:
    genome = col.split("|")[0]
    species = col.split("|", 1)[1] if "|" in col else "NA"

    values = abund[["Sample", col]].copy()
    values.columns = ["Sample", "Abundance"]

    values["Fraction"] = values["Sample"].apply(
        lambda x: "PA" if "G08" in x else ("FL" if "L08" in x else "Unknown")
    )

    pa = values[values["Fraction"] == "PA"]["Abundance"]
    fl = values[values["Fraction"] == "FL"]["Abundance"]

    mean_pa = pa.mean()
    mean_fl = fl.mean()

    prev_pa = (pa > 0).sum()
    prev_fl = (fl > 0).sum()
    total_prev = prev_pa + prev_fl

    total_mean = mean_pa + mean_fl
    log2_pa_fl = np.log2((mean_pa + 1e-9) / (mean_fl + 1e-9))

    if total_prev < 3:
        lifestyle = "Low_confidence_rare"
    elif log2_pa_fl > 1:
        lifestyle = "PA_associated"
    elif log2_pa_fl < -1:
        lifestyle = "FL_associated"
    else:
        lifestyle = "Generalist"

    records.append({
        "Genome": genome,
        "Species": species,
        "mean_PA": mean_pa,
        "mean_FL": mean_fl,
        "log2_PA_FL": log2_pa_fl,
        "prevalence_PA": prev_pa,
        "prevalence_FL": prev_fl,
        "total_prevalence": total_prev,
        "total_mean_abundance": total_mean,
        "genome_lifestyle": lifestyle
    })

class_df = pd.DataFrame(records)
class_df = class_df.sort_values("total_mean_abundance", ascending=False)

out_class = f"{OUTDIR}/genome_FL_PA_classification_1824.tsv"
class_df.to_csv(out_class, sep="\t", index=False)

summary = class_df["genome_lifestyle"].value_counts().reset_index()
summary.columns = ["genome_lifestyle", "n_genomes"]
summary.to_csv(f"{OUTDIR}/genome_FL_PA_classification_1824_summary.tsv", sep="\t", index=False)

print("\nClassification summary:")
print(summary.to_string(index=False))
print(f"\nSaved: {out_class}")
