#!/usr/bin/env python3
"""
merge_transporters_lifestyle.py
================================
Author:  Jojy John
Date:    June 2026
Purpose: Merge transporter counts (from DIAMOND TCDB annotation) with genome
         lifestyle assignments (PA_associated / FL_associated) produced by
         classify_genome_lifestyle.R

Input:
    transporters.tsv               -- transporter counts per genome (tab-delimited)
    genome_metadata_lifestyle.csv  -- lifestyle reference (output of R script)

Output:
    transporters_with_lifestyle.csv   -- merged table with all genomes
    transporters_matched.xlsx         -- matched genomes only (Excel)
    merge_transporters_report.txt     -- summary of matching and statistics

Expected transporter columns (flexible — script adapts to what is present):
    ABC_gene_count     -- ABC transporter genes
    TBDT_count         -- TonB-dependent transporters
    MFS_count          -- Major facilitator superfamily transporters
    Other_count        -- Other transporter families

Usage:
    python merge_transporters_lifestyle.py
    python merge_transporters_lifestyle.py --transporters my_transporters.tsv \\
                                            --lifestyle my_lifestyle.csv
"""

import argparse
import os
import pandas as pd

# ==============================================================================
# 0. Argument parsing
# ==============================================================================
parser = argparse.ArgumentParser(
    description="Merge transporter counts with lifestyle assignments."
)
parser.add_argument(
    "--transporters", default="input/transporters.tsv",
    help="Path to transporter counts TSV file (default: input/transporters.tsv)"
)
parser.add_argument(
    "--lifestyle", default="output/genome_metadata_lifestyle.csv",
    help="Path to lifestyle reference CSV (default: output/genome_metadata_lifestyle.csv)"
)
parser.add_argument(
    "--outdir", default="output",
    help="Output directory (default: output/)"
)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ==============================================================================
# 1. Load input files
# ==============================================================================
print(f"Loading transporter table: {args.transporters}")

# Auto-detect separator from file extension
if args.transporters.endswith(".csv"):
    transporters = pd.read_csv(args.transporters)
else:
    transporters = pd.read_csv(args.transporters, sep="\t")

print(f"Loading lifestyle reference: {args.lifestyle}")
lifestyle = pd.read_csv(args.lifestyle)

print(f"\nTransporter table: {transporters.shape[0]} rows × {transporters.shape[1]} columns")
print(f"Lifestyle ref:     {lifestyle.shape[0]} rows × {lifestyle.shape[1]} columns")
print(f"\nTransporter columns: {transporters.columns.tolist()}")

# ==============================================================================
# 2. Clean column names and Genome IDs
# ==============================================================================
transporters.columns = transporters.columns.str.strip()
lifestyle.columns    = lifestyle.columns.str.strip()

transporters["Genome"] = transporters["Genome"].astype(str).str.strip()
lifestyle["Genome"]    = lifestyle["Genome"].astype(str).str.strip()

# Remove trailing ".0" from Excel float conversion
transporters["Genome"] = transporters["Genome"].str.replace(r"\.0$", "", regex=True)
lifestyle["Genome"]    = lifestyle["Genome"].str.replace(r"\.0$", "", regex=True)

print(f"\nTransporter genome ID preview:\n{transporters['Genome'].head(5).tolist()}")
print(f"Lifestyle ref genome ID preview:\n{lifestyle['Genome'].head(5).tolist()}")

# ==============================================================================
# 3. Diagnose mismatches
# ==============================================================================
transport_ids  = set(transporters["Genome"])
lifestyle_ids  = set(lifestyle["Genome"])

in_both         = transport_ids & lifestyle_ids
transport_only  = transport_ids - lifestyle_ids
lifestyle_only  = lifestyle_ids - transport_ids

print(f"\nID matching summary:")
print(f"  In both tables:             {len(in_both)}")
print(f"  Transporter only:           {len(transport_only)}")
print(f"  Lifestyle only:             {len(lifestyle_only)}")

if len(transport_only) > 0:
    print(f"\nFirst 10 unmatched transporter IDs:")
    print(list(transport_only)[:10])

# ==============================================================================
# 4. Merge
# ==============================================================================
lifestyle_cols = [
    "Genome", "Genome_size_Mbp", "Lifestyle", "log2_PA_FL",
    "Mean_PA", "Mean_FL", "Prev_PA", "Prev_FL",
    "PA_Dominance_Percent", "FL_Dominance_Percent",
    "PA_Higher_Count", "FL_Higher_Count", "Total_Paired_Comparisons"
]
lifestyle_cols = [c for c in lifestyle_cols if c in lifestyle.columns]

merged_full = pd.merge(
    transporters,
    lifestyle[lifestyle_cols],
    on="Genome",
    how="left",
    indicator=True
)

matched = merged_full[merged_full["_merge"] == "both"].drop(columns="_merge")
full    = merged_full.drop(columns="_merge")

print(f"\nMerge results:")
print(f"  Matched genomes:  {len(matched)}")
print(f"  Unmatched:        {(merged_full['_merge'] == 'left_only').sum()}")

# ==============================================================================
# 5. Save outputs
# ==============================================================================
out_csv  = os.path.join(args.outdir, "transporters_with_lifestyle.csv")
out_xlsx = os.path.join(args.outdir, "transporters_matched.xlsx")
out_rep  = os.path.join(args.outdir, "merge_transporters_report.txt")

full.to_csv(out_csv, index=False)
matched.to_excel(out_xlsx, index=False)

print(f"\nOutputs saved:")
print(f"  {out_csv}   (all genomes)")
print(f"  {out_xlsx}  (matched only)")

# ==============================================================================
# 6. Summary statistics by lifestyle
# ==============================================================================
# Detect transporter count columns dynamically
known_transporter_cols = [
    "ABC_gene_count", "TBDT_count", "MFS_count",
    "Other_count", "Total_transporters"
]
transport_cols = [c for c in known_transporter_cols if c in matched.columns]

# Also catch any column ending in "_count" that isn't Genome-related
transport_cols += [
    c for c in matched.columns
    if c.endswith("_count") and c not in known_transporter_cols
    and c != "Genome"
]
transport_cols = list(dict.fromkeys(transport_cols))  # deduplicate, preserve order

if "Lifestyle" in matched.columns and transport_cols:
    print(f"\nTransporter counts by lifestyle (mean ± std):")
    summary = matched.groupby("Lifestyle")[transport_cols].agg(["mean", "std"]).round(2)
    print(summary.to_string())

    # Also compute genome-size-normalized counts if Genome_size_Mbp is available
    if "Genome_size_Mbp" in matched.columns:
        print(f"\nTransporter counts per Mbp (genome-size normalized, mean):")
        norm_df = matched.copy()
        for col in transport_cols:
            norm_df[f"{col}_per_Mbp"] = norm_df[col] / norm_df["Genome_size_Mbp"]
        norm_cols = [f"{c}_per_Mbp" for c in transport_cols]
        print(norm_df.groupby("Lifestyle")[norm_cols].mean().round(3).to_string())

elif not transport_cols:
    print("\nNo transporter count columns detected for summary.")
    print(f"Expected columns: {known_transporter_cols}")
    print(f"Actual columns:   {transporters.columns.tolist()}")

# ==============================================================================
# 7. Wilcoxon test (PA vs FL) for each transporter family
# ==============================================================================
if "Lifestyle" in matched.columns and transport_cols:
    try:
        from scipy import stats

        print("\nWilcoxon rank-sum tests: PA_associated vs FL_associated")
        print(f"{'Column':<30} {'PA mean':>10} {'FL mean':>10} {'p-value':>12} {'sig':>6}")
        print("-" * 72)

        pa_df = matched[matched["Lifestyle"] == "PA_associated"]
        fl_df = matched[matched["Lifestyle"] == "FL_associated"]

        for col in transport_cols:
            pa_vals = pa_df[col].dropna()
            fl_vals = fl_df[col].dropna()

            if len(pa_vals) < 3 or len(fl_vals) < 3:
                continue

            stat, p = stats.ranksums(pa_vals, fl_vals)

            sig = ""
            if p <= 0.001:   sig = "***"
            elif p <= 0.01:  sig = "**"
            elif p <= 0.05:  sig = "*"
            elif p <= 0.10:  sig = "."

            print(
                f"{col:<30} {pa_vals.mean():>10.4f} {fl_vals.mean():>10.4f} "
                f"{p:>12.4e} {sig:>6}"
            )

    except ImportError:
        print("\nInstall scipy for Wilcoxon tests: pip install scipy")

# ==============================================================================
# 8. Write report
# ==============================================================================
with open(out_rep, "w") as f:
    f.write("Transporter-Lifestyle Merge Report\n")
    f.write("=" * 40 + "\n")
    f.write(f"Transporter input: {args.transporters}\n")
    f.write(f"Lifestyle input:   {args.lifestyle}\n")
    f.write(f"Genomes in transporter table: {len(transport_ids)}\n")
    f.write(f"Genomes in lifestyle ref:     {len(lifestyle_ids)}\n")
    f.write(f"Matched:                      {len(in_both)}\n")
    f.write(f"Unmatched:                    {len(transport_only)}\n")
    if len(transport_only) > 0:
        f.write("\nUnmatched transporter IDs (first 20):\n")
        for g in list(transport_only)[:20]:
            f.write(f"  {g}\n")
    f.write("\nTransporter columns found:\n")
    for col in transport_cols:
        f.write(f"  {col}\n")

print(f"\nReport saved: {out_rep}")
print("\n✓ merge_transporters_lifestyle.py complete.")
