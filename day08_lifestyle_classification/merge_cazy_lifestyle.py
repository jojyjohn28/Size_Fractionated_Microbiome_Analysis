#!/usr/bin/env python3
"""
merge_cazy_lifestyle.py
=======================
Author:  Jojy John
Date:    June 2026
Purpose: Merge CAZyme counts (from dbCAN) with genome lifestyle assignments
         (PA_associated / FL_associated) produced by classify_genome_lifestyle.R

Input:
    cazy_genome.xlsx          -- CAZyme counts per genome (from dbCAN overview.txt)
    genome_metadata_lifestyle.csv  -- lifestyle reference (output of R script)

Output:
    CAZy_with_lifestyle.xlsx  -- merged table, matched genomes only
    CAZy_with_lifestyle_full.csv  -- merged table with all genomes (unmatched as NaN)
    merge_cazy_report.txt     -- summary of matching results

Usage:
    python merge_cazy_lifestyle.py
    python merge_cazy_lifestyle.py --cazy my_cazy.xlsx --lifestyle my_lifestyle.csv
"""

import argparse
import os
import pandas as pd

# ==============================================================================
# 0. Argument parsing (optional — defaults work for standard file layout)
# ==============================================================================
parser = argparse.ArgumentParser(description="Merge CAZyme counts with lifestyle assignments.")
parser.add_argument("--cazy",      default="input/cazy_genome.xlsx",
                    help="Path to CAZyme counts Excel file (default: input/cazy_genome.xlsx)")
parser.add_argument("--lifestyle", default="output/genome_metadata_lifestyle.csv",
                    help="Path to lifestyle reference CSV (default: output/genome_metadata_lifestyle.csv)")
parser.add_argument("--outdir",    default="output",
                    help="Output directory (default: output/)")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ==============================================================================
# 1. Load input files
# ==============================================================================
print(f"Loading CAZyme table: {args.cazy}")
cazy = pd.read_excel(args.cazy)

print(f"Loading lifestyle reference: {args.lifestyle}")
lifestyle = pd.read_csv(args.lifestyle)

print(f"\nCAZyme table: {cazy.shape[0]} rows × {cazy.shape[1]} columns")
print(f"Lifestyle ref: {lifestyle.shape[0]} rows × {lifestyle.shape[1]} columns")

# ==============================================================================
# 2. Clean column names and Genome IDs
# ==============================================================================
# Strip whitespace from all column names
cazy.columns      = cazy.columns.str.strip()
lifestyle.columns = lifestyle.columns.str.strip()

# Convert Genome IDs to string and strip whitespace
# (Excel sometimes reads accessions as floats — .astype(str) handles this)
cazy["Genome"]      = cazy["Genome"].astype(str).str.strip()
lifestyle["Genome"] = lifestyle["Genome"].astype(str).str.strip()

# Remove trailing ".0" caused by Excel float conversion (e.g. "GCA_12345.0" → "GCA_12345")
cazy["Genome"]      = cazy["Genome"].str.replace(r"\.0$", "", regex=True)
lifestyle["Genome"] = lifestyle["Genome"].str.replace(r"\.0$", "", regex=True)

print(f"\nCAZyme table genome ID preview:\n{cazy['Genome'].head(5).tolist()}")
print(f"Lifestyle ref genome ID preview:\n{lifestyle['Genome'].head(5).tolist()}")

# ==============================================================================
# 3. Diagnose mismatches before merging
# ==============================================================================
cazy_ids      = set(cazy["Genome"])
lifestyle_ids = set(lifestyle["Genome"])

in_both       = cazy_ids & lifestyle_ids
cazy_only     = cazy_ids - lifestyle_ids
lifestyle_only = lifestyle_ids - cazy_ids

print(f"\nID matching summary:")
print(f"  In both tables:         {len(in_both)}")
print(f"  CAZy only (unmatched):  {len(cazy_only)}")
print(f"  Lifestyle only:         {len(lifestyle_only)}")

if len(cazy_only) > 0:
    print(f"\nFirst 10 unmatched CAZy IDs:")
    print(list(cazy_only)[:10])
    print("\nCommon causes:")
    print("  1. Version suffix mismatch: GCA_000153445.1 vs GCA_000153445")
    print("  2. Trailing whitespace (should be fixed by .str.strip())")
    print("  3. Excel float: 1234567.0 instead of 1234567")

# ==============================================================================
# 4. Merge
# ==============================================================================
# Select only the most useful lifestyle columns for the merge
lifestyle_cols = [
    "Genome", "Genome_size_Mbp", "Lifestyle", "log2_PA_FL",
    "Mean_PA", "Mean_FL", "Prev_PA", "Prev_FL",
    "PA_Dominance_Percent", "FL_Dominance_Percent",
    "PA_Higher_Count", "FL_Higher_Count", "Total_Paired_Comparisons"
]

# Keep only columns that exist (graceful handling if some are missing)
lifestyle_cols = [c for c in lifestyle_cols if c in lifestyle.columns]

merged = pd.merge(
    cazy,
    lifestyle[lifestyle_cols],
    on="Genome",
    how="left",
    indicator=True
)

matched   = merged[merged["_merge"] == "both"].drop(columns="_merge")
full_left = merged.drop(columns="_merge")

print(f"\nMerge results:")
print(f"  Matched genomes:   {len(matched)}")
print(f"  Unmatched (left):  {(merged['_merge'] == 'left_only').sum()}")

# ==============================================================================
# 5. Save outputs
# ==============================================================================
out_xlsx = os.path.join(args.outdir, "CAZy_with_lifestyle.xlsx")
out_csv  = os.path.join(args.outdir, "CAZy_with_lifestyle_full.csv")
out_rep  = os.path.join(args.outdir, "merge_cazy_report.txt")

matched.to_excel(out_xlsx, index=False)
full_left.to_csv(out_csv, index=False)

print(f"\nOutputs saved:")
print(f"  {out_xlsx}  (matched genomes only)")
print(f"  {out_csv}   (all genomes, unmatched lifestyle = NaN)")

# ==============================================================================
# 6. Summary statistics by lifestyle
# ==============================================================================
if "Lifestyle" in matched.columns:
    # Identify CAZyme family columns (GH, GT, PL, CE, AA, CBM)
    cazyme_cols = [c for c in matched.columns
                   if c.startswith(("GH", "GT", "PL", "CE", "AA", "CBM"))
                   or c in ("Total_CAZymes",)]

    if cazyme_cols:
        print("\nCAZyme counts by lifestyle (mean):")
        print(matched.groupby("Lifestyle")[cazyme_cols].mean().round(2).to_string())
    else:
        print("\nNo CAZyme family columns detected for summary.")
        print("Expected columns starting with GH, GT, PL, CE, AA, or CBM.")

# ==============================================================================
# 7. Write report
# ==============================================================================
with open(out_rep, "w") as f:
    f.write("CAZy-Lifestyle Merge Report\n")
    f.write("=" * 40 + "\n")
    f.write(f"CAZyme input:      {args.cazy}\n")
    f.write(f"Lifestyle input:   {args.lifestyle}\n")
    f.write(f"Genomes in CAZy:   {len(cazy_ids)}\n")
    f.write(f"Genomes in ref:    {len(lifestyle_ids)}\n")
    f.write(f"Matched:           {len(in_both)}\n")
    f.write(f"CAZy only:         {len(cazy_only)}\n")
    if len(cazy_only) > 0:
        f.write("Unmatched CAZy IDs (first 20):\n")
        for g in list(cazy_only)[:20]:
            f.write(f"  {g}\n")

print(f"\nReport saved: {out_rep}")
print("\n✓ merge_cazy_lifestyle.py complete.")
