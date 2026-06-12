################################################################################
# FL vs PA Genome Lifestyle Classification
# Author:  Jojy John
# Date:    June 2026
# Purpose: Classify representative estuarine genomes as free-living (FL) or
#          particle-associated (PA) based on abundance across paired
#          size-fractionated metagenomes.
#
# Input:   Abd_updated_new.xlsx  — genome abundance table
#            Columns: Genome, Genome_size_Mbp, Num_contigs, *G08 (PA), *L08 (FL)
#
# Outputs:
#   Abd_updated_new_with_lifestyle.xlsx  — full table with all new columns
#   Abd_updated_new_with_lifestyle.csv
#   genome_metadata_lifestyle.csv        — reduced metadata for downstream merging
################################################################################

# ==============================================================================
# 0. Libraries
# ==============================================================================
library(readxl)
library(dplyr)
library(openxlsx)

# ==============================================================================
# 1. Set paths
# ==============================================================================
INPUT_FILE    <- "input/Abd_updated_new.xlsx"
OUT_XLSX      <- "output/Abd_updated_new_with_lifestyle.xlsx"
OUT_CSV       <- "output/Abd_updated_new_with_lifestyle.csv"
OUT_META      <- "output/genome_metadata_lifestyle.csv"

dir.create("output", showWarnings = FALSE)

# ==============================================================================
# 2. Read abundance table
# ==============================================================================
df <- read_excel(INPUT_FILE)

cat("Loaded abundance table:", nrow(df), "genomes,", ncol(df), "columns\n")
cat("Column names preview:\n")
print(head(names(df), 10))

# ==============================================================================
# 3. Identify PA and FL sample columns
# ==============================================================================
# Convention: G08 suffix = Particle-Associated (PA)
#             L08 suffix = Free-Living (FL)

pa_cols <- grep("G08$", names(df), value = TRUE)
fl_cols <- grep("L08$", names(df), value = TRUE)

cat("\nPA columns (G08):", length(pa_cols), "\n")
cat("FL columns (L08):", length(fl_cols), "\n")

if (length(pa_cols) == 0 || length(fl_cols) == 0) {
  stop("No PA or FL columns found. Check column naming convention (G08/L08 suffixes).")
}

# ==============================================================================
# 4. Calculate mean abundance, prevalence, and log2 fold-change
# ==============================================================================
df2 <- df %>%
  mutate(
    # Mean relative abundance across all samples in each fraction
    Mean_PA    = rowMeans(select(., all_of(pa_cols)), na.rm = TRUE),
    Mean_FL    = rowMeans(select(., all_of(fl_cols)), na.rm = TRUE),

    # Prevalence: number of samples where genome was detected (abundance > 0)
    Prev_PA    = rowSums(select(., all_of(pa_cols)) > 0, na.rm = TRUE),
    Prev_FL    = rowSums(select(., all_of(fl_cols)) > 0, na.rm = TRUE),

    Total_Prev = Prev_PA + Prev_FL,

    # Log2 fold-change: positive = PA-biased, negative = FL-biased
    # Pseudocount 1e-9 prevents log(0)
    log2_PA_FL = log2((Mean_PA + 1e-9) / (Mean_FL + 1e-9))
  )

cat("\nLog2 fold-change summary:\n")
print(summary(df2$log2_PA_FL))

# ==============================================================================
# 5. Lifestyle classification
# ==============================================================================
df2 <- df2 %>%
  mutate(
    Lifestyle = case_when(
      Mean_PA > Mean_FL ~ "PA_associated",
      Mean_FL > Mean_PA ~ "FL_associated",
      TRUE              ~ "Equal"
    )
  )

cat("\nLifestyle classification:\n")
print(table(df2$Lifestyle))

# ==============================================================================
# 6. Paired dominance analysis
# ==============================================================================
# Match PA and FL columns by sample base name (strip G08/L08 suffix)

pa_base     <- sub("G08$", "", pa_cols)
fl_base     <- sub("L08$", "", fl_cols)
common_base <- intersect(pa_base, fl_base)

cat("\nMatched PA-FL pairs:", length(common_base), "\n")

# Report unmatched columns
pa_only <- setdiff(pa_base, fl_base)
fl_only <- setdiff(fl_base, pa_base)

if (length(pa_only) > 0) {
  cat("PA samples without matching FL sample:\n")
  print(pa_only)
}
if (length(fl_only) > 0) {
  cat("FL samples without matching PA sample:\n")
  print(fl_only)
}

# Reconstruct matched column names
matched_pa_cols <- paste0(common_base, "G08")
matched_fl_cols <- paste0(common_base, "L08")

# Verify columns exist in dataframe
stopifnot(all(matched_pa_cols %in% names(df2)))
stopifnot(all(matched_fl_cols %in% names(df2)))

# Count paired dominance per genome
df2$PA_Higher_Count <- 0L
df2$FL_Higher_Count <- 0L
df2$Equal_Count     <- 0L

for (i in seq_along(matched_pa_cols)) {
  pa_vec <- df2[[matched_pa_cols[i]]]
  fl_vec <- df2[[matched_fl_cols[i]]]

  df2$PA_Higher_Count <- df2$PA_Higher_Count + as.integer(pa_vec > fl_vec)
  df2$FL_Higher_Count <- df2$FL_Higher_Count + as.integer(fl_vec > pa_vec)
  df2$Equal_Count     <- df2$Equal_Count     + as.integer(pa_vec == fl_vec)
}

df2$Total_Paired_Comparisons <-
  df2$PA_Higher_Count + df2$FL_Higher_Count + df2$Equal_Count

# Dominance percentage (avoid division by zero)
df2$PA_Dominance_Percent <- ifelse(
  df2$Total_Paired_Comparisons > 0,
  100 * df2$PA_Higher_Count / df2$Total_Paired_Comparisons,
  NA_real_
)

df2$FL_Dominance_Percent <- ifelse(
  df2$Total_Paired_Comparisons > 0,
  100 * df2$FL_Higher_Count / df2$Total_Paired_Comparisons,
  NA_real_
)

cat("\nPaired dominance summary (PA_Dominance_Percent):\n")
print(summary(df2$PA_Dominance_Percent))

# ==============================================================================
# 7. Save full output
# ==============================================================================
write.xlsx(df2, OUT_XLSX, overwrite = TRUE)
write.csv(df2,  OUT_CSV,  row.names = FALSE)

cat("\nFull output saved:\n  ", OUT_XLSX, "\n  ", OUT_CSV, "\n")

# ==============================================================================
# 8. Save reduced metadata file (for merging with functional annotations)
# ==============================================================================
metadata_df <- df2 %>%
  select(
    Genome,
    Genome_size_Mbp,
    Num_contigs,
    Mean_PA,
    Mean_FL,
    Prev_PA,
    Prev_FL,
    Total_Prev,
    log2_PA_FL,
    Lifestyle,
    PA_Higher_Count,
    FL_Higher_Count,
    Equal_Count,
    Total_Paired_Comparisons,
    PA_Dominance_Percent,
    FL_Dominance_Percent
  )

write.csv(metadata_df, OUT_META, row.names = FALSE)

cat("Metadata file saved:\n  ", OUT_META, "\n")
cat("\nFinal genome counts:\n")
print(table(metadata_df$Lifestyle))

# ==============================================================================
# 9. Session info
# ==============================================================================
sink("output/sessionInfo_lifestyle_classification.txt")
sessionInfo()
sink()

message("\n✓ Lifestyle classification complete. Outputs in output/")
