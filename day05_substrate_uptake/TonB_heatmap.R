############################################################
## TonB-only: Universal ordered heatmaps (MG + MT)
## - Top 50 TonB families (global across all samples)
## - Column order fixed to your custom CP -> DE blocks
## - Colored bars: Bay, Season, Salinity, Size_fraction
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(pheatmap)
})

TOP_N <- 50

############################################################
## 1) Helpers
############################################################

clean_sample <- function(x){
  x %>% as.character() %>% str_trim() %>% str_replace_all("\\s+", "")
}

# MG sample names: CPBay_* / DEBay_* -> CP_* / DE_* and remove extra underscores
normalize_mg_to_meta <- function(x){
  x <- clean_sample(x)
  x <- str_replace(x, "^CPBay_", "CP_")
  x <- str_replace(x, "^DEBay_", "DE_")
  x <- str_replace_all(x, "_", "")
  x <- str_replace(x, "^(CP)(Spr|Sum|Fall)", "CP_\\2")
  x <- str_replace(x, "^(DE)(Spr|Sum|Fall)", "DE_\\2")
  x
}

# MT sample names: remove suffixes like _RNA1/_RNA2/_R1/_R2/_singleton(s)
normalize_mt_to_meta <- function(x){
  x <- clean_sample(x)
  
  x <- str_replace(x, "(_RNA1|_RNA2|_R1|_R2)$", "")
  x <- str_replace(x, "(_singletons?|_singleton|_singl|_single|_unpaired|_orphan)$", "")
  
  x <- str_replace(x, "(\\.RNA1|\\.RNA2|\\.R1|\\.R2)$", "")
  x <- str_replace(x, "(\\.singletons?|\\.singleton|\\.singl|\\.single|\\.unpaired|\\.orphan)$", "")
  
  x
}

safe_row_z <- function(m) {
  mu  <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  m2 <- sweep(m, 1, mu, "-")
  m2 <- sweep(m2, 1, sdv, "/")
  m2[!is.finite(m2)] <- 0
  m2
}

# Make metadata consistent so ordering works reliably
standardize_meta <- function(meta){
  meta <- meta %>%
    mutate(
      Sample = clean_sample(Sample),
      
      Bay = toupper(str_trim(as.character(Bay))),
      Bay = case_when(
        grepl("^CP", Bay) ~ "CP",
        grepl("^DE", Bay) ~ "DE",
        TRUE ~ Bay
      ),
      
      Season = str_to_title(str_trim(as.character(Season))),
      Season = recode(Season,
                      "Spr"="Spring", "Spring"="Spring",
                      "Sum"="Summer", "Summer"="Summer",
                      "Fall"="Fall", "Autumn"="Fall"
      ),
      
      Size_fraction = toupper(str_trim(as.character(Size_fraction))),
      Size_fraction = case_when(
        Size_fraction %in% c("FL","FREE","FREE-LIVING","FREELIVING") ~ "FL",
        Size_fraction %in% c("PA","PARTICLE","PARTICLE-ATTACHED","PARTICLEATTACHED") ~ "PA",
        TRUE ~ Size_fraction
      ),
      
      Salinity = str_to_title(str_trim(as.character(Salinity))),
      Salinity = case_when(
        Salinity %in% c("Low","Medium","High") ~ Salinity,
        suppressWarnings(!is.na(as.numeric(Salinity))) ~ {
          x <- as.numeric(Salinity)
          ifelse(x <= 5, "Low", ifelse(x <= 18, "Medium", "High"))
        },
        TRUE ~ Salinity
      )
    )
  meta
}

# Your desired x-axis order:
# CP: Spring FL, Spring PA, Summer FL, Summer PA
# DE: Fall FL, Fall PA, Spring FL, Spring PA, Summer FL, Summer PA
make_group_key <- function(meta){
  meta$Season <- factor(meta$Season, levels = c("Fall","Spring","Summer"))
  meta$Size_fraction <- factor(meta$Size_fraction, levels = c("FL","PA"))
  meta$Salinity <- factor(meta$Salinity, levels = c("Low","Medium","High"))
  
  block <- paste(meta$Bay, meta$Season, meta$Size_fraction, sep="|")
  desired <- c(
    "CP|Spring|FL", "CP|Spring|PA", "CP|Summer|FL", "CP|Summer|PA",
    "DE|Fall|FL",   "DE|Fall|PA",   "DE|Spring|FL", "DE|Spring|PA",
    "DE|Summer|FL", "DE|Summer|PA"
  )
  
  meta$BlockOrder <- match(block, desired)
  meta
}

############################################################
## 2) TonB heatmap function (MG or MT)
############################################################

make_tonb_heatmap_all <- function(tonb_xlsx,
                                  meta_xlsx,
                                  is_mt,
                                  out_pdf,
                                  main_title = "Top 50 TonB-dependent transporters - ordered"){
  
  # ---- load TonB table ----
  df <- read_excel(tonb_xlsx, sheet = 1)
  
  # first column must be the TonB family/class label
  names(df)[1] <- "TonB"
  
  # normalize sample column names
  old_cols <- colnames(df)[-1]
  new_cols <- if (is_mt) normalize_mt_to_meta(old_cols) else normalize_mg_to_meta(old_cols)
  colnames(df) <- c("TonB", new_cols)
  
  # build matrix
  if (is_mt) {
    # sum duplicates after normalization (e.g., RNA1/RNA2 -> same sample)
    df_long <- df %>%
      pivot_longer(-TonB, names_to="Sample", values_to="val") %>%
      mutate(val = as.numeric(val)) %>%
      group_by(TonB, Sample) %>%
      summarise(val = sum(val, na.rm = TRUE), .groups="drop")
    
    mat <- df_long %>%
      pivot_wider(names_from = Sample, values_from = val, values_fill = 0) %>%
      as.data.frame() %>%
      column_to_rownames("TonB") %>%
      as.matrix()
    
  } else {
    mat <- df %>%
      as.data.frame() %>%
      column_to_rownames("TonB")
    mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
    mat <- as.matrix(mat)
    mat[is.na(mat)] <- 0
  }
  
  mat[is.na(mat)] <- 0
  
  # ---- load metadata ----
  meta <- read_excel(meta_xlsx, sheet = 1) %>% as.data.frame()
  meta <- standardize_meta(meta)
  rownames(meta) <- meta$Sample
  
  # ---- overlap ----
  common <- intersect(colnames(mat), meta$Sample)
  if (length(common) < 2) {
    cat("\nTonB-only examples:\n"); print(head(setdiff(colnames(mat), meta$Sample), 30))
    cat("\nMeta-only examples:\n"); print(head(setdiff(meta$Sample, colnames(mat)), 30))
    stop("Too few overlapping samples between TonB matrix and metadata.")
  }
  
  mat  <- mat[, common, drop=FALSE]
  meta2 <- meta[common, , drop=FALSE]
  
  # ---- keep non-zero rows, select Top N globally ----
  mat <- mat[rowSums(mat) > 0, , drop=FALSE]
  if (nrow(mat) < 2) stop("Too few non-zero TonB rows to plot.")
  
  rs <- rowSums(mat)
  top_ids <- names(sort(rs, decreasing = TRUE))[seq_len(min(TOP_N, length(rs)))]
  mat <- mat[top_ids, , drop=FALSE]
  
  # ---- transform (for visualization) ----
  mat_z <- safe_row_z(log1p(mat))
  
  # ---- enforce x-axis order ----
  meta2 <- make_group_key(meta2)
  
  # push samples that do not match desired blocks to the end
  if (any(is.na(meta2$BlockOrder))) {
    cat("\nWARNING: Some samples could not be assigned BlockOrder; pushed to end:\n")
    print(meta2[is.na(meta2$BlockOrder), c("Bay","Season","Size_fraction","Salinity")])
    meta2$BlockOrder[is.na(meta2$BlockOrder)] <- 999
  }
  
  ord <- order(meta2$BlockOrder, meta2$Salinity, rownames(meta2))
  mat_z <- mat_z[, ord, drop=FALSE]
  
  ann <- meta2[ord, c("Bay","Season","Salinity","Size_fraction"), drop=FALSE]
  
  # ---- plot ----
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  
  pdf(out_pdf, width = 14, height = 8)
  pheatmap(
    mat_z,
    annotation_col = ann,
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 8,
    cluster_cols = FALSE,     # keep your custom order!
    main = paste0(main_title, "\nlog1p + row Z-score (visualization)")
  )
  dev.off()
  
  message("Wrote: ", out_pdf)
}

############################################################
## 3) RUN TonB (MG + MT)
############################################################

meta_xlsx <- "~/Jojy_Research_Sync/Fl_vs_Pa/substrte_uptake_finalfiles_feb19/all_files/final/metadata.xlsx"
out_dir   <- "~/Jojy_Research_Sync/Fl_vs_Pa/substrte_uptake_finalfiles_feb19/all_files/final/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Metagenome TonB
make_tonb_heatmap_all(
  tonb_xlsx = "~/Jojy_Research_Sync/Fl_vs_Pa/substrte_uptake_finalfiles_feb19/all_files/final/MG_TonB_final.xlsx",
  meta_xlsx = meta_xlsx,
  is_mt     = FALSE,
  out_pdf   = file.path(out_dir, "MG_TonB_Top50_ALL_ordered.pdf"),
  main_title= "Top 50 TonB-dependent transporters (Metagenome) - ALL samples ordered"
)

# Metatranscriptome TonB (expression)
make_tonb_heatmap_all(
  tonb_xlsx = "~/Jojy_Research_Sync/Fl_vs_Pa/substrte_uptake_finalfiles_feb19/all_files/final/MT_TonB_final.xlsx",
  meta_xlsx = meta_xlsx,
  is_mt     = TRUE,
  out_pdf   = file.path(out_dir, "MT_TonB_Top50_ALL_ordered.pdf"),
  main_title= "Top 50 TonB-dependent transporters (MT expression) - ALL samples ordered"
)