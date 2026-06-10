################################################################################
# Co-occurrence Network Analysis — FL vs PA Estuarine Metagenomes
# Author: Jojy John
# Date:   April 14, 2026
# Description:
#   Genome-resolved co-occurrence network analysis from relative abundance data.
#   Builds Spearman + CLR networks for four ecological groups (CP_FL, CP_PA,
#   DE_FL, DE_PA), computes network statistics, detects hubs and keystone taxa
#   via Zi-Pi analysis, and tests module ecology against salinity and fraction.
################################################################################

# ==============================================================================
# 0. Libraries
# ==============================================================================
library(tidyverse)
library(microeco)      # microtable object
library(NetCoMi)       # network construction
library(igraph)        # graph operations and centrality
library(tidygraph)     # tidy interface to igraph
library(ggraph)        # ggplot2-style network plotting
library(ggrepel)       # non-overlapping text labels
library(readxl)        # read taxonomy from Excel
library(viridis)       # colour palettes
library(patchwork)     # combine ggplots

# ==============================================================================
# 1. Set working directory and create output folders
# ==============================================================================
setwd("~/Jojy_Research_Sync/Fl_vs_Pa/Co-occurance_network_analysis")

dir.create("results_figures", showWarnings = FALSE)
dir.create("results_tables",  showWarnings = FALSE)

# ==============================================================================
# 2. Load input files
# ==============================================================================
abund <- read.delim(
  "genome_relative_abundance.tsv",
  header = TRUE, row.names = 1, check.names = FALSE
)

meta <- read.delim(
  "metadata_updated.csv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Read taxonomy — adjust path/sheet as needed
tax_raw <- read_excel("gtdb_taxonomy.xlsx")   # or read.delim() for TSV

# ==============================================================================
# 3. Clean and align inputs
# ==============================================================================

# Strip hidden whitespace from IDs — silent killer of sample matching
rownames(abund) <- trimws(rownames(abund))
meta$sample     <- trimws(meta$sample)

# Keep only samples present in both objects
common_samples <- intersect(rownames(abund), meta$sample)
abund <- abund[common_samples, ]
meta  <- meta[match(common_samples, meta$sample), ]

# Non-negotiable alignment check
stopifnot(all(rownames(abund) == meta$sample))

# Transpose: genomes as rows, samples as columns
abund_t <- t(abund) %>% as.data.frame()

# ==============================================================================
# 4. Group labels
# ==============================================================================
meta <- meta %>%
  mutate(
    Bay2      = case_when(
      Bay == "Chesapeake" ~ "CP",
      Bay == "Delaware"   ~ "DE",
      TRUE                ~ Bay
    ),
    Fraction2 = case_when(
      size_fraction == "Free Living"        ~ "FL",
      size_fraction == "Particle Attached"  ~ "PA",
      TRUE                                  ~ size_fraction
    ),
    Group = paste(Bay2, Fraction2, sep = "_")
  )
rownames(meta) <- meta$sample

groups <- c("CP_FL", "CP_PA", "DE_FL", "DE_PA")

# ==============================================================================
# 5. Parse GTDB taxonomy
# ==============================================================================
tax_parsed <- tax_raw %>%
  separate(
    gtdb_taxonomy,
    into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
    sep = ";", fill = "right"
  ) %>%
  mutate(across(Domain:Species, ~gsub("^[a-z]__", "", .x)))

# ==============================================================================
# 6. Global prevalence and abundance filtering
# ==============================================================================
prev_global      <- rowSums(abund_t > 0) / ncol(abund_t)
abund_filt       <- abund_t[prev_global >= 0.30, ]                 # ≥30% of all samples
abund_filt       <- abund_filt[rowMeans(abund_filt) > 0.0005, ]   # mean abundance >0.05%

cat("Genomes after global filter:", nrow(abund_filt), "\n")

# ==============================================================================
# 7. Build microeco object
# ==============================================================================
dataset <- microtable$new(otu_table = abund_filt, sample_table = meta)
dataset$tidy_dataset()
dataset$otu_table <- as.matrix(dataset$otu_table)

# ==============================================================================
# 8. Network construction (Spearman + CLR) — per group
# ==============================================================================

# Helper: extract samples for a group and apply within-group filters
get_group_matrix <- function(group_name,
                              prev_cutoff  = 0.40,
                              abund_cutoff = 0.001) {
  samps <- meta$sample[meta$Group == group_name]
  mat   <- abund_filt[, samps, drop = FALSE]

  prev_within  <- rowSums(mat > 0) / ncol(mat)
  mean_within  <- rowMeans(mat)
  mat          <- mat[prev_within >= prev_cutoff & mean_within > abund_cutoff, ]

  cat(group_name, ":", ncol(mat), "samples,", nrow(mat), "genomes\n")
  t(mat)   # netConstruct expects samples × taxa
}

networks <- list()

for (grp in groups) {
  mat_filt <- get_group_matrix(grp)

  if (nrow(mat_filt) < 5 || ncol(mat_filt) < 10) {
    message("Skipping ", grp, " — too few samples or taxa")
    next
  }

  net_g <- netConstruct(
    data        = mat_filt,
    measure     = "spearman",
    normMethod  = "clr",         # centered log-ratio before ranking
    zeroMethod  = "pseudo",      # pseudo-count for zeros
    sparsMethod = "threshold",
    thresh      = 0.30,          # |ρ| ≥ 0.30
    verbose     = 1
  )

  networks[[grp]] <- net_g
}

# ==============================================================================
# 9. Convert to igraph and annotate
# ==============================================================================

make_igraph <- function(net_obj, grp_name) {
  adj <- net_obj$adjaMat1

  g <- graph_from_adjacency_matrix(
    adj, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  g <- delete_edges(g, E(g)[weight == 0])

  E(g)$sign       <- ifelse(E(g)$weight > 0, "positive", "negative")
  E(g)$abs_weight <- abs(E(g)$weight)
  V(g)$degree     <- degree(g)

  # Attach taxonomy
  V(g)$Phylum <- tax_parsed$Phylum[match(V(g)$name, tax_parsed$genome)]
  V(g)$Genus  <- tax_parsed$Genus[match(V(g)$name,  tax_parsed$genome)]
  V(g)$Phylum[is.na(V(g)$Phylum)] <- "Unknown"
  V(g)$Genus[is.na(V(g)$Genus)]   <- "Unknown"

  V(g)$group <- grp_name
  g
}

graphs <- lapply(groups, function(grp) {
  if (!is.null(networks[[grp]])) make_igraph(networks[[grp]], grp)
})
names(graphs) <- groups

# ==============================================================================
# 10. Network statistics
# ==============================================================================
calc_net_stats <- function(x, grp_name) {
  data.frame(
    Group           = grp_name,
    Nodes           = gorder(x),
    Edges           = gsize(x),
    Density         = round(edge_density(x), 4),
    AvgDegree       = round(mean(degree(x)), 2),
    Transitivity    = round(transitivity(x, type = "global"), 4),
    Modularity      = round(modularity(cluster_walktrap(x)), 4),
    Components      = components(x)$no,
    AvgPathLength   = round(mean_distance(x, directed = FALSE, unconnected = TRUE), 3),
    Diameter        = diameter(x, directed = FALSE, unconnected = TRUE),
    PositiveEdges   = sum(E(x)$weight > 0),
    NegativeEdges   = sum(E(x)$weight < 0),
    Pos_Neg_Ratio   = round(sum(E(x)$weight > 0) / (sum(E(x)$weight < 0) + 1), 3)
  )
}

network_stats <- bind_rows(
  lapply(groups, function(grp) calc_net_stats(graphs[[grp]], grp))
)

write.csv(network_stats, "results_tables/network_stats_final.csv", row.names = FALSE)
print(network_stats)

# ==============================================================================
# 11. Hub detection — centrality metrics
# ==============================================================================
calc_node_table <- function(x, grp_name) {
  data.frame(
    Group       = grp_name,
    Genome      = V(x)$name,
    Degree      = degree(x),
    Betweenness = round(betweenness(x, directed = FALSE), 2),
    Closeness   = round(closeness(x, normalized = TRUE), 4),
    Eigenvector = round(eigen_centrality(x)$vector, 4),
    Phylum      = V(x)$Phylum,
    Genus       = V(x)$Genus
  )
}

all_nodes <- bind_rows(
  lapply(groups, function(grp) calc_node_table(graphs[[grp]], grp))
)

write.csv(all_nodes, "results_tables/node_table_with_taxonomy.csv", row.names = FALSE)

top_hubs <- all_nodes %>%
  group_by(Group) %>%
  slice_max(Degree, n = 20) %>%
  ungroup()

write.csv(top_hubs, "results_tables/top_hubs.csv", row.names = FALSE)

# ==============================================================================
# 12. Zi-Pi keystone role analysis
# ==============================================================================
calc_zipi <- function(g_obj, grp_name) {
  wc   <- cluster_walktrap(g_obj)
  memb <- membership(wc)
  A    <- as.matrix(as_adjacency_matrix(g_obj, attr = "weight", sparse = FALSE))
  A_bin <- (A != 0) * 1

  n     <- vcount(g_obj)
  nodes <- V(g_obj)$name
  zi    <- numeric(n)
  pi    <- numeric(n)

  for (i in seq_len(n)) {
    m_i   <- memb[i]
    mod_m <- which(memb == m_i)

    k_i   <- sum(A_bin[i, mod_m]) - A_bin[i, i]  # within-module degree
    k_m   <- sapply(mod_m, function(j) sum(A_bin[j, mod_m]) - A_bin[j, j])
    zi[i] <- if (sd(k_m) == 0) 0 else (k_i - mean(k_m)) / sd(k_m)

    k_total <- sum(A_bin[i, ]) - A_bin[i, i]
    if (k_total == 0) {
      pi[i] <- 0
    } else {
      mod_ids  <- unique(memb)
      pi_terms <- sapply(mod_ids, function(m) {
        k_im <- sum(A_bin[i, which(memb == m)])
        (k_im / k_total)^2
      })
      pi[i] <- 1 - sum(pi_terms)
    }
  }

  data.frame(
    Group  = grp_name,
    Genome = nodes,
    Zi     = round(zi, 4),
    Pi     = round(pi, 4),
    Module = as.integer(memb),
    Phylum = V(g_obj)$Phylum,
    Genus  = V(g_obj)$Genus
  ) %>%
    mutate(
      Role = case_when(
        Zi > 2.5  & Pi > 0.62 ~ "Network hub",
        Zi > 2.5  & Pi <= 0.62 ~ "Module hub",
        Zi <= 2.5 & Pi > 0.62 ~ "Connector",
        TRUE                   ~ "Peripheral"
      )
    )
}

zipi <- bind_rows(lapply(groups, function(grp) calc_zipi(graphs[[grp]], grp)))

write.csv(zipi, "results_tables/zipi_roles.csv", row.names = FALSE)

keystone <- zipi %>%
  filter(Role != "Peripheral") %>%
  arrange(Group, desc(Zi))

write.csv(keystone, "results_tables/keystone_taxa.csv", row.names = FALSE)

cat("Keystone role summary:\n")
print(table(keystone$Role))

# ==============================================================================
# 13. Module membership and ecology tests
# ==============================================================================

# Module membership table
module_membership <- zipi %>%
  select(Group, Genome, Module, Phylum, Genus, Role)

write.csv(module_membership, "results_tables/module_membership.csv", row.names = FALSE)

# Phylum composition per module
module_phylum_summary <- module_membership %>%
  group_by(Group, Module, Phylum) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Group, Module) %>%
  mutate(Prop = n / sum(n)) %>%
  ungroup()

write.csv(
  module_phylum_summary,
  "results_tables/module_phylum_summary.csv",
  row.names = FALSE
)

# Module relative abundance per sample (for ecology tests)
otu_long <- abund_filt %>%
  rownames_to_column("Genome") %>%
  pivot_longer(-Genome, names_to = "sample", values_to = "Abundance")

module_sample_abund <- module_membership %>%
  select(Group, Genome, Module) %>%
  inner_join(otu_long, by = "Genome") %>%
  group_by(Group, Module, sample) %>%
  summarise(ModuleAbundance = sum(Abundance), .groups = "drop") %>%
  left_join(meta %>% select(sample, Salinity, Fraction2, season), by = "sample")

# Module–salinity correlations (Spearman)
module_salinity_stats <- module_sample_abund %>%
  group_by(Group, Module) %>%
  summarise(
    n            = n(),
    rho_salinity = cor(ModuleAbundance, Salinity, method = "spearman", use = "pairwise.complete.obs"),
    p_salinity   = cor.test(ModuleAbundance, Salinity, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(p_adj_salinity = p.adjust(p_salinity, method = "BH")) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      p_adj_salinity <= 0.001 ~ "***",
      p_adj_salinity <= 0.01  ~ "**",
      p_adj_salinity <= 0.05  ~ "*",
      p_adj_salinity <= 0.10  ~ ".",
      TRUE                    ~ ""
    )
  )

write.csv(
  module_salinity_stats,
  "results_tables/module_salinity_stats.csv",
  row.names = FALSE
)

# Module FL vs PA comparison (Wilcoxon per module per bay)
module_fraction_stats <- module_sample_abund %>%
  group_by(Group, Module) %>%
  summarise(
    p_fraction = tryCatch(
      wilcox.test(ModuleAbundance ~ Fraction2)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(p_adj_fraction = p.adjust(p_fraction, method = "BH")) %>%
  ungroup()

write.csv(
  module_fraction_stats,
  "results_tables/module_fraction_stats.csv",
  row.names = FALSE
)

# ==============================================================================
# 14. Figure 3A — Network visualisation (phylum-coloured, 2×2 panel)
# ==============================================================================

# Keep only top N edges for readability
subset_top_edges <- function(g_obj, top_n = 300) {
  edf <- igraph::as_data_frame(g_obj, what = "edges") %>%
    arrange(desc(abs_weight)) %>%
    slice_head(n = top_n)
  keep_nodes <- unique(c(edf$from, edf$to))
  vdf <- igraph::as_data_frame(g_obj, what = "vertices") %>%
    filter(name %in% keep_nodes)
  graph_from_data_frame(edf, vertices = vdf, directed = FALSE)
}

plot_net_phylum <- function(g_obj, title_text, top_labels = 15) {
  tbl <- as_tbl_graph(g_obj) %>%
    activate(nodes) %>%
    mutate(label_me = rank(-degree, ties.method = "first") <= top_labels)

  ggraph(tbl, layout = "fr") +
    geom_edge_link(aes(color = sign, width = abs_weight), alpha = 0.25) +
    geom_node_point(aes(size = degree, fill = Phylum),
                    shape = 21, color = "black", alpha = 0.9) +
    geom_node_text(aes(label = ifelse(label_me, Genus, "")),
                   repel = TRUE, size = 2.5, max.overlaps = 20) +
    scale_edge_color_manual(values = c(positive = "darkgreen", negative = "red3")) +
    scale_edge_width(range = c(0.2, 1.2), guide = "none") +
    scale_size(range = c(2, 8), guide = "none") +
    scale_fill_viridis_d(option = "turbo") +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.title = element_text(face = "bold")
    ) +
    ggtitle(title_text)
}

net_plots <- lapply(groups, function(grp) {
  g_sub <- subset_top_edges(graphs[[grp]], top_n = 300)
  plot_net_phylum(g_sub, title_text = grp)
})
names(net_plots) <- groups

p_network_panel <- (net_plots[["CP_PA"]] | net_plots[["CP_FL"]]) /
                   (net_plots[["DE_PA"]] | net_plots[["DE_FL"]])

ggsave(
  "results_figures/Figure3A_network_phylum_panel.pdf",
  p_network_panel,
  width = 14, height = 11
)
ggsave(
  "results_figures/Figure3A_network_phylum_panel.tiff",
  p_network_panel,
  width = 14, height = 11, dpi = 600, compression = "lzw"
)

# ==============================================================================
# 15. Figure 3B — Zi-Pi role plot (clean main-text version)
# ==============================================================================
zipi_plot <- zipi %>%
  mutate(
    Role  = factor(Role,  levels = c("Peripheral", "Connector", "Module hub", "Network hub")),
    Group = factor(Group, levels = c("CP_FL", "CP_PA", "DE_FL", "DE_PA"))
  )

# Save role count summary for supplementary table
zipi_summary <- zipi_plot %>%
  group_by(Group, Role) %>%
  summarise(n = n(), .groups = "drop")

write.csv(
  zipi_summary,
  "results_tables/Figure3B_ZiPi_role_summary.csv",
  row.names = FALSE
)

p3B <- ggplot(zipi_plot, aes(x = Pi, y = Zi, color = Role, shape = Role)) +
  geom_hline(yintercept = 2.5,  linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = 0.62, linetype = "dashed", color = "grey35") +
  geom_point(size = 2.1, alpha = 0.70) +
  facet_wrap(~Group, nrow = 2) +
  scale_color_manual(
    values = c(
      "Peripheral"  = "#4daf4a",
      "Connector"   = "#e41a1c",
      "Module hub"  = "#377eb8",
      "Network hub" = "#984ea3"
    ),
    na.value = "grey70"
  ) +
  labs(
    x      = "Participation coefficient (Pi)",
    y      = "Within-module connectivity (Zi)",
    color  = "Network role",
    shape  = "Network role",
    title  = "Zi-Pi classification of genome-resolved network roles"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text   = element_text(face = "bold", size = 13),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid   = element_blank()
  )

ggsave(
  "results_figures/Figure3B_ZiPi_role_plot.pdf",
  p3B, width = 8.5, height = 6.5
)
ggsave(
  "results_figures/Figure3B_ZiPi_role_plot.tiff",
  p3B, width = 8.5, height = 6.5, dpi = 600, compression = "lzw"
)

# ==============================================================================
# 16. Figure 3C — Module salinity correlation bar chart
# ==============================================================================
module_sal_plot <- module_salinity_stats %>%
  mutate(Group = factor(Group, levels = c("CP_FL", "CP_PA", "DE_FL", "DE_PA")))

p3C <- ggplot(
  module_sal_plot,
  aes(x = as.factor(Module), y = rho_salinity, fill = rho_salinity)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = sig),
    vjust = ifelse(module_sal_plot$rho_salinity >= 0, -0.4, 1.3),
    size = 5
  ) +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#b2182b",
    midpoint = 0,
    name     = "Spearman rho"
  ) +
  labs(
    x     = "Network module",
    y     = "Correlation with salinity (Spearman rho)",
    title = "Salinity-associated network modules"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text   = element_text(face = "bold", size = 13),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  "results_figures/Figure3C_module_salinity_correlations.pdf",
  p3C, width = 10, height = 4.8
)
ggsave(
  "results_figures/Figure3C_module_salinity_correlations.tiff",
  p3C, width = 10, height = 4.8, dpi = 600, compression = "lzw"
)

# ==============================================================================
# 17. Session info
# ==============================================================================
sink("results_tables/sessionInfo.txt")
sessionInfo()
sink()

message("✓ All outputs written to results_figures/ and results_tables/")
