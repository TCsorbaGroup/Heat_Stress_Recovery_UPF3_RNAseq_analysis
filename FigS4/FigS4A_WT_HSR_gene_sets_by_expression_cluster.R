# Distribution of Col-0 HSR gene sets across the six expression-trajectory clusters

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ------------------------------ Settings ------------------------------
alpha_cutoff <- 0.05
lfc_cutoff <- 1
cluster_file <- "NMD_expression_trajectory_top2000_by_variance_k6_gene_clusters.csv"
output_prefix <- "Col_HSR_gene_sets_by_expression_cluster"
cluster_levels <- paste0("Cluster ", 1:6)

# ------------------------- Helpers and checks -------------------------
if (!exists("dds", inherits = TRUE)) {
  stop("The DESeq2 object 'dds' is not available. Load it before running this script.")
}
if (!file.exists(cluster_file)) {
  stop("Cannot find cluster assignments: ", cluster_file,
       "\nRun the top-2000 k-means clustering script first in this working directory.")
}

clean_gene_id <- function(x) {
  x <- sub("^gene:", "", as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

# -------------------- 1. Read cluster assignments --------------------
cluster_assignments <- read.csv(cluster_file, stringsAsFactors = FALSE,
                                check.names = FALSE) %>%
  transmute(
    gene = clean_gene_id(gene),
    Cluster = suppressWarnings(as.integer(as.character(Cluster)))
  ) %>%
  filter(!is.na(gene), nzchar(gene), !is.na(Cluster), Cluster %in% 1:6) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Cluster = factor(Cluster, levels = 1:6, labels = cluster_levels))

if (nrow(cluster_assignments) == 0L) {
  stop("No valid cluster assignments found in ", cluster_file,
       ". Expected columns named 'gene' and 'Cluster'.")
}
if (length(unique(as.character(cluster_assignments$Cluster))) != 6L) {
  warning("The cluster assignment file does not contain all six clusters.")
}

message("Cluster membership read from ", cluster_file, ":")
print(cluster_assignments %>% dplyr::count(Cluster, name = "n_cluster_genes"))

# ---------------- 2. Define Col-0-induced gene groups ----------------
expected_coefficients <- c("treatment_LAT_vs_NT", "treatment_REC_vs_NT")
missing_coefficients <- setdiff(expected_coefficients, DESeq2::resultsNames(dds))
if (length(missing_coefficients) > 0L) {
  stop("Missing expected DESeq2 coefficient(s): ",
       paste(missing_coefficients, collapse = ", "),
       "\nAvailable coefficients: ", paste(DESeq2::resultsNames(dds), collapse = ", "))
}

res_lat_wt <- DESeq2::results(dds, name = "treatment_LAT_vs_NT")
res_rec_wt <- DESeq2::results(dds, name = "treatment_REC_vs_NT")

lat_results <- as.data.frame(res_lat_wt) %>%
  rownames_to_column("gene") %>%
  mutate(gene = clean_gene_id(gene)) %>%
  transmute(gene, log2FC_LAT = log2FoldChange, padj_LAT = padj)

rec_results <- as.data.frame(res_rec_wt) %>%
  rownames_to_column("gene") %>%
  mutate(gene = clean_gene_id(gene)) %>%
  transmute(gene, log2FC_REC = log2FoldChange, padj_REC = padj)

gene_classification <- full_join(lat_results, rec_results, by = "gene") %>%
  mutate(
    is_WT_LAT_induced = !is.na(padj_LAT) & padj_LAT < alpha_cutoff &
      !is.na(log2FC_LAT) & log2FC_LAT > lfc_cutoff,
    is_WT_REC_induced = !is.na(padj_REC) & padj_REC < alpha_cutoff &
      !is.na(log2FC_REC) & log2FC_REC > lfc_cutoff,
    hsr_pattern = case_when(
      is_WT_LAT_induced & is_WT_REC_induced ~ "Sustained HS responsive",
      is_WT_LAT_induced & !is_WT_REC_induced ~ "Acute HS responsive",
      !is_WT_LAT_induced & is_WT_REC_induced ~ "Delayed HS responsive",
      TRUE ~ "Other expressed"
    )
  )

# Define Col-0-expressed genes using normalized counts in any Col-0 sample,
# across NT, LAT and REC. This only defines the non-induced remainder set.
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
rownames(norm_counts) <- clean_gene_id(rownames(norm_counts))
sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  rownames_to_column("Sample")
if (!all(c("genotype", "treatment") %in% names(sample_info))) {
  stop("colData(dds) must contain 'genotype' and 'treatment' columns.")
}
wt_labels <- c("COL", "COL-0", "WT", "WILDTYPE", "WILD TYPE")
wt_samples <- sample_info$Sample[
  toupper(trimws(as.character(sample_info$genotype))) %in% wt_labels
]
wt_samples <- intersect(wt_samples, colnames(norm_counts))
if (length(wt_samples) == 0L) {
  stop("No Col-0 samples were identified from colData(dds)$genotype.")
}
wt_expressed_genes <- rownames(norm_counts)[
  rowSums(norm_counts[, wt_samples, drop = FALSE] > 0, na.rm = TRUE) > 0
] %>% unique()

# The four non-exclusive summary sets plus the remainder.
# Group 1 (HS responsive) is the union LAT-induced OR REC-induced.
gene_sets <- list(
  "HS responsive" = gene_classification %>%
    filter(is_WT_LAT_induced | is_WT_REC_induced) %>% pull(gene) %>% unique(),
  "Acute HS responsive" = gene_classification %>%
    filter(is_WT_LAT_induced & !is_WT_REC_induced) %>% pull(gene) %>% unique(),
  "Sustained HS responsive" = gene_classification %>%
    filter(is_WT_LAT_induced & is_WT_REC_induced) %>% pull(gene) %>% unique(),
  "Delayed HS responsive" = gene_classification %>%
    filter(!is_WT_LAT_induced & is_WT_REC_induced) %>% pull(gene) %>% unique(),
  "Other expressed" = gene_classification %>%
    filter(!is_WT_LAT_induced & !is_WT_REC_induced, gene %in% wt_expressed_genes) %>%
    pull(gene) %>% unique()
)

group_order <- names(gene_sets)
group_totals <- tibble(
  GeneGroup = group_order,
  n_total_group = vapply(gene_sets, function(x) length(unique(x)), integer(1))
)
message("Total genes in each HSR gene set:")
print(group_totals)

# ---------------- 3. Count memberships across clusters ----------------
cluster_gene_ids <- cluster_assignments %>%
  mutate(Cluster = as.character(Cluster)) %>%
  dplyr::select(gene, Cluster)

count_by_cluster <- bind_rows(lapply(group_order, function(grp) {
  genes <- unique(gene_sets[[grp]])
  assigned <- tibble(gene = genes) %>% left_join(cluster_gene_ids, by = "gene")
  in_clusters <- assigned %>%
    filter(!is.na(Cluster)) %>%
    dplyr::count(Cluster, name = "n_genes")
  bind_rows(
    tibble(Cluster = cluster_levels) %>%
      left_join(in_clusters, by = "Cluster") %>%
      mutate(n_genes = replace_na(n_genes, 0L), GeneGroup = grp),
    tibble(
      Cluster = "Outside top 2000",
      n_genes = sum(is.na(assigned$Cluster)),
      GeneGroup = grp
    )
  )
})) %>%
  left_join(group_totals, by = "GeneGroup") %>%
  mutate(
    percentage = ifelse(n_total_group > 0, 100 * n_genes / n_total_group, NA_real_),
    GeneGroup = factor(GeneGroup, levels = group_order),
    Cluster = factor(Cluster, levels = c(cluster_levels, "Outside top 2000")),
    cell_label = ifelse(is.na(percentage), "n = 0\nNA", sprintf("n = %d\n%.1f%%", n_genes, percentage))
  )

# Export all gene-set totals and the complete count/percentage table.
write.csv(gene_classification, paste0(output_prefix, "_DESeq2_classification.csv"), row.names = FALSE)
write.csv(group_totals, paste0(output_prefix, "_gene_group_totals.csv"), row.names = FALSE)
write.csv(count_by_cluster, paste0(output_prefix, "_cluster_counts_and_percentages.csv"), row.names = FALSE)
write.csv(
  bind_rows(lapply(group_order, function(grp) tibble(GeneGroup = grp, gene = unique(gene_sets[[grp]])))),
  paste0(output_prefix, "_gene_sets.csv"), row.names = FALSE
)

# -------------------------- 4. Table figure ---------------------------
group_label_map <- setNames(
  paste0(group_order, "\n(N = ", group_totals$n_total_group, ")"),
  group_order
)
count_by_cluster <- count_by_cluster %>%
  mutate(GeneGroupLabel = factor(
    unname(group_label_map[as.character(GeneGroup)]),
    levels = unname(group_label_map[group_order])
  ))

p_table <- ggplot(count_by_cluster,
                  aes(x = GeneGroupLabel, y = Cluster, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = cell_label), size = 5.0, lineheight = 0.95) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", na.value = "grey95",
                      name = "% of\ngene group", limits = c(0, 100)) +
  scale_y_discrete(limits = rev(c(cluster_levels, "Outside top 2000")), drop = FALSE) +
  labs(
    title = "Col-0 heat-stress-response gene groups across expression clusters",
    subtitle = paste0(
      "DESeq2 criteria in Col-0: padj < ", alpha_cutoff,
      " and log2FC > ", lfc_cutoff,
      "; genes clustered from the top-", nrow(cluster_assignments), " trajectory set"
    ),
    x = NULL, y = NULL,
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 13, colour = "black", lineheight = 0.95),
    axis.text.y = element_text(size = 15, colour = "black"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 13, hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = margin(8, 12, 8, 8)
  )

print(p_table)

ggsave(paste0(output_prefix, ".png"), p_table,
       width = 13, height = 5.8, units = "in", dpi = 300, bg = "white")
ggsave(paste0(output_prefix, ".pdf"), p_table,
       width = 13, height = 5.8, units = "in", bg = "white")

message("Done. Outputs saved with prefix: ", output_prefix)
