# K-means clustering of gene-expression trajectories (k = 6)
# Uses DESeq2-normalized counts.
# Per-gene profiles: log2(normalized count + 1), averaged across biological
# replicates within each genotype x treatment, then z-score standardized across
# the six genotype-treatment means for each gene.

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------- Settings -----------------------------
# Optional selection criterion for the genes:
#   "variance" = most variable across the six genotype-treatment means
#   "mean"     = highest average expression across the six means
top_selection <- "variance"
n_top_genes <- 2000
k_clusters <- 6
set.seed(123)

# Output names automatically reflect selection mode and k.
out_prefix <- paste0(
  "NMD_expression_trajectory_top", n_top_genes,
  "_by_", top_selection, "_k", k_clusters
)

# ----------------------- 1. Validate and clean IDs ------------------
if (!exists("dds", inherits = TRUE)) {
  stop("The DESeq2 object 'dds' was not found. Load it before running this script.")
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

norm_counts <- DESeq2::counts(dds, normalized = TRUE)
original_gene_ids <- rownames(norm_counts)
cleaned_gene_ids <- clean_gene_id(original_gene_ids)

if (anyNA(cleaned_gene_ids) || any(!nzchar(cleaned_gene_ids))) {
  stop("One or more gene IDs are empty or missing after cleaning.")
}

duplicate_ids <- unique(cleaned_gene_ids[duplicated(cleaned_gene_ids)])
if (length(duplicate_ids) > 0L) {
  stop(
    "Cleaning gene IDs created duplicate IDs. No rows were merged. Examples: ",
    paste(head(duplicate_ids, 10), collapse = ", "),
    ". Inspect the count matrix and resolve duplicates before clustering."
  )
}
rownames(norm_counts) <- cleaned_gene_ids

# Sample metadata from dds. Map Col/Col-0 to Col-0 and UPF3 to upf3.
sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  tibble::rownames_to_column("Sample")

if (!all(c("genotype", "treatment") %in% names(sample_info))) {
  stop("colData(dds) must contain 'genotype' and 'treatment' columns.")
}

sample_info <- sample_info %>%
  mutate(
    genotype_raw = toupper(trimws(as.character(genotype))),
    treatment_raw = toupper(trimws(as.character(treatment))),
    Genotype = case_when(
      genotype_raw %in% c("COL", "COL-0", "WT", "WILDTYPE", "WILD TYPE") ~ "Col-0",
      genotype_raw %in% c("UPF3", "UPF3 MUTANT", "UPF3-MUTANT", "UPF-3") ~ "upf3",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      treatment_raw %in% c("NT", "NON-TREATED", "UNTREATED") ~ "NT",
      treatment_raw %in% c("LAT", "P/T", "PT") ~ "LAT",
      treatment_raw %in% c("REC", "RECOVERY", "P/T +1D", "PT +1D") ~ "REC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Sample %in% colnames(norm_counts))

if (anyNA(sample_info$Genotype) || anyNA(sample_info$Treatment)) {
  bad <- sample_info %>%
    filter(is.na(Genotype) | is.na(Treatment)) %>%
    dplyr::select(Sample, genotype, treatment)
  print(bad)
  stop("Could not map one or more genotype/treatment labels. Check the printed sample metadata.")
}

if (anyDuplicated(sample_info$Sample)) {
  stop("Sample metadata contain duplicated sample names.")
}
if (!setequal(sample_info$Sample, colnames(norm_counts))) {
  missing_meta <- setdiff(colnames(norm_counts), sample_info$Sample)
  extra_meta <- setdiff(sample_info$Sample, colnames(norm_counts))
  if (length(missing_meta) > 0L) {
    stop("No usable metadata for count-matrix samples: ", paste(missing_meta, collapse = ", "))
  }
}

sample_info <- sample_info %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC"))
  )

# Ensure all six conditions contain at least one sample.
sample_counts <- sample_info %>% dplyr::count(Genotype, Treatment, .drop = FALSE)
print(sample_counts)
required_conditions <- expand_grid(
  Genotype = c("Col-0", "upf3"),
  Treatment = c("NT", "LAT", "REC")
)
missing_conditions <- anti_join(
  required_conditions,
  sample_info %>% distinct(Genotype, Treatment) %>% mutate(across(everything(), as.character)),
  by = c("Genotype", "Treatment")
)
if (nrow(missing_conditions) > 0) {
  print(missing_conditions)
  stop("At least one genotype-treatment condition has no samples.")
}

# -------------------- 2. Build six-condition profiles ---------------
log_expr <- log2(norm_counts + 1)

expr_long <- as.data.frame(log_expr) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "Sample",
    values_to = "log2_expr"
  ) %>%
  left_join(sample_info %>% dplyr::select(Sample, Genotype, Treatment), by = "Sample") %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC"))
  )

# Mean expression across biological replicates for each genotype-treatment.
gene_condition_means <- expr_long %>%
  group_by(gene, Genotype, Treatment) %>%
  summarise(mean_log2_expr = mean(log2_expr, na.rm = TRUE), .groups = "drop") %>%
  mutate(Condition = paste(as.character(Genotype), as.character(Treatment), sep = "_")) %>%
  dplyr::select(gene, Condition, mean_log2_expr) %>%
  pivot_wider(names_from = Condition, values_from = mean_log2_expr)

condition_cols <- c(
  "Col-0_NT", "Col-0_LAT", "Col-0_REC",
  "upf3_NT", "upf3_LAT", "upf3_REC"
)
missing_cols <- setdiff(condition_cols, colnames(gene_condition_means))
if (length(missing_cols) > 0) {
  stop("Missing expression-profile columns: ", paste(missing_cols, collapse = ", "))
}

# Require complete and finite six-condition profiles.
profile_df <- gene_condition_means %>%
  filter(if_all(all_of(condition_cols), ~ is.finite(.x)))
profile_mat_all <- as.matrix(profile_df[, condition_cols, drop = FALSE])
rownames(profile_mat_all) <- profile_df$gene

profile_variance <- apply(profile_mat_all, 1, var)
profile_mean <- rowMeans(profile_mat_all)
valid <- is.finite(profile_variance) & profile_variance > 0 & is.finite(profile_mean)
profile_mat_all <- profile_mat_all[valid, , drop = FALSE]
profile_variance <- profile_variance[valid]
profile_mean <- profile_mean[valid]

if (!top_selection %in% c("variance", "mean")) {
  stop('top_selection must be either "variance" or "mean".')
}
ranking_metric <- if (top_selection == "variance") profile_variance else profile_mean
ranking <- order(ranking_metric, decreasing = TRUE)
n_keep <- min(n_top_genes, length(ranking))
selected_genes <- rownames(profile_mat_all)[ranking[seq_len(n_keep)]]
profile_mat <- profile_mat_all[selected_genes, , drop = FALSE]

if (nrow(profile_mat) <= k_clusters) {
  stop("Too few genes remain for k-means clustering.")
}
if (n_keep < n_top_genes) {
  warning("Only ", n_keep, " genes had complete, non-constant profiles; clustering all of them.")
}

# Z-score each gene across its six genotype-treatment mean values.
z_mat <- t(scale(t(profile_mat)))
z_mat <- z_mat[complete.cases(z_mat), , drop = FALSE]
if (nrow(z_mat) <= k_clusters) stop("Too few genes remain after z-score standardization.")

cat("\nGene counts at each stage:\n")
cat("Genes in dds after DESeq2 filtering: ", nrow(norm_counts), "\n", sep = "")
cat("Genes with complete six-condition profiles: ", nrow(profile_df), "\n", sep = "")
cat("Genes eligible after non-constant profile filter: ", length(profile_variance), "\n", sep = "")
cat("Genes selected for clustering: ", nrow(z_mat), "\n", sep = "")

# ------------------------- 3. K-means (k=6) -------------------------
set.seed(123)
km <- kmeans(z_mat, centers = k_clusters, nstart = 50, iter.max = 100)

cluster_assignments <- tibble(
  gene = rownames(z_mat),
  Cluster = factor(km$cluster, levels = seq_len(k_clusters))
)

write.csv(cluster_assignments,
          paste0(out_prefix, "_gene_clusters.csv"),
          row.names = FALSE)
write.csv(
  data.frame(gene = rownames(z_mat), z_mat, check.names = FALSE),
  paste0(out_prefix, "_zscore_profiles.csv"),
  row.names = FALSE
)

# ------------------------ 4. Trajectory plot ------------------------
z_long <- as.data.frame(z_mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = all_of(condition_cols),
    names_to = "Condition",
    values_to = "z_score"
  ) %>%
  separate(Condition, into = c("Genotype", "Treatment"), sep = "_", remove = TRUE) %>%
  left_join(cluster_assignments %>% mutate(Cluster = as.character(Cluster)), by = "gene") %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC")),
    Cluster = factor(Cluster, levels = as.character(seq_len(k_clusters)))
  )

cluster_means <- z_long %>%
  group_by(Cluster, Genotype, Treatment) %>%
  summarise(mean_z = mean(z_score, na.rm = TRUE), .groups = "drop")

cluster_sizes <- cluster_assignments %>%
  dplyr::count(Cluster, name = "n_genes") %>%
  mutate(label = paste0("n = ", n_genes))

z_min <- min(z_long$z_score, na.rm = TRUE)
z_max <- max(z_long$z_score, na.rm = TRUE)
z_span <- z_max - z_min
if (!is.finite(z_span) || z_span == 0) z_span <- 1

cluster_sizes <- cluster_sizes %>%
  mutate(
    x = factor("REC", levels = c("NT", "LAT", "REC")),
    y = z_min + 0.06 * z_span
  )

trajectory_colors <- c("Col-0" = "#595959", "upf3" = "#E60012")
selection_label <- if (top_selection == "variance") "variance" else "mean expression"

p <- ggplot() +
  geom_line(
    data = z_long,
    aes(x = Treatment, y = z_score,
        group = interaction(gene, Genotype), color = Genotype),
    alpha = 0.10,
    linewidth = 0.25
  ) +
  geom_line(
    data = cluster_means,
    aes(x = Treatment, y = mean_z,
        group = Genotype, color = Genotype),
    linewidth = 1.15
  ) +
  geom_point(
    data = cluster_means,
    aes(x = Treatment, y = mean_z, color = Genotype),
    size = 1.6
  ) +
  geom_text(
    data = cluster_sizes,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    fontface = "italic",
    family = "serif",
    size = 7
  ) +
  facet_wrap(~Cluster, ncol = 3) +
  scale_color_manual(
    values = trajectory_colors,
    breaks = c("upf3", "Col-0"),
    labels = c(expression(italic(upf3)), "Col-0")
  ) +
  scale_x_discrete(
    limits = c("NT", "LAT", "REC"),
    labels = c("NT", "P/T", "P/T +1d"),
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = paste0("Expression trajectory clusters: top ", nrow(z_mat), " genes by ", selection_label),
    x = NULL,
    y = "Scaled expression (z-score)",
    color = "genotype"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 17),
    strip.background = element_rect(fill = "grey80", colour = "grey35", linewidth = 0.5),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major = element_line(colour = "grey87", linewidth = 0.4),
    panel.spacing = grid::unit(0.65, "lines"),
    plot.title = element_text(size = 17, hjust = 0, face = "bold"),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p)

# ------------------------------ 5. Save -----------------------------
ggsave(paste0(out_prefix, ".png"), p,
       width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave(paste0(out_prefix, ".pdf"), p,
       width = 10, height = 7, units = "in", bg = "white")

write.csv(cluster_means,
          paste0(out_prefix, "_cluster_mean_trajectories.csv"),
          row.names = FALSE)

print(cluster_sizes %>% dplyr::select(Cluster, n_genes))
