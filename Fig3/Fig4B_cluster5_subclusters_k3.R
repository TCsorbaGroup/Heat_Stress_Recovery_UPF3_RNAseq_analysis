# Sub-clustering of Cluster 5 (k = 3)
#
# Input files are produced by the main trajectory-clustering script:
#   NMD_expression_trajectory_top2000_by_variance_k6_gene_clusters.csv
#   NMD_expression_trajectory_top2000_by_variance_k6_zscore_profiles.csv
#
# This script takes the genes assigned to one parent cluster in the original
# k = 6 clustering and splits their existing six-condition z-score profiles
# into k = 3 subclusters. It does not rerun the original top-2000 gene selection
# or recalculate the z-scores.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ----------------------------- Settings -----------------------------
input_prefix <- "NMD_expression_trajectory_top2000_by_variance_k6"
parent_cluster <- 5
k_subclusters <- 3
set.seed(123)

out_prefix <- paste0("NMD_cluster", parent_cluster,
                     "_subclusters_k", k_subclusters)

cluster_file <- paste0(input_prefix, "_gene_clusters.csv")
profile_file <- paste0(input_prefix, "_zscore_profiles.csv")

condition_cols <- c(
  "Col-0_NT", "Col-0_LAT", "Col-0_REC",
  "upf3_NT", "upf3_LAT", "upf3_REC"
)

# ----------------------------- 1. Read inputs ------------------------
if (!file.exists(cluster_file)) {
  stop("Cannot find ", cluster_file,
       ". Run the original trajectory-clustering script first, or set input_prefix correctly.")
}
if (!file.exists(profile_file)) {
  stop("Cannot find ", profile_file,
       ". Run the original trajectory-clustering script first, or set input_prefix correctly.")
}

cluster_assignments <- read.csv(
  cluster_file, stringsAsFactors = FALSE, check.names = FALSE
)
z_profiles <- read.csv(
  profile_file, stringsAsFactors = FALSE, check.names = FALSE
)

required_cluster_cols <- c("gene", "Cluster")
if (!all(required_cluster_cols %in% names(cluster_assignments))) {
  stop("The cluster assignment file must contain columns: gene, Cluster.")
}
if (!all(c("gene", condition_cols) %in% names(z_profiles))) {
  stop("The z-score profile file must contain gene and these columns: ",
       paste(condition_cols, collapse = ", "))
}

# Validate gene identifiers instead of silently discarding duplicate rows.
if (anyNA(cluster_assignments$gene) || any(!nzchar(trimws(cluster_assignments$gene)))) {
  stop("The cluster assignment file contains missing or empty gene IDs.")
}
if (anyDuplicated(cluster_assignments$gene)) {
  examples <- unique(cluster_assignments$gene[duplicated(cluster_assignments$gene)])
  stop("Duplicated gene IDs in the cluster assignment file (examples: ",
       paste(head(examples, 10), collapse = ", "), ").")
}
if (anyNA(z_profiles$gene) || any(!nzchar(trimws(z_profiles$gene)))) {
  stop("The z-score profile file contains missing or empty gene IDs.")
}
if (anyDuplicated(z_profiles$gene)) {
  examples <- unique(z_profiles$gene[duplicated(z_profiles$gene)])
  stop("Duplicated gene IDs in the z-score profile file (examples: ",
       paste(head(examples, 10), collapse = ", "), ").")
}

cluster_assignments <- cluster_assignments %>%
  mutate(Cluster = suppressWarnings(as.integer(as.character(Cluster))))
if (anyNA(cluster_assignments$Cluster)) {
  stop("Some Cluster values are missing or not integers.")
}

# Ensure the six profile columns are numeric. Stop if a non-missing value
# cannot be interpreted as a number, rather than silently converting it to NA.
for (col in condition_cols) {
  original <- z_profiles[[col]]
  converted <- suppressWarnings(as.numeric(as.character(original)))
  bad_numeric <- !is.na(original) & is.na(converted)
  if (any(bad_numeric)) {
    stop("Non-numeric values found in profile column '", col, "'.")
  }
  z_profiles[[col]] <- converted
}

parent_genes <- cluster_assignments %>%
  filter(Cluster == parent_cluster) %>%
  pull(gene)

if (length(parent_genes) == 0L) {
  stop("No genes found in parent Cluster ", parent_cluster,
       ". Check parent_cluster and the input cluster assignment file.")
}

missing_profile_genes <- setdiff(parent_genes, z_profiles$gene)
parent_profiles_all <- z_profiles %>%
  filter(gene %in% parent_genes)
parent_profiles <- parent_profiles_all %>%
  filter(if_all(all_of(condition_cols), ~ is.finite(.x)))

if (nrow(parent_profiles) < k_subclusters) {
  stop("Only ", nrow(parent_profiles),
       " genes in parent Cluster ", parent_cluster,
       " have complete, finite profiles; this is too few for k = ",
       k_subclusters, ".")
}

profile_mat <- as.matrix(parent_profiles[, condition_cols, drop = FALSE])
storage.mode(profile_mat) <- "numeric"
rownames(profile_mat) <- parent_profiles$gene

if (nrow(unique(as.data.frame(profile_mat))) < k_subclusters) {
  stop("Fewer than ", k_subclusters,
       " distinct profiles are available; k-means cannot form all subclusters.")
}

message("Parent cluster: ", parent_cluster)
message("Genes assigned to parent cluster: ", length(parent_genes))
message("Genes without a matching profile row: ", length(missing_profile_genes))
message("Genes with a profile row: ", nrow(parent_profiles_all))
message("Genes with complete, finite profiles used for clustering: ", nrow(profile_mat))
if (length(missing_profile_genes) > 0L) {
  message("Examples of genes without profiles: ",
          paste(head(missing_profile_genes, 10), collapse = ", "))
}
if (nrow(parent_profiles_all) > nrow(parent_profiles)) {
  message("Profiles omitted because at least one condition was missing or non-finite: ",
          nrow(parent_profiles_all) - nrow(parent_profiles))
}

# ----------------------------- 2. K-means ----------------------------
set.seed(123)
km_sub <- kmeans(
  profile_mat,
  centers = k_subclusters,
  nstart = 50,
  iter.max = 100
)

subcluster_assignments <- tibble(
  gene = rownames(profile_mat),
  ParentCluster = parent_cluster,
  Subcluster = factor(km_sub$cluster, levels = seq_len(k_subclusters))
)

print(subcluster_assignments %>% dplyr::count(Subcluster, name = "n_genes"))

write.csv(
  subcluster_assignments,
  paste0(out_prefix, "_gene_assignments.csv"),
  row.names = FALSE
)

sub_profiles_df <- data.frame(
  gene = rownames(profile_mat),
  profile_mat,
  check.names = FALSE
) %>%
  left_join(subcluster_assignments, by = "gene")

write.csv(
  sub_profiles_df,
  paste0(out_prefix, "_zscore_profiles.csv"),
  row.names = FALSE
)

# ----------------------------- 3. Long data --------------------------
z_long <- sub_profiles_df %>%
  pivot_longer(
    cols = all_of(condition_cols),
    names_to = "Condition",
    values_to = "z_score"
  ) %>%
  separate(Condition, into = c("Genotype", "Treatment"),
           sep = "_", remove = TRUE) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC")),
    Subcluster = factor(Subcluster, levels = seq_len(k_subclusters))
  )

subcluster_means <- z_long %>%
  group_by(Subcluster, Genotype, Treatment) %>%
  summarise(mean_z = mean(z_score), .groups = "drop")

subcluster_sizes <- subcluster_assignments %>%
  dplyr::count(Subcluster, name = "n_genes") %>%
  mutate(label = paste0("n = ", n_genes))

z_min <- min(z_long$z_score, na.rm = TRUE)
z_max <- max(z_long$z_score, na.rm = TRUE)
z_span <- z_max - z_min
if (!is.finite(z_span) || z_span == 0) z_span <- 1
subcluster_sizes <- subcluster_sizes %>%
  mutate(
    x = factor("REC", levels = c("NT", "LAT", "REC")),
    y = z_min + 0.06 * z_span
  )

# ----------------------------- 4. Plot -------------------------------
trajectory_colors <- c("Col-0" = "#595959", "upf3" = "#E60012")

p_sub <- ggplot() +
  geom_line(
    data = z_long,
    aes(
      x = Treatment, y = z_score,
      group = interaction(gene, Genotype), color = Genotype
    ),
    alpha = 0.12,
    linewidth = 0.3
  ) +
  geom_line(
    data = subcluster_means,
    aes(x = Treatment, y = mean_z, group = Genotype, color = Genotype),
    linewidth = 1.35
  ) +
  geom_point(
    data = subcluster_means,
    aes(x = Treatment, y = mean_z, color = Genotype),
    size = 2.1
  ) +
  geom_text(
    data = subcluster_sizes,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    fontface = "plain",
    family = "serif",
    size = 5.5
  ) +
  facet_wrap(~Subcluster, ncol = 3) +
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
    title = paste0("Cluster ", parent_cluster,
                   " subclusters (k = ", k_subclusters, ")"),
    x = NULL,
    y = "Scaled expression (z-score)",
    color = "Genotype"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    strip.background = element_rect(
      fill = "grey80", colour = "grey35", linewidth = 0.5
    ),
    strip.text = element_text(size = 16, face = "plain"),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.major = element_line(colour = "grey87", linewidth = 0.35),
    panel.spacing = grid::unit(0.6, "lines"),
    plot.title = element_text(size = 20, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_sub)

# ----------------------------- 5. Save -------------------------------
ggsave(
  paste0(out_prefix, ".png"), p_sub,
  width = 10.5, height = 4.7, units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(out_prefix, ".pdf"), p_sub,
  width = 10.5, height = 4.7, units = "in", bg = "white"
)

# Save the genotype x treatment mean trajectory for each subcluster too.
write.csv(
  subcluster_means,
  paste0(out_prefix, "_mean_trajectories.csv"),
  row.names = FALSE
)
