# NMD target gene expression comparison
# Descriptive values: mean log2(normalized count + 1) across genes, per biological sample.
# Statistical tests: Welch two-sample t-tests across biological replicates, separately for each
# gene group and treatment; BH correction across the six group-by-treatment tests.
# Requires an existing DESeq2 object named `dds` and NMD_target_genes.txt.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# -------------------------------------------------------------------
# 1. Input and gene-ID cleaning
# -------------------------------------------------------------------
nmd_file <- "NMD_target_genes.txt"
if (!exists("dds", inherits = TRUE)) {
  stop("The DESeq2 object `dds` was not found. Run/load the DESeq2 analysis first.")
}
if (!file.exists(nmd_file)) {
  stop("Cannot find ", nmd_file, " in working directory: ", getwd())
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^gene:", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

# Extract Arabidopsis gene IDs robustly, including IDs embedded in other text.
nmd_lines <- readLines(nmd_file, warn = FALSE)
nmd_lines <- sub("#.*$", "", nmd_lines)
nmd_matches <- regmatches(
  nmd_lines,
  gregexpr("AT[1-5CM]G[0-9]{5}(\\.[0-9]+)?", nmd_lines, ignore.case = TRUE)
)
nmd_genes <- unique(clean_gene_id(unlist(nmd_matches, use.names = FALSE)))
nmd_genes <- nmd_genes[!is.na(nmd_genes) & nzchar(nmd_genes)]
if (length(nmd_genes) == 0L) {
  stop("No Arabidopsis TAIR gene IDs were identified in ", nmd_file, ".")
}

# -------------------------------------------------------------------
# 2. Normalized counts and sample metadata
# -------------------------------------------------------------------
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
clean_ids <- clean_gene_id(rownames(norm_counts))
if (anyDuplicated(clean_ids)) {
  dup_ids <- unique(clean_ids[duplicated(clean_ids)])
  stop(
    "Gene-ID cleaning created duplicate IDs (e.g. ", paste(head(dup_ids, 5), collapse = ", "),
    "). Resolve duplicate rows before averaging expression; do not silently count a gene twice."
  )
}
rownames(norm_counts) <- clean_ids

sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  rownames_to_column("Sample")
if (!all(c("genotype", "treatment") %in% names(sample_info))) {
  stop("colData(dds) must contain `genotype` and `treatment` columns.")
}

sample_info <- sample_info %>%
  mutate(
    Genotype = case_when(
      tolower(as.character(genotype)) %in% c("col", "wt", "col-0", "wildtype", "wild_type") ~ "Col-0",
      tolower(as.character(genotype)) %in% c("upf3", "upf-3") ~ "upf3",
      TRUE ~ NA_character_
    ),
    treatment_code = toupper(as.character(treatment)),
    Treatment = case_when(
      treatment_code == "NT" ~ "NT",
      treatment_code == "LAT" ~ "P/T",
      treatment_code == "REC" ~ "P/T +1d",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Sample %in% colnames(norm_counts))

if (anyNA(sample_info$Genotype) || anyNA(sample_info$Treatment)) {
  bad <- sample_info %>% filter(is.na(Genotype) | is.na(Treatment)) %>% pull(Sample)
  stop("Unrecognized genotype/treatment labels in colData(dds) for: ", paste(bad, collapse = ", "))
}
if (nrow(sample_info) != ncol(norm_counts) || anyDuplicated(sample_info$Sample)) {
  stop("Sample metadata and normalized-count columns do not match one-to-one.")
}

sample_info <- sample_info %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "P/T", "P/T +1d"))
  ) %>%
  dplyr::select(Sample, Genotype, Treatment)

sample_design_check <- sample_info %>% count(Genotype, Treatment, name = "n_biological_samples")
print(sample_design_check)
if (any(sample_design_check$n_biological_samples < 2L)) {
  warning("At least one genotype-by-treatment group has fewer than two samples; some tests may be NA.")
}

matched_nmd_genes <- intersect(nmd_genes, rownames(norm_counts))
other_genes <- setdiff(rownames(norm_counts), matched_nmd_genes)
if (length(matched_nmd_genes) == 0L) {
  stop("None of the IDs in NMD_target_genes.txt matched genes in dds after ID cleaning.")
}
if (length(other_genes) == 0L) {
  stop("No 'All other genes' remain after defining the NMD target set.")
}

cat("Unique NMD IDs in file: ", length(nmd_genes), "\n", sep = "")
cat("NMD target genes matched to dds: ", length(matched_nmd_genes), "\n", sep = "")
cat("Other genes in dds: ", length(other_genes), "\n", sep = "")
cat("DESeq2 genes total: ", nrow(norm_counts), "\n", sep = "")

# Save the list of target IDs and matching status for transparent auditing.
gene_membership <- tibble(gene = rownames(norm_counts)) %>%
  mutate(
    GeneGroup = if_else(gene %in% matched_nmd_genes, "NMD target genes", "All other genes"),
    listed_in_NMD_file = gene %in% nmd_genes
  )
write.csv(gene_membership, "NMD_target_expression_gene_membership.csv", row.names = FALSE)

# -------------------------------------------------------------------
# 3. Aggregate expression across genes within each biological sample
# -------------------------------------------------------------------
log_norm <- log2(norm_counts + 1)
expr_long <- as.data.frame(log_norm) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Sample", values_to = "log2_expr") %>%
  left_join(sample_info, by = "Sample") %>%
  mutate(
    GeneGroup = if_else(gene %in% matched_nmd_genes, "NMD target genes", "All other genes")
  ) %>%
  filter(!is.na(Genotype), !is.na(Treatment))

sample_group_means <- expr_long %>%
  group_by(GeneGroup, Genotype, Treatment, Sample) %>%
  summarise(
    sample_mean_expr = mean(log2_expr, na.rm = TRUE),
    n_genes = sum(is.finite(log2_expr)),
    .groups = "drop"
  )

summary_expr <- sample_group_means %>%
  group_by(GeneGroup, Genotype, Treatment) %>%
  summarise(
    mean_expr = mean(sample_mean_expr, na.rm = TRUE),
    sd_expr = sd(sample_mean_expr, na.rm = TRUE),
    n_samples = sum(is.finite(sample_mean_expr)),
    se_expr = if_else(n_samples > 1L, sd_expr / sqrt(n_samples), NA_real_),
    n_genes = max(n_genes, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------
# 4. Genotype tests within each gene group and treatment
# -------------------------------------------------------------------
# Independent biological samples: Welch's two-sample t-test (not paired by gene).
run_replicate_test <- function(dat) {
  x <- dat$sample_mean_expr[dat$Genotype == "Col-0" & is.finite(dat$sample_mean_expr)]
  y <- dat$sample_mean_expr[dat$Genotype == "upf3" & is.finite(dat$sample_mean_expr)]
  tt <- NULL
  if (length(x) >= 2L && length(y) >= 2L) {
    tt <- tryCatch(stats::t.test(y, x, var.equal = FALSE), error = function(e) NULL)
  }
  tibble(
    n_Col0 = length(x),
    n_upf3 = length(y),
    mean_Col0 = if (length(x)) mean(x) else NA_real_,
    mean_upf3 = if (length(y)) mean(y) else NA_real_,
    difference_upf3_minus_Col0 = if (length(x) && length(y)) mean(y) - mean(x) else NA_real_,
    ci_low = if (is.null(tt)) NA_real_ else unname(tt$conf.int[1]),
    ci_high = if (is.null(tt)) NA_real_ else unname(tt$conf.int[2]),
    pvalue = if (is.null(tt)) NA_real_ else tt$p.value
  )
}

test_keys <- sample_group_means %>% distinct(GeneGroup, Treatment)
test_results <- bind_rows(lapply(seq_len(nrow(test_keys)), function(i) {
  gg <- as.character(test_keys$GeneGroup[i])
  tr <- as.character(test_keys$Treatment[i])
  dat <- sample_group_means %>% filter(GeneGroup == gg, Treatment == tr)
  bind_cols(tibble(GeneGroup = gg, Treatment = tr), run_replicate_test(dat))
})) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    significance = case_when(
      is.na(padj) ~ "NA",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("\nWelch t-tests across biological replicate-level gene-set means; BH correction across all six tests:\n")
print(test_results)

# -------------------------------------------------------------------
# 5. Significance label positions and shared y-axis limits
# -------------------------------------------------------------------
sig_positions <- summary_expr %>%
  group_by(GeneGroup, Treatment) %>%
  summarise(y_position = max(mean_expr + ifelse(is.na(se_expr), 0, se_expr), na.rm = TRUE) + 0.12,
            .groups = "drop") %>%
  left_join(test_results %>% dplyr::select(GeneGroup, Treatment, significance), by = c("GeneGroup", "Treatment"))

range_max <- max(c(summary_expr$mean_expr + ifelse(is.na(summary_expr$se_expr), 0, summary_expr$se_expr),
                   sig_positions$y_position + 0.15), na.rm = TRUE)
y_lower <- 5  # Requested fixed lower limit; values below 5 will be clipped from the display.
y_upper <- max(ceiling((range_max + 0.10) * 2) / 2, y_lower + 0.5)
if (!is.finite(y_upper) || y_upper <= y_lower) {
  stop("Could not determine a valid shared y-axis range from the sample-level expression summaries.")
}

# -------------------------------------------------------------------
# 6. Plot
# -------------------------------------------------------------------
genotype_colors <- c("Col-0" = "grey65", "upf3" = "#D73027")

theme_common <- theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.text = element_text(size = 15),
    legend.key.width = grid::unit(0.70, "cm"),
    legend.spacing.x = grid::unit(0.25, "cm"),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_text(size = 12.5),
    axis.ticks.y = element_line()
  )

make_mean_expr_plot <- function(group_name, show_y = TRUE) {
  plot_data <- summary_expr %>% filter(GeneGroup == group_name)
  labels <- sig_positions %>% filter(GeneGroup == group_name)

  ggplot(plot_data, aes(x = Treatment, y = mean_expr, color = Genotype, group = Genotype)) +
    geom_line(linewidth = 1.0) +
    geom_errorbar(
      aes(
        ymin = pmax(0, mean_expr - se_expr),
        ymax = mean_expr + se_expr
      ),
      width = 0.12,
      linewidth = 0.65,
      na.rm = TRUE
    ) +
    geom_point(size = 2.7) +
    geom_text(
      data = labels,
      aes(x = Treatment, y = y_position, label = significance),
      inherit.aes = FALSE,
      color = "black",
      size = 5
    ) +
    scale_color_manual(
      values = genotype_colors,
      breaks = c("Col-0", "upf3"),
      labels = c("Col-0", expression(italic(upf3)))
    ) +
    scale_x_discrete(limits = c("NT", "P/T", "P/T +1d"), drop = FALSE) +
    scale_y_continuous(breaks = pretty(c(y_lower, y_upper), n = 6), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(y_lower, y_upper)) +
    labs(
      title = group_name,
      x = NULL,
      y = if (show_y) expression(log[2]("normalized count" + 1)) else NULL,
      color = NULL
    ) +
    theme_common +
    theme(
      axis.title.y = if (show_y) element_text() else element_blank(),
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank()
    )
}

p_targets <- make_mean_expr_plot("NMD target genes", show_y = TRUE)
p_others <- make_mean_expr_plot("All other genes", show_y = FALSE)

expression_mean_plot_SE <- (p_targets | p_others) +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.key.width = grid::unit(0.70, "cm"),
    legend.spacing.x = grid::unit(0.25, "cm")
  )

print(expression_mean_plot_SE)

# -------------------------------------------------------------------
# 7. Save figure and auditable data tables
# -------------------------------------------------------------------
ggsave(
  "NMD_target_expression_mean_SE.png",
  expression_mean_plot_SE,
  width = 6, height = 5.5, units = "in", dpi = 300, bg = "white"
)
ggsave(
  "NMD_target_expression_mean_SE.pdf",
  expression_mean_plot_SE,
  width = 6, height = 5.5, units = "in", bg = "white"
)

write.csv(sample_group_means, "NMD_target_expression_sample_level_means.csv", row.names = FALSE)
write.csv(summary_expr, "NMD_target_expression_summary_biological_replicates_SE.csv", row.names = FALSE)
write.csv(test_results, "NMD_target_expression_Col0_vs_upf3_biological_replicate_tests.csv", row.names = FALSE)

message("Done. Main plot summarized expression across genes within each biological sample.")
message("The inferential unit is the biological sample; p-values are BH-adjusted across six tests.")
