############################################################################
# Five-panel Col-0 heat-stress-response gene-set expression figure
# Uses an existing DESeq2 object named `dds`.
#
# Gene sets defined from Col-0 treatment effects:
#   1) HSR-responsive: induced in Col-0 at LAT OR REC vs NT
#   2) Acute: induced at LAT, but not at REC
#   3) Sustained: induced at both LAT and REC
#   4) Delayed: induced at REC, but not at LAT
#   5) Other Col-0-expressed genes: expressed in at least one Col-0 sample and
#      not meeting the induction criterion.
#
# Induction criterion: BH-adjusted p < 0.05 and log2FC > 1.
# Plotted values: log2(DESeq2 normalized count + 1), averaged over biological
# replicates per gene. Genotype comparisons are paired Wilcoxon tests
# across gene-level means and should be interpreted descriptively.
#
# The HSR-responsive set intentionally overlaps with the Acute, Sustained, and
# Delayed subsets. The three subsets are mutually exclusive by definition.
############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(SummarizedExperiment)
})

# -------------------------------------------------------------------
# 1. Settings
# -------------------------------------------------------------------
alpha_cutoff <- 0.05
lfc_cutoff <- 1

y_axis_min <- 3
y_axis_max <- NA_real_  # NA = determine upper limit from the data

output_prefix <- "Col_HSR_five_gene_sets_expression"

group_order <- c(
  "HS responsive",
  "Acute HS responsive",
  "Sustained HS responsive",
  "Delayed HS responsive",
  "Other expressed"
)

genotype_levels_expected <- c("Col-0", "upf3")
treatment_raw_levels <- c("NT", "LAT", "REC")
treatment_plot_levels <- c("NT", "P/T", "P/T +1d")

# -------------------------------------------------------------------
# 2. Helpers and checks
# -------------------------------------------------------------------
if (!exists("dds", inherits = TRUE)) {
  stop("Required object 'dds' was not found. Run the DESeq2 analysis first.")
}
if (!inherits(dds, "DESeqDataSet")) {
  stop("Object 'dds' must be a DESeqDataSet.")
}
if (!all(c("genotype", "treatment") %in% names(SummarizedExperiment::colData(dds)))) {
  stop("colData(dds) must contain 'genotype' and 'treatment' columns.")
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

normalize_genotype <- function(x) {
  x <- tolower(gsub("[[:space:]]+", "_", trimws(as.character(x))))
  dplyr::case_when(
    x %in% c("col", "wt", "col-0", "col_0", "wildtype", "wild_type") ~ "Col-0",
    x %in% c("upf3", "upf-3", "upf_3") ~ "upf3",
    TRUE ~ NA_character_
  )
}

normalize_treatment <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("NT", "NO TREATMENT", "UNTREATED") ~ "NT",
    x %in% c("LAT", "P/T", "PT") ~ "LAT",
    x %in% c("REC", "P/T +1D", "PT +1D", "PT+1D") ~ "REC",
    TRUE ~ NA_character_
  )
}

required_coef_names <- c("treatment_LAT_vs_NT", "treatment_REC_vs_NT")
observed_coef_names <- DESeq2::resultsNames(dds)
missing_coef_names <- setdiff(required_coef_names, observed_coef_names)
if (length(missing_coef_names) > 0L) {
  stop(
    "Missing expected DESeq2 coefficient(s): ",
    paste(missing_coef_names, collapse = ", "),
    "\nAvailable names: ", paste(observed_coef_names, collapse = ", ")
  )
}

# These main treatment coefficients must represent Col-0 vs NT effects.
meta <- as.data.frame(SummarizedExperiment::colData(dds))
genotype_ref <- if (is.factor(meta$genotype)) levels(meta$genotype)[1] else NA_character_
treatment_ref <- if (is.factor(meta$treatment)) levels(meta$treatment)[1] else NA_character_
if (is.na(genotype_ref) || is.na(treatment_ref)) {
  stop("The 'genotype' and 'treatment' columns in colData(dds) must be factors with explicit reference levels.")
}
if (is.na(normalize_genotype(genotype_ref)) || normalize_genotype(genotype_ref) != "Col-0") {
  stop("The genotype reference level must be Col/Col-0/WT for treatment_LAT_vs_NT and treatment_REC_vs_NT to represent Col-0 effects. Current reference: ", genotype_ref)
}
if (normalize_treatment(treatment_ref) != "NT") {
  stop("The treatment reference level must be NT. Current reference: ", treatment_ref)
}

# -------------------------------------------------------------------
# 3. Define Col-0 treatment-responsive genes using DESeq2 results
# -------------------------------------------------------------------
res_lat_wt <- DESeq2::results(dds, name = "treatment_LAT_vs_NT", alpha = alpha_cutoff)
res_rec_wt <- DESeq2::results(dds, name = "treatment_REC_vs_NT", alpha = alpha_cutoff)

result_to_table <- function(res, suffix) {
  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene)) %>%
    filter(!is.na(gene), nzchar(gene))
  if (anyDuplicated(out$gene)) {
    dup <- unique(out$gene[duplicated(out$gene)])
    stop("Duplicate gene IDs after normalization in the ", suffix,
         " DESeq2 result (examples: ", paste(head(dup, 10), collapse = ", "), ").")
  }
  out
}

lat_results <- result_to_table(res_lat_wt, "LAT") %>%
  transmute(gene, log2FC_LAT = log2FoldChange, pvalue_LAT = pvalue, padj_LAT = padj)
rec_results <- result_to_table(res_rec_wt, "REC") %>%
  transmute(gene, log2FC_REC = log2FoldChange, pvalue_REC = pvalue, padj_REC = padj)

gene_classification <- full_join(lat_results, rec_results, by = "gene") %>%
  mutate(
    is_WT_LAT_induced = !is.na(padj_LAT) & padj_LAT < alpha_cutoff &
      !is.na(log2FC_LAT) & log2FC_LAT > lfc_cutoff,
    is_WT_REC_induced = !is.na(padj_REC) & padj_REC < alpha_cutoff &
      !is.na(log2FC_REC) & log2FC_REC > lfc_cutoff,
    is_HSR_responsive = is_WT_LAT_induced | is_WT_REC_induced,
    is_acute = is_WT_LAT_induced & !is_WT_REC_induced,
    is_sustained = is_WT_LAT_induced & is_WT_REC_induced,
    is_delayed = !is_WT_LAT_induced & is_WT_REC_induced
  )

# -------------------------------------------------------------------
# 4. Normalized counts and Col-0-expressed genes
# -------------------------------------------------------------------
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
count_gene_ids <- clean_gene_id(rownames(norm_counts))
if (anyNA(count_gene_ids) || any(!nzchar(count_gene_ids))) {
  stop("The dds object contains missing or empty gene IDs after normalization.")
}
if (anyDuplicated(count_gene_ids)) {
  dup <- unique(count_gene_ids[duplicated(count_gene_ids)])
  stop("Gene IDs in dds are duplicated after ID normalization (examples: ",
       paste(head(dup, 10), collapse = ", "), ").")
}
rownames(norm_counts) <- count_gene_ids

sample_info_raw <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  tibble::rownames_to_column("Sample")
if (anyDuplicated(sample_info_raw$Sample)) {
  stop("Duplicate sample names found in colData(dds).")
}
if (!all(sample_info_raw$Sample %in% colnames(norm_counts))) {
  stop("Sample names in colData(dds) do not all match normalized-count columns.")
}

sample_mapping <- sample_info_raw %>%
  transmute(
    Sample,
    Genotype = normalize_genotype(genotype),
    Treatment_raw = normalize_treatment(treatment)
  )

unknown_samples <- sample_mapping %>%
  filter(is.na(Genotype) | is.na(Treatment_raw)) %>%
  distinct(Sample)
if (nrow(unknown_samples) > 0L) {
  stop("Unrecognized genotype/treatment labels in colData(dds): ",
       paste(unknown_samples$Sample, collapse = ", "),
       ". Check the mapping helpers before proceeding.")
}

sample_mapping <- sample_mapping %>%
  mutate(
    Genotype = factor(Genotype, levels = genotype_levels_expected),
    Treatment = dplyr::recode(
      Treatment_raw,
      NT = "NT",
      LAT = "P/T",
      REC = "P/T +1d"
    ),
    Treatment = factor(Treatment, levels = treatment_plot_levels)
  )

required_combinations <- tidyr::expand_grid(
  Genotype = factor(genotype_levels_expected, levels = genotype_levels_expected),
  Treatment = factor(treatment_plot_levels, levels = treatment_plot_levels)
)
missing_combinations <- required_combinations %>%
  anti_join(sample_mapping %>% distinct(Genotype, Treatment), by = c("Genotype", "Treatment"))
if (nrow(missing_combinations) > 0L) {
  stop("At least one genotype x treatment combination has no samples. Missing combinations: ",
       paste(paste(missing_combinations$Genotype, missing_combinations$Treatment, sep = "/"), collapse = ", "))
}

wt_samples <- sample_mapping %>% filter(Genotype == "Col-0") %>% pull(Sample)
wt_samples <- intersect(wt_samples, colnames(norm_counts))
if (length(wt_samples) == 0L) stop("No Col-0 samples found in dds.")

wt_expressed_genes <- rownames(norm_counts)[
  rowSums(norm_counts[, wt_samples, drop = FALSE] > 0, na.rm = TRUE) > 0
] %>% unique()

# -------------------------------------------------------------------
# 5. Build the five gene sets
# -------------------------------------------------------------------
# HSR-responsive intentionally overlaps with Acute/Sustained/Delayed.
# The three response subclasses do not overlap with one another.
responsive_genes <- unique(gene_classification$gene[gene_classification$is_HSR_responsive])

other_wt_expressed_genes <- setdiff(wt_expressed_genes, responsive_genes)

gene_sets <- bind_rows(
  gene_classification %>%
    filter(is_HSR_responsive) %>%
    transmute(gene, GeneGroup = group_order[1]),
  gene_classification %>%
    filter(is_acute) %>%
    transmute(gene, GeneGroup = group_order[2]),
  gene_classification %>%
    filter(is_sustained) %>%
    transmute(gene, GeneGroup = group_order[3]),
  gene_classification %>%
    filter(is_delayed) %>%
    transmute(gene, GeneGroup = group_order[4]),
  tibble(gene = other_wt_expressed_genes, GeneGroup = group_order[5])
) %>%
  filter(!is.na(gene), nzchar(gene)) %>%
  distinct(gene, GeneGroup) %>%
  mutate(GeneGroup = factor(GeneGroup, levels = group_order))

if (nrow(gene_sets) == 0L) stop("No genes were assigned to any expression gene set.")

gene_set_sizes <- gene_sets %>%
  dplyr::count(GeneGroup, name = "n_genes") %>%
  right_join(tibble(GeneGroup = factor(group_order, levels = group_order)), by = "GeneGroup") %>%
  mutate(n_genes = tidyr::replace_na(n_genes, 0L)) %>%
  arrange(GeneGroup)
print(gene_set_sizes)

message("Gene sets: HSR-responsive overlaps its three response subclasses by design.")
message("Other expressed genes are Col-0-expressed genes not classified as HSR-induced.")

# -------------------------------------------------------------------
# 6. Expression profiles: replicate-average per gene, then mean across genes
# -------------------------------------------------------------------
expr_long <- as.data.frame(log2(norm_counts + 1)) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Sample", values_to = "log2_expr") %>%
  inner_join(sample_mapping %>% dplyr::select(Sample, Genotype, Treatment), by = "Sample") %>%
  inner_join(gene_sets, by = "gene")

if (nrow(expr_long) == 0L) {
  stop("No expression rows matched the gene sets and dds samples. Check ID normalization.")
}

gene_means <- expr_long %>%
  group_by(gene, GeneGroup, Genotype, Treatment) %>%
  summarise(gene_mean = mean(log2_expr, na.rm = TRUE), .groups = "drop")

summary_expr <- gene_means %>%
  group_by(GeneGroup, Genotype, Treatment) %>%
  summarise(
    mean_expr = mean(gene_mean, na.rm = TRUE),
    sd_expr = sd(gene_mean, na.rm = TRUE),
    n_genes = sum(!is.na(gene_mean)),
    se_expr = if_else(n_genes > 1L, sd_expr / sqrt(n_genes), NA_real_),
    .groups = "drop"
  )

if (nrow(summary_expr) == 0L || !any(is.finite(summary_expr$mean_expr))) {
  stop("No finite summary expression values were calculated.")
}

# -------------------------------------------------------------------
# 7. Descriptive paired Col-0 vs upf3 comparison across genes
# -------------------------------------------------------------------
gene_wide <- gene_means %>%
  dplyr::select(gene, GeneGroup, Treatment, Genotype, gene_mean) %>%
  pivot_wider(names_from = Genotype, values_from = gene_mean)

safe_paired_wilcox <- function(wt, upf3_values) {
  ok <- !is.na(wt) & !is.na(upf3_values)
  if (sum(ok) < 2L) return(NA_real_)
  if (all(wt[ok] == upf3_values[ok])) return(1)
  tryCatch(
    suppressWarnings(stats::wilcox.test(wt[ok], upf3_values[ok],
                                        paired = TRUE, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

test_results <- gene_wide %>%
  group_by(GeneGroup, Treatment) %>%
  summarise(
    n_paired = sum(!is.na(`Col-0`) & !is.na(upf3)),
    pvalue = safe_paired_wilcox(`Col-0`, upf3),
    .groups = "drop"
  ) %>%
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

# -------------------------------------------------------------------
# 8. Plot preparation and five-panel figure
# -------------------------------------------------------------------
if (length(y_axis_min) != 1L || !is.finite(y_axis_min)) {
  stop("y_axis_min must be one finite number.")
}
if (!is.na(y_axis_max) && (!is.finite(y_axis_max) || y_axis_max <= y_axis_min)) {
  stop("y_axis_max must be NA or a finite value greater than y_axis_min.")
}

upper_by_group <- summary_expr %>%
  mutate(mean_plus_se = mean_expr + if_else(is.na(se_expr), 0, se_expr)) %>%
  group_by(GeneGroup, Treatment) %>%
  summarise(y_max = max(mean_plus_se, na.rm = TRUE), .groups = "drop")

max_summary <- max(summary_expr$mean_expr + ifelse(is.na(summary_expr$se_expr), 0, summary_expr$se_expr),
                   na.rm = TRUE)
if (!is.finite(max_summary)) max_summary <- max(summary_expr$mean_expr, na.rm = TRUE)
label_offset <- max(0.25, 0.04 * max(1, max_summary))

sig_positions <- upper_by_group %>%
  left_join(test_results, by = c("GeneGroup", "Treatment")) %>%
  mutate(y_position = y_max + label_offset)

auto_y_max <- max(c(summary_expr$mean_expr + ifelse(is.na(summary_expr$se_expr), 0, summary_expr$se_expr),
                    sig_positions$y_position), na.rm = TRUE)
if (!is.finite(auto_y_max)) auto_y_max <- max_summary + 1
y_upper <- if (is.na(y_axis_max)) {
  max(y_axis_min + 1, auto_y_max + max(0.3, 0.06 * abs(auto_y_max)))
} else {
  y_axis_max
}
if (y_upper <= y_axis_min) y_upper <- y_axis_min + 1

genotype_colors <- c("Col-0" = "grey35", "upf3" = "#D73027")

make_hsr_plot <- function(group_name, show_y = FALSE) {
  plot_data <- summary_expr %>% filter(as.character(GeneGroup) == group_name)
  labels <- sig_positions %>%
    filter(as.character(GeneGroup) == group_name, !is.na(significance), significance != "NA")
  group_n <- gene_set_sizes$n_genes[match(group_name, as.character(gene_set_sizes$GeneGroup))]
  if (length(group_n) == 0L || is.na(group_n)) group_n <- 0L

  p <- ggplot(plot_data,
              aes(x = Treatment, y = mean_expr, color = Genotype, group = Genotype)) +
    geom_line(linewidth = 2.0, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = pmax(y_axis_min, mean_expr - se_expr), ymax = mean_expr + se_expr),
      width = 0.15, linewidth = 0.75, na.rm = TRUE
    ) +
    geom_point(size = 4.0, na.rm = TRUE) +
    geom_text(
      data = labels,
      aes(x = Treatment, y = y_position, label = significance),
      inherit.aes = FALSE, color = "black", size = 7.0, fontface = "bold", na.rm = TRUE
    ) +
    scale_color_manual(
      values = genotype_colors,
      breaks = c("Col-0", "upf3"),
      labels = c("Col-0", expression(italic(upf3))),
      drop = FALSE
    ) +
    scale_x_discrete(limits = treatment_plot_levels, drop = FALSE) +
    scale_y_continuous(
      breaks = pretty(c(y_axis_min, y_upper), n = 6),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(y_axis_min, y_upper), clip = "off") +
    labs(
      title = group_name,
      subtitle = paste0("n = ", group_n ),
      x = NULL,
      y = if (show_y) expression(log[2]("normalized count" + 1)) else NULL,
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 1)),
      plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(t = 0, b = 4)),
      axis.title.y = if (show_y) element_text(size = 17, face = "bold") else element_blank(),
      axis.text.y = if (show_y) element_text(size = 13.5, face = "bold", colour = "black") else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.text.x = element_text(size = 13.5, face = "bold", colour = "black"),
      plot.margin = margin(5, 7, 5, 5)
    )
  p
}

plots <- lapply(seq_along(group_order), function(i) {
  make_hsr_plot(group_order[i], show_y = (i == 1L))
})
names(plots) <- group_order

hsr_expression_plot <- patchwork::wrap_plots(
  plots,
  ncol = 5,
  guides = "collect"
) & theme(legend.position = "top", legend.text = element_text(size = 10, face = "bold"))

print(hsr_expression_plot)

# -------------------------------------------------------------------
# 9. Save figure and result tables
# -------------------------------------------------------------------
ggsave(paste0(output_prefix, ".png"), hsr_expression_plot,
       width = 15, height = 5.5, units = "in", dpi = 300, bg = "white")
ggsave(paste0(output_prefix, ".pdf"), hsr_expression_plot,
       width = 15, height = 8.5, units = "in", bg = "white")

write.csv(gene_classification,
          paste0(output_prefix, "_gene_classification.csv"), row.names = FALSE)
write.csv(gene_sets,
          paste0(output_prefix, "_gene_set_membership.csv"), row.names = FALSE)
write.csv(gene_set_sizes,
          paste0(output_prefix, "_gene_set_sizes.csv"), row.names = FALSE)
write.csv(summary_expr,
          paste0(output_prefix, "_expression_summary_SE_across_genes.csv"), row.names = FALSE)
write.csv(gene_means,
          paste0(output_prefix, "_gene_level_expression_means.csv"), row.names = FALSE)
write.csv(test_results,
          paste0(output_prefix, "_Col0_vs_upf3_paired_Wilcoxon_across_genes.csv"), row.names = FALSE)

message("Finished. Output prefix: ", output_prefix)
message("Gene set membership: HSR-responsive overlaps Acute/Sustained/Delayed by design.")
message("Expression is log2(normalized count + 1); SE is calculated across genes, not biological replicates.")
message("Paired Wilcoxon tests compare per-gene means between genotypes and are descriptive, not biological-replicate tests.")
message("Panel y-axis range: ", y_axis_min, " to ", y_upper)
