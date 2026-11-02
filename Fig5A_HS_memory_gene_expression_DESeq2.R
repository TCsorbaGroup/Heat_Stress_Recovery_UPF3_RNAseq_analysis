############################################################################
# Heat-stress memory gene expression plot
# HSFA2, HSFA3, HSP18.2, HSP22, plus HSP101 as a control panel.
# Significance stars use the DESeq2 genotype contrast for each treatment;
# padj is the genome-wide BH-adjusted p-value from the fitted DESeq2 model.
# Requires an existing DESeq2 object named `dds` with genotype and treatment
# columns in colData(dds), and the expected interaction-model coefficients.
############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# 0. Check DESeq2 object and model
# -------------------------------------------------------------------------
if (!exists("dds", inherits = TRUE)) {
  stop("Required object 'dds' was not found. Run the DESeq2 analysis first.")
}
if (!inherits(dds, "DESeqDataSet")) {
  stop("Object 'dds' must be a DESeqDataSet.")
}

required_metadata <- c("genotype", "treatment")
metadata_names <- names(SummarizedExperiment::colData(dds))
if (!all(required_metadata %in% metadata_names)) {
  stop("colData(dds) must contain columns: ", paste(required_metadata, collapse = ", "))
}

required_coefficients <- c(
  "genotype_UPF3_vs_Col",
  "genotypeUPF3.treatmentLAT",
  "genotypeUPF3.treatmentREC"
)
missing_coefficients <- setdiff(required_coefficients, DESeq2::resultsNames(dds))
if (length(missing_coefficients)) {
  stop(
    "The expected DESeq2 coefficients were not found: ",
    paste(missing_coefficients, collapse = ", "),
    ". Check the design and factor reference levels in dds."
  )
}

# -------------------------------------------------------------------------
# 1. Target genes
# -------------------------------------------------------------------------
hs_memory_genes <- tibble::tribble(
  ~gene,       ~Gene,
  "AT2G26150", "HSFA2",
  "AT5G03720", "HSFA3",
  "AT5G59720", "HSP18.2",
  "AT4G10250", "HSP22",
  "AT1G74310", "HSP101 (control)"
)

gene_order <- hs_memory_genes$Gene

gene_colors <- c("Col-0" = "grey35", "upf3" = "#E00016")

treatment_order <- c("NT", "LAT", "REC")
treatment_labels <- c("NT" = "NT", "LAT" = "P/T", "REC" = "P/T +1d")

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

# -------------------------------------------------------------------------
# 2. Extract normalized counts and validate gene IDs
# -------------------------------------------------------------------------
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
count_gene_ids <- clean_gene_id(rownames(norm_counts))

if (anyNA(count_gene_ids) || any(!nzchar(count_gene_ids))) {
  stop("The row names of normalized counts contain missing/empty gene IDs after cleaning.")
}
if (anyDuplicated(count_gene_ids)) {
  duplicate_ids <- unique(count_gene_ids[duplicated(count_gene_ids)])
  stop(
    "Gene IDs are duplicated after normalization; refusing to merge them silently. Examples: ",
    paste(head(duplicate_ids, 10), collapse = ", ")
  )
}
rownames(norm_counts) <- count_gene_ids

found_genes <- intersect(hs_memory_genes$gene, rownames(norm_counts))
missing_genes <- setdiff(hs_memory_genes$gene, rownames(norm_counts))
if (length(missing_genes)) {
  warning(
    "These heat-stress memory genes were not found in dds and will be omitted: ",
    paste(missing_genes, collapse = ", ")
  )
}
if (!length(found_genes)) {
  stop("None of the target genes were found in rownames(dds). Check the gene-ID format.")
}
gene_info <- hs_memory_genes %>% filter(gene %in% found_genes)

# -------------------------------------------------------------------------
# 3. Use colData(dds) for sample annotation rather than parsing names only
# -------------------------------------------------------------------------
sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  tibble::rownames_to_column("Sample")

if (!all(colnames(norm_counts) %in% sample_info$Sample)) {
  stop("Some normalized-count columns do not match sample names in colData(dds).")
}
sample_info <- sample_info[match(colnames(norm_counts), sample_info$Sample), , drop = FALSE]
if (anyNA(sample_info$Sample) || !identical(sample_info$Sample, colnames(norm_counts))) {
  stop("Could not align colData(dds) with normalized-count columns.")
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
      treatment_raw %in% c("NT", "UNTREATED", "NON-TREATED", "CONTROL") ~ "NT",
      treatment_raw %in% c("LAT", "P/T", "PT") ~ "LAT",
      treatment_raw %in% c("REC", "RECOVERY", "P/T +1D", "PT+1D", "PT +1D") ~ "REC",
      TRUE ~ NA_character_
    )
  )

if (anyNA(sample_info$Genotype) || anyNA(sample_info$Treatment)) {
  bad <- sample_info %>% filter(is.na(Genotype) | is.na(Treatment)) %>% pull(Sample)
  stop(
    "Could not map genotype/treatment metadata for sample(s): ",
    paste(bad, collapse = ", "),
    ". Check colData(dds)$genotype and $treatment."
  )
}

sample_info <- sample_info %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = treatment_order)
  )

sample_grid <- sample_info %>% dplyr::count(Genotype, Treatment, name = "n_samples")
required_grid <- expand_grid(
  Genotype = factor(c("Col-0", "upf3"), levels = c("Col-0", "upf3")),
  Treatment = factor(treatment_order, levels = treatment_order)
)
missing_combinations <- anti_join(required_grid, sample_grid, by = c("Genotype", "Treatment"))
if (nrow(missing_combinations)) {
  stop("At least one genotype x treatment combination has no samples in colData(dds).")
}
if (any(sample_grid$n_samples < 2L)) {
  warning("At least one genotype x treatment group has fewer than two samples; its SE may be undefined.")
}
message("Sample counts by genotype and treatment:")
print(sample_grid)

# Add 1 to each biological sample's normalized count (without log transformation),
# then calculate mean and SE across biological samples.
expr_data <- as.data.frame(norm_counts[found_genes, , drop = FALSE] + 1) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Sample", values_to = "expr") %>%
  left_join(sample_info %>% dplyr::select(Sample, Genotype, Treatment), by = "Sample") %>%
  left_join(gene_info, by = "gene") %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = treatment_order)
  )

if (!nrow(expr_data) || anyNA(expr_data$Genotype) || anyNA(expr_data$Treatment)) {
  stop("Expression values could not be assigned to valid genotype/treatment groups.")
}

summary_expr <- expr_data %>%
  group_by(gene, Gene, Genotype, Treatment) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    sd_expr = sd(expr, na.rm = TRUE),
    n_samples = sum(!is.na(expr)),
    se_expr = if_else(n_samples > 1L, sd_expr / sqrt(n_samples), NA_real_),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# 4. DESeq2 genotype contrasts and significance stars
# -------------------------------------------------------------------------
# These are treatment-specific upf3-vs-Col contrasts from the interaction model.
# alpha=0.05 aligns independent filtering with the adjusted-p threshold shown below.
get_genotype_result <- function(treatment) {
  res <- switch(
    treatment,
    NT = DESeq2::results(dds, name = "genotype_UPF3_vs_Col", alpha = 0.05),
    LAT = DESeq2::results(
      dds,
      contrast = list(c("genotype_UPF3_vs_Col", "genotypeUPF3.treatmentLAT")),
      alpha = 0.05
    ),
    REC = DESeq2::results(
      dds,
      contrast = list(c("genotype_UPF3_vs_Col", "genotypeUPF3.treatmentREC")),
      alpha = 0.05
    ),
    stop("Unknown treatment: ", treatment)
  )

  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene), Treatment = treatment) %>%
    filter(gene %in% gene_info$gene) %>%
    dplyr::select(gene, Treatment, log2FoldChange, pvalue, padj)

  if (anyDuplicated(out$gene)) {
    stop("Duplicate target-gene rows after ID cleaning in DESeq2 results for treatment ", treatment, ".")
  }
  out
}

sig_data <- bind_rows(lapply(treatment_order, get_genotype_result)) %>%
  left_join(gene_info, by = "gene") %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01 ~ "**",
      !is.na(padj) & padj < 0.05 ~ "*",
      !is.na(padj) ~ "ns",
      TRUE ~ "NA"
    ),
    Treatment = factor(Treatment, levels = treatment_order),
    Gene = factor(Gene, levels = gene_order)
  )

# Put significant labels above the higher genotype mean +/- SE within each facet.
factor_ranges <- summary_expr %>%
  group_by(Gene) %>%
  summarise(
    y_min = min(mean_expr - tidyr::replace_na(se_expr, 0), na.rm = TRUE),
    y_max = max(mean_expr + tidyr::replace_na(se_expr, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(offset = pmax(0.15, 0.08 * pmax(y_max - y_min, 1)))

sig_positions <- summary_expr %>%
  group_by(gene, Gene, Treatment) %>%
  summarise(
    y_position = max(mean_expr + tidyr::replace_na(se_expr, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(factor_ranges %>% dplyr::select(Gene, offset), by = "Gene") %>%
  mutate(y_position = y_position + offset) %>%
  left_join(sig_data %>% dplyr::select(gene, Treatment, significance), by = c("gene", "Treatment"))

# -------------------------------------------------------------------------
# 5. Plot
# -------------------------------------------------------------------------
hs_memory_plot <- ggplot(
  summary_expr,
  aes(x = Treatment, y = mean_expr, color = Genotype, group = Genotype)
) +
  geom_line(linewidth = 1.15) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_expr - se_expr),
      ymax = mean_expr + se_expr
    ),
    width = 0.10,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(size = 2.8) +
  geom_text(
    data = sig_positions %>% filter(significance %in% c("*", "**", "***")),
    aes(x = Treatment, y = y_position, label = significance),
    inherit.aes = FALSE,
    color = "black",
    fontface = "bold",
    size = 7
  ) +
  facet_wrap(~Gene, nrow = 1, scales = "free_y", drop = TRUE) +
  scale_color_manual(
    values = gene_colors,
    breaks = c("Col-0", "upf3"),
    labels = c("Col-0", expression(italic(upf3)))
  ) +
  scale_x_discrete(
    limits = treatment_order,
    labels = treatment_labels,
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.24))) +
  labs(
    title = "Heat-stress memory gene expression",
    subtitle = "Mean ± SE of normalized count + 1; stars: DESeq2 genotype contrast, BH-adjusted p-values",
    x = NULL,
    y = "Normalized count + 1",
    color = "Genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.background = element_rect(fill = "grey85", colour = "grey30", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(colour = "black", face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.spacing = grid::unit(0.8, "lines"),
    plot.margin = margin(7, 10, 7, 7)
  )

print(hs_memory_plot)

# -------------------------------------------------------------------------
# 6. Save plot and tables
# -------------------------------------------------------------------------
ggsave(
  "HS_memory_gene_expression.png",
  hs_memory_plot,
  width = 13, height = 4.8, units = "in", dpi = 300, bg = "white"
)
ggsave(
  "HS_memory_gene_expression.pdf",
  hs_memory_plot,
  width = 13, height = 4.8, units = "in", bg = "white"
)
write.csv(summary_expr, "HS_memory_gene_expression_summary.csv", row.names = FALSE)
write.csv(sig_data, "HS_memory_gene_expression_DESeq2_tests.csv", row.names = FALSE)
write.csv(sample_grid, "HS_memory_gene_expression_sample_counts.csv", row.names = FALSE)

cat("\nHeat-stress memory genes found in dds:\n")
print(gene_info)
cat("\nSample counts in each genotype/treatment combination:\n")
print(sample_grid)
cat("\nDESeq2 genotype contrasts for the target genes (BH-adjusted p-values):\n")
print(sig_data %>% dplyr::select(Gene, Treatment, log2FoldChange, pvalue, padj, significance))
