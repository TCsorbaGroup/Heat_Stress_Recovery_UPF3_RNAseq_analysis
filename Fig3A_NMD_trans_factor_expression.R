# NMD trans-factor expression plot
# UPF3, SMG7, UPF1 and UPF2 across NT, P/T and P/T +1d treatments.
# Mean log2(DESeq2-normalized count + 1) +/- SE across biological replicates.
# Significance stars are based on DESeq2 genome-wide BH-adjusted padj values.
# Requires an existing DESeqDataSet named `dds`.

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# -------------------------------------------------------------------
# 1. Target genes
# -------------------------------------------------------------------
trans_factors <- tibble::tribble(
  ~gene,        ~Factor,
  "AT1G33980", "UPF3",
  "AT5G19400", "SMG7",
  "AT5G47010", "UPF1",
  "AT2G39260", "UPF2"
) %>%
  mutate(Factor = factor(Factor, levels = c("UPF3", "SMG7", "UPF1", "UPF2")))

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^gene:", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

if (!exists("dds", inherits = TRUE)) {
  stop("The DESeq2 object 'dds' was not found. Load it before running this script.")
}
if (!all(c("genotype", "treatment") %in% names(as.data.frame(SummarizedExperiment::colData(dds))))) {
  stop("colData(dds) must contain 'genotype' and 'treatment' columns.")
}

# -------------------------------------------------------------------
# 2. Extract DESeq2-normalized counts and clean identifiers
# -------------------------------------------------------------------
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
rownames(norm_counts) <- clean_gene_id(rownames(norm_counts))

if (anyDuplicated(rownames(norm_counts))) {
  stop("Cleaning gene IDs created duplicate row names in normalized counts. Resolve duplicates before plotting.")
}

found_genes <- intersect(as.character(trans_factors$gene), rownames(norm_counts))
missing_genes <- setdiff(as.character(trans_factors$gene), rownames(norm_counts))

if (length(missing_genes) > 0L) {
  warning("Target genes not found in dds and omitted: ", paste(missing_genes, collapse = ", "))
}
if (length(found_genes) == 0L) {
  stop("None of the four target genes were found in rownames(dds). Check gene ID format.")
}

factor_info <- trans_factors %>%
  mutate(gene = as.character(gene)) %>%
  filter(gene %in% found_genes) %>%
  mutate(Factor = factor(as.character(Factor), levels = c("UPF3", "SMG7", "UPF1", "UPF2")))

# Use the sample metadata in dds as the source of truth instead of inferring
# genotype and treatment from the sample names.
sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  tibble::rownames_to_column("Sample") %>%
  transmute(
    Sample,
    Genotype = case_when(
      tolower(as.character(genotype)) %in% c("col", "wt", "col-0", "wildtype", "wild_type") ~ "Col-0",
      tolower(as.character(genotype)) %in% c("upf3", "upf-3") ~ "upf3",
      TRUE ~ NA_character_
    ),
    Treatment = toupper(as.character(treatment))
  ) %>%
  filter(Treatment %in% c("NT", "LAT", "REC"), !is.na(Genotype)) %>%
  mutate(
    Genotype = factor(Genotype, levels = c("Col-0", "upf3")),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC"))
  )

if (!all(colnames(norm_counts) %in% sample_info$Sample)) {
  missing_samples <- setdiff(colnames(norm_counts), sample_info$Sample)
  if (length(missing_samples)) {
    stop("Some count-matrix samples could not be matched to dds genotype/treatment metadata: ",
         paste(missing_samples, collapse = ", "))
  }
}

expr_data <- as.data.frame(log2(norm_counts[found_genes, , drop = FALSE] + 1)) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Sample", values_to = "log2_expr") %>%
  left_join(sample_info, by = "Sample") %>%
  filter(!is.na(Genotype), !is.na(Treatment)) %>%
  left_join(factor_info, by = "gene") %>%
  mutate(
    Factor = factor(as.character(Factor), levels = c("UPF3", "SMG7", "UPF1", "UPF2"))
  )

if (nrow(expr_data) == 0L) {
  stop("No samples were classified. Check colData(dds) genotype/treatment values.")
}

# Mean and SE across the biological replicates for each gene / genotype / treatment.
summary_expr <- expr_data %>%
  group_by(gene, Factor, Genotype, Treatment) %>%
  summarise(
    mean_expr = mean(log2_expr, na.rm = TRUE),
    sd_expr = sd(log2_expr, na.rm = TRUE),
    n_samples = sum(!is.na(log2_expr)),
    se_expr = if_else(n_samples > 1L, sd_expr / sqrt(n_samples), NA_real_),
    .groups = "drop"
  )

# -------------------------------------------------------------------
# 3. Genotype contrasts from the existing DESeq2 model for significance stars
#    These are model-based gene-level results, with genome-wide BH adjustment.
# -------------------------------------------------------------------
required_results <- c("genotype_UPF3_vs_Col", "genotypeUPF3.treatmentLAT", "genotypeUPF3.treatmentREC")
missing_results <- setdiff(required_results, DESeq2::resultsNames(dds))
if (length(missing_results) > 0L) {
  stop("Expected DESeq2 coefficient(s) not found: ", paste(missing_results, collapse = ", "),
       ". Available: ", paste(DESeq2::resultsNames(dds), collapse = ", "))
}

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

  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene), Treatment = treatment) %>%
    filter(gene %in% as.character(factor_info$gene)) %>%
    dplyr::select(gene, Treatment, log2FoldChange, pvalue, padj)
}

sig_data <- bind_rows(
  get_genotype_result("NT"),
  get_genotype_result("LAT"),
  get_genotype_result("REC")
) %>%
  left_join(factor_info %>% mutate(Factor = as.character(Factor)), by = "gene") %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01  ~ "**",
      !is.na(padj) & padj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Treatment = factor(Treatment, levels = c("NT", "LAT", "REC")),
    Factor = factor(Factor, levels = c("UPF3", "SMG7", "UPF1", "UPF2"))
  )

factor_ranges <- summary_expr %>%
  group_by(Factor) %>%
  summarise(
    y_min = min(pmax(0, mean_expr - se_expr), na.rm = TRUE),
    y_max = max(mean_expr + se_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(offset = pmax(0.15, 0.08 * pmax(y_max - y_min, 1)))

sig_positions <- summary_expr %>%
  group_by(gene, Factor, Treatment) %>%
  summarise(y_position = max(mean_expr + se_expr, na.rm = TRUE), .groups = "drop") %>%
  left_join(factor_ranges %>% dplyr::select(Factor, offset), by = "Factor") %>%
  mutate(y_position = y_position + offset) %>%
  left_join(sig_data %>% dplyr::select(gene, Treatment, significance), by = c("gene", "Treatment"))

# -------------------------------------------------------------------
# 4. Plot
# -------------------------------------------------------------------
genotype_colors <- c("Col-0" = "grey35", "upf3" = "#E00016")

trans_factor_plot <- ggplot(
  summary_expr,
  aes(x = Treatment, y = mean_expr, color = Genotype, group = Genotype)
) +
  geom_line(linewidth = 1.0) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_expr - se_expr), ymax = mean_expr + se_expr),
    width = 0.10, linewidth = 0.55, na.rm = TRUE
  ) +
  geom_point(size = 2.5) +
  geom_text(
    data = sig_positions %>% filter(!is.na(significance), significance != ""),
    aes(x = Treatment, y = y_position, label = significance),
    inherit.aes = FALSE, color = "black", size = 5
  ) +
  facet_wrap(~Factor, nrow = 1, scales = "free_y", drop = FALSE) +
  scale_color_manual(
    values = genotype_colors,
    breaks = c("Col-0", "upf3"),
    labels = c("Col-0", expression(italic(upf3)))
  ) +
  scale_x_discrete(
    limits = c("NT", "LAT", "REC"),
    labels = c("NT", "P/T", "P/T +1d"),
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  labs(
    title = "NMD trans-factor expression",
    x = NULL,
    y = expression(log[2]("normalized count" + 1)),
    color = "Genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Főcím
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    
    # Panelcímek: UPF3, SMG7, UPF1, UPF2
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey30",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 14,
      face = "plain"
    ),
    
    # Jelmagyarázat
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.key.width = grid::unit(0.70, "cm"),
    legend.spacing.x = grid::unit(0.25, "cm"),
    
    # Tengelycímek
    axis.title.y = element_text(size = 14),
    
    # Tengelyfeliratok és kezelések
    axis.text.x = element_text(
      size = 12.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12.5,
      colour = "black"
    ),
    
    # Panelek közötti távolság
    panel.spacing = grid::unit(0.7, "lines"),
    
    plot.margin = margin(5, 8, 5, 5)
  )

print(trans_factor_plot)

# -------------------------------------------------------------------
# 5. Save plot and summary tables
# -------------------------------------------------------------------
ggplot2::ggsave(
  "NMD_trans_factor_expression.png", trans_factor_plot,
  width = 10.5, height = 3.8, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  "NMD_trans_factor_expression.pdf", trans_factor_plot,
  width = 10.5, height = 3.8, units = "in", bg = "white"
)
write.csv(summary_expr, "NMD_trans_factor_expression_summary.csv", row.names = FALSE)
write.csv(sig_data, "NMD_trans_factor_expression_DESeq2_tests.csv", row.names = FALSE)

cat("\nGenes found (display order):\n")
print(factor_info)
cat("\nSample counts in each genotype/treatment combination:\n")
print(expr_data %>% distinct(Sample, Genotype, Treatment) %>% count(Genotype, Treatment))
