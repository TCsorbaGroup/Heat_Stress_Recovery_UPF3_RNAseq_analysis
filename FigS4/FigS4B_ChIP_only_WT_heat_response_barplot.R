###############################################################################
# WT heat response of HsfA2/HsfA3 ChIP-only genes
# Outputs: one horizontal bar plot (PNG + PDF), category counts and gene table.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readxl)
})

# ------------------------------ Settings -------------------------------------
chip_file  <- "13059_2023_2970_MOESM3_ESM.xlsx"
chip_sheet <- "gene-clusters-with-short-descri"
alpha_cutoff <- 0.05
lfc_class_cutoff <- 1
output_prefix <- "ChIP_only_WT_heat_response"

# ------------------------------ Helpers --------------------------------------
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

extract_agi <- function(x) {
  x <- as.character(x)
  hits <- regmatches(
    x,
    gregexpr("AT[1-5CM]G[0-9]{5}(?:\\.[0-9]+)?", x,
             ignore.case = TRUE, perl = TRUE)
  )
  vapply(hits, function(z) {
    if (!length(z)) return(NA_character_)
    clean_gene_id(z[[1]])
  }, character(1), USE.NAMES = FALSE)
}

# ----------------------- 1. Validate DESeq2 object ---------------------------
if (!exists("dds", inherits = TRUE) || !inherits(get("dds", inherits = TRUE), "DESeqDataSet")) {
  stop("Load the current DESeq2 object `dds` before running this script.")
}
dds <- get("dds", inherits = TRUE)
cd <- as.data.frame(SummarizedExperiment::colData(dds)) %>% rownames_to_column("Sample")
if (!all(c("genotype", "treatment") %in% names(cd))) {
  stop("colData(dds) must contain `genotype` and `treatment` columns.")
}
coef_needed <- c("treatment_LAT_vs_NT", "treatment_REC_vs_NT")
coef_missing <- setdiff(coef_needed, DESeq2::resultsNames(dds))
if (length(coef_missing)) {
  stop("Missing expected DESeq2 coefficient(s): ", paste(coef_missing, collapse = ", "),
       ". Available: ", paste(DESeq2::resultsNames(dds), collapse = ", "))
}

# The coefficient interpretation below requires Col-0 and NT as reference levels.
if (is.factor(cd$genotype) && !tolower(as.character(levels(cd$genotype)[1])) %in%
    c("col", "col-0", "wt", "wildtype", "wild_type")) {
  stop("The first genotype factor level must be Col-0/WT for these treatment coefficients.")
}
if (is.factor(cd$treatment) && as.character(levels(cd$treatment)[1]) != "NT") {
  stop("The first treatment factor level must be NT for these treatment coefficients.")
}

cd <- cd %>% mutate(
  genotype_chr = tolower(as.character(genotype)),
  treatment_chr = toupper(as.character(treatment)),
  is_wt = genotype_chr %in% c("col", "col-0", "wt", "wildtype", "wild_type")
)
wt_samples <- cd$Sample[cd$is_wt]
if (!length(wt_samples)) stop("No Col-0/WT samples found in colData(dds).")

# Normalized counts and cleaned IDs are used to define expression detectability.
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
ids <- clean_gene_id(rownames(norm_counts))
if (anyNA(ids) || any(!nzchar(ids))) stop("Missing/empty gene IDs in dds.")
if (anyDuplicated(ids)) stop("Gene IDs are duplicated after cleaning in dds; resolve IDs before proceeding.")
rownames(norm_counts) <- ids
wt_samples <- intersect(wt_samples, colnames(norm_counts))
if (!length(wt_samples)) stop("No Col-0/WT sample names matched the normalized count matrix.")
expr_ok <- rownames(norm_counts)[rowSums(norm_counts[, wt_samples, drop = FALSE] > 0, na.rm = TRUE) > 0]

# ----------------------- 2. Define the three HSR classes ---------------------
# In the current model, these treatment coefficients represent Col-0 LAT/REC vs NT.
res_lat <- as.data.frame(DESeq2::results(dds, name = "treatment_LAT_vs_NT", alpha = alpha_cutoff)) %>%
  rownames_to_column("gene_id") %>% mutate(gene_id = clean_gene_id(gene_id)) %>%
  transmute(gene_id, lat_log2FC = log2FoldChange, lat_padj = padj)
res_rec <- as.data.frame(DESeq2::results(dds, name = "treatment_REC_vs_NT", alpha = alpha_cutoff)) %>%
  rownames_to_column("gene_id") %>% mutate(gene_id = clean_gene_id(gene_id)) %>%
  transmute(gene_id, rec_log2FC = log2FoldChange, rec_padj = padj)

wt_effects <- full_join(res_lat, res_rec, by = "gene_id") %>% mutate(
  lat_induced_class = !is.na(lat_padj) & lat_padj < alpha_cutoff &
    !is.na(lat_log2FC) & lat_log2FC > lfc_class_cutoff,
  rec_induced_class = !is.na(rec_padj) & rec_padj < alpha_cutoff &
    !is.na(rec_log2FC) & rec_log2FC > lfc_class_cutoff,
  is_acute = lat_induced_class & !rec_induced_class,
  is_sustained = lat_induced_class & rec_induced_class,
  is_delayed = !lat_induced_class & rec_induced_class
)

acute_genes <- wt_effects$gene_id[wt_effects$is_acute]
sustained_genes <- wt_effects$gene_id[wt_effects$is_sustained]
delayed_genes <- wt_effects$gene_id[wt_effects$is_delayed]
hsr_union <- unique(c(acute_genes, sustained_genes, delayed_genes))

# ----------------------- 3. Read HsfA2/HsfA3 ChIP clusters -------------------
if (!file.exists(chip_file)) {
  stop("ChIP Excel file not found: ", chip_file, " (working directory: ", getwd(), ")")
}
chip_raw <- as.data.frame(readxl::read_excel(chip_file, sheet = chip_sheet), stringsAsFactors = FALSE)
if (ncol(chip_raw) < 2L) stop("The ChIP worksheet must contain gene ID and cluster columns.")
chip_gene <- extract_agi(chip_raw[[1]])
cluster_raw <- as.character(chip_raw[[2]])
cluster_num <- suppressWarnings(as.integer(gsub("^[[:space:]]*[cC]", "", cluster_raw)))
need_parse <- is.na(cluster_num) & !is.na(cluster_raw)
if (any(need_parse)) {
  parsed <- regmatches(cluster_raw[need_parse], regexpr("[0-9]+", cluster_raw[need_parse]))
  cluster_num[need_parse] <- suppressWarnings(as.integer(parsed))
}
chip_table <- tibble(gene_id = chip_gene, cluster = cluster_num) %>%
  filter(!is.na(gene_id), nzchar(gene_id), cluster %in% c(7L, 11L, 12L, 14L)) %>%
  distinct(gene_id, cluster)
if (nrow(chip_table) == 0L) stop("No genes found in ChIP clusters c7/c11/c12/c14; check the worksheet and column order.")
chip_genes <- unique(chip_table$gene_id)
q_chip_only <- setdiff(chip_genes, hsr_union)
if (!length(q_chip_only)) stop("The ChIP-only gene set is empty after excluding the three HSR classes.")

# Keep a representative cluster for each gene (lowest cluster number, matching reference logic).
chip_cl <- chip_table %>% group_by(gene_id) %>% slice_min(cluster, n = 1, with_ties = FALSE) %>% ungroup()

# ----------------------- 4. Classify WT heat response ------------------------
left_tab <- tibble(gene_id = q_chip_only) %>%
  left_join(chip_cl, by = "gene_id") %>%
  left_join(wt_effects, by = "gene_id") %>%
  mutate(
    expressed = gene_id %in% expr_ok,
    wt_response = case_when(
      !expressed ~ "not / lowly expressed",
      (!is.na(lat_padj) & lat_padj < alpha_cutoff & !is.na(lat_log2FC) & lat_log2FC > 0) |
        (!is.na(rec_padj) & rec_padj < alpha_cutoff & !is.na(rec_log2FC) & rec_log2FC > 0) ~
        "weakly induced (below the class cut-off)",
      (!is.na(lat_padj) & lat_padj < alpha_cutoff & !is.na(lat_log2FC) & lat_log2FC < 0) |
        (!is.na(rec_padj) & rec_padj < alpha_cutoff & !is.na(rec_log2FC) & rec_log2FC < 0) ~
        "repressed by heat",
      is.na(lat_padj) & is.na(rec_padj) ~ "not tested",
      TRUE ~ "no significant change"
    )
  )

# Preserve reference-style category ordering in the plot; 'not tested' is retained if present.
response_order <- c(
  "no significant change",
  "weakly induced (below the class cut-off)",
  "not / lowly expressed",
  "repressed by heat",
  "not tested"
)
response_counts <- left_tab %>% dplyr::count(wt_response, name = "n") %>%
  mutate(wt_response = factor(wt_response, levels = response_order)) %>%
  arrange(wt_response)

cat("\nHSR class sizes (strict induction criterion):\n")
print(c(Acute = length(unique(acute_genes)), Sustained = length(unique(sustained_genes)),
        Delayed = length(unique(delayed_genes))))
cat("\nChIP target genes in c7/c11/c12/c14: ", length(chip_genes), "\n", sep = "")
cat("ChIP-only genes: ", length(q_chip_only), "\n", sep = "")
cat("WT heat-response categories:\n")
print(response_counts)

# ----------------------- 5. Horizontal bar plot ------------------------------
plot_data <- response_counts %>% filter(!is.na(wt_response))
wt_plot <- ggplot(plot_data, aes(x = reorder(as.character(wt_response), n), y = n)) +
  geom_col(fill = "#7B3294", width = 0.65) +
  geom_text(aes(label = n), hjust = -0.2, size = 6) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = scales::pretty_breaks()) +
  labs(
    x = NULL,
    y = "Gene number",
    title = paste0("WT heat response of the ", length(q_chip_only), " ChIP-only genes")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 19),
    axis.text.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18)
  )

print(wt_plot)
ggsave(paste0(output_prefix, ".png"), wt_plot, width = 10, height = 5.8, units = "in", dpi = 300, bg = "white")
ggsave(paste0(output_prefix, ".pdf"), wt_plot, width = 10, height = 5.8, units = "in", bg = "white")
write.csv2(left_tab, paste0(output_prefix, "_gene_table.csv"), row.names = FALSE)
write.csv2(response_counts, paste0(output_prefix, "_category_counts.csv"), row.names = FALSE)
message("Finished. Output prefix: ", output_prefix)
