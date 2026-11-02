############################################################################
# NMD cis-feature analysis for the current DESeq2 workflow
# Primary annotation: Araport11 GFF3
# Primary uORF evidence: RiboTaper calls + RNA/Ribo-seq detectability filter
#
# Panels:
#   A: maximum annotated 3'UTR length across isoforms (bp, log scale)
#   B: 3'UTR-intron prevalence in a representative transcript
#   C: RiboTaper-supported uORF prevalence among Ribo-seq-detectable genes
#
# Key methodological choices:
#   * DE thresholds remain padj < 0.05 and |log2FC| >= 1.
#   * Down / Unchanged / Up are mutually exclusive among testable genes.
#   * padj/log2FC NA genes are reported as Not testable, NOT silently classified
#     as Unchanged and are excluded from the three plotted status groups.
#   * 3'UTR length uses the maximum annotated isoform length; representative
#     transcript length and any-isoform intron status are also saved as checks.
#   * Primary uORF-negative calls are only made where the gene is detectable in
#     the Ribo-seq expression table (TPM >= ribo_tpm_min). Positive experimental
#     calls are retained even if that TPM filter does not flag the gene.
#   * P-values are BH-adjusted separately within each feature family across the
#     six planned comparisons
############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readxl)
  library(patchwork)
})

# -------------------------------------------------------------------------
# 0. Settings
# -------------------------------------------------------------------------
annotation_file <- "Araport11_GFF3_genes_transposons.20241001.gff.gz"
ribotaper_file  <- "Dataset S1_RiboTaper.xlsx"
cips_file       <- "Dataset S5_CiPS.xlsx"
ribo_expr_file  <- "GSE183264_Ribo.genes.results"
ribo_tpm_min    <- 1

padj_cutoff <- 0.05
lfc_cutoff  <- 1

treatment_order <- c("NT", "LAT", "REC")
status_order <- c("Down", "Unch.", "Up")
status_colors <- c(Down = "#08519C", Unch. = "grey80", Up = "#E3120B")

# -------------------------------------------------------------------------
# 0a. Required object / package / input checks
# -------------------------------------------------------------------------
required_objects <- c("dds", "res_NT", "res_LAT", "res_REC")
for (obj in required_objects) {
  if (!exists(obj, inherits = TRUE)) {
    stop("Required object '", obj, "' was not found. Run the first DESeq2 script first.")
  }
}
if (!inherits(dds, "DESeqDataSet")) stop("'dds' must be a DESeqDataSet object.")

required_files <- c(annotation_file, ribotaper_file, cips_file, ribo_expr_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Required input file(s) not found in the working directory:\n - ",
    paste(missing_files, collapse = "\n - "),
    "\nCurrent working directory: ", getwd(),
    "\nThe Ribo-seq expression file is required for a defensible uORF-negative denominator."
  )
}

# -------------------------------------------------------------------------
# 1. Identifier and input helpers
# -------------------------------------------------------------------------
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

clean_tx_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x
}

# Expand GFF3 Parent attributes so each returned row has exactly one Parent.
expand_parent_rows <- function(df) {
  if (nrow(df) == 0L) return(df[0, , drop = FALSE])
  if (!"Parent" %in% names(df)) stop("The imported GFF has no 'Parent' column.")
  parent_values <- lapply(seq_len(nrow(df)), function(i) {
    value <- df$Parent[[i]]
    if (is.null(value) || length(value) == 0L || all(is.na(value))) return(character(0))
    value <- paste(as.character(value), collapse = ",")
    vals <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
    vals[nzchar(vals) & !is.na(vals)]
  })
  keep <- lengths(parent_values) > 0L
  if (!any(keep)) return(df[0, , drop = FALSE])
  source_rows <- which(keep)
  out <- df[rep(source_rows, lengths(parent_values)), , drop = FALSE]
  out$Parent <- unlist(parent_values, use.names = FALSE)
  rownames(out) <- NULL
  out
}

extract_agi <- function(df, col = NULL) {
  if (!is.null(col)) {
    if (!col %in% names(df)) {
      stop("Expected column '", col, "' not found. Available columns: ", paste(names(df), collapse = ", "))
    }
    df <- df[, col, drop = FALSE]
  }
  ids <- unlist(lapply(df, function(x) {
    if (is.character(x) || is.factor(x)) as.character(x) else character(0)
  }), use.names = FALSE)
  matches <- regmatches(ids, gregexpr("AT[1-5CM]G[0-9]{5}", ids, ignore.case = TRUE))
  unique(toupper(unlist(matches, use.names = FALSE)))
}

read_any_table <- function(path, sheet = 1, skip = 0) {
  as.data.frame(readxl::read_excel(path, sheet = sheet, skip = skip))
}

# Stop rather than silently merge different rows if ID-cleaning collapses IDs.
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
raw_count_ids <- rownames(norm_counts)
count_gene_ids <- clean_gene_id(raw_count_ids)
if (anyDuplicated(count_gene_ids)) {
  duplicate_ids <- unique(count_gene_ids[duplicated(count_gene_ids)])
  stop("Gene IDs are duplicated after normalization; investigate before analysis: ",
       paste(head(duplicate_ids, 10), collapse = ", "))
}
rownames(norm_counts) <- count_gene_ids
all_genes <- unique(count_gene_ids[!is.na(count_gene_ids) & nzchar(count_gene_ids)])

sample_info <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
  tibble::rownames_to_column("Sample")
if (!all(c("genotype", "treatment") %in% names(sample_info))) {
  stop("colData(dds) must contain 'genotype' and 'treatment'.")
}
if (!all(sample_info$Sample %in% colnames(norm_counts))) {
  stop("Sample names in colData(dds) do not all match normalized-count columns.")
}

# -------------------------------------------------------------------------
# 2. Treatment-specific expression filter and testable DE status
# -------------------------------------------------------------------------
expressed_genes_by_treatment <- setNames(lapply(treatment_order, function(trt) {
  samples <- sample_info$Sample[as.character(sample_info$treatment) == trt]
  if (!length(samples)) stop("No samples found for treatment: ", trt)
  samples <- intersect(samples, colnames(norm_counts))
  if (!length(samples)) stop("No matching count columns found for treatment: ", trt)
  rownames(norm_counts)[rowSums(norm_counts[, samples, drop = FALSE] > 0, na.rm = TRUE) > 0] |>
    unique()
}), treatment_order)

result_to_table <- function(res, result_name) {
  out <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene)) %>%
    filter(!is.na(gene), nzchar(gene))
  if (anyDuplicated(out$gene)) {
    dup <- unique(out$gene[duplicated(out$gene)])
    stop("Duplicate gene IDs in ", result_name, " after normalization: ", paste(head(dup, 10), collapse = ", "))
  }
  out %>% dplyr::select(gene, log2FoldChange, pvalue, padj)
}

results_tbl <- list(
  NT = result_to_table(res_NT, "res_NT"),
  LAT = result_to_table(res_LAT, "res_LAT"),
  REC = result_to_table(res_REC, "res_REC")
)

prepare_status <- function(result_tbl, treatment_name) {
  tibble(gene = all_genes) %>%
    left_join(result_tbl, by = "gene") %>%
    mutate(
      Testability = if_else(!is.na(padj) & !is.na(log2FoldChange), "Testable", "Not testable"),
      Status = case_when(
        is.na(padj) | is.na(log2FoldChange) ~ "Not testable",
        padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down",
        padj < padj_cutoff & log2FoldChange >=  lfc_cutoff ~ "Up",
        TRUE ~ "Unch."
      ),
      Treatment = treatment_name
    )
}

status_by_treatment <- setNames(lapply(treatment_order, function(trt) {
  prepare_status(results_tbl[[trt]], trt)
}), treatment_order)

# Count testability separately. NA padj/log2FC genes are not called Unchanged.
testability_summary <- bind_rows(lapply(treatment_order, function(trt) {
  status_tbl <- status_by_treatment[[trt]]
  expressed <- expressed_genes_by_treatment[[trt]]
  d <- status_tbl %>% filter(gene %in% expressed)
  tibble(
    Treatment = trt,
    n_expressed = n_distinct(d$gene),
    n_testable = sum(d$Status %in% status_order),
    n_not_testable = sum(d$Status == "Not testable"),
    n_Down = sum(d$Status == "Down"),
    n_Unch. = sum(d$Status == "Unch."),
    n_Up = sum(d$Status == "Up")
  )
}))
print(testability_summary)

# The plotted status groups only contain expressed genes with a usable adjusted p-value and log2FC.
groups <- setNames(lapply(treatment_order, function(trt) {
  status_tbl <- status_by_treatment[[trt]]
  expressed <- expressed_genes_by_treatment[[trt]]
  setNames(lapply(status_order, function(st) {
    unique(status_tbl$gene[status_tbl$Status == st & status_tbl$gene %in% expressed])
  }), status_order)
}), treatment_order)

status_sizes <- bind_rows(lapply(treatment_order, function(trt) {
  tibble(
    Treatment = trt,
    Status = factor(status_order, levels = status_order),
    n_expressed_testable = vapply(status_order, function(st) length(groups[[trt]][[st]]), integer(1))
  )
}))
print(status_sizes)

# -------------------------------------------------------------------------
# 3. 3'UTR length and intron status
# -------------------------------------------------------------------------
gff <- rtracklayer::import(annotation_file)
gff_df <- as.data.frame(gff)
required_gff_cols <- c("type", "start", "end", "ID", "Parent")
if (!all(required_gff_cols %in% names(gff_df))) {
  stop("Imported GFF lacks required columns: ", paste(setdiff(required_gff_cols, names(gff_df)), collapse = ", "))
}

mrna_raw <- gff_df %>%
  filter(type %in% c("mRNA", "transcript")) %>%
  dplyr::select(ID, Parent)
mrna_expanded <- expand_parent_rows(mrna_raw) %>%
  transmute(transcript = clean_tx_id(ID), gene = clean_gene_id(Parent)) %>%
  filter(!is.na(transcript), nzchar(transcript), !is.na(gene), nzchar(gene)) %>%
  distinct(transcript, gene)

map_tx_gene <- function(tx) {
  mapped <- mrna_expanded$gene[match(tx, mrna_expanded$transcript)]
  # Only use transcript-to-gene fallback when the transcript ID itself clearly embeds an AGI ID.
  is_agi_tx <- grepl("^AT[1-5CM]G[0-9]{5}(\\.[0-9]+)?$", tx, ignore.case = TRUE)
  fallback <- ifelse(is_agi_tx, clean_gene_id(tx), NA_character_)
  ifelse(!is.na(mapped) & nzchar(mapped), mapped, fallback)
}

utr_raw <- gff_df %>%
  filter(type == "three_prime_UTR") %>%
  dplyr::select(Parent, start, end)
if (nrow(utr_raw) == 0L) stop("No 'three_prime_UTR' features were found in the GFF file.")

utr_expanded <- expand_parent_rows(utr_raw) %>%
  transmute(
    transcript = clean_tx_id(Parent),
    start = as.numeric(start),
    end = as.numeric(end),
    width = as.numeric(end) - as.numeric(start) + 1
  ) %>%
  filter(!is.na(transcript), nzchar(transcript), is.finite(width), width > 0) %>%
  distinct(transcript, start, end, .keep_all = TRUE) %>%
  mutate(gene = map_tx_gene(transcript)) %>%
  filter(!is.na(gene), nzchar(gene))

utr_by_tx <- utr_expanded %>%
  group_by(gene, transcript) %>%
  summarise(
    utr3_length = sum(width, na.rm = TRUE),
    n_utr_segments = n_distinct(paste(start, end, sep = ":")),
    .groups = "drop"
  )

cds_raw <- gff_df %>%
  filter(type == "CDS") %>%
  dplyr::select(Parent, start, end)
cds_expanded <- expand_parent_rows(cds_raw) %>%
  transmute(
    transcript = clean_tx_id(Parent),
    start = as.numeric(start),
    end = as.numeric(end),
    width = as.numeric(end) - as.numeric(start) + 1
  ) %>%
  filter(!is.na(transcript), nzchar(transcript), is.finite(width), width > 0) %>%
  distinct(transcript, start, end, .keep_all = TRUE)

cds_by_tx <- cds_expanded %>%
  group_by(transcript) %>%
  summarise(cds_length = sum(width, na.rm = TRUE), .groups = "drop")

tx_table <- utr_by_tx %>%
  left_join(cds_by_tx, by = "transcript") %>%
  mutate(cds_length = tidyr::replace_na(cds_length, 0)) %>%
  arrange(gene, desc(cds_length), desc(utr3_length), transcript)

representative_tx <- tx_table %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(
    gene,
    representative_transcript = transcript,
    representative_utr3_length = utr3_length,
    representative_utr_segments = n_utr_segments,
    representative_cds_length = cds_length
  )

# The maximum annotated 3'UTR across isoforms.
gene_utr_length <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    max_utr3_length = max(utr3_length, na.rm = TRUE),
    n_annotated_utr_isoforms = n_distinct(transcript),
    .groups = "drop"
  ) %>%
  left_join(representative_tx, by = "gene")

utr_intron_any <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    has_3utr_intron_any = any(n_utr_segments > 1, na.rm = TRUE),
    n_annotated_utr_isoforms = n_distinct(transcript),
    .groups = "drop"
  )

gene_features_3utr <- representative_tx %>%
  transmute(
    gene,
    representative_transcript,
    has_3utr_intron = representative_utr_segments > 1,
    representative_utr3_length,
    representative_utr_segments
  ) %>%
  left_join(utr_intron_any, by = "gene") %>%
  mutate(has_3utr_intron = as.logical(has_3utr_intron))

# -------------------------------------------------------------------------
# 4. uORF evidence and Ribo-seq detectability
# -------------------------------------------------------------------------
ribotaper_tab <- read_any_table(ribotaper_file, sheet = "S1A_uORFs", skip = 0)
cips_tab      <- read_any_table(cips_file, sheet = "A.ALL", skip = 2)

ribotaper_genes <- extract_agi(ribotaper_tab, "gene_id")
cips_genes      <- extract_agi(cips_tab, "gene_id")
if (!length(ribotaper_genes)) stop("No AGI gene IDs extracted from the RiboTaper table.")
if (!length(cips_genes)) warning("No AGI gene IDs extracted from the CiPS table; CiPS robustness analysis may be empty.")

if (!"AA" %in% names(cips_tab)) {
  warning("CiPS table has no 'AA' column; the >=10 aa robustness set will be empty.")
  cips_10aa_genes <- character(0)
} else {
  cips_10aa_tab <- cips_tab %>%
    mutate(.AA_numeric = suppressWarnings(as.numeric(as.character(AA)))) %>%
    filter(!is.na(.AA_numeric), .AA_numeric >= 10)
  cips_10aa_genes <- extract_agi(cips_10aa_tab, "gene_id")
}

# This table defines which genes could reasonably be called uORF-negative in the assay.
ribo_expr <- read.delim(ribo_expr_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("gene_id", "TPM") %in% names(ribo_expr))) {
  stop("Ribo-seq expression file must contain columns 'gene_id' and 'TPM'.")
}
ribo_expr_clean <- tibble(
  gene = clean_gene_id(ribo_expr$gene_id),
  TPM = suppressWarnings(as.numeric(as.character(ribo_expr$TPM)))
) %>%
  filter(!is.na(gene), nzchar(gene), is.finite(TPM)) %>%
  group_by(gene) %>%
  summarise(TPM = max(TPM, na.rm = TRUE), .groups = "drop")
ribo_detectable_genes <- ribo_expr_clean %>%
  filter(TPM >= ribo_tpm_min) %>%
  pull(gene) %>%
  unique()
if (!length(ribo_detectable_genes)) stop("No genes pass the Ribo-seq detectability threshold; check the TPM file and threshold.")

uorf_external_sets <- list(
  "Ribo-seq uORF (RiboTaper)" = ribotaper_genes,
  "Ribo-seq uORF (CiPS)" = cips_genes,
  "Ribo-seq uORF (CiPS, >=10 aa)" = cips_10aa_genes
)
if (length(ribotaper_genes) && length(cips_genes)) {
  uorf_external_sets[["Ribo-seq uORF (RiboTaper or CiPS)"]] <- union(ribotaper_genes, cips_genes)
}
main_uorf_definition <- "Ribo-seq uORF (RiboTaper)"
main_uorf_col <- make.names(main_uorf_definition)

# Positive calls are retained as positive experimental evidence. A non-positive
# gene is a valid negative only if it was detectable at the transcript level.
uorf_df <- tibble(GENEID = all_genes)
for (nm in names(uorf_external_sets)) {
  target_set <- unique(clean_gene_id(uorf_external_sets[[nm]]))
  flag <- rep(NA, nrow(uorf_df))
  if (length(target_set) > 0L) {
    positive <- uorf_df$GENEID %in% target_set
    detectable <- uorf_df$GENEID %in% ribo_detectable_genes
    flag[detectable & !positive] <- FALSE
    # An experimental positive call is retained even if RNA TPM is below threshold;
    # the direct RiboTaper/CiPS evidence itself demonstrates a positive observation.
    flag[positive] <- TRUE
  } else {
    message("No IDs were parsed for uORF definition '", nm,
            "'; its robustness flags are left NA rather than treating all genes as negatives.")
  }
  uorf_df[[make.names(nm)]] <- flag
}
uorf_df$has_uORF <- uorf_df[[main_uorf_col]]

ribo_positive_not_detectable <- intersect(ribotaper_genes, setdiff(all_genes, ribo_detectable_genes))
message("RiboTaper target genes: ", length(ribotaper_genes),
        "; detected genes (TPM >= ", ribo_tpm_min, "): ", length(ribo_detectable_genes),
        "; RiboTaper positives below TPM threshold (retained as known positives): ",
        length(ribo_positive_not_detectable))
message("RiboTaper calls overlapping DESeq2 universe: ", length(intersect(ribotaper_genes, all_genes)))
message("CiPS target genes: ", length(cips_genes),
        "; CiPS >=10 aa: ", length(cips_10aa_genes))

# -------------------------------------------------------------------------
# 5. Assemble data for the three panels
# -------------------------------------------------------------------------
make_group_gene_table <- function() {
  bind_rows(lapply(treatment_order, function(trt) {
    bind_rows(lapply(status_order, function(st) {
      tibble(gene = groups[[trt]][[st]], Treatment = trt, Status = st)
    }))
  })) %>%
    distinct(Treatment, Status, gene) %>%
    mutate(
      Treatment = factor(Treatment, levels = treatment_order),
      Status = factor(Status, levels = status_order)
    )
}

group_gene_table <- make_group_gene_table()

# Maximum annotated 3'UTR length across isoforms; only positive finite lengths are evaluable.
len_data <- group_gene_table %>%
  left_join(gene_utr_length, by = "gene") %>%
  filter(is.finite(max_utr3_length), max_utr3_length > 0) %>%
  mutate(Status = factor(Status, levels = status_order))

# Sensitivity view for length: representative-transcript 3'UTR length.
rep_len_data <- group_gene_table %>%
  left_join(gene_utr_length, by = "gene") %>%
  filter(is.finite(representative_utr3_length), representative_utr3_length > 0) %>%
  mutate(Status = factor(Status, levels = status_order))

# Main status from the representative transcript; any-isoform flag is saved as robustness check.
intron_gene_status <- gene_features_3utr %>%
  dplyr::select(gene, has_3utr_intron, has_3utr_intron_any) %>%
  filter(!is.na(has_3utr_intron))

# No 3'UTR filter; only Ribo-detectable negatives or experimental positives are evaluable.
uorf_gene_status <- uorf_df %>%
  transmute(gene = GENEID, has_uORF = as.logical(has_uORF))

calculate_percentage_data <- function(feature_status, flag_col, feature_name) {
  bind_rows(lapply(treatment_order, function(trt) {
    bind_rows(lapply(status_order, function(st) {
      genes <- unique(groups[[trt]][[st]])
      elig <- tibble(gene = genes) %>%
        left_join(feature_status, by = "gene") %>%
        filter(!is.na(.data[[flag_col]]))
      n_total <- nrow(elig)
      n_positive <- if (n_total > 0L) sum(elig[[flag_col]] %in% TRUE) else 0L
      tibble(
        Treatment = trt, Status = st, Feature = feature_name,
        n_total = n_total, n_positive = n_positive,
        n_negative = n_total - n_positive,
        percentage = if (n_total > 0L) 100 * n_positive / n_total else NA_real_
      )
    }))
  })) %>%
    mutate(Treatment = factor(Treatment, levels = treatment_order),
           Status = factor(Status, levels = status_order))
}

intron_summary <- calculate_percentage_data(intron_gene_status, "has_3utr_intron", "3′UTR intron (representative transcript)")
intron_any_summary <- calculate_percentage_data(intron_gene_status, "has_3utr_intron_any", "3′UTR intron (any isoform)")
uorf_summary <- calculate_percentage_data(uorf_gene_status, "has_uORF", "RiboTaper uORF")

# Detailed coverage table shows how many group members can contribute to each panel.
annotation_coverage <- bind_rows(lapply(treatment_order, function(trt) {
  bind_rows(lapply(status_order, function(st) {
    g <- unique(groups[[trt]][[st]])
    tibble(
      Treatment = trt, Status = st, n_status_genes = length(g),
      n_with_annotated_3UTR = sum(g %in% gene_utr_length$gene),
      n_with_representative_3UTR_intron_call = sum(g %in% gene_features_3utr$gene),
      n_RiboTaper_evaluable = sum(!is.na(uorf_df$has_uORF[match(g, uorf_df$GENEID)])),
      n_RiboTaper_positive = sum(uorf_df$has_uORF[match(g, uorf_df$GENEID)] %in% TRUE, na.rm = TRUE)
    )
  }))
}))

# -------------------------------------------------------------------------
# 6. Statistical tests
#    Contrasts: Down vs Unchanged and Unchanged vs Up, separately per treatment.
#    BH adjustment across these six contrasts within each feature family.
# -------------------------------------------------------------------------
p_to_stars <- function(p) {
  if (length(p) == 0L || is.na(p) || !is.finite(p)) return("n/a")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

comparisons <- list(c("Down", "Unch."), c("Unch.", "Up"))

get_wilcox_tests <- function(dat, value_col, metric_label) {
  out <- bind_rows(lapply(treatment_order, function(trt) {
    bind_rows(lapply(seq_along(comparisons), function(i) {
      comp <- comparisons[[i]]
      x <- dat[[value_col]][dat$Treatment == trt & dat$Status == comp[1]]
      y <- dat[[value_col]][dat$Treatment == trt & dat$Status == comp[2]]
      x <- x[is.finite(x) & x > 0]
      y <- y[is.finite(y) & y > 0]
      p <- if (length(x) >= 2L && length(y) >= 2L) {
        tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else NA_real_
      tibble(
        feature = metric_label, Treatment = trt, group1 = comp[1], group2 = comp[2],
        comparison_index = i, n1 = length(x), n2 = length(y),
        median1 = if (length(x)) median(x) else NA_real_,
        median2 = if (length(y)) median(y) else NA_real_,
        p = p
      )
    }))
  }))
  out %>% mutate(p_adj_BH = p.adjust(p, method = "BH"), label = vapply(p_adj_BH, p_to_stars, character(1)),
                Treatment = factor(Treatment, levels = treatment_order))
}

get_fisher_tests <- function(dat, metric_label) {
  out <- bind_rows(lapply(treatment_order, function(trt) {
    bind_rows(lapply(seq_along(comparisons), function(i) {
      comp <- comparisons[[i]]
      d1 <- dat %>% filter(Treatment == trt, Status == comp[1])
      d2 <- dat %>% filter(Treatment == trt, Status == comp[2])
      n1 <- if (nrow(d1)) sum(d1$n_total, na.rm = TRUE) else 0L
      p1 <- if (nrow(d1)) sum(d1$n_positive, na.rm = TRUE) else 0L
      n2 <- if (nrow(d2)) sum(d2$n_total, na.rm = TRUE) else 0L
      p2 <- if (nrow(d2)) sum(d2$n_positive, na.rm = TRUE) else 0L
      p <- or <- ci_low <- ci_high <- NA_real_
      if (n1 > 0L && n2 > 0L) {
        tab <- matrix(c(p1, n1 - p1, p2, n2 - p2), nrow = 2, byrow = TRUE,
                      dimnames = list(c(comp[1], comp[2]), c("positive", "negative")))
        ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
        if (!is.null(ft)) {
          p <- ft$p.value
          or <- unname(ft$estimate)
          if (length(ft$conf.int) == 2L) {
            ci_low <- unname(ft$conf.int[1])
            ci_high <- unname(ft$conf.int[2])
          }
        }
      }
      tibble(
        feature = metric_label, Treatment = trt, group1 = comp[1], group2 = comp[2],
        comparison_index = i, n1 = n1, positive1 = p1, n2 = n2, positive2 = p2,
        odds_ratio = or, CI_low = ci_low, CI_high = ci_high, p = p
      )
    }))
  }))
  out %>% mutate(p_adj_BH = p.adjust(p, method = "BH"), label = vapply(p_adj_BH, p_to_stars, character(1)),
                Treatment = factor(Treatment, levels = treatment_order))
}

len_tests <- get_wilcox_tests(len_data, "max_utr3_length", "Maximum 3′UTR length")
rep_len_tests <- get_wilcox_tests(rep_len_data, "representative_utr3_length", "Representative-transcript 3′UTR length")
intron_tests <- get_fisher_tests(intron_summary, "3′UTR intron (representative transcript)")
intron_any_tests <- get_fisher_tests(intron_any_summary, "3′UTR intron (any isoform)")
uorf_tests <- get_fisher_tests(uorf_summary, "RiboTaper-supported uORF")

# -------------------------------------------------------------------------
# 7. Plot
# -------------------------------------------------------------------------
slide_width_cm <- 19.05
slide_height_cm <- 25.40
individual_panel_fraction <- 1 / 6
slide_area_in2 <- (slide_width_cm / 2.54) * (slide_height_cm / 2.54)
individual_panel_area_in2 <- slide_area_in2 * individual_panel_fraction
individual_plot_aspect_ratio <- 4.50 / 2.78
individual_plot_width_in <- sqrt(individual_panel_area_in2 * individual_plot_aspect_ratio)
individual_plot_height_in <- sqrt(individual_panel_area_in2 / individual_plot_aspect_ratio)

facet_treatment_labels <- c(
  NT = "NT",
  LAT = "P/T",
  REC = "P/T +1d"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", colour = "grey70", linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 9.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour = "black", size = 8.0),
    axis.text.y = element_text(colour = "black", size = 8.0),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 4)),
    plot.title = element_blank(),
    plot.margin = margin(3, 4, 3, 3, unit = "pt")
  )

if (nrow(len_data) == 0L) stop("No genes with positive annotated 3′UTR lengths remain for Panel A.")

len_max_by_treatment <- len_data %>%
  group_by(Treatment) %>%
  summarise(max_value = max(max_utr3_length, na.rm = TRUE), .groups = "drop")

sig_len <- len_tests %>%
  filter(!is.na(p_adj_BH)) %>%
  left_join(len_max_by_treatment, by = "Treatment") %>%
  mutate(
    x1 = match(group1, status_order), x2 = match(group2, status_order),
    y_position = max_value * if_else(comparison_index == 1L, 1.15, 1.45),
    label_y = y_position * 1.5
  )
len_y_upper <- max(c(len_data$max_utr3_length, sig_len$label_y), na.rm = TRUE) * 1.08
len_y_upper <- max(len_y_upper, 100)

pA <- ggplot(len_data, aes(x = Status, y = max_utr3_length, fill = Status)) +
  geom_boxplot(width = 0.68, outlier.alpha = 0.12, na.rm = TRUE) +
  geom_segment(data = sig_len, aes(x = x1, xend = x2, y = y_position, yend = y_position),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = sig_len, aes(x = x1, xend = x1, y = y_position * 0.96, yend = y_position),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = sig_len, aes(x = x2, xend = x2, y = y_position * 0.96, yend = y_position),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_len, aes(x = (x1 + x2) / 2, y = label_y, label = label),
            inherit.aes = FALSE, size = 4.0,) +
  facet_wrap(~Treatment, nrow = 1,
             labeller = labeller(Treatment = facet_treatment_labels)) +
  scale_x_discrete(limits = status_order, drop = FALSE) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(1, len_y_upper), clip = "off") +
  labs(y = "Maximum 3′UTR length (bp, log scale)", x = NULL) +
  base_theme

make_percentage_plot <- function(dat, tests, y_title) {
  if (!nrow(dat)) stop("No data available for: ", y_title)
  ymax_by_treatment <- dat %>%
    group_by(Treatment) %>%
    summarise(max_pct = if (all(is.na(percentage))) 0 else max(percentage, na.rm = TRUE), .groups = "drop")
  sig_pos <- tests %>%
    filter(!is.na(p_adj_BH)) %>%
    left_join(ymax_by_treatment, by = "Treatment") %>%
    mutate(
      x1 = match(group1, status_order), x2 = match(group2, status_order),
      y_position = max_pct + if_else(comparison_index == 1L, 5, 12),
      label_y = y_position + 1.5
    )
  candidates <- c(dat$percentage, sig_pos$label_y)
  candidates <- candidates[is.finite(candidates)]
  bar_y_max <- if (length(candidates)) min(125, max(25, max(candidates) + 4)) else 25

  ggplot(dat, aes(x = Status, y = percentage, fill = Status)) +
    geom_col(width = 0.7, na.rm = TRUE) +
    # Do not print n on this figure series: each feature has a different eligible denominator.
    geom_text(aes(label = ifelse(is.na(percentage), "", sprintf("%.1f%%", percentage))),
              vjust = -0.30, size = 3.5, na.rm = TRUE) +
    geom_segment(data = sig_pos, aes(x = x1, xend = x2, y = y_position, yend = y_position),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = sig_pos, aes(x = x1, xend = x1, y = y_position - 1.2, yend = y_position),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = sig_pos, aes(x = x2, xend = x2, y = y_position - 1.2, yend = y_position),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = sig_pos, aes(x = (x1 + x2) / 2, y = label_y, label = label),
              inherit.aes = FALSE, size = 4.0, na.rm = TRUE) +
    facet_wrap(~Treatment, nrow = 1,
               labeller = labeller(Treatment = facet_treatment_labels)) +
    scale_x_discrete(limits = status_order, drop = FALSE) +
    scale_fill_manual(values = status_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, bar_y_max),
                       breaks = seq(0, min(100, floor(bar_y_max / 10) * 10), by = 10),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(y = y_title, x = NULL) +
    base_theme
}

pB <- make_percentage_plot(
  intron_summary, intron_tests,
  "Genes with a 3′UTR intron (%)"
)
pC <- make_percentage_plot(
  uorf_summary, uorf_tests,
  "Genes with RiboTaper-supported uORF (%)"
)

combined_plot <- pA + pB + pC +
  plot_layout(ncol = 3, widths = c(1.15, 1, 1))
print(combined_plot)

# -------------------------------------------------------------------------
# 8. Export tables and figures
# -------------------------------------------------------------------------
write.csv(status_sizes, "NMD_cis_feature_status_gene_counts.csv", row.names = FALSE)
write.csv(testability_summary, "NMD_cis_feature_status_testability.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(status_by_treatment), "NMD_cis_feature_gene_status_by_treatment.csv", row.names = FALSE)
write.csv(gene_utr_length, "NMD_cis_feature_GFF_UTR_lengths_by_gene.csv", row.names = FALSE)
write.csv(gene_features_3utr, "NMD_cis_feature_GFF_UTR_intron_status_by_gene.csv", row.names = FALSE)
write.csv(uorf_df, "NMD_cis_feature_RiboTaper_uORF_status_by_gene.csv", row.names = FALSE)
write.csv(group_gene_table, "NMD_cis_feature_gene_groups_testable_expressed.csv", row.names = FALSE)
write.csv(annotation_coverage, "NMD_cis_feature_annotation_coverage.csv", row.names = FALSE)
write.csv(len_data, "NMD_cis_feature_panelA_UTR_length_data.csv", row.names = FALSE)
write.csv(rep_len_data, "NMD_cis_feature_panelA_representative_UTR_length_robustness_data.csv", row.names = FALSE)
write.csv(intron_summary, "NMD_cis_feature_panelB_intron_summary.csv", row.names = FALSE)
write.csv(intron_any_summary, "NMD_cis_feature_panelB_any_isoform_intron_robustness.csv", row.names = FALSE)
write.csv(uorf_summary, "NMD_cis_feature_panelC_uORF_summary.csv", row.names = FALSE)
write.csv(len_tests, "NMD_cis_feature_panelA_Wilcoxon_tests_BH.csv", row.names = FALSE)
write.csv(rep_len_tests, "NMD_cis_feature_panelA_representative_length_Wilcoxon_robustness_BH.csv", row.names = FALSE)
write.csv(intron_tests, "NMD_cis_feature_panelB_Fisher_tests_BH.csv", row.names = FALSE)
write.csv(intron_any_tests, "NMD_cis_feature_panelB_any_isoform_Fisher_robustness_BH.csv", row.names = FALSE)
write.csv(uorf_tests, "NMD_cis_feature_panelC_Fisher_tests_BH.csv", row.names = FALSE)

# Feature-denominator transparency: status-gene counts vs evaluable counts per panel.
panel_denominators <- bind_rows(
  len_data %>% count(Treatment, Status, name = "n_evaluable") %>% mutate(Feature = "Maximum 3′UTR length"),
  rep_len_data %>% count(Treatment, Status, name = "n_evaluable") %>% mutate(Feature = "Representative 3′UTR length"),
  intron_summary %>% transmute(Treatment, Status, n_evaluable = n_total, Feature = "Representative 3′UTR intron"),
  intron_any_summary %>% transmute(Treatment, Status, n_evaluable = n_total, Feature = "Any-isoform 3′UTR intron"),
  uorf_summary %>% transmute(Treatment, Status, n_evaluable = n_total, Feature = "RiboTaper uORF")
)
write.csv(panel_denominators, "NMD_cis_feature_panel_evaluable_denominators.csv", row.names = FALSE)

# Sensitivity/robustness table across external uORF definitions using the same assay-evaluable denominator.
robustness_rows <- bind_rows(lapply(names(uorf_external_sets), function(nm) {
  flag_col <- make.names(nm)
  temp_df <- uorf_df %>% transmute(gene = GENEID, flag = .data[[flag_col]])
  bind_rows(lapply(treatment_order, function(trt) {
    bind_rows(lapply(status_order, function(st) {
      g <- unique(groups[[trt]][[st]])
      x <- tibble(gene = g) %>% left_join(temp_df, by = "gene") %>% filter(!is.na(flag))
      tibble(
        definition = nm, Treatment = trt, Status = st,
        n_status_genes = length(g), n_evaluable = nrow(x),
        n_positive = sum(x$flag %in% TRUE),
        percentage = if (nrow(x)) 100 * mean(x$flag %in% TRUE) else NA_real_
      )
    }))
  }))
}))
write.csv(robustness_rows, "NMD_cis_feature_uORF_definition_robustness.csv", row.names = FALSE)

# Explicit audit of external calls and ID overlap with the DESeq2-filtered universe.
external_set_audit <- bind_rows(lapply(names(uorf_external_sets), function(nm) {
  s <- unique(clean_gene_id(uorf_external_sets[[nm]]))
  tibble(
    definition = nm,
    n_unique_external_genes = length(s),
    n_in_DESeq2_universe = length(intersect(s, all_genes)),
    n_in_Ribo_detectable_set = length(intersect(s, ribo_detectable_genes)),
    n_in_both_DESeq2_and_Ribo_detectable = length(intersect(intersect(s, all_genes), ribo_detectable_genes))
  )
}))
write.csv(external_set_audit, "NMD_cis_feature_external_uORF_set_audit.csv", row.names = FALSE)

# Combined figure for overview
ggsave("Figure_NMD_cis_features_DESeq2_GFF.png", combined_plot,
       width = 13.5, height = 4.7, dpi = 300, bg = "white")
ggsave("Figure_NMD_cis_features_DESeq2_GFF.pdf", combined_plot,
       width = 13.5, height = 4.7, bg = "white")

ggsave("Figure_A_UTR_length_DESeq2_GFF.png", pA,
       width = individual_plot_width_in, height = individual_plot_height_in,
       dpi = 300, bg = "white")
ggsave("Figure_A_UTR_length_DESeq2_GFF.pdf", pA,
       width = individual_plot_width_in, height = individual_plot_height_in,
       bg = "white")
ggsave("Figure_B_3UTR_intron_DESeq2_GFF.png", pB,
       width = individual_plot_width_in, height = individual_plot_height_in,
       dpi = 300, bg = "white")
ggsave("Figure_B_3UTR_intron_DESeq2_GFF.pdf", pB,
       width = individual_plot_width_in, height = individual_plot_height_in,
       bg = "white")
ggsave("Figure_C_uORF_RiboTaper_no_UTR_filter.png", pC,
       width = individual_plot_width_in, height = individual_plot_height_in,
       dpi = 300, bg = "white")
ggsave("Figure_C_uORF_RiboTaper_no_UTR_filter.pdf", pC,
       width = individual_plot_width_in, height = individual_plot_height_in,
       bg = "white")

message("Completed NMD cis-feature analysis.")
message("Primary uORF definition: ", main_uorf_definition)
message("Panel C denominator: genes with positive RiboTaper evidence or transcript-level TPM >= ", ribo_tpm_min, ".")
message("DE status NA results are reported as Not testable, not as Unchanged.")
message("BH-adjusted p-values are used for the displayed significance labels.")
