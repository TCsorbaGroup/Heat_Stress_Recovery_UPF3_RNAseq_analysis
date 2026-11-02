###############################################################################
# GO Biological Process enrichment for HSR-responsive genes elevated in upf3
# at NT.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(AnnotationDbi)
  library(stringr)
  library(rrvgo)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
alpha_cutoff <- 0.05
lfc_cutoff <- 1
fdr_cutoff <- 0.05
rrvgo_threshold <- 0.7
max_terms_to_plot <- 10
min_go_set_size <- 10
output_prefix <- "HSR_upf3_high_in_NT_GO_BP_full_Arabidopsis"

# -----------------------------------------------------------------------------
# 2. Helpers and input validation
# -----------------------------------------------------------------------------
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

clean_gene_set <- function(x) {
  x <- clean_gene_id(x)
  unique(x[!is.na(x) & nzchar(x)])
}

normalize_genotype <- function(x) {
  x <- tolower(gsub("[[:space:]]+", "", trimws(as.character(x))))
  dplyr::case_when(
    x %in% c("wt", "col", "col-0", "col0", "col_0", "wildtype", "wild_type", "wild-type") ~ "WT",
    x %in% c("upf3", "upf-3", "upf_3") ~ "upf3",
    TRUE ~ NA_character_
  )
}

normalize_condition <- function(x) {
  x <- toupper(gsub("[[:space:]]+", "", trimws(as.character(x))))
  dplyr::case_when(
    x %in% c("NT", "NOTREATMENT", "UNTREATED") ~ "NT",
    x %in% c("LAT", "P/T", "PT") ~ "LAT",
    x %in% c("REC", "P/T+1D", "PT+1D") ~ "REC",
    TRUE ~ NA_character_
  )
}

results_to_table <- function(res, label = "DESeq2 result") {
  x <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")
  x$gene_id <- clean_gene_id(x$gene_id)
  if (anyNA(x$gene_id) || any(!nzchar(x$gene_id))) {
    stop("Missing gene IDs in ", label, " after ID cleanup.")
  }
  if (anyDuplicated(x$gene_id)) {
    dups <- unique(x$gene_id[duplicated(x$gene_id)])
    stop("Duplicate gene IDs after cleanup in ", label, ": ",
         paste(head(dups, 10), collapse = ", "))
  }
  if (!all(c("log2FoldChange", "pvalue", "padj") %in% names(x))) {
    stop("The ", label, " table lacks log2FoldChange, pvalue, and/or padj. Available: ",
         paste(names(x), collapse = ", "))
  }
  tibble(
    gene_id = x$gene_id,
    log2FC = as.numeric(x$log2FoldChange),
    pvalue = as.numeric(x$pvalue),
    padj = as.numeric(x$padj)
  )
}

fit_wt_hsr_results <- function(raw_counts, sample_info) {
  wt_meta <- sample_info %>%
    filter(genotype == "WT") %>%
    mutate(condition = factor(condition, levels = c("NT", "LAT", "REC"))) %>%
    as.data.frame()
  if (nrow(wt_meta) < 2L || any(table(wt_meta$condition) < 2L)) {
    stop("The WT-only model needs at least two Col-0 samples in each of NT, LAT, and REC.")
  }
  rownames(wt_meta) <- wt_meta$sample
  col_data <- data.frame(
    treatment = factor(wt_meta$condition, levels = c("NT", "LAT", "REC")),
    row.names = wt_meta$sample
  )
  dds_wt <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw_counts[, wt_meta$sample, drop = FALSE]),
    colData = col_data,
    design = ~ treatment
  )
  dds_wt <- DESeq2::DESeq(dds_wt, quiet = TRUE)
  list(
    lat = DESeq2::results(dds_wt, contrast = c("treatment", "LAT", "NT"), alpha = alpha_cutoff),
    rec = DESeq2::results(dds_wt, contrast = c("treatment", "REC", "NT"), alpha = alpha_cutoff)
  )
}

fit_nt_genotype_result <- function(raw_counts, sample_info) {
  nt_meta <- sample_info %>%
    filter(condition == "NT") %>%
    mutate(genotype = factor(genotype, levels = c("WT", "upf3"))) %>%
    as.data.frame()
  if (nrow(nt_meta) < 4L || any(table(nt_meta$genotype) < 2L)) {
    stop("The NT-only model needs at least two NT samples for each genotype.")
  }
  rownames(nt_meta) <- nt_meta$sample
  col_data <- data.frame(
    genotype = factor(nt_meta$genotype, levels = c("WT", "upf3")),
    row.names = nt_meta$sample
  )
  dds_nt <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(raw_counts[, nt_meta$sample, drop = FALSE]),
    colData = col_data,
    design = ~ genotype
  )
  dds_nt <- DESeq2::DESeq(dds_nt, quiet = TRUE)
  DESeq2::results(dds_nt, contrast = c("genotype", "upf3", "WT"), alpha = alpha_cutoff)
}

if (!exists("dds", inherits = TRUE) || !inherits(get("dds", inherits = TRUE), "DESeqDataSet")) {
  stop("A DESeq2 object named 'dds' is required. Run the DESeq2 analysis first.")
}
dds <- get("dds", inherits = TRUE)
meta <- as.data.frame(SummarizedExperiment::colData(dds))
if (!all(c("genotype", "treatment") %in% names(meta))) {
  stop("colData(dds) must contain columns named 'genotype' and 'treatment'.")
}

raw_counts <- DESeq2::counts(dds, normalized = FALSE)
if (is.null(rownames(raw_counts)) || is.null(colnames(raw_counts))) {
  stop("The raw-count matrix from dds must have gene and sample names.")
}
if (any(!is.finite(raw_counts)) || any(raw_counts < 0)) {
  stop("Raw counts in dds contain non-finite or negative values.")
}
if (is.null(rownames(meta)) || any(!nzchar(rownames(meta))) ||
    !setequal(rownames(meta), colnames(raw_counts))) {
  stop("Sample names in colData(dds) do not match the raw-count columns.")
}
meta <- meta[colnames(raw_counts), , drop = FALSE]

clean_ids <- clean_gene_id(rownames(raw_counts))
if (anyNA(clean_ids) || any(!nzchar(clean_ids))) stop("Some dds gene IDs are empty after cleanup.")
if (anyDuplicated(clean_ids)) {
  dups <- unique(clean_ids[duplicated(clean_ids)])
  stop("Gene IDs are duplicated after cleanup: ", paste(head(dups, 10), collapse = ", "))
}
rownames(raw_counts) <- clean_ids

sample_info <- tibble(
  sample = colnames(raw_counts),
  genotype = normalize_genotype(meta$genotype),
  condition = normalize_condition(meta$treatment)
)
if (anyNA(sample_info$genotype) || anyNA(sample_info$condition)) {
  bad <- sample_info %>% filter(is.na(genotype) | is.na(condition)) %>% pull(sample)
  stop("Unrecognized genotype/treatment labels: ", paste(bad, collapse = ", "))
}
sample_info <- sample_info %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "upf3")),
    condition = factor(condition, levels = c("NT", "LAT", "REC"))
  )

sample_distribution <- sample_info %>% count(genotype, condition, .drop = FALSE)
message("Sample distribution (expected three replicates per cell in this session):")
print(as.data.frame(sample_distribution), row.names = FALSE)

needed_combos <- tidyr::expand_grid(
  genotype = factor(c("WT", "upf3"), levels = c("WT", "upf3")),
  condition = factor(c("NT", "LAT", "REC"), levels = c("NT", "LAT", "REC"))
)
missing_combos <- needed_combos %>%
  anti_join(sample_info %>% distinct(genotype, condition), by = c("genotype", "condition"))
if (nrow(missing_combos) > 0L) {
  stop("Missing genotype x condition combinations: ",
       paste(paste(missing_combos$genotype, missing_combos$condition, sep = "/"), collapse = ", "))
}
if (any(table(sample_info$genotype)[c("WT", "upf3")] < 2L)) {
  stop("At least two samples per genotype are required.")
}
if (any(table(sample_info$condition)[c("NT", "LAT", "REC")] < 2L)) {
  stop("At least two samples per treatment are required.")
}

# -----------------------------------------------------------------------------
# 3. Define HSR-responsive genes 
# -----------------------------------------------------------------------------
design_terms <- attr(stats::terms(DESeq2::design(dds)), "term.labels")
has_interaction <- any(grepl(":", design_terms, fixed = TRUE) &
                         grepl("genotype", design_terms, fixed = TRUE) &
                         grepl("treatment", design_terms, fixed = TRUE))
genotype_ref <- if (is.factor(meta$genotype)) levels(meta$genotype)[1] else NA_character_
treatment_ref <- if (is.factor(meta$treatment)) levels(meta$treatment)[1] else NA_character_
genotype_ref_norm <- normalize_genotype(genotype_ref)
treatment_ref_norm <- normalize_condition(treatment_ref)
coef_names <- DESeq2::resultsNames(dds)
required_hsr_coef <- c("treatment_LAT_vs_NT", "treatment_REC_vs_NT")
can_use_full_for_hsr <- has_interaction &&
  identical(genotype_ref_norm, "WT") && identical(treatment_ref_norm, "NT") &&
  all(required_hsr_coef %in% coef_names)

hsr_res_list <- NULL
if (can_use_full_for_hsr) {
  hsr_res_list <- tryCatch(
    list(
      lat = DESeq2::results(dds, name = "treatment_LAT_vs_NT", alpha = alpha_cutoff),
      rec = DESeq2::results(dds, name = "treatment_REC_vs_NT", alpha = alpha_cutoff)
    ),
    error = function(e) {
      message("Could not extract WT treatment coefficients; refitting WT-only model: ", conditionMessage(e))
      NULL
    }
  )
}
if (is.null(hsr_res_list)) {
  message("Defining WT HSR genes with a WT-only DESeq2 model: LAT/REC vs NT.")
  hsr_res_list <- fit_wt_hsr_results(raw_counts, sample_info)
  hsr_method <- "WT-only DESeq2 model, LAT/REC vs NT"
} else {
  message("Using dds treatment coefficients at WT/Col-0 reference and NT baseline.")
  hsr_method <- "existing dds treatment coefficients at WT/Col-0 reference"
}

lat_res <- results_to_table(hsr_res_list$lat, "WT LAT vs NT")
names(lat_res)[names(lat_res) == "log2FC"] <- "log2FC_LAT"
names(lat_res)[names(lat_res) == "pvalue"] <- "pvalue_LAT"
names(lat_res)[names(lat_res) == "padj"] <- "padj_LAT"
rec_res <- results_to_table(hsr_res_list$rec, "WT REC vs NT")
names(rec_res)[names(rec_res) == "log2FC"] <- "log2FC_REC"
names(rec_res)[names(rec_res) == "pvalue"] <- "pvalue_REC"
names(rec_res)[names(rec_res) == "padj"] <- "padj_REC"

hsr_classification <- full_join(lat_res, rec_res, by = "gene_id") %>%
  mutate(
    induced_LAT = !is.na(padj_LAT) & padj_LAT < alpha_cutoff &
      !is.na(log2FC_LAT) & log2FC_LAT > lfc_cutoff,
    induced_REC = !is.na(padj_REC) & padj_REC < alpha_cutoff &
      !is.na(log2FC_REC) & log2FC_REC > lfc_cutoff,
    HSR_responsive = induced_LAT | induced_REC,
    HSR_class = case_when(
      induced_LAT & induced_REC ~ "Sustained",
      induced_LAT & !induced_REC ~ "Acute/transient",
      !induced_LAT & induced_REC ~ "Delayed",
      TRUE ~ "Not HSR-responsive"
    )
  )

hsr_genes <- clean_gene_set(hsr_classification$gene_id[hsr_classification$HSR_responsive])
if (length(hsr_genes) == 0L) {
  stop("No WT/Col-0 HSR-responsive genes passed padj < 0.05 and log2FC > 1 at LAT or REC vs NT.")
}
message("HSR definition method: ", hsr_method)
message("Total WT/Col-0 HSR-responsive genes (LAT OR REC): ", length(hsr_genes))

# -----------------------------------------------------------------------------
# 4. Select HSR genes significantly >=2-fold higher in upf3 at NT
# -----------------------------------------------------------------------------
raw_genotype_values <- unique(trimws(as.character(meta$genotype)))
raw_wt_values <- raw_genotype_values[normalize_genotype(raw_genotype_values) == "WT"]
raw_upf3_values <- raw_genotype_values[normalize_genotype(raw_genotype_values) == "upf3"]
if (length(raw_wt_values) != 1L || length(raw_upf3_values) != 1L) {
  stop("Could not map exactly one raw genotype label to WT and upf3 in colData(dds).")
}
can_use_full_for_nt <- has_interaction &&
  identical(genotype_ref_norm, "WT") && identical(treatment_ref_norm, "NT")
nt_res_raw <- NULL
if (can_use_full_for_nt) {
  nt_res_raw <- tryCatch(
    DESeq2::results(
      dds,
      contrast = c("genotype", raw_upf3_values, raw_wt_values),
      alpha = alpha_cutoff
    ),
    error = function(e) {
      message("Could not extract NT genotype contrast; refitting NT-only model: ", conditionMessage(e))
      NULL
    }
  )
}
if (is.null(nt_res_raw)) {
  message("Estimating the NT-specific upf3-vs-WT contrast from NT samples only.")
  nt_res_raw <- fit_nt_genotype_result(raw_counts, sample_info)
  nt_method <- "NT-only DESeq2 model, upf3 vs WT"
} else {
  nt_method <- "existing dds genotype contrast at WT/Col-0 and NT reference"
}

nt_de <- results_to_table(nt_res_raw, "NT upf3 vs WT")
names(nt_de)[names(nt_de) == "log2FC"] <- "log2FC_upf3_vs_WT"
names(nt_de)[names(nt_de) == "pvalue"] <- "pvalue_NT"
names(nt_de)[names(nt_de) == "padj"] <- "padj_NT"

nt_upf3_high <- nt_de %>%
  filter(
    !is.na(padj_NT), padj_NT < alpha_cutoff,
    !is.na(log2FC_upf3_vs_WT), log2FC_upf3_vs_WT >= lfc_cutoff
  ) %>%
  pull(gene_id) %>% unique()
candidate_genes <- intersect(hsr_genes, nt_upf3_high)

message("NT contrast method: ", nt_method)
message("Genes >=2-fold and significant in upf3 vs WT at NT: ", length(nt_upf3_high))
message("Intersection with total WT HSR-responsive genes: ", length(candidate_genes))

# Save complete DESeq2 tables and explicit gene lists with fold-change statistics.
write.csv(hsr_classification,
          paste0(output_prefix, "_all_gene_HSR_classification.csv"), row.names = FALSE)
write.csv(nt_de,
          paste0(output_prefix, "_NT_upf3_vs_WT_DESeq2_all_genes.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = sort(hsr_genes)),
          paste0(output_prefix, "_all_HSR_responsive_gene_ids.csv"), row.names = FALSE)

# All HSR-responsive genes, including their NT upf3-vs-WT effect size and test.
hsr_genes_with_changes <- hsr_classification %>%
  filter(HSR_responsive) %>%
  left_join(nt_de, by = "gene_id") %>%
  mutate(
    NT_upf3_high_significant_2fold = !is.na(padj_NT) & padj_NT < alpha_cutoff &
      !is.na(log2FC_upf3_vs_WT) & log2FC_upf3_vs_WT >= lfc_cutoff,
    selected_for_GO_candidate_set = gene_id %in% candidate_genes
  ) %>%
  arrange(padj_NT)
write.csv(
  hsr_genes_with_changes,
  paste0(output_prefix, "_all_HSR_genes_with_NT_change_statistics.csv"),
  row.names = FALSE
)

# All genes meeting the NT upf3-vs-WT statistical/effect-size cutoffs,
# whether or not they belong to the WT-defined HSR-responsive set.
nt_upf3_high_table <- nt_de %>%
  filter(gene_id %in% nt_upf3_high) %>%
  left_join(
    hsr_classification %>% dplyr::select(gene_id, HSR_responsive, HSR_class,
      log2FC_LAT, pvalue_LAT, padj_LAT, log2FC_REC, pvalue_REC, padj_REC),
    by = "gene_id"
  ) %>%
  mutate(in_WT_HSR_responsive_set = gene_id %in% hsr_genes) %>%
  arrange(padj_NT)
write.csv(
  nt_upf3_high_table,
  paste0(output_prefix, "_all_NT_upf3_high_genes_with_change_statistics.csv"),
  row.names = FALSE
)

# Exact foreground before excluding genes that lack GO BP annotation.
candidate_table <- hsr_classification %>%
  filter(gene_id %in% candidate_genes) %>%
  left_join(nt_de, by = "gene_id") %>%
  arrange(padj_NT)
write.csv(
  candidate_table,
  paste0(output_prefix, "_selected_HSR_genes_upf3_high_at_NT.csv"),
  row.names = FALSE
)

if (length(candidate_genes) == 0L) {
  stop("No genes met both definitions: WT HSR-responsive and NT upf3 vs WT log2FC >= 1 with BH padj < 0.05. Gene tables were saved.")
}

# -----------------------------------------------------------------------------
# 5. Full Arabidopsis TAIR gene universe, restricted only to GO BP annotation
# -----------------------------------------------------------------------------
valid_tair_ids <- clean_gene_set(AnnotationDbi::keys(org.At.tair.db, keytype = "TAIR"))
if (length(valid_tair_ids) == 0L) stop("No TAIR identifiers were retrieved from org.At.tair.db.")
message("Genome-wide TAIR identifiers in org.At.tair.db: ", length(valid_tair_ids))

go_annot <- suppressMessages(suppressWarnings(
  AnnotationDbi::select(
    org.At.tair.db,
    keys = valid_tair_ids,
    keytype = "TAIR",
    columns = c("GO", "ONTOLOGY")
  )
))

bp_annotated_genes <- unique(clean_gene_set(
  go_annot$TAIR[
    !is.na(go_annot$GO) & !is.na(go_annot$ONTOLOGY) & go_annot$ONTOLOGY == "BP"
  ]
))
go_background <- intersect(valid_tair_ids, bp_annotated_genes)
if (length(go_background) == 0L) {
  stop("No genome-wide TAIR genes with GO BP annotations were found in org.At.tair.db.")
}
go_foreground <- intersect(candidate_genes, go_background)
if (length(go_foreground) == 0L) {
  stop("None of the selected HSR/upf3-high genes map to the TAIR GO BP background. Gene tables were saved.")
}

n_not_valid_tair <- length(setdiff(candidate_genes, valid_tair_ids))
n_no_bp_annotation <- length(setdiff(intersect(candidate_genes, valid_tair_ids), go_background))
summary_table <- tibble(
  item = c(
    "WT/Col-0 HSR-responsive genes (LAT OR REC)",
    "Genes significantly >=2-fold higher in upf3 than WT at NT",
    "Selected genes intersecting HSR-responsive set",
    "Selected genes with valid TAIR IDs",
    "Selected genes without GO BP annotation",
    "Selected genes used as GO BP foreground",
    "All TAIR IDs in org.At.tair.db",
    "Full Arabidopsis GO BP-annotated background"
  ),
  n_genes = c(
    length(hsr_genes), length(nt_upf3_high), length(candidate_genes),
    length(intersect(candidate_genes, valid_tair_ids)), n_no_bp_annotation,
    length(go_foreground), length(valid_tair_ids), length(go_background)
  )
)
message("")
message("Foreground/background summary:")
print(as.data.frame(summary_table), row.names = FALSE)
write.csv(summary_table, paste0(output_prefix, "_gene_counts_and_GO_universe.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = sort(valid_tair_ids)), paste0(output_prefix, "_full_TAIR_gene_universe.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = sort(go_background)), paste0(output_prefix, "_GO_BP_background.csv"), row.names = FALSE)
write.csv(data.frame(gene_id = sort(go_foreground)), paste0(output_prefix, "_selected_genes_with_GO_BP_annotation.csv"), row.names = FALSE)

go_foreground_with_changes <- candidate_table %>%
  filter(gene_id %in% go_foreground) %>%
  mutate(used_in_GO_BP_test = TRUE) %>%
  arrange(padj_NT)
write.csv(
  go_foreground_with_changes,
  paste0(output_prefix, "_GO_BP_foreground_gene_list_with_change_statistics.csv"),
  row.names = FALSE
)

candidate_without_bp <- candidate_table %>%
  filter(gene_id %in% setdiff(candidate_genes, go_foreground)) %>%
  mutate(used_in_GO_BP_test = FALSE,
         exclusion_reason = ifelse(gene_id %in% valid_tair_ids,
                                   "No GO BP annotation in org.At.tair.db",
                                   "Gene ID not found in org.At.tair.db TAIR universe")) %>%
  arrange(padj_NT)
write.csv(
  candidate_without_bp,
  paste0(output_prefix, "_selected_candidates_excluded_from_GO_BP_with_change_statistics.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. GO BP enrichment
# -----------------------------------------------------------------------------
message("")
message("Running GO Biological Process over-representation test...")
message("Foreground: selected WT HSR-responsive genes elevated >=2-fold in upf3 at NT.")
message("Background: all TAIR genes with GO BP annotation in org.At.tair.db.")

ego <- clusterProfiler::enrichGO(
  gene = go_foreground,
  universe = go_background,
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = min_go_set_size,
  readable = FALSE
)
go_results <- as.data.frame(ego)
if (nrow(go_results) > 0L) go_results <- go_results %>% arrange(p.adjust, pvalue)
write.csv(go_results, paste0(output_prefix, "_GO_BP_enrichment_all_terms.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. rrvgo semantic reduction (Wang similarity, threshold 0.7) and dotplot
# -----------------------------------------------------------------------------
format_unreduced_terms <- function(terms) {
  if (nrow(terms) == 0L) return(tibble())
  terms %>% transmute(
    GO_ID = ID,
    Description,
    p.adjust,
    Count,
    GeneRatio,
    n_terms_collapsed = 1L,
    log10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
  )
}

reduce_go_terms <- function(terms) {
  terms <- terms %>% distinct(ID, .keep_all = TRUE)
  if (nrow(terms) <= 1L) return(format_unreduced_terms(terms))

  sim_matrix <- tryCatch(
    rrvgo::calculateSimMatrix(
      unique(terms$ID),
      orgdb = "org.At.tair.db",
      keytype = "TAIR",
      ont = "BP",
      method = "Wang"
    ),
    error = function(e) {
      message("rrvgo similarity calculation failed; retaining all significant terms. Reason: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(sim_matrix) || nrow(sim_matrix) == 0L) return(format_unreduced_terms(terms))

  terms_with_similarity <- terms %>% filter(ID %in% rownames(sim_matrix))
  terms_without_similarity <- terms %>% filter(!ID %in% rownames(sim_matrix))
  if (nrow(terms_with_similarity) <= 1L) {
    return(bind_rows(format_unreduced_terms(terms_with_similarity), format_unreduced_terms(terms_without_similarity)))
  }

  sim_matrix <- sim_matrix[terms_with_similarity$ID, terms_with_similarity$ID, drop = FALSE]
  scores <- setNames(
    -log10(pmax(terms_with_similarity$p.adjust, .Machine$double.xmin)),
    terms_with_similarity$ID
  )
  reduced <- tryCatch(
    as.data.frame(rrvgo::reduceSimMatrix(
      sim_matrix,
      scores = scores,
      threshold = rrvgo_threshold,
      orgdb = "org.At.tair.db",
      keytype = "TAIR"
    )),
    error = function(e) {
      message("rrvgo reduction failed; retaining all significant terms. Reason: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(reduced) || nrow(reduced) == 0L || !all(c("go", "parent") %in% names(reduced))) {
    return(format_unreduced_terms(terms))
  }

  collapsed_counts <- reduced %>%
    group_by(parent) %>%
    summarise(n_terms_collapsed = n_distinct(go), .groups = "drop")
  representatives <- reduced %>%
    filter(go == parent) %>%
    distinct(parent, .keep_all = TRUE) %>%
    dplyr::select(parent)
  representative_terms <- representatives %>%
    left_join(collapsed_counts, by = "parent") %>%
    left_join(terms_with_similarity, by = c("parent" = "ID")) %>%
    transmute(
      GO_ID = parent,
      Description,
      p.adjust,
      Count,
      GeneRatio,
      n_terms_collapsed = as.integer(n_terms_collapsed),
      log10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
    )

  bind_rows(representative_terms, format_unreduced_terms(terms_without_similarity)) %>%
    distinct(GO_ID, .keep_all = TRUE)
}

if (nrow(go_results) == 0L) {
  message("No GO BP terms were returned; full enrichment table has been saved and no dotplot will be produced.")
  write.csv(tibble(), paste0(output_prefix, "_GO_BP_reduced_terms.csv"), row.names = FALSE)
} else {
  sig_go_results <- go_results %>% filter(!is.na(p.adjust), p.adjust < fdr_cutoff)
  message("GO BP terms passing BH-adjusted p < ", fdr_cutoff, ": ", nrow(sig_go_results))

  if (nrow(sig_go_results) == 0L) {
    message("No GO BP terms passed the FDR cutoff; no significance-only dotplot will be generated.")
    write.csv(tibble(), paste0(output_prefix, "_GO_BP_reduced_terms.csv"), row.names = FALSE)
  } else {
    reduced_go <- reduce_go_terms(sig_go_results) %>% arrange(desc(log10FDR), p.adjust)
    write.csv(reduced_go, paste0(output_prefix, "_GO_BP_reduced_terms.csv"), row.names = FALSE)

    plot_data <- reduced_go %>%
      filter(is.finite(log10FDR), !is.na(Description)) %>%
      arrange(desc(log10FDR)) %>%
      slice_head(n = max_terms_to_plot)

    if (nrow(plot_data) == 0L) {
      message("No representative terms remained for plotting; see reduced GO table.")
    } else {
      plot_data <- plot_data %>%
        mutate(TermKey = paste0(Description, " [GO:", GO_ID, "]"))
      plot_data$TermKey <- factor(plot_data$TermKey, levels = rev(unique(plot_data$TermKey)))

      p_go <- ggplot(plot_data, aes(x = log10FDR, y = TermKey)) +
        geom_point(aes(size = Count, colour = log10FDR), alpha = 0.95) +
        scale_colour_gradient(
          low = "#3B0066", high = "#E60026",
          name = expression(-log[10](FDR))
        ) +
        scale_size_continuous(range = c(2.5, 8), breaks = pretty(plot_data$Count, n = 4), name = "Gene count") +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        scale_y_discrete(labels = function(x) stringr::str_wrap(sub(" \\[GO:.*$", "", x), width = 40)) +
        labs(
          title = "GO BP enrichment of HS-responsive genes elevated in upf3 at NT",
          subtitle = paste0(
            ""
          ),
          x = expression(-log[10](FDR)),
          y = NULL,
          caption = paste0(
            "HSR: WT/Col-0 padj < ", alpha_cutoff,
            " and log2FC > ", lfc_cutoff,
            " at LAT or REC vs NT. Foreground: NT upf3 vs WT padj < ", alpha_cutoff,
            " and log2FC >= ", lfc_cutoff, ". Background: all TAIR genes annotated to GO BP."
          )
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title = element_text(size = 14, hjust = 0, face = "bold"),
          plot.subtitle = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 14),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          plot.caption = element_blank(),
          plot.margin = margin(8, 10, 8, 8)
        )

      print(p_go)
      plot_height <- max(5.0, 2.0 + 0.42 * nrow(plot_data))
      ggsave(
        filename = paste0(output_prefix, "_dotplot.png"), plot = p_go,
        width = 8.5, height = plot_height, units = "in", dpi = 300, bg = "white"
      )
      ggsave(
        filename = paste0(output_prefix, "_dotplot.pdf"), plot = p_go,
        width = 8.5, height = plot_height, units = "in", bg = "white"
      )
      write.csv(plot_data, paste0(output_prefix, "_GO_BP_dotplot_terms.csv"), row.names = FALSE)
      message("Saved PNG dotplot: ", paste0(output_prefix, "_dotplot.png"))
      message("Saved PDF dotplot: ", paste0(output_prefix, "_dotplot.pdf"))
    }
  }
}

# -----------------------------------------------------------------------------
# 8. Final report
# -----------------------------------------------------------------------------
message("")
message("HSR upf3-high-at-NT GO BP analysis completed.")
message("HSR definition method: ", hsr_method)
message("NT genotype contrast method: ", nt_method)
message("WT/Col-0 HSR-responsive genes: ", length(hsr_genes))
message("Selected HSR genes >=2-fold higher in upf3 at NT: ", length(candidate_genes))
message("Selected genes with GO BP annotation: ", length(go_foreground))
message("Saved all HSR genes with NT change statistics: ", paste0(output_prefix, "_all_HSR_genes_with_NT_change_statistics.csv"))
message("Saved all NT upf3-high genes with change statistics: ", paste0(output_prefix, "_all_NT_upf3_high_genes_with_change_statistics.csv"))
message("Saved exact GO BP foreground with change statistics: ", paste0(output_prefix, "_GO_BP_foreground_gene_list_with_change_statistics.csv"))
message("Saved candidates excluded from GO BP test: ", paste0(output_prefix, "_selected_candidates_excluded_from_GO_BP_with_change_statistics.csv"))
message("Full TAIR universe size: ", length(valid_tair_ids))
message("Effective full-genome GO BP background size: ", length(go_background))
message("HSR cutoff: padj < ", alpha_cutoff, " and log2FC > ", lfc_cutoff)
message("NT cutoff: padj < ", alpha_cutoff, " and log2FC >= ", lfc_cutoff)
message("GO cutoff: BH-adjusted p < ", fdr_cutoff, "; rrvgo Wang threshold = ", rrvgo_threshold)
message("Output prefix: ", output_prefix)
