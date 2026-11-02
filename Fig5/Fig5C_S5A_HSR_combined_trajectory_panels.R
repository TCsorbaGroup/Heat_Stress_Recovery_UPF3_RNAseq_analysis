###############################################################################
# Self-contained combined HSR trajectory figure: Transient / Sustained / Delayed
#
# HSR class definition (same cutoffs as the supplied Col-0 HSR script):
#   Induced in WT/Col-0 if BH padj < 0.05 and log2FC > 1 versus NT.
#   Transient (acute): LAT-induced, not REC-induced
#   Sustained:         LAT-induced and REC-induced
#   Delayed:           REC-induced, not LAT-induced
#
# Correlations:
#   Pearson correlation across the three WT expression means (NT, LAT, REC)
#   for each HSR gene against the matching WT HSFA2 or HSFA3 expression profile.
#   r > 0.7 -> included in that TF-correlated panel.
#   If r_HSFA2 > 0.7 and r_HSFA3 > 0.7, a gene appears in BOTH panels.
#   If both correlations are <= 0.7, it appears in the Neither panel.
#   Genes with undefined correlations cannot be assigned and are reported.
#
# Trajectories:
#   Expression is log2(FPKM + 1) if fpkm_wide exists and expression_source='auto';
#   otherwise log2(DESeq2 normalized count + 1) from dds.
#   Each gene is z-scored across all six genotype x condition means.
#   Grey lines = individual WT trajectories; thick lines = genotype means;
#   dashed blue line = scaled WT reference trajectory for the matching TF.
# Outputs: one combined PNG + PDF, plus CSV tables for correlations, membership,
#          panel counts and genes that could not be assigned to a correlation group.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

# ------------------------------- Settings ------------------------------------
expression_source <- "auto"  # "auto" or "dds_norm"
alpha_cutoff <- 0.05
lfc_cutoff <- 1
correlation_cutoff <- 0.7
condition_levels <- c("NT", "LAT", "REC")
genotype_levels <- c("WT", "upf3")
class_names <- c("Transient (acute)", "Sustained", "Delayed")
correlation_group_levels <- c("HSFA2-correlated", "HSFA3-correlated", "Neither")

# Preserve the two exclusions of HSFA2 and HSFA2's adjacent gene. HSFA3 itself is also
# excluded from the target set because its correlation with its own profile is
# trivially 1 and would not be an informative target-gene correlation.
target_exclusions <- c("AT2G26150", "AT2G26140", "AT5G03720")

output_dir <- getwd()
output_prefix <- "HSR_all_genes_HSFA2_HSFA3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# --------------------------- Package checks ----------------------------------
required_pkgs <- c("ggplot2", "dplyr", "tidyr", "tibble", "patchwork")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("Install required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

# --------------------------- Helper functions --------------------------------
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

normalise_label <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", trimws(as.character(x))))
}

normalise_genotype <- function(x) {
  z <- tolower(gsub("[[:space:]_]", "", trimws(as.character(x))))
  dplyr::case_when(
    z %in% c("col", "col0", "col-0", "wt", "wildtype") ~ "WT",
    z %in% c("upf3", "upf-3") ~ "upf3",
    TRUE ~ NA_character_
  )
}

normalise_condition <- function(x) {
  z <- toupper(gsub("[^[:alnum:]]+", "", trimws(as.character(x))))
  dplyr::case_when(
    z %in% c("NT", "CONTROL", "NOTREATMENT", "UNTREATED", "NOHS", "NOHEAT") ~ "NT",
    z %in% c("LAT", "PT", "HEAT", "HS", "ACUTE") ~ "LAT",
    z %in% c("REC", "RECOVERY", "RECOVERED", "PT1D", "PTPLUS1D") ~ "REC",
    TRUE ~ NA_character_
  )
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || length(y) < 3L ||
      !is.finite(stats::sd(x)) || !is.finite(stats::sd(y)) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = "pearson"))
}

get_symbols <- function(ids) {
  ids <- unique(clean_gene_set(ids))
  out <- tibble(gene_id = ids, SYMBOL = NA_character_)

  if (requireNamespace("org.At.tair.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    ann <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        org.At.tair.db::org.At.tair.db,
        keys = ids,
        keytype = "TAIR",
        columns = "SYMBOL"
      )),
      error = function(e) NULL
    )

    if (!is.null(ann) && nrow(ann) > 0L && all(c("TAIR", "SYMBOL") %in% names(ann))) {
      ann <- ann %>%
        mutate(TAIR = clean_gene_id(TAIR)) %>%
        group_by(TAIR) %>%
        summarise(
          SYMBOL = paste(unique(SYMBOL[!is.na(SYMBOL) & nzchar(SYMBOL)]), collapse = " / "),
          .groups = "drop"
        )
      out <- out %>%
        left_join(ann, by = c("gene_id" = "TAIR"), suffix = c("", ".ann")) %>%
        mutate(SYMBOL = ifelse(!is.na(SYMBOL.ann) & nzchar(SYMBOL.ann), SYMBOL.ann, SYMBOL)) %>%
        select(gene_id, SYMBOL)
    }
  }

  out %>% mutate(label = ifelse(is.na(SYMBOL) | SYMBOL == "", gene_id, SYMBOL))
}

# ------------------------- 1. Recover HSR classes -----------------------------
hsr_classes <- NULL
hsr_source <- NULL

if (is.null(hsr_classes) && exists("gene_classification", inherits = TRUE)) {
  gc <- as.data.frame(get("gene_classification", inherits = TRUE))
  gene_col <- intersect(tolower(names(gc)), c("gene", "gene_id", "geneid", "tair"))
  if (nrow(gc) > 0L && length(gene_col) > 0L &&
      all(c("is_acute", "is_sustained", "is_delayed") %in% names(gc))) {
    gi <- match(gene_col[1], tolower(names(gc)))
    hsr_classes <- list(
      "Transient (acute)" = clean_gene_set(gc[[gi]][as_flag(gc$is_acute)]),
      "Sustained" = clean_gene_set(gc[[gi]][as_flag(gc$is_sustained)]),
      "Delayed" = clean_gene_set(gc[[gi]][as_flag(gc$is_delayed)])
    )
    hsr_source <- "existing gene_classification table"
  }
}

if (is.null(hsr_classes) && exists("gene_sets", inherits = TRUE)) {
  gs <- as.data.frame(get("gene_sets", inherits = TRUE))
  lower_names <- tolower(names(gs))
  gene_idx <- match(TRUE, lower_names %in% c("gene", "gene_id", "geneid", "tair"), nomatch = 0)
  class_idx <- match(TRUE, lower_names %in% c("genegroup", "gene_group", "group", "class"), nomatch = 0)

  if (gene_idx > 0L && class_idx > 0L) {
    raw_label <- trimws(as.character(gs[[class_idx]]))
    norm_label <- normalise_label(raw_label)
    acute_labels <- unique(raw_label[grepl("acute|transient", norm_label)])
    sustained_labels <- unique(raw_label[grepl("sustain", norm_label)])
    delayed_labels <- unique(raw_label[grepl("delayed|reconly", norm_label)])

    if (length(acute_labels) == 1L && length(sustained_labels) == 1L && length(delayed_labels) == 1L) {
      hsr_classes <- list(
        "Transient (acute)" = clean_gene_set(gs[[gene_idx]][raw_label == acute_labels[1]]),
        "Sustained" = clean_gene_set(gs[[gene_idx]][raw_label == sustained_labels[1]]),
        "Delayed" = clean_gene_set(gs[[gene_idx]][raw_label == delayed_labels[1]])
      )
      hsr_source <- "existing gene_sets table"
    }
  }
}

if (is.null(hsr_classes)) {
  class_vector_objects <- c("transient_genes", "sustained_genes", "delayed_genes")
  if (all(vapply(class_vector_objects, exists, logical(1), inherits = TRUE))) {
    hsr_classes <- list(
      "Transient (acute)" = clean_gene_set(get("transient_genes", inherits = TRUE)),
      "Sustained" = clean_gene_set(get("sustained_genes", inherits = TRUE)),
      "Delayed" = clean_gene_set(get("delayed_genes", inherits = TRUE))
    )
    hsr_source <- "existing transient_genes / sustained_genes / delayed_genes objects"
  }
}

if (is.null(hsr_classes) && exists("dds", inherits = TRUE)) {
  if (!requireNamespace("DESeq2", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("DESeq2 and SummarizedExperiment are required to rebuild HSR classes from dds.")
  }

  dds_obj <- get("dds", inherits = TRUE)
  if (!inherits(dds_obj, "DESeqDataSet")) stop("Object 'dds' is not a DESeqDataSet.")
  cd <- as.data.frame(SummarizedExperiment::colData(dds_obj))
  if (!all(c("genotype", "treatment") %in% names(cd))) {
    stop("colData(dds) must contain 'genotype' and 'treatment' columns to rebuild HSR classes.")
  }

  geno <- normalise_genotype(cd$genotype)
  cond <- normalise_condition(cd$treatment)
  wt_idx <- which(geno == "WT" & cond %in% condition_levels)
  if (!length(wt_idx)) stop("No WT/Col-0 samples were found in dds.")

  dds_wt <- dds_obj[, wt_idx]
  wt_cond <- factor(cond[wt_idx], levels = condition_levels)
  SummarizedExperiment::colData(dds_wt)$treatment <- wt_cond
  DESeq2::design(dds_wt) <- ~ treatment
  # Re-estimate size factors/dispersions after restricting to WT samples.
  DESeq2::sizeFactors(dds_wt) <- NULL
  DESeq2::normalizationFactors(dds_wt) <- NULL
  dds_wt <- DESeq2::DESeq(dds_wt, quiet = TRUE)

  res_lat <- as.data.frame(DESeq2::results(dds_wt, name = "treatment_LAT_vs_NT", alpha = alpha_cutoff)) %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_id = clean_gene_id(gene_id)) %>%
    transmute(gene_id, lfc_lat = log2FoldChange, padj_lat = padj)
  res_rec <- as.data.frame(DESeq2::results(dds_wt, name = "treatment_REC_vs_NT", alpha = alpha_cutoff)) %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_id = clean_gene_id(gene_id)) %>%
    transmute(gene_id, lfc_rec = log2FoldChange, padj_rec = padj)

  rebuilt_cls <- full_join(res_lat, res_rec, by = "gene_id") %>%
    mutate(
      induced_lat = !is.na(padj_lat) & padj_lat < alpha_cutoff &
        !is.na(lfc_lat) & lfc_lat > lfc_cutoff,
      induced_rec = !is.na(padj_rec) & padj_rec < alpha_cutoff &
        !is.na(lfc_rec) & lfc_rec > lfc_cutoff
    )

  hsr_classes <- list(
    "Transient (acute)" = clean_gene_set(rebuilt_cls$gene_id[rebuilt_cls$induced_lat & !rebuilt_cls$induced_rec]),
    "Sustained" = clean_gene_set(rebuilt_cls$gene_id[rebuilt_cls$induced_lat & rebuilt_cls$induced_rec]),
    "Delayed" = clean_gene_set(rebuilt_cls$gene_id[!rebuilt_cls$induced_lat & rebuilt_cls$induced_rec])
  )
  hsr_source <- paste0("rebuilt with WT-only DESeq2 contrasts (padj < ", alpha_cutoff,
                       "; log2FC > ", lfc_cutoff, ")")
}

if (is.null(hsr_classes)) {
  stop("Could not recover HSR classes. Keep transient_genes/sustained_genes/delayed_genes, gene_classification, gene_sets, or dds in the session.")
}
if (!all(class_names %in% names(hsr_classes))) {
  stop("The recovered HSR classes do not contain all required class names: ", paste(class_names, collapse = ", "))
}
hsr_classes <- hsr_classes[class_names]
hsr_classes <- lapply(hsr_classes, clean_gene_set)

# Apply exclusions to target genes only. TF reference expression remains available.
hsr_classes <- lapply(hsr_classes, function(x) setdiff(x, clean_gene_set(target_exclusions)))
if (any(vapply(hsr_classes, length, integer(1)) == 0L)) {
  empty_classes <- names(hsr_classes)[vapply(hsr_classes, length, integer(1)) == 0L]
  stop("At least one HSR class is empty after exclusions: ", paste(empty_classes, collapse = ", "))
}

cat("HSR classes recovered from: ", hsr_source, "\n", sep = "")
cat("Genes per HSR class (all HSR-responsive genes; no ChIP filter): ",
    paste(names(hsr_classes), lengths(hsr_classes), sep = " = ", collapse = "; "), "\n", sep = "")
cat("Target genes excluded from correlation analysis: ",
    paste(target_exclusions, collapse = ", "), "\n", sep = "")

# Save the gene membership actually used in the correlation analysis.
class_membership_all <- bind_rows(lapply(names(hsr_classes), function(cl) {
  tibble(HSR_class = cl, gene_id = hsr_classes[[cl]])
}))

# ------------------------- 2. Prepare expression matrix -----------------------
source_used <- expression_source
if (identical(expression_source, "auto")) {
  source_used <- if (exists("fpkm_wide", inherits = TRUE)) "fpkm_wide" else "dds_norm"
}
if (!source_used %in% c("fpkm_wide", "dds_norm")) {
  stop("expression_source must be 'auto', 'fpkm_wide', or 'dds_norm'.")
}

if (source_used == "fpkm_wide") {
  if (!exists("fpkm_wide", inherits = TRUE)) stop("expression_source='fpkm_wide' but object 'fpkm_wide' is not available.")
  expr_mat <- as.matrix(get("fpkm_wide", inherits = TRUE))
  if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) {
    stop("fpkm_wide must have gene IDs as row names and sample IDs as column names.")
  }
  suppressWarnings(storage.mode(expr_mat) <- "numeric")
  expr_mat[!is.finite(expr_mat) | expr_mat < 0] <- NA_real_
  expr_mat <- log2(expr_mat + 1)
  expression_label <- "log2(FPKM + 1)"
} else {
  if (!exists("dds", inherits = TRUE)) stop("expression_source='dds_norm' but object 'dds' is not available.")
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 to use normalized counts.")
  expr_mat <- DESeq2::counts(get("dds", inherits = TRUE), normalized = TRUE)
  suppressWarnings(storage.mode(expr_mat) <- "numeric")
  if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) stop("The normalized count matrix has no gene/sample names.")
  expr_mat[!is.finite(expr_mat) | expr_mat < 0] <- NA_real_
  expr_mat <- log2(expr_mat + 1)
  expression_label <- "log2(normalized count + 1)"
}

rownames(expr_mat) <- clean_gene_id(rownames(expr_mat))
if (anyNA(rownames(expr_mat)) || any(!nzchar(rownames(expr_mat)))) {
  stop("Expression matrix contains empty gene IDs after cleanup.")
}
if (anyDuplicated(rownames(expr_mat))) {
  expr_df <- as.data.frame(expr_mat, check.names = FALSE) %>%
    rownames_to_column("gene_id") %>%
    group_by(gene_id) %>%
    summarise(across(where(is.numeric), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)), .groups = "drop")
  expr_mat <- as.matrix(expr_df %>% column_to_rownames("gene_id"))
}

# Recover sample metadata from dds where possible; otherwise parse sample names.
get_sample_meta <- function(sample_names) {
  ans <- NULL
  if (exists("dds", inherits = TRUE) && requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    dds_obj <- get("dds", inherits = TRUE)
    cd <- as.data.frame(SummarizedExperiment::colData(dds_obj))
    cd$sample <- rownames(cd)
    if (all(sample_names %in% cd$sample) && all(c("genotype", "treatment") %in% names(cd))) {
      ans <- cd %>%
        filter(sample %in% sample_names) %>%
        transmute(
          sample,
          genotype = normalise_genotype(genotype),
          condition = normalise_condition(treatment)
        )
      ans <- ans[match(sample_names, ans$sample), , drop = FALSE]
      rownames(ans) <- NULL
    }
  }

  if (is.null(ans) || nrow(ans) != length(sample_names) ||
      anyNA(ans$genotype) || anyNA(ans$condition)) {
    if (!requireNamespace("stringr", quietly = TRUE)) {
      stop("Install stringr to infer genotype/treatment labels from sample names.")
    }
    mm <- stringr::str_match(
      sample_names,
      "(?i)(Col(?:-0)?|WT|UPF3|UPF-3)[_.-]+(NT|LAT|REC|P/T(?:[_. -]*\\+?[ ]*1[dD])?)[_.-]*([0-9]+)?"
    )
    if (is.null(mm) || anyNA(mm[, 2:3])) {
      stop("Could not recover sample metadata. Make sure dds colData has 'genotype'/'treatment', or use sample names like Col_NT_1 or UPF3_REC_3.")
    }
    ans <- data.frame(
      sample = sample_names,
      genotype = normalise_genotype(mm[, 2]),
      condition = normalise_condition(mm[, 3]),
      stringsAsFactors = FALSE
    )
  }

  ans %>% filter(!is.na(genotype), !is.na(condition), condition %in% condition_levels)
}

meta <- get_sample_meta(colnames(expr_mat))
if (!nrow(meta)) stop("No samples with recognized genotype and condition labels were found.")
if (!all(genotype_levels %in% unique(meta$genotype))) stop("Expression data must include both WT/Col-0 and upf3 samples.")
if (!all(condition_levels %in% unique(meta$condition))) stop("Expression data must include NT, LAT and REC samples.")
if (anyDuplicated(meta$sample)) stop("Duplicate sample names found in sample metadata.")
expr_mat <- expr_mat[, meta$sample, drop = FALSE]
meta <- meta[match(colnames(expr_mat), meta$sample), , drop = FALSE]

expr_long <- as.data.frame(expr_mat, check.names = FALSE) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "expr") %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(genotype), !is.na(condition), is.finite(expr))

traj_long <- expr_long %>%
  group_by(gene_id, genotype, condition) %>%
  summarise(value = mean(expr, na.rm = TRUE), n_rep = n(), .groups = "drop") %>%
  mutate(
    genotype = factor(genotype, levels = genotype_levels),
    condition = factor(condition, levels = condition_levels)
  )

cat("Expression source: ", source_used, " (", expression_label, ")\n", sep = "")
cat("Number of genes in expression matrix: ", nrow(expr_mat), "\n", sep = "")
cat("Sample distribution:\n")
print(meta %>% count(genotype, condition, name = "n_samples"))

# -------------------------- 3. WT TF reference profiles -----------------------
ref_ids <- c(HSFA2 = "AT2G26150", HSFA3 = "AT5G03720")
ref_table <- traj_long %>%
  filter(gene_id %in% unname(ref_ids), as.character(genotype) == "WT") %>%
  group_by(gene_id, condition) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(TF = ifelse(gene_id == ref_ids[["HSFA2"]], "HSFA2", "HSFA3")) %>%
  select(TF, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value)

get_ref_vec <- function(tf) {
  x <- ref_table %>% filter(TF == tf)
  if (nrow(x) != 1L || !all(condition_levels %in% names(x))) return(rep(NA_real_, 3L))
  as.numeric(unlist(x[1, condition_levels], use.names = FALSE))
}

hsfa2_vec <- get_ref_vec("HSFA2")
hsfa3_vec <- get_ref_vec("HSFA3")
if (any(!is.finite(hsfa2_vec))) stop("Cannot obtain a complete WT expression profile for HSFA2 (AT2G26150) across NT/LAT/REC.")
if (any(!is.finite(hsfa3_vec))) stop("Cannot obtain a complete WT expression profile for HSFA3 (AT5G03720) across NT/LAT/REC.")
if (stats::sd(hsfa2_vec) == 0) stop("HSFA2 WT profile has zero variance; Pearson correlation cannot be calculated.")
if (stats::sd(hsfa3_vec) == 0) stop("HSFA3 WT profile has zero variance; Pearson correlation cannot be calculated.")

ref_long <- bind_rows(
  tibble(
    TF = "HSFA2",
    correlation_group = "HSFA2-correlated",
    condition = factor(condition_levels, levels = condition_levels),
    ref_z = as.numeric(scale(hsfa2_vec))
  ),
  tibble(
    TF = "HSFA3",
    correlation_group = "HSFA3-correlated",
    condition = factor(condition_levels, levels = condition_levels),
    ref_z = as.numeric(scale(hsfa3_vec))
  )
) %>% filter(is.finite(ref_z))

# ------------------------- 4. Correlations for ALL HSR genes ------------------
hsr_ids <- unique(unlist(hsr_classes, use.names = FALSE))

wt_profiles <- traj_long %>%
  filter(gene_id %in% hsr_ids, as.character(genotype) == "WT") %>%
  select(gene_id, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value)
for (cond_name in condition_levels) {
  if (!cond_name %in% names(wt_profiles)) wt_profiles[[cond_name]] <- NA_real_
}

wt_profiles_complete <- wt_profiles %>%
  filter(if_all(all_of(condition_levels), ~ is.finite(.x)))

if (!nrow(wt_profiles_complete)) {
  stop("No HSR genes have a complete WT profile across NT, LAT and REC.")
}

profile_matrix <- as.matrix(wt_profiles_complete[, condition_levels, drop = FALSE])
wt_profiles_complete$r_HSFA2 <- apply(profile_matrix, 1, function(x) safe_cor(as.numeric(x), hsfa2_vec))
wt_profiles_complete$r_HSFA3 <- apply(profile_matrix, 1, function(x) safe_cor(as.numeric(x), hsfa3_vec))

# Preserve all class members in a table, including genes whose correlations
# are undefined or which lack a complete WT expression profile.
all_class_profiles <- class_membership_all %>%
  left_join(wt_profiles_complete %>% select(gene_id, all_of(condition_levels), r_HSFA2, r_HSFA3), by = "gene_id") %>%
  left_join(get_symbols(class_membership_all$gene_id), by = "gene_id") %>%
  mutate(
    corr_status = case_when(
      is.finite(r_HSFA2) & is.finite(r_HSFA3) & r_HSFA2 > correlation_cutoff & r_HSFA3 > correlation_cutoff ~ "Both > cutoff",
      is.finite(r_HSFA2) & is.finite(r_HSFA3) & r_HSFA2 > correlation_cutoff ~ "HSFA2 only",
      is.finite(r_HSFA2) & is.finite(r_HSFA3) & r_HSFA3 > correlation_cutoff ~ "HSFA3 only",
      is.finite(r_HSFA2) & is.finite(r_HSFA3) ~ "Neither",
      TRUE ~ "Correlation unavailable"
    )
  )

class_profile_table <- all_class_profiles %>%
  filter(is.finite(r_HSFA2), is.finite(r_HSFA3)) %>%
  mutate(HSR_class = factor(HSR_class, levels = class_names))

if (!nrow(class_profile_table)) stop("No class genes had finite correlations with both HSFA2 and HSFA3.")

cat("Genes per class available for correlation assignment:\n")
print(class_profile_table %>% count(HSR_class, name = "n_genes_with_correlations"))
cat("Correlation group membership (dual-correlated genes count in both TF groups):\n")
print(class_profile_table %>% group_by(HSR_class) %>% summarise(
  n_HSFA2_group = sum(r_HSFA2 > correlation_cutoff),
  n_HSFA3_group = sum(r_HSFA3 > correlation_cutoff),
  n_both = sum(r_HSFA2 > correlation_cutoff & r_HSFA3 > correlation_cutoff),
  n_neither = sum(r_HSFA2 <= correlation_cutoff & r_HSFA3 <= correlation_cutoff),
  .groups = "drop"
))

# -------------------------- 5. Correlation-group membership -------------------
trajectory_membership <- bind_rows(
  class_profile_table %>%
    filter(r_HSFA2 > correlation_cutoff) %>%
    transmute(HSR_class = as.character(HSR_class), gene_id, label,
              correlation_group = "HSFA2-correlated"),
  class_profile_table %>%
    filter(r_HSFA3 > correlation_cutoff) %>%
    transmute(HSR_class = as.character(HSR_class), gene_id, label,
              correlation_group = "HSFA3-correlated"),
  class_profile_table %>%
    filter(r_HSFA2 <= correlation_cutoff, r_HSFA3 <= correlation_cutoff) %>%
    transmute(HSR_class = as.character(HSR_class), gene_id, label,
              correlation_group = "Neither")
) %>%
  distinct(HSR_class, gene_id, correlation_group, .keep_all = TRUE) %>%
  mutate(
    HSR_class = factor(HSR_class, levels = class_names),
    correlation_group = factor(correlation_group, levels = correlation_group_levels)
  )

if (!nrow(trajectory_membership)) stop("No genes were available for the combined trajectory figure.")

# Six-point trajectory filter: each gene must have WT/upf3 x NT/LAT/REC.
trajectory_data <- traj_long %>%
  filter(gene_id %in% trajectory_membership$gene_id) %>%
  group_by(gene_id) %>%
  filter(n_distinct(paste(as.character(genotype), as.character(condition), sep = "_")) == 6L) %>%
  mutate(
    value_sd = stats::sd(value, na.rm = TRUE),
    # value_sd is repeated across the six rows in each grouped gene.
    z = dplyr::if_else(
      is.finite(value_sd) & value_sd > 0,
      as.numeric(scale(value)),
      rep(NA_real_, length(value))
    )
  ) %>%
  ungroup() %>%
  filter(is.finite(z)) %>%
  select(-value_sd) %>%
  inner_join(trajectory_membership, by = "gene_id") %>%
  mutate(
    HSR_class = factor(as.character(HSR_class), levels = class_names),
    correlation_group = factor(as.character(correlation_group), levels = correlation_group_levels),
    genotype = factor(as.character(genotype), levels = genotype_levels),
    condition = factor(as.character(condition), levels = condition_levels)
  )

if (!nrow(trajectory_data)) {
  stop("No correlation-selected genes have complete six-point WT/upf3 expression profiles.")
}

# Count unique genes per panel after all trajectory filters.
panel_counts <- trajectory_data %>%
  distinct(HSR_class, correlation_group, gene_id) %>%
  count(HSR_class, correlation_group, name = "n")

panel_grid <- tidyr::expand_grid(
  HSR_class = factor(class_names, levels = class_names),
  correlation_group = factor(correlation_group_levels, levels = correlation_group_levels)
) %>%
  left_join(panel_counts, by = c("HSR_class", "correlation_group")) %>%
  mutate(
    n = tidyr::replace_na(n, 0L),
    panel = paste0(as.character(HSR_class), "\n",
                   as.character(correlation_group), " (n = ", n, ")")
  )

panel_levels <- panel_grid$panel
panel_labels <- setNames(panel_grid$panel, panel_grid$panel)

trajectory_data <- trajectory_data %>%
  mutate(
    HSR_class_key = as.character(HSR_class),
    correlation_group_key = as.character(correlation_group)
  ) %>%
  left_join(
    panel_grid %>% transmute(
      HSR_class_key = as.character(HSR_class),
      correlation_group_key = as.character(correlation_group),
      panel
    ),
    by = c("HSR_class_key", "correlation_group_key")
  ) %>%
  mutate(panel = factor(panel, levels = panel_levels)) %>%
  select(-HSR_class_key, -correlation_group_key)

# Match TF reference profiles to the corresponding HSFA2/HSFA3 facets only.
ref_plot <- ref_long %>%
  filter(correlation_group != "Neither") %>%
  mutate(correlation_group_key = as.character(correlation_group)) %>%
  inner_join(
    panel_grid %>%
      filter(as.character(correlation_group) != "Neither") %>%
      transmute(
        HSR_class_key = as.character(HSR_class),
        correlation_group_key = as.character(correlation_group),
        panel
      ),
    by = "correlation_group_key"
  ) %>%
  mutate(panel = factor(panel, levels = panel_levels))

# -------------------------- 6. Combined 3 x 3 plot ----------------------------
p_trajectory_all <- ggplot(trajectory_data, aes(x = condition, y = z)) +
  geom_line(
    data = trajectory_data %>% filter(as.character(genotype) == "WT"),
    aes(group = interaction(gene_id, panel, drop = TRUE)),
    inherit.aes = TRUE,
    color = "grey75", alpha = 0.60, linewidth = 0.45,
    na.rm = TRUE
  ) +
  stat_summary(
    aes(color = genotype, group = genotype),
    fun = mean, geom = "line", linewidth = 1.25, na.rm = TRUE
  ) +
  stat_summary(
    aes(color = genotype, group = genotype),
    fun = mean, geom = "point", size = 2.3, na.rm = TRUE
  ) +
  geom_line(
    data = ref_plot,
    aes(x = condition, y = ref_z, group = 1),
    inherit.aes = FALSE,
    linetype = "dashed", color = "#08519C", linewidth = 1,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~panel, ncol = 3, drop = FALSE,
    labeller = labeller(panel = as_labeller(panel_labels))
  ) +
  scale_color_manual(values = c("WT" = "grey30", "upf3" = "#E3120B"), drop = FALSE) +
  scale_x_discrete(drop = FALSE, limits = condition_levels,
                   labels = c("NT", "P/T", "P/T +1d")) +
  labs(
    x = NULL,
    y = "Scaled expression (z-score)",
    color = NULL,
    title = "Expression trajectories across all WT-defined HSR genes",
    subtitle = paste0(
      "All HSR genes (no ChIP filter); Pearson r > ", correlation_cutoff,
      ". Dual-correlated genes appear in both TF panels. Grey = individual WT genes; ",
      "solid lines = genotype means; dashed blue = matching WT TF profile."
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    panel.spacing = grid::unit(0.8, "lines")
  )

print(p_trajectory_all)

# -------------------------- 7. Save plot and tables ---------------------------
out_png <- file.path(output_dir, paste0(output_prefix, "_all_HSR_trajectories.png"))
out_pdf <- file.path(output_dir, paste0(output_prefix, "_all_HSR_trajectories.pdf"))
ggsave(out_png, p_trajectory_all, width = 15, height = 11, dpi = 300, bg = "white")
ggsave(out_pdf, p_trajectory_all, width = 15, height = 11, bg = "white")

# All class genes and their correlations/status (including unavailable values).
correlation_csv <- file.path(output_dir, paste0(output_prefix, "_all_HSR_gene_correlations.csv"))
write.csv(all_class_profiles, correlation_csv, row.names = FALSE, na = "")

membership_csv <- file.path(output_dir, paste0(output_prefix, "_trajectory_membership_by_correlation_group.csv"))
write.csv(
  trajectory_membership %>% mutate(across(where(is.factor), as.character)),
  membership_csv, row.names = FALSE, na = ""
)

panel_csv <- file.path(output_dir, paste0(output_prefix, "_panel_gene_counts.csv"))
write.csv(panel_grid %>% mutate(across(where(is.factor), as.character)), panel_csv, row.names = FALSE)

class_sizes_csv <- file.path(output_dir, paste0(output_prefix, "_HSR_class_gene_counts.csv"))
complete_trajectory_counts <- trajectory_data %>%
  distinct(HSR_class, gene_id) %>%
  count(HSR_class, name = "n_with_complete_six_point_trajectory") %>%
  mutate(HSR_class = as.character(HSR_class))

class_sizes <- tibble(HSR_class = class_names) %>%
  mutate(
    n_HSR_targets_after_exclusions = vapply(hsr_classes[class_names], length, integer(1)),
    n_with_finite_HSFA2_and_HSFA3_correlations = vapply(class_names, function(cl) {
      sum(as.character(all_class_profiles$HSR_class) == cl &
            is.finite(all_class_profiles$r_HSFA2) & is.finite(all_class_profiles$r_HSFA3))
    }, integer(1))
  ) %>%
  left_join(complete_trajectory_counts, by = "HSR_class") %>%
  mutate(n_with_complete_six_point_trajectory = tidyr::replace_na(n_with_complete_six_point_trajectory, 0L))
write.csv(class_sizes, class_sizes_csv, row.names = FALSE)

cat("\nCombined trajectory figure saved:\n",
    normalizePath(out_png, mustWork = TRUE), "\n",
    normalizePath(out_pdf, mustWork = TRUE), "\n\n", sep = "")
cat("Correlation table saved:\n", normalizePath(correlation_csv, mustWork = TRUE), "\n", sep = "")
cat("Membership table saved:\n", normalizePath(membership_csv, mustWork = TRUE), "\n", sep = "")
cat("Panel count table saved:\n", normalizePath(panel_csv, mustWork = TRUE), "\n", sep = "")
cat("Class sizes saved:\n", normalizePath(class_sizes_csv, mustWork = TRUE), "\n\n", sep = "")
cat("HSR class source: ", hsr_source, "\n", sep = "")
cat("Expression source: ", source_used, " (", expression_label, ")\n", sep = "")
cat("Correlation threshold: r > ", correlation_cutoff, "\n", sep = "")
cat("Important: all WT-defined HSR class genes were eligible; no ChIP-positive filter was applied.\n")
