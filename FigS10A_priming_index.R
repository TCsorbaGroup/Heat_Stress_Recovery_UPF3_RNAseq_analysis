############################################################################
# Robustness of the NT primed state:
# Priming index with and without NMD targets
#
# Priming index:
#   log2FC(upf3 / Col-0 at NT) /
#   max(log2FC(Col-0 LAT / NT), log2FC(Col-0 REC / NT))
#
# Variants:
#   1. All HSR genes
#   2. Without empirical NMD targets
#   3. Without empirical NMD targets, 3'UTR-intron genes and uORF genes
#
# Groups:
#   All HSR classes
#   Transient (acute)
#   Transient (acute), ChIP-bound
#   Transient (acute), not bound
#   Sustained
#   Sustained, ChIP-bound
#   Sustained, not bound
#   Delayed
#   Delayed, ChIP-bound
#   Delayed, not bound
############################################################################


# ----------------------------- 1. Packages -----------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readxl)
})


# ----------------------------- 2. Settings -----------------------------

alpha_cutoff <- 0.05
lfc_cutoff <- 1

empirical_nmd_target_file <- "NMD_target_genes.txt"
chip_file <- "13059_2023_2970_MOESM3_ESM.xlsx"
chip_sheet <- "gene-clusters-with-short-descri"

output_prefix <- "Figure_NT_priming_without_NMD_targets_index"

priming_clip_min <- -0.5
priming_clip_max <- 1.5

# Plot group order.
nt_grp <- c(
  "All HSR classes",
  "Transient (acute)",
  "Transient (acute), ChIP-bound",
  "Transient (acute), not bound",
  "Sustained",
  "Sustained, ChIP-bound",
  "Sustained, not bound",
  "Delayed",
  "Delayed, ChIP-bound",
  "Delayed, not bound"
)

variant_order <- c(
  "All genes",
  "Without empirical NMD targets",
  "Without empirical targets,\n3'UTR-intron and uORF genes"
)

nt_x_lab <- c(
  "All HSR classes" = "All HSR classes",
  "Transient (acute)" = "Transient (acute)",
  "Transient (acute), ChIP-bound" = "Transient (acute)\nChIP-bound",
  "Transient (acute), not bound" = "Transient (acute)\nnot bound",
  "Sustained" = "Sustained",
  "Sustained, ChIP-bound" = "Sustained\nChIP-bound",
  "Sustained, not bound" = "Sustained\nnot bound",
  "Delayed" = "Delayed",
  "Delayed, ChIP-bound" = "Delayed\nChIP-bound",
  "Delayed, not bound" = "Delayed\nnot bound"
)

nt_cols <- c(
  "All HSR classes" = "grey65",
  "Transient (acute)" = "#E69F00",
  "Sustained" = "#D55E00",
  "Delayed" = "#2166AC"
)


# ----------------------------- 3. Helper functions ---------------------

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub(
    "^(gene|transcript|mRNA|rna):",
    "",
    x,
    ignore.case = TRUE
  )
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}


clean_gene_set <- function(x) {
  x <- clean_gene_id(x)
  unique(x[!is.na(x) & nzchar(x)])
}


extract_agi_ids <- function(x) {
  x <- as.character(unlist(x, use.names = FALSE))
  
  matches <- regmatches(
    x,
    gregexpr(
      "AT[1-5CM]G[0-9]{5}(?:\\.[0-9]+)?",
      x,
      ignore.case = TRUE,
      perl = TRUE
    )
  )
  
  clean_gene_set(
    unlist(matches, use.names = FALSE)
  )
}


normalise_genotype <- function(x) {
  z <- tolower(
    gsub("[^[:alnum:]]+", "", trimws(as.character(x)))
  )
  
  dplyr::case_when(
    z %in% c("col", "col0", "wt", "wildtype") ~ "Col-0",
    z == "upf3" ~ "upf3",
    TRUE ~ NA_character_
  )
}


normalise_treatment <- function(x) {
  z <- toupper(
    gsub("[[:space:]]+", "", trimws(as.character(x)))
  )
  
  dplyr::case_when(
    z %in% c("NT", "NOTREATMENT", "UNTREATED") ~ "NT",
    z %in% c("LAT", "PT", "P/T", "HEAT", "HEATSTRESS") ~ "LAT",
    z %in% c("REC", "PT+1D", "P/T+1D", "RECOVERY", "RECOVERED") ~ "REC",
    TRUE ~ NA_character_
  )
}


as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA, x != 0))
  }
  
  z <- tolower(trimws(as.character(x)))
  
  out <- rep(NA, length(z))
  
  out[z %in% c("true", "t", "1", "yes", "y", "positive", "pos")] <- TRUE
  
  out[z %in% c("false", "f", "0", "no", "n", "negative", "neg")] <- FALSE
  
  out
}

extract_lfc_table <- function(x, source_label = "input table") {
  
  if (is.null(x)) {
    return(NULL)
  }
  
  d <- tryCatch(
    as.data.frame(x),
    error = function(e) NULL
  )
  
  if (is.null(d) || nrow(d) == 0L) {
    return(NULL)
  }
  
  norm_names <- tolower(
    gsub("[^[:alnum:]]+", "", names(d))
  )
  
  id_idx <- which(
    norm_names %in% c("geneid", "gene", "tair", "agi", "geneidentifier")
  )
  
  if (length(id_idx) > 0L) {
    gene_ids <- clean_gene_id(d[[id_idx[1]]])
  } else {
    gene_ids <- clean_gene_id(rownames(d))
  }
  
  
  lfc_idx <- which(
    norm_names %in% c(
      "log2fc",
      "log2foldchange",
      "lfc",
      "logfc"
    )
  )
  
  
  if (length(lfc_idx) == 0L) {
    return(NULL)
  }
  
  
  padj_idx <- which(
    norm_names %in% c(
      "padj",
      "adjpvalue",
      "padjusted",
      "p adj",
      "padjbh",
      "fdr"
    )
  )
  
  
  out <- tibble(
    gene_id = gene_ids,
    log2FC = suppressWarnings(
      as.numeric(d[[lfc_idx[1]]])
    ),
    padj = if (length(padj_idx) > 0L) {
      suppressWarnings(
        as.numeric(d[[padj_idx[1]]])
      )
    } else {
      NA_real_
    }
  ) %>%
    filter(
      !is.na(gene_id),
      nzchar(gene_id)
    )
  
  
  if (anyDuplicated(out$gene_id)) {
    stop(
      "Duplicate gene IDs in ", source_label,
      " after gene-ID cleanup."
    )
  }
  
  out
}


canonical_class <- function(x) {
  z <- tolower(
    gsub("[^[:alnum:]]+", "", trimws(as.character(x)))
  )
  
  dplyr::case_when(
    grepl("acute|transient", z) ~ "Acute",
    grepl("sustain", z) ~ "Sustained",
    grepl("delayed", z) ~ "Delayed",
    TRUE ~ NA_character_
  )
}


find_gene_id_column <- function(df) {
  norm_names <- tolower(
    gsub("[^[:alnum:]]+", "", names(df))
  )
  
  idx <- which(
    norm_names %in% c("geneid", "gene", "tair", "agi", "geneidentifier")
  )
  
  if (length(idx) == 0L) {
    return(NA_integer_)
  }
  
  idx[1]
}


find_flag_column <- function(df, aliases) {
  norm_names <- tolower(
    gsub("[^[:alnum:]]+", "", names(df))
  )
  
  aliases <- tolower(
    gsub("[^[:alnum:]]+", "", aliases)
  )
  
  idx <- which(norm_names %in% aliases)
  
  if (length(idx) == 0L) {
    return(NA_integer_)
  }
  
  idx[1]
}


extract_flag_table <- function(
    df,
    aliases,
    output_flag,
    source_label = "feature table"
) {
  
  if (is.null(df)) {
    return(NULL)
  }
  
  df <- tryCatch(
    as.data.frame(df),
    error = function(e) NULL
  )
  
  if (is.null(df) || nrow(df) == 0L) {
    return(NULL)
  }
  
  
  gene_idx <- find_gene_id_column(df)
  flag_idx <- find_flag_column(df, aliases)
  
  
  if (is.na(gene_idx) || is.na(flag_idx)) {
    return(NULL)
  }
  
  
  out <- tibble(
    gene_id = clean_gene_id(df[[gene_idx]]),
    flag = as_logical_flag(df[[flag_idx]])
  ) %>%
    filter(
      !is.na(gene_id),
      nzchar(gene_id)
    )
  
  
  if (anyDuplicated(out$gene_id)) {
    
    # If the source is not gene-unique, collapse to a conservative,
    # gene-level flag. TRUE wins; otherwise retain FALSE if known.
    out <- out %>%
      group_by(gene_id) %>%
      summarise(
        flag = if (any(flag %in% TRUE)) {
          TRUE
        } else if (any(!is.na(flag))) {
          FALSE
        } else {
          NA
        },
        .groups = "drop"
      )
  }
  
  names(out)[2] <- output_flag
  
  out
}


# ----------------------------- 4. Validate dds -------------------------

if (!exists("dds", inherits = TRUE)) {
  stop(
    "The DESeq2 object 'dds' is not available. ",
    "Load the current DESeq2 dataset first."
  )
}

if (!inherits(dds, "DESeqDataSet")) {
  stop("Object 'dds' must be a DESeqDataSet.")
}


raw_counts <- DESeq2::counts(
  dds,
  normalized = FALSE
)


count_gene_ids <- clean_gene_id(
  rownames(raw_counts)
)


if (anyNA(count_gene_ids) || any(!nzchar(count_gene_ids))) {
  stop("The dds count matrix contains missing or empty gene IDs.")
}


if (anyDuplicated(count_gene_ids)) {
  stop("Duplicate gene IDs remain in the dds count matrix after cleanup.")
}


rownames(raw_counts) <- count_gene_ids


sample_info <- as.data.frame(
  SummarizedExperiment::colData(dds)
) %>%
  tibble::rownames_to_column("Sample")


if (!all(c("genotype", "treatment") %in% names(sample_info))) {
  stop(
    "colData(dds) must contain 'genotype' and 'treatment'."
  )
}


sample_info <- sample_info %>%
  mutate(
    Genotype = normalise_genotype(genotype),
    Condition = normalise_treatment(treatment)
  )


if (anyDuplicated(sample_info$Sample)) {
  stop("Duplicate sample names in colData(dds).")
}


if (anyNA(sample_info$Genotype) || anyNA(sample_info$Condition)) {
  bad <- sample_info %>%
    filter(is.na(Genotype) | is.na(Condition)) %>%
    pull(Sample)
  
  stop(
    "Unrecognised genotype or treatment labels in colData(dds): ",
    paste(bad, collapse = ", "),
    ". Update normalise_genotype() or normalise_treatment()."
  )
}


if (!all(sample_info$Sample %in% colnames(raw_counts))) {
  stop("Sample names in colData(dds) do not match raw-count columns.")
}


sample_info <- sample_info[
  match(colnames(raw_counts), sample_info$Sample),
  ,
  drop = FALSE
]


if (anyNA(sample_info$Sample)) {
  stop("Could not match all dds count columns to sample metadata.")
}


# ----------------------------- 5. Recover genotype and WT effects -------

# Use the existing NT/LAT/REC contrasts when they are available and
# contain the expected log2FC columns. Otherwise refit DESeq2 within
# each condition, with genotype as the sole predictor.

fit_genotype_effect_for_condition <- function(condition_name) {
  
  selected <- sample_info %>%
    filter(
      Condition == condition_name,
      Genotype %in% c("Col-0", "upf3")
    )
  
  
  if (nrow(selected) == 0L) {
    stop("No samples found for condition ", condition_name, ".")
  }
  
  
  selected$Genotype <- factor(
    selected$Genotype,
    levels = c("Col-0", "upf3")
  )
  
  
  sample_counts <- table(selected$Genotype)
  
  if (any(sample_counts < 2L)) {
    stop(
      "At least two biological replicates per genotype are required in ",
      condition_name, ". Found: ",
      paste(names(sample_counts), sample_counts, sep = "=", collapse = ", ")
    )
  }
  
  
  counts_sub <- raw_counts[
    ,
    selected$Sample,
    drop = FALSE
  ]
  
  
  coldata_sub <- data.frame(
    genotype_fit = factor(
      selected$Genotype,
      levels = c("Col-0", "upf3")
    ),
    row.names = selected$Sample
  )
  
  
  dds_sub <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = coldata_sub,
    design = ~ genotype_fit
  )
  
  
  dds_sub <- DESeq2::DESeq(
    dds_sub,
    quiet = TRUE
  )
  
  
  res <- DESeq2::results(
    dds_sub,
    contrast = c("genotype_fit", "upf3", "Col-0"),
    alpha = alpha_cutoff
  )
  
  
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    transmute(
      gene_id = clean_gene_id(gene_id),
      log2FC = log2FoldChange,
      padj = padj
    )
}


get_genotype_effect <- function(condition_name) {
  
  if (exists(condition_name, inherits = TRUE)) {
    
    candidate <- extract_lfc_table(
      get(condition_name, inherits = TRUE),
      paste0("existing ", condition_name, " table")
    )
    
    
    if (!is.null(candidate)) {
      message(
        "Using existing genotype-effect table for ",
        condition_name, "."
      )
      return(candidate)
    }
  }
  
  
  message(
    "Refitting genotype effect for ",
    condition_name,
    " from dds raw counts."
  )
  
  
  fit_genotype_effect_for_condition(condition_name)
}


nt_effect <- get_genotype_effect("NT") %>%
  transmute(
    gene_id,
    nt_log2FC = log2FC,
    nt_padj = padj
  )


# Prefer existing WT heat-response results when available.
wt_lat <- NULL
wt_rec <- NULL


if (
  exists("wt_tabs", inherits = TRUE) &&
  is.list(get("wt_tabs", inherits = TRUE))
) {
  
  wt_tabs_obj <- get("wt_tabs", inherits = TRUE)
  
  
  if ("LAT" %in% names(wt_tabs_obj)) {
    wt_lat <- extract_lfc_table(
      wt_tabs_obj$LAT,
      "wt_tabs$LAT"
    )
  }
  
  
  if ("REC" %in% names(wt_tabs_obj)) {
    wt_rec <- extract_lfc_table(
      wt_tabs_obj$REC,
      "wt_tabs$REC"
    )
  }
}


# If either table is missing, fit treatment effects using Col-0 samples only.
if (is.null(wt_lat) || is.null(wt_rec)) {
  
  message(
    "Existing wt_tabs are incomplete; fitting Col-0 treatment effects ",
    "from the current dds object."
  )
  
  
  wt_meta <- sample_info %>%
    filter(
      Genotype == "Col-0",
      Condition %in% c("NT", "LAT", "REC")
    )
  
  
  treatment_counts <- table(
    factor(
      wt_meta$Condition,
      levels = c("NT", "LAT", "REC")
    )
  )
  
  
  if (any(treatment_counts < 2L)) {
    stop(
      "At least two Col-0 biological replicates are required in NT, LAT ",
      "and REC to calculate WT induction. Found: ",
      paste(
        names(treatment_counts),
        treatment_counts,
        sep = "=",
        collapse = ", "
      )
    )
  }
  
  
  counts_wt <- raw_counts[
    ,
    wt_meta$Sample,
    drop = FALSE
  ]
  
  
  coldata_wt <- data.frame(
    treatment_fit = factor(
      wt_meta$Condition,
      levels = c("NT", "LAT", "REC")
    ),
    row.names = wt_meta$Sample
  )
  
  
  dds_wt <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_wt,
    colData = coldata_wt,
    design = ~ treatment_fit
  )
  
  
  dds_wt <- DESeq2::DESeq(
    dds_wt,
    quiet = TRUE
  )
  
  
  res_lat <- as.data.frame(
    DESeq2::results(
      dds_wt,
      contrast = c("treatment_fit", "LAT", "NT"),
      alpha = alpha_cutoff
    )
  ) %>%
    tibble::rownames_to_column("gene_id") %>%
    transmute(
      gene_id = clean_gene_id(gene_id),
      log2FC = log2FoldChange,
      padj = padj
    )
  
  
  res_rec <- as.data.frame(
    DESeq2::results(
      dds_wt,
      contrast = c("treatment_fit", "REC", "NT"),
      alpha = alpha_cutoff
    )
  ) %>%
    tibble::rownames_to_column("gene_id") %>%
    transmute(
      gene_id = clean_gene_id(gene_id),
      log2FC = log2FoldChange,
      padj = padj
    )
  
  
  if (is.null(wt_lat)) wt_lat <- res_lat
  if (is.null(wt_rec)) wt_rec <- res_rec
}


if (
  anyDuplicated(wt_lat$gene_id) ||
  anyDuplicated(wt_rec$gene_id)
) {
  stop("Duplicate gene IDs in the WT treatment-effect tables.")
}


wt_effect <- full_join(
  wt_lat %>%
    transmute(
      gene_id,
      wt_log2FC_LAT = log2FC,
      wt_padj_LAT = padj
    ),
  
  wt_rec %>%
    transmute(
      gene_id,
      wt_log2FC_REC = log2FC,
      wt_padj_REC = padj
    ),
  
  by = "gene_id"
)


wt_effect$wt_induction_max <- apply(
  as.matrix(
    wt_effect[, c("wt_log2FC_LAT", "wt_log2FC_REC")]
  ),
  1,
  function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }
)


# The denominator is the maximum of the two WT induction estimates.
# A non-positive or unavailable denominator does not define a meaningful
# priming ratio and is excluded from the boxplots.
priming_table <- nt_effect %>%
  inner_join(
    wt_effect,
    by = "gene_id"
  ) %>%
  mutate(
    priming_index = ifelse(
      is.finite(nt_log2FC) &
        is.finite(wt_induction_max) &
        wt_induction_max > 0,
      nt_log2FC / wt_induction_max,
      NA_real_
    )
  )


# ----------------------------- 6. Recover the HSR classes ---------------

hsr_classes <- NULL
hsr_source <- NULL

if (exists("gene_sets", inherits = TRUE)) {
  
  gs <- get("gene_sets", inherits = TRUE)
  
  
  if (
    is.data.frame(gs) &&
    all(c("gene", "GeneGroup") %in% names(gs))
  ) {
    
    group_class <- vapply(
      as.character(gs$GeneGroup),
      canonical_class,
      character(1)
    )
    
    
    candidate_classes <- list(
      
      Acute = clean_gene_set(
        gs$gene[group_class == "Acute"]
      ),
      
      Sustained = clean_gene_set(
        gs$gene[group_class == "Sustained"]
      ),
      
      Delayed = clean_gene_set(
        gs$gene[group_class == "Delayed"]
      )
    )
    
    
    if (all(lengths(candidate_classes) > 0L)) {
      hsr_classes <- candidate_classes
      hsr_source <- "existing gene_sets data.frame"
    }
  }
}

if (
  is.null(hsr_classes) &&
  exists("gene_classification", inherits = TRUE)
) {
  
  gc <- get(
    "gene_classification",
    inherits = TRUE
  )
  
  
  if (
    is.data.frame(gc) &&
    all(
      c(
        "gene",
        "is_acute",
        "is_sustained",
        "is_delayed"
      ) %in% names(gc)
    )
  ) {
    
    hsr_classes <- list(
      
      Acute = clean_gene_set(
        gc$gene[!is.na(gc$is_acute) & gc$is_acute]
      ),
      
      Sustained = clean_gene_set(
        gc$gene[!is.na(gc$is_sustained) & gc$is_sustained]
      ),
      
      Delayed = clean_gene_set(
        gc$gene[!is.na(gc$is_delayed) & gc$is_delayed]
      )
    )
    
    hsr_source <- "existing gene_classification flags"
    
    
  } else if (
    is.data.frame(gc) &&
    all(
      c(
        "gene",
        "is_WT_LAT_induced",
        "is_WT_REC_induced"
      ) %in% names(gc)
    )
  ) {
    
    lat_induced <- !is.na(gc$is_WT_LAT_induced) &
      gc$is_WT_LAT_induced
    
    rec_induced <- !is.na(gc$is_WT_REC_induced) &
      gc$is_WT_REC_induced
    
    
    hsr_classes <- list(
      
      Acute = clean_gene_set(
        gc$gene[lat_induced & !rec_induced]
      ),
      
      Sustained = clean_gene_set(
        gc$gene[lat_induced & rec_induced]
      ),
      
      Delayed = clean_gene_set(
        gc$gene[!lat_induced & rec_induced]
      )
    )
    
    hsr_source <- "existing LAT/REC induction flags"
  }
}

if (is.null(hsr_classes)) {
  
  lat_induced <- !is.na(wt_lat$padj) &
    wt_lat$padj < alpha_cutoff &
    !is.na(wt_lat$log2FC) &
    wt_lat$log2FC > lfc_cutoff
  
  
  rec_induced <- !is.na(wt_rec$padj) &
    wt_rec$padj < alpha_cutoff &
    !is.na(wt_rec$log2FC) &
    wt_rec$log2FC > lfc_cutoff
  
  
  lat_genes <- clean_gene_set(
    wt_lat$gene_id[lat_induced]
  )
  
  rec_genes <- clean_gene_set(
    wt_rec$gene_id[rec_induced]
  )
  
  
  hsr_classes <- list(
    
    Acute = setdiff(lat_genes, rec_genes),
    
    Sustained = intersect(lat_genes, rec_genes),
    
    Delayed = setdiff(rec_genes, lat_genes)
  )
  
  
  hsr_source <- "reconstructed from Col-0 treatment effects"
}


hsr_classes <- lapply(
  hsr_classes[c("Acute", "Sustained", "Delayed")],
  clean_gene_set
)


if (any(lengths(hsr_classes) == 0L)) {
  print(lengths(hsr_classes))
  stop(
    "At least one of the three HSR classes is empty. ",
    "Check the current HSR gene-set definitions."
  )
}


hsr_pairwise <- combn(
  names(hsr_classes),
  2,
  simplify = FALSE
)


hsr_overlap <- vapply(
  hsr_pairwise,
  function(pair) {
    length(
      intersect(
        hsr_classes[[pair[1]]],
        hsr_classes[[pair[2]]]
      )
    )
  },
  integer(1)
)


if (any(hsr_overlap > 0L)) {
  stop(
    "The Acute, Sustained and Delayed classes are not mutually exclusive ",
    "after gene-ID cleanup."
  )
}


message("HSR class source: ", hsr_source)


# ----------------------------- 7. Recover ChIP target genes -------------

chip_targets <- NULL
chip_source <- NULL


if (exists("hsfa23_all_relevant", inherits = TRUE)) {
  
  chip_targets <- extract_agi_ids(
    get("hsfa23_all_relevant", inherits = TRUE)
  )
  
  if (length(chip_targets) > 0L) {
    chip_source <- "existing hsfa23_all_relevant object"
  }
}


if (
  (is.null(chip_targets) || length(chip_targets) == 0L) &&
  exists("chip_all", inherits = TRUE)
) {
  
  chip_targets <- extract_agi_ids(
    get("chip_all", inherits = TRUE)
  )
  
  if (length(chip_targets) > 0L) {
    chip_source <- "existing chip_all object"
  }
}

if (
  (is.null(chip_targets) || length(chip_targets) == 0L) &&
  exists("q_group_sets", inherits = TRUE)
) {
  
  q_sets <- get(
    "q_group_sets",
    inherits = TRUE
  )
  
  
  if (
    is.list(q_sets) &&
    !is.null(names(q_sets))
  ) {
    
    bound_names <- names(q_sets)[
      grepl("ChIP-bound", names(q_sets), ignore.case = TRUE) &
        !grepl("not bound", names(q_sets), ignore.case = TRUE)
    ]
    
    
    if (length(bound_names) > 0L) {
      chip_targets <- clean_gene_set(
        unlist(
          q_sets[bound_names],
          use.names = FALSE
        )
      )
      
      chip_source <- "existing ChIP-bound q_group_sets"
    }
  }
}

if (
  (is.null(chip_targets) || length(chip_targets) == 0L)
) {
  
  if (!file.exists(chip_file)) {
    stop(
      "Could not recover the broad HsfA2/HsfA3 ChIP target set. ",
      "Load hsfa23_all_relevant/chip_all, or make this file available: ",
      chip_file
    )
  }
  
  
  chip_df <- as.data.frame(
    readxl::read_excel(
      chip_file,
      sheet = chip_sheet
    )
  )
  
  
  if (ncol(chip_df) < 2L) {
    stop("ChIP worksheet must contain gene IDs and cluster labels.")
  }
  
  
  cluster_txt <- as.character(
    chip_df[[2]]
  )
  
  
  cluster_num <- ifelse(
    grepl("[0-9]+", cluster_txt),
    sub(
      "^[^0-9]*([0-9]+).*",
      "\\1",
      cluster_txt
    ),
    NA_character_
  )
  
  
  cluster_num <- suppressWarnings(
    as.integer(cluster_num)
  )
  
  
  chip_targets <- extract_agi_ids(
    chip_df[[1]][
      cluster_num %in% c(7, 11, 12, 14)
    ]
  )
  
  
  chip_source <- paste0(
    "Excel: ", chip_file,
    " [", chip_sheet, "]"
  )
}


if (length(chip_targets) == 0L) {
  stop("The recovered ChIP target set is empty.")
}


message("ChIP target source: ", chip_source)
message("Unique ChIP target genes: ", length(chip_targets))


# ----------------------------- 8. Recover NMD and structural features ----

de_features <- NULL


if (exists("de_features_unique", inherits = TRUE)) {
  
  de_features <- as.data.frame(
    get("de_features_unique", inherits = TRUE)
  )
  
  
  gene_idx <- find_gene_id_column(
    de_features
  )
  
  
  if (is.na(gene_idx)) {
    stop(
      "Object de_features_unique exists but has no recognisable gene-ID column."
    )
  }
  
  
  de_features$gene_id <- clean_gene_id(
    de_features[[gene_idx]]
  )
  
  
  if (anyNA(de_features$gene_id) || any(!nzchar(de_features$gene_id))) {
    stop("Missing/empty gene IDs in de_features_unique.")
  }
  
  
  if (anyDuplicated(de_features$gene_id)) {
    stop(
      "de_features_unique contains duplicated gene IDs after cleanup. ",
      "Resolve these before applying gene exclusions."
    )
  }
}

if (!file.exists(empirical_nmd_target_file)) {
  stop(
    "Cannot find the empirical NMD-target list: ",
    empirical_nmd_target_file,
    "\nWorking directory: ",
    getwd()
  )
}

nmd_lines <- readLines(
  empirical_nmd_target_file,
  warn = FALSE
)

empirical_nmd_targets <- extract_agi_ids(nmd_lines)

if (length(empirical_nmd_targets) == 0L) {
  stop(
    "No Arabidopsis TAIR gene IDs could be extracted from ",
    empirical_nmd_target_file,
    ". Check the file contents and identifier format."
  )
}

empirical_nmd_targets <- clean_gene_set(
  empirical_nmd_targets
)

nmd_target_source <- normalizePath(
  empirical_nmd_target_file,
  winslash = "/",
  mustWork = TRUE
)

message(
  "Empirical NMD-target source: ",
  nmd_target_source
)

message(
  "Unique empirical NMD-target genes: ",
  length(empirical_nmd_targets)
)

utr_flag_table <- NULL


if (!is.null(de_features)) {
  
  utr_flag_table <- extract_flag_table(
    de_features,
    aliases = c(
      "has_3utr_intron",
      "has_3prime_utr_intron",
      "utr_intron"
    ),
    output_flag = "has_3utr_intron",
    source_label = "de_features_unique"
  )
}


if (
  is.null(utr_flag_table) &&
  exists("gene_features_3utr", inherits = TRUE)
) {
  
  utr_flag_table <- extract_flag_table(
    as.data.frame(
      get("gene_features_3utr", inherits = TRUE)
    ),
    aliases = c(
      "has_3utr_intron",
      "has_3prime_utr_intron",
      "utr_intron"
    ),
    output_flag = "has_3utr_intron",
    source_label = "gene_features_3utr"
  )
}


utr_candidate_files <- c(
  "Col_HSR_five_gene_sets_cis_feature_comparison_representative_utr_intron_by_gene.csv",
  "Col_HSR_five_gene_sets_cis_feature_comparison_utr_intron_by_gene.csv"
)


if (is.null(utr_flag_table)) {
  
  existing_file <- utr_candidate_files[
    file.exists(utr_candidate_files)
  ]
  
  
  if (length(existing_file) > 0L) {
    
    utr_file_df <- read.csv(
      existing_file[1],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    
    utr_flag_table <- extract_flag_table(
      utr_file_df,
      aliases = c(
        "has_3utr_intron",
        "has_3prime_utr_intron",
        "utr_intron"
      ),
      output_flag = "has_3utr_intron",
      source_label = existing_file[1]
    )
  }
}


if (is.null(utr_flag_table)) {
  stop(
    "The has_3utr_intron annotation could not be recovered. ",
    "Load de_features_unique or gene_features_3utr, or provide the ",
    "per-gene structural-feature output from the preceding cis-feature analysis."
  )
}

uorf_flag_table <- NULL


if (!is.null(de_features)) {
  
  uorf_flag_table <- extract_flag_table(
    de_features,
    aliases = c(
      "has_uORF",
      "RiboTaper",
      "has_ribotaper_uorf"
    ),
    output_flag = "has_uORF",
    source_label = "de_features_unique"
  )
}


if (
  is.null(uorf_flag_table) &&
  exists("uorf_df", inherits = TRUE)
) {
  
  uorf_flag_table <- extract_flag_table(
    as.data.frame(
      get("uorf_df", inherits = TRUE)
    ),
    aliases = c(
      "has_uORF",
      "RiboTaper",
      "has_ribotaper_uorf"
    ),
    output_flag = "has_uORF",
    source_label = "uorf_df"
  )
}


uorf_candidate_files <- c(
  "Col_HSR_five_gene_sets_cis_feature_comparison_uorf_status_by_gene.csv",
  "Col_HSR_five_gene_sets_cis_feature_comparison_uORF_status_by_gene.csv"
)


if (is.null(uorf_flag_table)) {
  
  existing_file <- uorf_candidate_files[
    file.exists(uorf_candidate_files)
  ]
  
  
  if (length(existing_file) > 0L) {
    
    uorf_file_df <- read.csv(
      existing_file[1],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    
    uorf_flag_table <- extract_flag_table(
      uorf_file_df,
      aliases = c(
        "has_uORF",
        "RiboTaper",
        "has_ribotaper_uorf"
      ),
      output_flag = "has_uORF",
      source_label = existing_file[1]
    )
  }
}


if (is.null(uorf_flag_table)) {
  stop(
    "The has_uORF annotation could not be recovered. ",
    "Load de_features_unique or uorf_df, or provide the per-gene ",
    "uORF status table from the preceding cis-feature analysis."
  )
}


# Merge structural annotations. Missing entries remain NA rather than
# being treated as confirmed negatives.
feature_table <- full_join(
  utr_flag_table,
  uorf_flag_table,
  by = "gene_id"
) %>%
  mutate(
    has_3utr_intron = as_logical_flag(has_3utr_intron),
    has_uORF = as_logical_flag(has_uORF)
  )

if (!is.null(de_features)) {
  
  utr_from_de <- extract_flag_table(
    de_features,
    aliases = c(
      "has_3utr_intron",
      "has_3prime_utr_intron",
      "utr_intron"
    ),
    output_flag = "has_3utr_intron_de",
    source_label = "de_features_unique"
  )
  
  
  uorf_from_de <- extract_flag_table(
    de_features,
    aliases = c(
      "has_uORF",
      "RiboTaper",
      "has_ribotaper_uorf"
    ),
    output_flag = "has_uORF_de",
    source_label = "de_features_unique"
  )
  
  
  if (!is.null(utr_from_de)) {
    feature_table <- full_join(
      feature_table,
      utr_from_de,
      by = "gene_id"
    ) %>%
      mutate(
        has_3utr_intron = coalesce(
          has_3utr_intron,
          has_3utr_intron_de
        )
      ) %>%
      dplyr::select(-has_3utr_intron_de)
  }
  
  
  if (!is.null(uorf_from_de)) {
    feature_table <- full_join(
      feature_table,
      uorf_from_de,
      by = "gene_id"
    ) %>%
      mutate(
        has_uORF = coalesce(
          has_uORF,
          has_uORF_de
        )
      ) %>%
      dplyr::select(-has_uORF_de)
  }
}


# Collapse duplicate entries conservatively if merging sources introduced them.
feature_table <- feature_table %>%
  group_by(gene_id) %>%
  summarise(
    has_3utr_intron = if (any(has_3utr_intron %in% TRUE)) {
      TRUE
    } else if (any(!is.na(has_3utr_intron))) {
      FALSE
    } else {
      NA
    },
    
    has_uORF = if (any(has_uORF %in% TRUE)) {
      TRUE
    } else if (any(!is.na(has_uORF))) {
      FALSE
    } else {
      NA
    },
    
    .groups = "drop"
  )


if (
  !any(feature_table$has_3utr_intron %in% TRUE) ||
  !any(feature_table$has_uORF %in% TRUE)
) {
  warning(
    "At least one structural feature has no TRUE-positive genes in the ",
    "loaded annotation. Verify that the intended current feature tables ",
    "were loaded."
  )
}


# ----------------------------- 9. Define the HSR gene table -------------

class_name_map <- c(
  "Acute" = "Transient (acute)",
  "Sustained" = "Sustained",
  "Delayed" = "Delayed"
)


hsr_class_table <- bind_rows(
  lapply(
    names(hsr_classes),
    function(nm) {
      tibble(
        gene_id = hsr_classes[[nm]],
        hsr_class = class_name_map[[nm]]
      )
    }
  )
) %>%
  filter(
    gene_id %in% rownames(raw_counts)
  ) %>%
  distinct(
    gene_id,
    hsr_class
  )


if (nrow(hsr_class_table) == 0L) {
  stop("No HSR genes match the gene IDs in the current dds object.")
}


if (anyDuplicated(hsr_class_table$gene_id)) {
  stop(
    "At least one gene belongs to multiple HSR classes after ID cleanup."
  )
}


# Retain genes with both the NT genotype effect and a finite positive
# WT heat-induction denominator, so a priming index can be calculated.
nt_hsr <- hsr_class_table %>%
  left_join(
    priming_table,
    by = "gene_id"
  ) %>%
  filter(
    is.finite(nt_log2FC),
    is.finite(wt_induction_max),
    wt_induction_max > 0,
    is.finite(priming_index)
  ) %>%
  left_join(
    feature_table,
    by = "gene_id"
  ) %>%
  mutate(
    status = if_else(
      gene_id %in% chip_targets,
      "ChIP-bound",
      "not bound"
    )
  )


if (nrow(nt_hsr) == 0L) {
  stop(
    "No HSR genes have a calculable finite priming index. ",
    "Check the NT and WT treatment-effect tables."
  )
}


if (anyDuplicated(nt_hsr$gene_id)) {
  stop("Duplicate HSR gene IDs in the final priming-index table.")
}


message("")
message("HSR genes with a valid priming index: ", nrow(nt_hsr))
message("Class counts:")
print(
  nt_hsr %>%
    dplyr::count(hsr_class, name = "n_genes")
)


# ----------------------------- 10. Define NMD exclusion variants --------

nt_nmd_emp <- empirical_nmd_targets


# Strict exclusion includes all genes with any of:
#   - empirical NMD target annotation
#   - TRUE has_3utr_intron flag
#   - TRUE has_uORF flag
#
# NA feature values are not treated as TRUE and therefore are not removed.
nt_nmd_strict <- unique(
  c(
    nt_nmd_emp,
    feature_table$gene_id[
      feature_table$has_3utr_intron %in% TRUE
    ],
    feature_table$gene_id[
      feature_table$has_uORF %in% TRUE
    ]
  )
)


nt_variants <- list(
  "All genes" = character(0),
  "Without empirical NMD targets" = nt_nmd_emp,
  "Without empirical targets,\n3'UTR-intron and uORF genes" = nt_nmd_strict
)


# ----------------------------- 11. Build group-specific data ------------

make_group_list <- function(hs) {
  
  grp <- list(
    "All HSR classes" = hs
  )
  
  
  for (cn in c("Transient (acute)", "Sustained", "Delayed")) {
    
    grp[[cn]] <- hs %>%
      filter(
        hsr_class == cn
      )
    
    
    grp[[paste0(cn, ", ChIP-bound")]] <- hs %>%
      filter(
        hsr_class == cn,
        status == "ChIP-bound"
      )
    
    
    grp[[paste0(cn, ", not bound")]] <- hs %>%
      filter(
        hsr_class == cn,
        status == "not bound"
      )
  }
  
  
  grp[nt_grp]
}


nt_sens_res <- lapply(
  names(nt_variants),
  function(vn) {
    
    excluded <- nt_variants[[vn]]
    
    
    # The same exclusion set is applied to all displayed HSR groups.
    hs <- nt_hsr %>%
      filter(
        !gene_id %in% excluded
      )
    
    
    removed_n <- length(
      intersect(
        unique(nt_hsr$gene_id),
        excluded
      )
    )
    
    
    vlab <- paste0(
      vn,
      "\n(",
      removed_n,
      " of ",
      n_distinct(nt_hsr$gene_id),
      " HSR genes removed)"
    )
    
    
    grp <- make_group_list(
      hs
    )
    
    
    # Summary table: medians are computed from the un-clipped index.
    # Percent labels in the figure are these medians multiplied by 100.
    summary_dat <- bind_rows(
      lapply(
        names(grp),
        function(gn) {
          
          d <- grp[[gn]] %>%
            filter(
              is.finite(priming_index)
            )
          
          
          tibble(
            variant = vlab,
            group = gn,
            n = nrow(d),
            median_priming_index = if (nrow(d) > 0L) {
              median(d$priming_index, na.rm = TRUE)
            } else {
              NA_real_
            }
          )
        }
      )
    )
    
    
    plot_dat <- bind_rows(
      lapply(
        names(grp),
        function(gn) {
          
          d <- grp[[gn]]
          
          if (nrow(d) == 0L) {
            return(NULL)
          }
          
          d %>%
            mutate(
              group = gn,
              variant = vlab
            )
        }
      )
    )
    
    
    list(
      tab = summary_dat,
      dat = plot_dat
    )
  }
)


nt_sens <- bind_rows(
  lapply(
    nt_sens_res,
    `[[`,
    "tab"
  )
)


nt_sens_dat <- bind_rows(
  lapply(
    nt_sens_res,
    `[[`,
    "dat"
  )
)


if (nrow(nt_sens_dat) == 0L) {
  stop("No grouped priming-index values remain for plotting.")
}


# Preserve the requested variant and group ordering.
variant_labels <- vapply(
  names(nt_variants),
  function(vn) {
    
    excluded <- nt_variants[[vn]]
    
    removed_n <- length(
      intersect(
        unique(nt_hsr$gene_id),
        excluded
      )
    )
    
    paste0(
      vn,
      "\n(",
      removed_n,
      " of ",
      n_distinct(nt_hsr$gene_id),
      " HSR genes removed)"
    )
  },
  character(1)
)


nt_sens$variant <- factor(
  nt_sens$variant,
  levels = unname(variant_labels)
)


nt_sens_dat$variant <- factor(
  nt_sens_dat$variant,
  levels = unname(variant_labels)
)


nt_sens$group <- factor(
  nt_sens$group,
  levels = nt_grp
)


nt_sens_dat$group <- factor(
  nt_sens_dat$group,
  levels = nt_grp
)


nt_sens <- nt_sens %>%
  mutate(
    lab = ifelse(
      is.finite(median_priming_index),
      paste0(
        round(100 * median_priming_index),
        "%"
      ),
      "n/a"
    )
  )


nt_sens_dat <- nt_sens_dat %>%
  mutate(
    cls = case_when(
      group == "All HSR classes" ~ "All HSR classes",
      grepl("^Transient", as.character(group)) ~ "Transient (acute)",
      grepl("^Sustained", as.character(group)) ~ "Sustained",
      grepl("^Delayed", as.character(group)) ~ "Delayed",
      TRUE ~ NA_character_
    ),
    
    priming_index_plot = pmax(
      pmin(
        priming_index,
        priming_clip_max
      ),
      priming_clip_min
    )
  )


# ----------------------------- 12. Export tables ------------------------

write.csv(
  nt_sens,
  paste0(
    output_prefix,
    "_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  nt_sens_dat,
  paste0(
    output_prefix,
    "_gene_level_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  feature_table,
  paste0(
    output_prefix,
    "_gene_feature_flags.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    gene_id = sort(unique(nt_nmd_emp))
  ),
  paste0(
    output_prefix,
    "_empirical_NMD_target_set.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    gene_id = sort(unique(nt_nmd_strict))
  ),
  paste0(
    output_prefix,
    "_strict_NMD_structural_exclusion_set.csv"
  ),
  row.names = FALSE
)

variant_audit <- tibble(
  variant = unname(variant_labels),
  n_HSR_genes_before = n_distinct(nt_hsr$gene_id),
  n_HSR_genes_removed = vapply(
    names(nt_variants),
    function(vn) {
      length(
        intersect(
          unique(nt_hsr$gene_id),
          nt_variants[[vn]]
        )
      )
    },
    integer(1)
  ),
  n_HSR_genes_remaining = n_HSR_genes_before - n_HSR_genes_removed
)


write.csv(
  variant_audit,
  paste0(
    output_prefix,
    "_variant_gene_counts.csv"
  ),
  row.names = FALSE
)


message("")
message("Priming-index sensitivity summary:")
print(
  nt_sens %>%
    mutate(
      variant = gsub("\n", " ", as.character(variant))
    ) %>%
    dplyr::select(
      variant,
      group,
      n,
      median_priming_index
    )
)


message("")
message("Variant gene counts:")
print(
  variant_audit %>%
    mutate(
      variant = gsub("\n", " ", variant)
    )
)


# ----------------------------- 13. Plot priming index -------------------

p_sens_b <- ggplot(
  nt_sens_dat,
  aes(
    x = group,
    y = priming_index_plot,
    fill = cls
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    colour = "grey25",
    linewidth = 0.65
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_text(
    data = nt_sens,
    aes(
      x = group,
      y = 1.40,
      label = lab
    ),
    inherit.aes = FALSE,
    size = 4.8,
    colour = "black",
    na.rm = TRUE
  ) +
  facet_grid(
    variant ~ .,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = nt_cols,
    guide = "none",
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = nt_grp,
    labels = nt_x_lab,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(
      priming_clip_min,
      priming_clip_max
    ),
    breaks = seq(
      priming_clip_min,
      priming_clip_max,
      by = 0.5
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  labs(
    x = NULL,
    y = "Priming index (clipped to -0.5 ... 1.5)",
    title = "Priming index with and without NMD targets",
    subtitle = paste0(
      "Index = upf3/Col-0 log2FC at NT / WT induction ",
      "(maximum of P/T and P/T +1d); ",
      "dotted line = full heat level; numbers = median"
    )
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 10.5,
      colour = "black",
      lineheight = 0.95,
      hjust = 0.5
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      margin = margin(r = 8)
    ),
    
    strip.text.y.right = element_text(
      size = 11,
      face = "bold",
      angle = 270
    ),
    
    strip.background = element_rect(
      fill = "grey90",
      colour = "grey60",
      linewidth = 0.5
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.3
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      margin = margin(b = 4)
    ),
    
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 7)
    ),
    
    panel.spacing.y = grid::unit(
      0.18,
      "lines"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


print(p_sens_b)


# ----------------------------- 14. Save figure --------------------------

pdf_file <- paste0(
  output_prefix,
  ".pdf"
)

png_file <- paste0(
  output_prefix,
  ".png"
)


ggsave(
  filename = pdf_file,
  plot = p_sens_b,
  width = 15.5,
  height = 9,
  units = "in",
  bg = "white"
)


ggsave(
  filename = png_file,
  plot = p_sens_b,
  width = 15.5,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)


message("")
message("NT priming-index robustness analysis completed.")
message("HSR class source: ", hsr_source)
message("ChIP target source: ", chip_source)
message("Empirical NMD target source: ", nmd_target_source)
message("Valid HSR genes with priming index: ", nrow(nt_hsr))
message("")
message("Saved figure: ", pdf_file)
message("Saved figure: ", png_file)
message("Saved summary table: ",
        paste0(output_prefix, "_summary.csv"))
message("Saved gene-level data: ",
        paste0(output_prefix, "_gene_level_data.csv"))
message("Saved variant counts: ",
        paste0(output_prefix, "_variant_gene_counts.csv"))