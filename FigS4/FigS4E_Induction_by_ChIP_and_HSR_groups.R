# =============================================================================
# Induction of HSR classes and ChIP groups in WT and upf3
#
# Main question:
#   For each gene group, how does expression change after heat treatment
#   in WT and upf3 relative to the WT NT baseline?
#
# Groups:
#   1. Transient (acute), ChIP-bound
#   2. Transient (acute), not bound
#   3. Sustained, ChIP-bound
#   4. Sustained, not bound
#   5. Delayed, ChIP-bound
#   6. Delayed, not bound
#   7. ChIP only (no HSR class)
#
# Conditions:
#   LAT = P/T
#   REC = recovery
#
# Statistics:
#   Paired Wilcoxon test across gene-level WT and upf3 expression changes.
#   BH correction across all available group-by-condition tests.
#
# Expression source:
#   Preferred: existing q_smp, expr_rel, q_group_sets and q_grp_levels.
#   Fallback: current DESeq2 object dds, using log2(normalized count + 1),
#             centered per gene relative to its WT NT mean.
# =============================================================================


# ----------------------------- 0. Packages -----------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})


# ----------------------------- 1. Settings -----------------------------

alpha_cutoff <- 0.05

output_prefix <- "HSR_groups_induction_Col-0_vs_upf3"

conditions_to_plot <- c("LAT", "REC")

condition_labels <- c(
  "LAT" = "P/T",
  "REC" = "REC"
)

genotype_levels <- c("Col-0", "upf3")


# ----------------------------- 2. Helper functions --------------------

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
  
  unique(
    x[!is.na(x) & nzchar(x)]
  )
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
    z %in% c("wt", "col", "col0", "wildtype") ~ "Col-0",
    z == "upf3" ~ "upf3",
    TRUE ~ NA_character_
  )
}


normalise_condition <- function(x) {
  z <- toupper(
    gsub("[[:space:]]+", "", trimws(as.character(x)))
  )
  
  dplyr::case_when(
    z %in% c("NT", "NOTREATMENT", "UNTREATED") ~ "NT",
    
    z %in% c("LAT", "PT", "P/T", "HEAT", "HEATSTRESS") ~ "LAT",
    
    z %in% c(
      "REC", "PT+1D", "P/T+1D", "RECOVERY", "RECOVERED"
    ) ~ "REC",
    
    TRUE ~ NA_character_
  )
}


safe_median <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  median(x)
}


safe_paired_wilcox <- function(wt, upf3_values) {
  ok <- is.finite(wt) & is.finite(upf3_values)
  
  wt <- wt[ok]
  upf3_values <- upf3_values[ok]
  
  if (length(wt) < 2L) {
    return(NA_real_)
  }
  
  if (all(wt == upf3_values)) {
    return(1)
  }
  
  tryCatch(
    suppressWarnings(
      stats::wilcox.test(
        upf3_values,
        wt,
        paired = TRUE,
        exact = FALSE
      )$p.value
    ),
    error = function(e) NA_real_
  )
}


p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}


format_group_axis <- function(x) {
  x <- sub(
    ", not bound$",
    "\n(not bound)",
    x
  )
  
  x <- sub(
    ", ChIP-bound$",
    "\n(ChIP-bound)",
    x
  )
  
  x <- sub(
    " \\(no HSR class\\)$",
    "",
    x
  )
  
  x
}


# ----------------------------- 3. Choose expression source --------------

required_q_objects <- c(
  "q_smp",
  "expr_rel",
  "q_group_sets",
  "q_grp_levels"
)

q_objects_available <- all(
  vapply(
    required_q_objects,
    function(nm) exists(nm, inherits = TRUE),
    logical(1)
  )
)


if (q_objects_available) {
  
  q_smp_obj <- get("q_smp", inherits = TRUE)
  expr_mat <- as.matrix(
    get("expr_rel", inherits = TRUE)
  )
  
  q_group_sets_obj <- get(
    "q_group_sets",
    inherits = TRUE
  )
  
  q_grp_levels_obj <- as.character(
    get("q_grp_levels", inherits = TRUE)
  )
  
  
  if (!is.data.frame(q_smp_obj)) {
    stop("'q_smp' must be a data frame.")
  }
  
  
  if (
    !all(
      c("sample", "genotype", "condition") %in%
      names(q_smp_obj)
    )
  ) {
    stop(
      "'q_smp' must contain sample, genotype and condition columns."
    )
  }
  
  
  if (
    is.null(rownames(expr_mat)) ||
    is.null(colnames(expr_mat))
  ) {
    stop(
      "'expr_rel' must have gene IDs as row names and sample IDs ",
      "as column names."
    )
  }
  
  
  if (!is.numeric(expr_mat)) {
    suppressWarnings(
      storage.mode(expr_mat) <- "numeric"
    )
  }
  
  
  if (anyDuplicated(clean_gene_id(rownames(expr_mat)))) {
    stop(
      "Duplicate gene IDs after identifier cleanup in expr_rel."
    )
  }
  
  
  rownames(expr_mat) <- clean_gene_id(
    rownames(expr_mat)
  )
  
  
  expression_source <- "existing expr_rel expression layer"
  
  
  expression_layer_label <- if (
    exists("q_layer", inherits = TRUE)
  ) {
    as.character(get("q_layer", inherits = TRUE))[1]
  } else {
    "existing relative-expression layer"
  }
  
  
  sample_info <- q_smp_obj %>%
    transmute(
      sample = as.character(sample),
      genotype = normalise_genotype(genotype),
      condition = normalise_condition(condition)
    )
  
  
  # q_smp should already describe the samples represented in expr_rel.
  if (anyDuplicated(sample_info$sample)) {
    stop("Duplicate sample IDs in q_smp.")
  }
  
  
  if (anyNA(sample_info$genotype) || anyNA(sample_info$condition)) {
    stop(
      "Some q_smp genotype or condition labels were not recognised. ",
      "Check normalise_genotype() and normalise_condition()."
    )
  }
  
  
  missing_expr_samples <- setdiff(
    sample_info$sample,
    colnames(expr_mat)
  )
  
  
  if (length(missing_expr_samples) > 0L) {
    stop(
      "Samples in q_smp are absent from expr_rel: ",
      paste(head(missing_expr_samples, 10), collapse = ", ")
    )
  }
  
  
  # Use the existing group definitions from the current workflow.
  ind_groups <- c(
    q_grp_levels_obj[
      grepl("bound", q_grp_levels_obj, ignore.case = TRUE)
    ],
    "ChIP only (no HSR class)"
  )
  
  ind_groups <- unique(ind_groups)
  
  
  if (!is.list(q_group_sets_obj) ||
      is.null(names(q_group_sets_obj))) {
    stop(
      "'q_group_sets' must be a named list of gene vectors."
    )
  }
  
  
  missing_group_sets <- setdiff(
    ind_groups,
    names(q_group_sets_obj)
  )
  
  
  if (length(missing_group_sets) > 0L) {
    stop(
      "Some requested groups are absent from q_group_sets: ",
      paste(missing_group_sets, collapse = ", ")
    )
  }
  
  
  group_sets <- lapply(
    q_group_sets_obj[ind_groups],
    clean_gene_set
  )
  
  
  names(group_sets) <- ind_groups
  
  
  message("Expression source: ", expression_source)
  message("Expression layer: ", expression_layer_label)
  
  
} else {
  
  # ---------------------------------------------------------------------
  # Fallback path: build the same data from the current DESeq2 object.
  # ---------------------------------------------------------------------
  
  if (!exists("dds", inherits = TRUE)) {
    stop(
      "Neither the complete q_* object set nor 'dds' is available. ",
      "Run the expression/HSR/ChIP workflow first."
    )
  }
  
  
  if (!requireNamespace("DESeq2", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "The fallback expression path requires DESeq2 and ",
      "SummarizedExperiment."
    )
  }
  
  
  dds_obj <- get("dds", inherits = TRUE)
  
  
  if (!inherits(dds_obj, "DESeqDataSet")) {
    stop("'dds' must be a DESeqDataSet.")
  }
  
  
  cd <- as.data.frame(
    SummarizedExperiment::colData(dds_obj)
  ) %>%
    rownames_to_column("sample")
  
  
  if (!all(c("genotype", "treatment") %in% names(cd))) {
    stop(
      "colData(dds) must contain genotype and treatment columns."
    )
  }
  
  
  norm_counts <- DESeq2::counts(
    dds_obj,
    normalized = TRUE
  )
  
  
  rownames(norm_counts) <- clean_gene_id(
    rownames(norm_counts)
  )
  
  
  if (anyDuplicated(rownames(norm_counts))) {
    stop(
      "Duplicate gene IDs in the normalised count matrix after cleanup."
    )
  }
  
  
  if (!all(cd$sample %in% colnames(norm_counts))) {
    stop(
      "Sample names in colData(dds) do not match the count matrix."
    )
  }
  
  
  # Use log2(normalized count + 1), consistently with the current
  # five-gene-set expression workflow.
  expr_mat <- log2(
    norm_counts + 1
  )
  
  
  sample_info <- cd %>%
    transmute(
      sample = as.character(sample),
      genotype = normalise_genotype(genotype),
      condition = normalise_condition(treatment)
    )
  
  
  # This plot requires WT and upf3 samples at NT, LAT and REC.
  relevant_samples <- sample_info %>%
    filter(
      genotype %in% c("Col-0", "upf3"),
      condition %in% c("NT", "LAT", "REC")
    )
  
  
  if (nrow(relevant_samples) == 0L) {
    stop(
      "No usable Col-0/upf3 samples with NT, LAT or REC treatments found."
    )
  }
  
  
  sample_info <- relevant_samples
  
  
  if (anyDuplicated(sample_info$sample)) {
    stop("Duplicate sample IDs in colData(dds).")
  }
  
  
  if (
    anyNA(sample_info$genotype) ||
    anyNA(sample_info$condition)
  ) {
    stop(
      "Unrecognised genotype/treatment labels remain in the selected ",
      "dds sample metadata."
    )
  }
  
  
  missing_expr_samples <- setdiff(
    sample_info$sample,
    colnames(expr_mat)
  )
  
  
  if (length(missing_expr_samples) > 0L) {
    stop(
      "Some selected samples are absent from the expression matrix."
    )
  }
  
  
  expression_source <- "dds normalized counts"
  
  expression_layer_label <- "log2(normalized count + 1)"
  
  
  # Try the HSR classes in the current workflow first.
  hsr_classes <- NULL
  
  
  if (exists("gene_sets", inherits = TRUE)) {
    
    gs <- get("gene_sets", inherits = TRUE)
    
    
    if (
      is.data.frame(gs) &&
      all(c("gene", "GeneGroup") %in% names(gs))
    ) {
      
      gs_labels <- as.character(gs$GeneGroup)
      gs_norm <- tolower(
        gsub("[^[:alnum:]]+", "", gs_labels)
      )
      
      
      hsr_classes <- list(
        Acute = clean_gene_set(
          gs$gene[
            !is.na(gs_norm) &
              grepl("acute|transient", gs_norm)
          ]
        ),
        
        Sustained = clean_gene_set(
          gs$gene[
            !is.na(gs_norm) &
              grepl("sustain", gs_norm)
          ]
        ),
        
        Delayed = clean_gene_set(
          gs$gene[
            !is.na(gs_norm) &
              grepl("delayed", gs_norm)
          ]
        )
      )
      
      
      if (any(lengths(hsr_classes) == 0L)) {
        hsr_classes <- NULL
      }
    }
  }
  
  
  # Alternative source: gene_classification flags.
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
        c("gene", "is_acute", "is_sustained", "is_delayed") %in%
        names(gc)
      )
    ) {
      
      hsr_classes <- list(
        Acute = clean_gene_set(
          gc$gene[
            !is.na(gc$is_acute) & gc$is_acute
          ]
        ),
        
        Sustained = clean_gene_set(
          gc$gene[
            !is.na(gc$is_sustained) & gc$is_sustained
          ]
        ),
        
        Delayed = clean_gene_set(
          gc$gene[
            !is.na(gc$is_delayed) & gc$is_delayed
          ]
        )
      )
      
    } else if (
      is.data.frame(gc) &&
      all(
        c(
          "gene",
          "is_Col-0_LAT_induced",
          "is_Col-0_REC_induced"
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
    }
  }
  
  
  # Alternative source: named hsr_classes list.
  if (
    is.null(hsr_classes) &&
    exists("hsr_classes", inherits = TRUE)
  ) {
    
    hs_tmp <- get(
      "hsr_classes",
      inherits = TRUE
    )
    
    
    if (
      is.list(hs_tmp) &&
      all(c("Acute", "Sustained", "Delayed") %in% names(hs_tmp))
    ) {
      hsr_classes <- lapply(
        hs_tmp[c("Acute", "Sustained", "Delayed")],
        clean_gene_set
      )
    }
  }
  
  
  if (is.null(hsr_classes)) {
    stop(
      "Could not recover the Acute, Sustained and Delayed HSR gene sets. ",
      "Run the current HSR gene-set workflow first."
    )
  }
  
  
  # Recover the broad HsfA2/HsfA3 ChIP targets.
  chip_targets <- NULL
  
  
  if (exists("hsfa23_all_relevant", inherits = TRUE)) {
    
    chip_targets <- extract_agi_ids(
      get("hsfa23_all_relevant", inherits = TRUE)
    )
  }
  
  
  if (
    (is.null(chip_targets) || length(chip_targets) == 0L) &&
    exists("chip_all", inherits = TRUE)
  ) {
    
    chip_object <- get("chip_all", inherits = TRUE)
    
    if (is.atomic(chip_object) && is.null(dim(chip_object))) {
      chip_targets <- extract_agi_ids(chip_object)
    }
  }
  
  
  if (
    (is.null(chip_targets) || length(chip_targets) == 0L) &&
    exists("chip_clusters", inherits = TRUE)
  ) {
    
    chip_df <- as.data.frame(
      get("chip_clusters", inherits = TRUE)
    )
    
    
    if (ncol(chip_df) < 2L) {
      stop(
        "The existing chip_clusters object needs gene and cluster columns."
      )
    }
    
    
    name_norm <- tolower(
      gsub("[^[:alnum:]]+", "", names(chip_df))
    )
    
    
    gene_col <- which(
      name_norm %in% c("geneid", "gene", "tair", "agi")
    )[1]
    
    
    cluster_col <- which(
      name_norm %in% c("cluster", "chipcluster", "clusterid")
    )[1]
    
    
    if (is.na(gene_col)) gene_col <- 1L
    if (is.na(cluster_col)) cluster_col <- 2L
    
    
    cluster_text <- as.character(
      chip_df[[cluster_col]]
    )
    
    
    cluster_num <- ifelse(
      grepl("[0-9]+", cluster_text),
      sub(
        "^[^0-9]*([0-9]+).*",
        "\\1",
        cluster_text
      ),
      NA_character_
    )
    
    
    cluster_num <- suppressWarnings(
      as.integer(cluster_num)
    )
    
    
    chip_targets <- extract_agi_ids(
      chip_df[[gene_col]][
        cluster_num %in% c(7, 11, 12, 14)
      ]
    )
  }
  
  
  if (
    (is.null(chip_targets) || length(chip_targets) == 0L) &&
    file.exists("13059_2023_2970_MOESM3_ESM.xlsx")
  ) {
    
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(
        "Package readxl is needed to read the ChIP Excel file. ",
        "Install it using install.packages('readxl')."
      )
    }
    
    
    chip_df <- as.data.frame(
      readxl::read_excel(
        "13059_2023_2970_MOESM3_ESM.xlsx",
        sheet = "gene-clusters-with-short-descri"
      )
    )
    
    
    if (ncol(chip_df) < 2L) {
      stop("The ChIP Excel sheet must have at least two columns.")
    }
    
    
    cluster_text <- as.character(chip_df[[2]])
    
    
    cluster_num <- ifelse(
      grepl("[0-9]+", cluster_text),
      sub(
        "^[^0-9]*([0-9]+).*",
        "\\1",
        cluster_text
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
  }
  
  
  if (
    is.null(chip_targets) ||
    length(chip_targets) == 0L
  ) {
    stop(
      "Could not recover the broad HsfA2/HsfA3 ChIP targets. ",
      "Load hsfa23_all_relevant/chip_all or make the ChIP Excel file available."
    )
  }
  
  
  chip_targets <- clean_gene_set(chip_targets)
  
  
  # Build the seven mutually exclusive groups.
  hsr_union <- unique(
    unlist(hsr_classes, use.names = FALSE)
  )
  
  
  group_sets <- list(
    
    "Transient (acute), ChIP-bound" = intersect(
      hsr_classes$Acute,
      chip_targets
    ),
    
    "Transient (acute), not bound" = setdiff(
      hsr_classes$Acute,
      chip_targets
    ),
    
    "Sustained, ChIP-bound" = intersect(
      hsr_classes$Sustained,
      chip_targets
    ),
    
    "Sustained, not bound" = setdiff(
      hsr_classes$Sustained,
      chip_targets
    ),
    
    "Delayed, ChIP-bound" = intersect(
      hsr_classes$Delayed,
      chip_targets
    ),
    
    "Delayed, not bound" = setdiff(
      hsr_classes$Delayed,
      chip_targets
    ),
    
    "ChIP only (no HSR class)" = setdiff(
      chip_targets,
      hsr_union
    )
  )
  
  
  group_sets <- lapply(
    group_sets,
    clean_gene_set
  )
  
  
  ind_groups <- names(group_sets)
  
  
  message("Expression source: ", expression_source)
  message("Expression layer: ", expression_layer_label)
}


# ----------------------------- 4. Validate group definitions ------------

if (length(ind_groups) != 7L) {
  stop(
    "Expected seven groups (six HSR x ChIP-bound/not-bound groups ",
    "plus ChIP-only). Found: ",
    length(ind_groups)
  )
}


if (anyDuplicated(ind_groups)) {
  stop("Duplicate group names were detected.")
}


group_sets <- lapply(
  group_sets[ind_groups],
  clean_gene_set
)


names(group_sets) <- ind_groups


group_sizes_raw <- tibble(
  group = ind_groups,
  n_genes_defined = vapply(
    group_sets,
    length,
    integer(1)
  )
)


message("")
message("Gene counts before expression filtering:")
print(group_sizes_raw)


# Check group disjointness. No gene should appear in more than one
# of the seven plotted groups.
membership_for_check <- bind_rows(
  lapply(
    ind_groups,
    function(gn) {
      tibble(
        group = gn,
        gene_id = group_sets[[gn]]
      )
    }
  )
)


overlapping_genes <- membership_for_check %>%
  dplyr::count(gene_id, name = "n_groups") %>%
  filter(n_groups > 1L)


if (nrow(overlapping_genes) > 0L) {
  stop(
    "The seven plotted gene groups are not mutually exclusive. ",
    "Examples: ",
    paste(
      head(overlapping_genes$gene_id, 10),
      collapse = ", "
    ),
    ". Check the HSR and ChIP gene-set definitions."
  )
}


# ----------------------------- 5. Prepare expression means ---------------

# The required per-gene group means are:
#   WT|NT    baseline
#   WT|LAT   P/T
#   upf3|LAT P/T
#   WT|REC   recovery
#   upf3|REC recovery

required_keys <- c(
  "Col-0|NT",
  "Col-0|LAT",
  "upf3|LAT",
  "Col-0|REC",
  "upf3|REC"
)


missing_samples_by_group <- bind_rows(
  lapply(
    required_keys,
    function(key) {
      
      parts <- strsplit(
        key,
        split = "|",
        fixed = TRUE
      )[[1]]
      
      
      n_samples <- sum(
        sample_info$genotype == parts[1] &
          sample_info$condition == parts[2]
      )
      
      
      tibble(
        sample_group = key,
        n_samples = n_samples
      )
    }
  )
)


message("")
message("Samples per required genotype-condition combination:")
print(missing_samples_by_group)


if (any(missing_samples_by_group$n_samples == 0L)) {
  stop(
    "At least one required genotype-condition combination has no samples. ",
    "Check the current sample metadata."
  )
}


# Average the supplied expression layer over biological replicates for
# each genotype-condition combination, gene by gene.
group_means <- matrix(
  NA_real_,
  nrow = nrow(expr_mat),
  ncol = length(required_keys),
  dimnames = list(
    rownames(expr_mat),
    required_keys
  )
)


for (key in required_keys) {
  
  parts <- strsplit(
    key,
    split = "|",
    fixed = TRUE
  )[[1]]
  
  
  selected_samples <- sample_info$sample[
    sample_info$genotype == parts[1] &
      sample_info$condition == parts[2]
  ]
  
  
  group_means[, key] <- rowMeans(
    expr_mat[
      ,
      selected_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  )
}


# If expr_rel is supplied, the reference workflow assumes it is already
# expressed relative to WT NT. Therefore, do not subtract the baseline
# a second time.
#
# If the dds fallback is used, expr_mat is log2(normalized count + 1);
# subtract the WT NT mean separately for each gene.
if (!q_objects_available) {
  
  baseline <- group_means[, "Col-0|NT"]
  
  change_mat <- sweep(
    group_means,
    MARGIN = 1,
    STATS = baseline,
    FUN = "-"
  )
  
} else {
  
  change_mat <- group_means
}


# Determine which genes are eligible for plotting.
if (
  q_objects_available &&
  exists("expr_ok", inherits = TRUE)
) {
  
  eligible_genes <- clean_gene_set(
    get("expr_ok", inherits = TRUE)
  )
  
} else if (!q_objects_available) {
  
  # Fallback expression filter:
  # retain genes with a positive DESeq2-normalized count in at least
  # one available sample.
  eligible_genes <- rownames(expr_mat)[
    rowSums(expr_mat > 0, na.rm = TRUE) > 0
  ] %>%
    clean_gene_set()
  
} else {
  
  # If expr_ok is unavailable, rely on genes represented in expr_rel.
  eligible_genes <- clean_gene_set(
    rownames(expr_mat)
  )
}


eligible_genes <- intersect(
  eligible_genes,
  rownames(change_mat)
)


if (length(eligible_genes) == 0L) {
  stop(
    "No expression-eligible genes match the expression matrix."
  )
}


# ----------------------------- 6. Build per-gene data --------------------

ind_long <- bind_rows(
  lapply(
    ind_groups,
    function(gn) {
      
      group_genes <- intersect(
        group_sets[[gn]],
        eligible_genes
      )
      
      
      if (length(group_genes) == 0L) {
        warning(
          "Group '", gn,
          "' contains no expression-eligible genes."
        )
        
        return(tibble())
      }
      
      
      bind_rows(
        lapply(
          conditions_to_plot,
          function(cc) {
            
            tibble(
              group = gn,
              gene_id = group_genes,
              condition = cc,
              WT = change_mat[
                group_genes,
                paste0("Col-0|", cc)
              ],
              upf3 = change_mat[
                group_genes,
                paste0("upf3|", cc)
              ]
            )
          }
        )
      )
    }
  )
)


if (nrow(ind_long) == 0L) {
  stop("No per-gene expression changes could be calculated.")
}


ind_long <- ind_long %>%
  mutate(
    group = factor(
      group,
      levels = ind_groups
    ),
    
    condition = factor(
      condition,
      levels = conditions_to_plot
    )
  )


# Export the per-gene values used in the figure.
write.csv(
  ind_long,
  paste0(
    output_prefix,
    "_gene_level_changes.csv"
  ),
  row.names = FALSE
)


# Exact memberships after the expression-eligibility filter.
group_membership_used <- bind_rows(
  lapply(
    ind_groups,
    function(gn) {
      tibble(
        group = gn,
        gene_id = intersect(
          group_sets[[gn]],
          eligible_genes
        )
      )
    }
  )
)


write.csv(
  group_membership_used,
  paste0(
    output_prefix,
    "_group_membership_used.csv"
  ),
  row.names = FALSE
)


group_sizes_used <- group_membership_used %>%
  dplyr::count(
    group,
    name = "n_genes_used"
  ) %>%
  right_join(
    tibble(
      group = factor(
        ind_groups,
        levels = ind_groups
      )
    ),
    by = "group"
  ) %>%
  mutate(
    n_genes_used = replace_na(
      n_genes_used,
      0L
    )
  ) %>%
  arrange(group)


write.csv(
  group_sizes_used,
  paste0(
    output_prefix,
    "_group_sizes.csv"
  ),
  row.names = FALSE
)


message("")
message("Expression-eligible gene counts by group:")
print(group_sizes_used)


# ----------------------------- 7. Paired Wilcoxon tests ----------------

ind_tests <- ind_long %>%
  group_by(
    group,
    condition
  ) %>%
  summarise(
    
    n = sum(
      is.finite(WT) &
        is.finite(upf3)
    ),
    
    median_WT = safe_median(
      WT[
        is.finite(WT) &
          is.finite(upf3)
      ]
    ),
    
    median_upf3 = safe_median(
      upf3[
        is.finite(WT) &
          is.finite(upf3)
      ]
    ),
    
    median_diff_upf3_minus_WT = safe_median(
      upf3[
        is.finite(WT) &
          is.finite(upf3)
      ] -
        WT[
          is.finite(WT) &
            is.finite(upf3)
        ]
    ),
    
    p = safe_paired_wilcox(
      WT,
      upf3
    ),
  
    .groups = "drop"
  ) %>%
  mutate(
    # BH correction across all group x condition comparisons in this panel.
    p_adj = p.adjust(
      p,
      method = "BH"
    ),
    
    sig = p_to_stars(p_adj),
    
    group = factor(
      group,
      levels = ind_groups
    ),
    
    condition = factor(
      condition,
      levels = conditions_to_plot
    )
  )


write.csv(
  ind_tests,
  paste0(
    output_prefix,
    "_tests.csv"
  ),
  row.names = FALSE
)


cat(
  "\nInduction relative to Col-0 NT: Col-0 versus upf3\n",
  "Paired Wilcoxon test across gene-level values; BH adjustment\n",
  sep = ""
)


print(
  as.data.frame(
    ind_tests %>%
      mutate(
        across(
          where(is.numeric),
          ~ signif(.x, 3)
        )
      )
  ),
  row.names = FALSE
)


# ----------------------------- 8. Prepare the plot ----------------------

ind_plot <- ind_long %>%
  pivot_longer(
    cols = c(WT, upf3),
    names_to = "genotype",
    values_to = "change"
  ) %>%
  filter(
    is.finite(change)
  ) %>%
  mutate(
    group = factor(
      group,
      levels = ind_groups
    ),
    
    genotype = factor(
      genotype,
      levels = genotype_levels
    ),
    
    condition = factor(
      condition,
      levels = conditions_to_plot
    )
  )


if (nrow(ind_plot) == 0L) {
  stop("No finite per-gene values are available for the plot.")
}


# Use the 99th percentile to keep extreme gene-level values from
# excessively compressing the boxes, as in the reference code.
q99 <- as.numeric(
  stats::quantile(
    ind_plot$change,
    probs = 0.99,
    na.rm = TRUE
  )
)


if (!is.finite(q99)) {
  stop("Could not calculate a finite upper plotting limit.")
}


ind_ylim <- if (q99 > 0) {
  q99 * 1.15
} else {
  max(ind_plot$change, na.rm = TRUE) + 1
}


# Make sure annotations remain above the central part of the data.
ind_ylim <- max(
  ind_ylim,
  max(ind_plot$change, na.rm = TRUE) * 1.05,
  1
)


# If the original q_layer exists, preserve its label.
y_axis_layer <- if (q_objects_available) {
  expression_layer_label
} else {
  expression_layer_label
}


# ----------------------------- 9. Draw figure ---------------------------

genotype_colours <- c(
  "Col-0" = "grey65",
  "upf3" = "#E3120B"
)


p_ind <- ggplot(
  ind_plot,
  aes(
    x = group,
    y = change,
    fill = genotype
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    position = position_dodge(0.75)
  ) +
  geom_text(
    data = ind_tests,
    aes(
      x = group,
      y = ind_ylim * 0.97,
      label = sig
    ),
    inherit.aes = FALSE,
    size = 10,
    fontface = "bold"
  ) +
  facet_wrap(
    ~condition,
    ncol = 1,
    labeller = as_labeller(condition_labels)
  ) +
  scale_fill_manual(
    values = genotype_colours,
    breaks = genotype_levels,
    name = NULL,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = ind_groups,
    labels = format_group_axis,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(NA, ind_ylim)
  ) +
  labs(
    x = NULL,
    
    y = paste0(
      "log2 change vs Col-0 NT per gene (",
      y_axis_layer,
      ")"
    ),
    
    title = "Induction of the HSR classes and ChIP groups in Col-0 and in upf3",
    
    subtitle = paste0(
      "Boxes = genes; stars: upf3 vs Col-0 per group, gene-level paired ",
      "Wilcoxon (genes as units, liberal; replicate-level test is in ",
      "the 3-vs-3 figure), BH"
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.text = element_text( size = 15),
    axis.text.x = element_text(
      size = 15,
      colour = "black",
      lineheight = 0.95
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 15
    ),
    axis.title.y = element_text( size = 18),
    strip.background = element_rect(
      fill = "grey90"
    ),
    strip.text = element_text( size = 15, face = "bold"),
    plot.title = element_text(
      face = "bold",
      size = 20
    ),
    
    plot.subtitle = element_blank(),
    
    panel.grid.minor = element_blank()
  )


print(p_ind)


# ----------------------------- 10. Save outputs -------------------------

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
  plot = p_ind,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)


ggsave(
  filename = png_file,
  plot = p_ind,
  width = 17.5,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)


message("")
message("Q3d induction figure completed.")
message("Expression source: ", expression_source)
message("Expression layer: ", expression_layer_label)
message("")
message("Saved figure: ", pdf_file)
message("Saved figure: ", png_file)
message(
  "Saved per-gene values: ",
  paste0(output_prefix, "_gene_level_changes.csv")
)
message(
  "Saved tests: ",
  paste0(output_prefix, "_tests.csv")
)
message(
  "Saved group membership: ",
  paste0(output_prefix, "_group_membership_used.csv")
)
message(
  "Saved group sizes: ",
  paste0(output_prefix, "_group_sizes.csv")
)
message("")
message(
  "Interpretation: the paired Wilcoxon tests treat genes as the paired units, ",
  "not biological replicates. These tests describe gene-set-wide shifts but ",
  "do not replace genotype comparisons based on biological replicate samples."
)