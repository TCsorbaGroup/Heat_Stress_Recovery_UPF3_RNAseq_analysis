############################################################################
# Structural cis-feature comparison for the five Col-0 HSR gene sets
# Panels:
#   A) Maximum annotated 3'UTR length across isoforms
#   B) Percentage with a 3'UTR intron in the representative transcript.
#   C) Percentage with a RiboTaper-supported uORF.
#      No 3'UTR annotation filter is applied to this panel.
#
# Gene-set definitions:
#   - Existing current-session gene_sets or gene_classification objects
#     are used preferentially.
#   - Fallback reconstruction uses the current DESeq2 object dds.
#   - HS responsive genes = union of Acute, Sustained and Delayed.
#   - Other expressed genes = all genes in the analysis universe outside
#     the three HSR classes.
# Statistics:
#   - Panel A: Wilcoxon rank-sum test, each foreground set vs reference.
#   - Panels B/C: two-sided Fisher's exact test, each foreground set
#     vs reference.
#   - BH adjustment is applied separately within each feature family
#     across the four foreground-vs-reference comparisons.
#   - All plot significance labels use BH-adjusted p-values.
#
# uORF evaluability:
#   - RiboTaper/CiPS-positive genes stay positive even if TPM < threshold.
#   - A negative call is assigned only to genes with TPM >= ribo_tpm_min.
#   - Genes with no positive call and insufficient TPM remain NA/unassessed.
############################################################################


# ----------------------------- 0. Packages -----------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})


# ----------------------------- 1. Settings -----------------------------

alpha_cutoff <- 0.05
lfc_cutoff <- 1

promoter_length <- 1000L  # retained for provenance; not used in this analysis

annotation_candidates <- c(
  "Araport11_GFF3_genes_transposons.20241001.gff.gz",
  "Araport11_GFF3_genes_transposons.20230105.gff.gz"
)

ribotaper_file <- "Dataset S1_RiboTaper.xlsx"
cips_file <- "Dataset S5_CiPS.xlsx"
ribo_expr_file <- "GSE183264_Ribo.genes.results"

# Positive RiboTaper/CiPS evidence remains positive even below this threshold.
# The threshold determines only whether a non-positive gene can be called
# evaluable negative.
ribo_tpm_min <- 1

# Significance stars MUST use BH-adjusted p-values.
significance_basis <- "BH"

# Odds/range and output settings.
out_prefix <- "Col_HSR_five_gene_sets_cis_feature_comparison"

group_order <- c(
  "Other expressed genes",
  "HS responsive genes",
  "Acute HS responsive genes",
  "Sustained HS responsive genes",
  "Delayed HS responsive genes"
)

reference_group <- "Other expressed genes"
foreground_groups <- setdiff(group_order, reference_group)

group_colors <- c(
  "Other expressed genes" = "grey75",
  "HS responsive genes" = "#6A3D9A",
  "Acute HS responsive genes" = "#E60012",
  "Sustained HS responsive genes" = "#164A91",
  "Delayed HS responsive genes" = "#26934F"
)

group_labels <- c(
  "Other expressed genes" = "Other expressed genes",
  "HS responsive genes" = "HS responsive genes",
  "Acute HS responsive genes" = "Acute HS responsive genes",
  "Sustained HS responsive genes" = "Sustained HS responsive genes",
  "Delayed HS responsive genes" = "Delayed HS responsive genes"
)


# ----------------------------- 2. Helpers ------------------------------

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


clean_tx_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub(
    "^(transcript|mRNA|rna):",
    "",
    x,
    ignore.case = TRUE
  )
  trimws(x)
}


clean_gene_set <- function(x) {
  x <- clean_gene_id(x)
  unique(x[!is.na(x) & nzchar(x)])
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


normalise_group_label <- function(x) {
  tolower(
    gsub("[^[:alnum:]]+", "", trimws(as.character(x)))
  )
}


# Recognise class labels
canonical_group_label <- function(x) {
  z <- normalise_group_label(x)
  
  if (is.na(z) || !nzchar(z)) {
    return(NA_character_)
  }
  
  if (grepl("otherexpress", z)) {
    return("Other expressed genes")
  }
  
  if (grepl("acute|transient", z)) {
    return("Acute HS responsive genes")
  }
  
  if (grepl("sustain", z)) {
    return("Sustained HS responsive genes")
  }
  
  if (grepl("delayed", z)) {
    return("Delayed HS responsive genes")
  }
  
  if (grepl("hsresponsive", z)) {
    return("HS responsive genes")
  }
  
  NA_character_
}


p_to_stars <- function(p) {
  if (length(p) == 0L || is.na(p) || !is.finite(p)) {
    return("n/a")
  }
  
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  
  "ns"
}


# Expand comma-separated Parent attributes and list-like rtracklayer fields
# so that every returned row has exactly one Parent value.
expand_parent_rows <- function(df) {
  if (nrow(df) == 0L) {
    return(df[0, , drop = FALSE])
  }
  
  if (!"Parent" %in% names(df)) {
    stop("The imported GFF3 has no Parent column.")
  }
  
  parent_values <- lapply(seq_len(nrow(df)), function(i) {
    
    value <- df$Parent[[i]]
    
    if (
      is.null(value) ||
      length(value) == 0L ||
      all(is.na(value))
    ) {
      return(character(0))
    }
    
    value <- paste(as.character(value), collapse = ",")
    
    vals <- trimws(
      unlist(
        strsplit(value, ",", fixed = TRUE),
        use.names = FALSE
      )
    )
    
    vals[nzchar(vals) & !is.na(vals)]
  })
  
  keep <- lengths(parent_values) > 0L
  
  if (!any(keep)) {
    return(df[0, , drop = FALSE])
  }
  
  source_rows <- which(keep)
  
  out <- df[
    rep(source_rows, lengths(parent_values)),
    ,
    drop = FALSE
  ]
  
  out$Parent <- unlist(
    parent_values[keep],
    use.names = FALSE
  )
  
  rownames(out) <- NULL
  out
}


extract_agi <- function(df, col = NULL) {
  
  if (!is.null(col)) {
    
    if (!col %in% names(df)) {
      stop(
        "Expected column '", col,
        "' not found. Available columns: ",
        paste(names(df), collapse = ", ")
      )
    }
    
    values <- as.character(df[[col]])
    
  } else {
    
    values <- unlist(
      lapply(df, function(x) {
        if (is.character(x) || is.factor(x)) {
          as.character(x)
        } else {
          character(0)
        }
      }),
      use.names = FALSE
    )
  }
  
  matches <- regmatches(
    values,
    gregexpr(
      "AT[1-5CM]G[0-9]{5}(?:\\.[0-9]+)?",
      values,
      ignore.case = TRUE,
      perl = TRUE
    )
  )
  
  clean_gene_set(
    unlist(matches, use.names = FALSE)
  )
}


read_excel_table <- function(path, sheet, skip = 0) {
  as.data.frame(
    readxl::read_excel(
      path,
      sheet = sheet,
      skip = skip
    )
  )
}


# ----------------------------- 3. Check inputs -------------------------

if (!exists("dds", inherits = TRUE)) {
  stop(
    "The DESeq2 object 'dds' is unavailable. ",
    "Load the current DESeq2 dataset first."
  )
}

if (!inherits(dds, "DESeqDataSet")) {
  stop("Object 'dds' must be a DESeqDataSet.")
}

if (
  length(significance_basis) != 1L ||
  significance_basis != "BH"
) {
  stop(
    "This workflow uses BH-adjusted p-values for every significance label. ",
    "Set significance_basis = 'BH'."
  )
}

if (
  length(ribo_tpm_min) != 1L ||
  !is.finite(ribo_tpm_min) ||
  ribo_tpm_min < 0
) {
  stop("ribo_tpm_min must be a non-negative finite number.")
}

found_annotations <- annotation_candidates[
  file.exists(annotation_candidates)
]

if (length(found_annotations) == 0L) {
  stop(
    "Could not find a GFF3 annotation file. Expected one of: ",
    paste(annotation_candidates, collapse = ", "),
    "\nWorking directory: ",
    getwd()
  )
}

annotation_file <- found_annotations[1]

required_files <- c(
  ribotaper_file,
  cips_file,
  ribo_expr_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required uORF input file(s) not found:\n - ",
    paste(missing_files, collapse = "\n - "),
    "\nWorking directory: ",
    getwd(),
    "\nThe Ribo-seq TPM table is required to identify evaluable negatives."
  )
}


# ----------------------------- 4. Load dds and analysis universe --------

norm_counts <- DESeq2::counts(
  dds,
  normalized = TRUE
)

count_gene_ids <- clean_gene_id(
  rownames(norm_counts)
)

if (anyNA(count_gene_ids) || any(!nzchar(count_gene_ids))) {
  stop(
    "The dds object contains missing or empty gene IDs after cleanup."
  )
}

if (anyDuplicated(count_gene_ids)) {
  dup <- unique(
    count_gene_ids[duplicated(count_gene_ids)]
  )
  
  stop(
    "Gene IDs in dds are duplicated after cleanup. Examples: ",
    paste(head(dup, 10), collapse = ", ")
  )
}

rownames(norm_counts) <- count_gene_ids


coldata_df <- as.data.frame(
  SummarizedExperiment::colData(dds)
) %>%
  tibble::rownames_to_column("Sample")


if (!all(c("genotype", "treatment") %in% names(coldata_df))) {
  stop(
    "colData(dds) must contain genotype and treatment columns."
  )
}


coldata_df <- coldata_df %>%
  mutate(
    Genotype = normalise_genotype(genotype)
  )


if (anyNA(coldata_df$Genotype)) {
  unknown <- unique(
    as.character(
      coldata_df$genotype[is.na(coldata_df$Genotype)]
    )
  )
  
  stop(
    "Unrecognised genotype label(s) in colData(dds): ",
    paste(unknown, collapse = ", "),
    ". Update normalise_genotype() if necessary."
  )
}


wt_samples <- intersect(
  coldata_df$Sample[
    coldata_df$Genotype == "Col-0"
  ],
  colnames(norm_counts)
)


if (length(wt_samples) == 0L) {
  stop("No Col-0 samples found in dds.")
}


wt_expressed_genes <- rownames(norm_counts)[
  rowSums(
    norm_counts[, wt_samples, drop = FALSE] > 0,
    na.rm = TRUE
  ) > 0
] %>%
  unique()


message("Col-0 samples: ", length(wt_samples))
message(
  "Genes expressed in at least one Col-0 sample: ",
  length(wt_expressed_genes)
)


# ----------------------------- 5. Import GFF and identify genes ----------

message("Importing GFF3 annotation: ", annotation_file)

gff <- rtracklayer::import(annotation_file)
gff_df <- as.data.frame(gff)


if (!all(c("type", "ID", "start", "end") %in% names(gff_df))) {
  stop(
    "GFF3 lacks expected type/ID/start/end fields."
  )
}


gff_gene_table <- gff_df %>%
  filter(
    tolower(type) == "gene"
  ) %>%
  transmute(
    gene = clean_gene_id(ID)
  ) %>%
  filter(
    !is.na(gene),
    nzchar(gene)
  ) %>%
  distinct()


if (nrow(gff_gene_table) == 0L) {
  stop("No gene features were found in the GFF3.")
}


if (anyDuplicated(gff_gene_table$gene)) {
  stop(
    "Duplicate gene identifiers remain in the GFF3 gene table after cleanup."
  )
}


gff_genes <- gff_gene_table$gene

dds_genes <- unique(
  clean_gene_id(rownames(norm_counts))
)


# This is the working universe for the five gene groups.
gene_universe <- Reduce(
  intersect,
  list(
    gff_genes,
    dds_genes,
    wt_expressed_genes
  )
)


if (length(gene_universe) == 0L) {
  stop(
    "No genes overlap between GFF3 annotation, dds, and the ",
    "Col-0-expressed gene set."
  )
}


message(
  "Analysis universe: GFF3-annotated, dds-represented, ",
  "Col-0-expressed genes = ",
  length(gene_universe)
)


# ----------------------------- 6. Recover HSR classes -------------------

hsr_classes_raw <- NULL
hsr_source <- NULL

if (exists("gene_sets", inherits = TRUE)) {
  
  gs_candidate <- get(
    "gene_sets",
    inherits = TRUE
  )
  
  
  if (
    is.data.frame(gs_candidate) &&
    all(c("gene", "GeneGroup") %in% names(gs_candidate))
  ) {
    
    mapped_groups <- vapply(
      as.character(gs_candidate$GeneGroup),
      canonical_group_label,
      character(1)
    )
    
    
    raw_classes <- list(
      
      Acute = clean_gene_set(
        gs_candidate$gene[
          mapped_groups == "Acute HS responsive genes"
        ]
      ),
      
      Sustained = clean_gene_set(
        gs_candidate$gene[
          mapped_groups == "Sustained HS responsive genes"
        ]
      ),
      
      Delayed = clean_gene_set(
        gs_candidate$gene[
          mapped_groups == "Delayed HS responsive genes"
        ]
      )
    )
    
    
    if (all(lengths(raw_classes) > 0L)) {
      hsr_classes_raw <- raw_classes
      hsr_source <- "existing gene_sets data.frame"
    }
  }
  
  if (
    is.null(hsr_classes_raw) &&
    is.list(gs_candidate) &&
    !is.data.frame(gs_candidate) &&
    !is.null(names(gs_candidate))
  ) {
    
    mapped_names <- vapply(
      names(gs_candidate),
      canonical_group_label,
      character(1)
    )
    
    
    idx_acute <- which(
      mapped_names == "Acute HS responsive genes"
    )
    
    idx_sustained <- which(
      mapped_names == "Sustained HS responsive genes"
    )
    
    idx_delayed <- which(
      mapped_names == "Delayed HS responsive genes"
    )
    
    
    if (
      length(idx_acute) == 1L &&
      length(idx_sustained) == 1L &&
      length(idx_delayed) == 1L
    ) {
      
      raw_classes <- list(
        Acute = clean_gene_set(
          gs_candidate[[idx_acute]]
        ),
        
        Sustained = clean_gene_set(
          gs_candidate[[idx_sustained]]
        ),
        
        Delayed = clean_gene_set(
          gs_candidate[[idx_delayed]]
        )
      )
      
      
      if (all(lengths(raw_classes) > 0L)) {
        hsr_classes_raw <- raw_classes
        hsr_source <- "existing named-list gene_sets object"
      }
    }
  }
}

if (
  is.null(hsr_classes_raw) &&
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
    
    hsr_classes_raw <- list(
      
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
    
    
    hsr_classes_raw <- list(
      
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

if (is.null(hsr_classes_raw)) {
  
  required_coef_names <- c(
    "treatment_LAT_vs_NT",
    "treatment_REC_vs_NT"
  )
  
  
  observed_coef_names <- DESeq2::resultsNames(dds)
  
  
  missing_coef_names <- setdiff(
    required_coef_names,
    observed_coef_names
  )
  
  
  if (length(missing_coef_names) > 0L) {
    stop(
      "Could not recover existing HSR gene classes, and expected DESeq2 ",
      "coefficients are missing: ",
      paste(missing_coef_names, collapse = ", ")
    )
  }
  
  
  design_terms <- attr(
    stats::terms(DESeq2::design(dds)),
    "term.labels"
  )
  
  
  has_interaction <- any(
    vapply(
      strsplit(design_terms, ":", fixed = TRUE),
      function(parts) {
        all(c("genotype", "treatment") %in% parts)
      },
      logical(1)
    )
  )
  
  
  observed_genotypes <- unique(
    coldata_df$Genotype
  )
  
  
  if (
    length(observed_genotypes) > 1L &&
    !has_interaction
  ) {
    stop(
      "Cannot reconstruct Col-0-specific HSR classes: dds contains ",
      "multiple genotypes but the design has no genotype:treatment ",
      "interaction. Use the existing Col-0 HSR gene sets or a suitable model."
    )
  }
  
  
  if (
    length(observed_genotypes) > 1L &&
    !is.factor(coldata_df$genotype)
  ) {
    stop(
      "A factor-valued genotype column is required to validate the ",
      "reference genotype for the fallback DESeq2 analysis."
    )
  }
  
  
  if (
    length(observed_genotypes) > 1L &&
    normalise_genotype(levels(coldata_df$genotype)[1]) != "Col-0"
  ) {
    stop(
      "The dds reference genotype is not Col-0; the treatment coefficients ",
      "cannot be assumed to represent Col-0 treatment effects."
    )
  }
  
  
  res_lat <- as.data.frame(
    DESeq2::results(
      dds,
      name = "treatment_LAT_vs_NT",
      alpha = alpha_cutoff
    )
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      gene = clean_gene_id(gene)
    )
  
  
  res_rec <- as.data.frame(
    DESeq2::results(
      dds,
      name = "treatment_REC_vs_NT",
      alpha = alpha_cutoff
    )
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      gene = clean_gene_id(gene)
    )
  
  
  lat_induced <- !is.na(res_lat$padj) &
    res_lat$padj < alpha_cutoff &
    !is.na(res_lat$log2FoldChange) &
    res_lat$log2FoldChange > lfc_cutoff
  
  
  rec_induced <- !is.na(res_rec$padj) &
    res_rec$padj < alpha_cutoff &
    !is.na(res_rec$log2FoldChange) &
    res_rec$log2FoldChange > lfc_cutoff
  
  
  lat_genes <- clean_gene_set(
    res_lat$gene[lat_induced]
  )
  
  rec_genes <- clean_gene_set(
    res_rec$gene[rec_induced]
  )
  
  
  hsr_classes_raw <- list(
    
    Acute = clean_gene_set(
      setdiff(lat_genes, rec_genes)
    ),
    
    Sustained = clean_gene_set(
      intersect(lat_genes, rec_genes)
    ),
    
    Delayed = clean_gene_set(
      setdiff(rec_genes, lat_genes)
    )
  )
  
  
  hsr_source <- "reconstructed from dds"
}


if (
  is.null(hsr_classes_raw) ||
  !all(c("Acute", "Sustained", "Delayed") %in% names(hsr_classes_raw))
) {
  stop(
    "Could not recover the Acute, Sustained and Delayed HSR classes."
  )
}


message("HSR class source: ", hsr_source)

hsr_classes <- lapply(
  hsr_classes_raw[c("Acute", "Sustained", "Delayed")],
  function(x) intersect(
    clean_gene_set(x),
    gene_universe
  )
)


# Check that Acute, Sustained and Delayed are mutually exclusive.
hsr_pairs <- combn(
  names(hsr_classes),
  2,
  simplify = FALSE
)


hsr_overlaps <- vapply(
  hsr_pairs,
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


if (any(hsr_overlaps > 0L)) {
  stop(
    "The HSR classes overlap after ID cleanup and universe filtering. ",
    "Check the upstream classification."
  )
}


hsr_union <- unique(
  unlist(hsr_classes, use.names = FALSE)
)


genes_by_group <- list(
  "Other expressed genes" = setdiff(
    gene_universe,
    hsr_union
  ),
  
  "HS responsive genes" = hsr_union,
  
  "Acute HS responsive genes" = hsr_classes$Acute,
  
  "Sustained HS responsive genes" = hsr_classes$Sustained,
  
  "Delayed HS responsive genes" = hsr_classes$Delayed
)


group_sizes <- tibble(
  Group = group_order,
  n_genes = vapply(
    genes_by_group[group_order],
    length,
    integer(1)
  )
)


message("Gene counts in the analysis universe:")
print(group_sizes)


# Export the actual group assignments.
gene_groups <- bind_rows(
  lapply(
    group_order,
    function(grp) {
      tibble(
        gene = genes_by_group[[grp]],
        Group = grp
      )
    }
  )
) %>%
  distinct(
    gene,
    Group
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = group_order
    )
  )


write.csv(
  gene_groups,
  paste0(out_prefix, "_gene_group_assignments.csv"),
  row.names = FALSE
)


write.csv(
  group_sizes,
  paste0(out_prefix, "_gene_group_sizes.csv"),
  row.names = FALSE
)


write.csv(
  data.frame(gene = gene_universe),
  paste0(out_prefix, "_analysis_universe.csv"),
  row.names = FALSE
)


# ----------------------------- 7. Map GFF transcripts to genes -----------

if (!"Parent" %in% names(gff_df)) {
  stop(
    "Imported GFF3 has no Parent field, which is required for transcript mapping."
  )
}


mrna_raw <- gff_df %>%
  filter(
    tolower(type) %in% c("mrna", "transcript")
  ) %>%
  dplyr::select(
    ID,
    Parent
  )


mrna_expanded <- expand_parent_rows(mrna_raw) %>%
  transmute(
    transcript = clean_tx_id(ID),
    gene = clean_gene_id(Parent)
  ) %>%
  filter(
    !is.na(transcript),
    nzchar(transcript),
    !is.na(gene),
    nzchar(gene)
  ) %>%
  distinct(
    transcript,
    gene
  )


map_tx_gene <- function(tx) {
  
  tx <- clean_tx_id(tx)
  
  mapped <- mrna_expanded$gene[
    match(tx, mrna_expanded$transcript)
  ]
  
  
  # Fallback only if a transcript ID directly embeds a TAIR gene ID.
  direct_agi <- grepl(
    "^AT[1-5CM]G[0-9]{5}(\\.[0-9]+)?$",
    tx,
    ignore.case = TRUE
  )
  
  
  fallback <- ifelse(
    direct_agi,
    clean_gene_id(tx),
    NA_character_
  )
  
  
  use_fallback <- is.na(mapped) | !nzchar(mapped)
  
  mapped[use_fallback] <- fallback[use_fallback]
  
  mapped
}


# ----------------------------- 8. 3'UTR length and representative isoform -

utr_raw <- gff_df %>%
  filter(
    tolower(type) == "three_prime_utr"
  ) %>%
  dplyr::select(
    Parent,
    start,
    end
  )


if (nrow(utr_raw) == 0L) {
  stop(
    "No three_prime_UTR features were found in the GFF3 annotation."
  )
}


utr_expanded <- expand_parent_rows(utr_raw) %>%
  transmute(
    transcript = clean_tx_id(Parent),
    start = as.numeric(start),
    end = as.numeric(end),
    width = as.numeric(end) - as.numeric(start) + 1
  ) %>%
  filter(
    !is.na(transcript),
    nzchar(transcript),
    is.finite(width),
    width > 0
  ) %>%
  distinct(
    transcript,
    start,
    end,
    .keep_all = TRUE
  ) %>%
  mutate(
    gene = map_tx_gene(transcript)
  ) %>%
  filter(
    !is.na(gene),
    nzchar(gene)
  )


if (nrow(utr_expanded) == 0L) {
  stop(
    "No usable 3'UTR segments could be mapped to genes."
  )
}


# Sum all annotated 3'UTR segments per transcript.
utr_by_tx <- utr_expanded %>%
  group_by(
    gene,
    transcript
  ) %>%
  summarise(
    utr3_length = sum(width, na.rm = TRUE),
    n_utr_segments = n_distinct(
      paste(start, end, sep = ":")
    ),
    .groups = "drop"
  )


# Panel A: the maximum total 3'UTR length across all isoforms per gene.
gene_utr_max <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    max_utr3_length = max(
      utr3_length,
      na.rm = TRUE
    ),
    n_annotated_utr_isoforms = n_distinct(transcript),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(max_utr3_length),
    max_utr3_length > 0
  )

cds_raw <- gff_df %>%
  filter(
    tolower(type) == "cds"
  ) %>%
  dplyr::select(
    Parent,
    start,
    end
  )


cds_expanded <- expand_parent_rows(cds_raw) %>%
  transmute(
    transcript = clean_tx_id(Parent),
    start = as.numeric(start),
    end = as.numeric(end),
    width = as.numeric(end) - as.numeric(start) + 1
  ) %>%
  filter(
    !is.na(transcript),
    nzchar(transcript),
    is.finite(width),
    width > 0
  ) %>%
  distinct(
    transcript,
    start,
    end,
    .keep_all = TRUE
  )


cds_by_tx <- cds_expanded %>%
  group_by(transcript) %>%
  summarise(
    cds_length = sum(width, na.rm = TRUE),
    .groups = "drop"
  )

representative_tx <- utr_by_tx %>%
  left_join(
    cds_by_tx,
    by = "transcript"
  ) %>%
  mutate(
    cds_length = tidyr::replace_na(
      cds_length,
      0
    )
  ) %>%
  arrange(
    gene,
    desc(cds_length),
    desc(utr3_length),
    transcript
  ) %>%
  group_by(gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  transmute(
    gene,
    representative_transcript = transcript,
    representative_utr3_length = utr3_length,
    representative_utr_segments = n_utr_segments,
    representative_cds_length = cds_length,
    has_3utr_intron = as.logical(n_utr_segments > 1)
  )


gene_features_3utr <- representative_tx %>%
  filter(
    gene %in% gene_universe
  )


gene_utr_length <- gene_utr_max %>%
  filter(
    gene %in% gene_universe
  ) %>%
  left_join(
    representative_tx,
    by = "gene"
  )

utr_intron_any <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    has_3utr_intron_any = any(
      n_utr_segments > 1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


gene_features_3utr <- gene_features_3utr %>%
  left_join(
    utr_intron_any,
    by = "gene"
  )


# ----------------------------- 9. uORF evidence and TPM evaluability ------

ribotaper_tab <- read_excel_table(
  ribotaper_file,
  sheet = "S1A_uORFs",
  skip = 0
)


cips_tab <- read_excel_table(
  cips_file,
  sheet = "A.ALL",
  skip = 2
)


ribotaper_genes <- extract_agi(
  ribotaper_tab,
  "gene_id"
)


cips_genes <- extract_agi(
  cips_tab,
  "gene_id"
)


if (length(ribotaper_genes) == 0L) {
  stop("No TAIR gene IDs could be extracted from the RiboTaper table.")
}


if (length(cips_genes) == 0L) {
  warning(
    "No TAIR gene IDs were extracted from the CiPS table. ",
    "CiPS sensitivity analyses may not be informative."
  )
}


if (!"AA" %in% names(cips_tab)) {
  
  warning(
    "CiPS input has no AA column; the CiPS >=10 aa sensitivity set will be empty."
  )
  
  cips_10aa_genes <- character(0)
  
} else {
  
  cips_10aa_tab <- cips_tab %>%
    mutate(
      .AA_numeric = suppressWarnings(
        as.numeric(as.character(AA))
      )
    ) %>%
    filter(
      !is.na(.AA_numeric),
      .AA_numeric >= 10
    )
  
  
  cips_10aa_genes <- extract_agi(
    cips_10aa_tab,
    "gene_id"
  )
}


uorf_external_sets <- list(
  "RiboTaper" = clean_gene_set(ribotaper_genes),
  "CiPS" = clean_gene_set(cips_genes),
  "CiPS >=10 aa" = clean_gene_set(cips_10aa_genes)
)


if (
  length(ribotaper_genes) > 0L &&
  length(cips_genes) > 0L
) {
  uorf_external_sets[["RiboTaper or CiPS"]] <- union(
    clean_gene_set(ribotaper_genes),
    clean_gene_set(cips_genes)
  )
}


# Read the Ribo-seq expression table used to define evaluable negatives.
ribo_expr <- read.delim(
  ribo_expr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (!all(c("gene_id", "TPM") %in% names(ribo_expr))) {
  stop(
    "The Ribo-seq expression file must contain columns 'gene_id' and 'TPM'."
  )
}


ribo_expr_clean <- tibble(
  gene = clean_gene_id(ribo_expr$gene_id),
  TPM = suppressWarnings(
    as.numeric(as.character(ribo_expr$TPM))
  )
) %>%
  filter(
    !is.na(gene),
    nzchar(gene),
    is.finite(TPM)
  ) %>%
  group_by(gene) %>%
  summarise(
    TPM = max(TPM, na.rm = TRUE),
    .groups = "drop"
  )


ribo_detectable_genes <- ribo_expr_clean %>%
  filter(
    TPM >= ribo_tpm_min
  ) %>%
  pull(gene) %>%
  unique()


if (length(ribo_detectable_genes) == 0L) {
  stop(
    "No genes passed the Ribo-seq TPM threshold. ",
    "Check the TPM file and ribo_tpm_min."
  )
}


# For each external uORF definition:
#   TRUE = positive uORF evidence, regardless of TPM;
#   FALSE = no positive call AND TPM is sufficiently high;
#   NA = no positive call and insufficient/no TPM evidence.
uorf_df <- tibble(
  gene = gene_universe
)


for (nm in names(uorf_external_sets)) {
  
  target_set <- uorf_external_sets[[nm]]
  
  positive <- uorf_df$gene %in% target_set
  detectable <- uorf_df$gene %in% ribo_detectable_genes
  
  flag <- rep(
    NA,
    nrow(uorf_df)
  )
  
  if (length(target_set) > 0L) {
    flag[detectable & !positive] <- FALSE
    flag[positive] <- TRUE
  } else {
    message(
      "No IDs were parsed for uORF definition '", nm,
      "'; its status remains NA rather than assigning unsupported negatives."
    )
  }
  
  uorf_df[[make.names(nm)]] <- flag
}


main_uorf_definition <- "RiboTaper"
main_uorf_col <- make.names(main_uorf_definition)

uorf_df$has_uORF <- uorf_df[[main_uorf_col]]


ribo_positive_not_detectable <- intersect(
  ribotaper_genes,
  setdiff(
    gene_universe,
    ribo_detectable_genes
  )
)


message(
  "RiboTaper genes: ", length(ribotaper_genes),
  "; Ribo-seq detectable genes (TPM >= ", ribo_tpm_min, "): ",
  length(ribo_detectable_genes),
  "; RiboTaper-positive genes retained below threshold: ",
  length(ribo_positive_not_detectable)
)


# ----------------------------- 10. Assemble feature data ----------------

all_group_genes <- bind_rows(
  lapply(
    group_order,
    function(grp) {
      tibble(
        gene = genes_by_group[[grp]],
        Group = grp
      )
    }
  )
) %>%
  distinct(
    Group,
    gene
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = group_order
    )
  )


# Panel A: genes with an annotated positive 3'UTR length.
utr_data <- all_group_genes %>%
  inner_join(
    gene_utr_max,
    by = "gene"
  ) %>%
  filter(
    is.finite(max_utr3_length),
    max_utr3_length > 0
  )


# Representative-transcript 3'UTR length sensitivity analysis.
rep_utr_data <- all_group_genes %>%
  inner_join(
    representative_tx,
    by = "gene"
  ) %>%
  filter(
    is.finite(representative_utr3_length),
    representative_utr3_length > 0
  )


# Panel B: intron status in the selected representative transcript.
intron_gene_status <- gene_features_3utr %>%
  transmute(
    gene,
    has_3utr_intron = as.logical(has_3utr_intron),
    has_3utr_intron_any = as.logical(has_3utr_intron_any)
  )


# General function for panel-specific positive percentages.
# NA calls are excluded from the denominator.
calculate_percentage <- function(
    feature_status,
    flag_col,
    feature_name
) {
  
  bind_rows(
    lapply(
      group_order,
      function(grp) {
        
        group_gene_ids <- unique(
          genes_by_group[[grp]]
        )
        
        
        eval_df <- tibble(
          gene = group_gene_ids
        ) %>%
          left_join(
            feature_status %>%
              dplyr::select(
                gene,
                flag = all_of(flag_col)
              ),
            by = "gene"
          ) %>%
          filter(
            !is.na(flag)
          )
        
        
        n_total <- nrow(eval_df)
        
        n_positive <- if (n_total > 0L) {
          sum(eval_df$flag %in% TRUE)
        } else {
          0L
        }
        
        
        tibble(
          Feature = feature_name,
          Group = grp,
          n_total = n_total,
          n_positive = n_positive,
          n_negative = n_total - n_positive,
          percentage = if (n_total > 0L) {
            100 * n_positive / n_total
          } else {
            NA_real_
          }
        )
      }
    )
  ) %>%
    mutate(
      Group = factor(
        Group,
        levels = group_order
      )
    )
}


intron_data <- calculate_percentage(
  intron_gene_status,
  "has_3utr_intron",
  "3'UTR intron (representative transcript)"
)


intron_any_data <- calculate_percentage(
  intron_gene_status,
  "has_3utr_intron_any",
  "3'UTR intron (any isoform)"
)


# Panel C: no 3'UTR annotation filter; evaluability is determined by
# RiboTaper status or Ribo-seq TPM detectability.
uorf_gene_status <- uorf_df %>%
  transmute(
    gene,
    has_uORF = as.logical(has_uORF)
  )


uorf_data <- calculate_percentage(
  uorf_gene_status,
  "has_uORF",
  "RiboTaper-supported uORF"
)


# Panel-specific denominators for transparency.
coverage_table <- bind_rows(
  lapply(
    group_order,
    function(grp) {
      
      genes <- unique(
        genes_by_group[[grp]]
      )
      
      
      matched_uorf_status <- uorf_df$has_uORF[
        match(genes, uorf_df$gene)
      ]
      
      
      tibble(
        Group = grp,
        n_genes_in_group = length(genes),
        n_with_annotated_3UTR = sum(
          genes %in% gene_utr_max$gene
        ),
        n_with_representative_transcript = sum(
          genes %in% gene_features_3utr$gene
        ),
        n_RiboTaper_evaluable = sum(
          !is.na(matched_uorf_status)
        ),
        n_RiboTaper_positive = sum(
          matched_uorf_status %in% TRUE,
          na.rm = TRUE
        )
      )
    }
  )
)


# ----------------------------- 11. Statistical tests --------------------

# Panel A: Wilcoxon rank-sum test, each foreground against reference.
get_wilcox_tests <- function(dat, value_col, metric_label) {
  
  ref <- dat %>%
    filter(
      as.character(Group) == reference_group
    ) %>%
    pull(all_of(value_col))
  
  
  bind_rows(
    lapply(
      foreground_groups,
      function(grp) {
        
        sg <- dat %>%
          filter(
            as.character(Group) == grp
          ) %>%
          pull(all_of(value_col))
        
        
        ref_use <- ref[
          is.finite(ref) & ref > 0
        ]
        
        sg_use <- sg[
          is.finite(sg) & sg > 0
        ]
        
        
        p <- if (
          length(ref_use) >= 2L &&
          length(sg_use) >= 2L
        ) {
          
          tryCatch(
            stats::wilcox.test(
              ref_use,
              sg_use,
              exact = FALSE
            )$p.value,
            error = function(e) NA_real_
          )
          
        } else {
          NA_real_
        }
        
        
        tibble(
          Feature = metric_label,
          Group = grp,
          reference_group = reference_group,
          n_reference = length(ref_use),
          n_group = length(sg_use),
          median_reference = if (length(ref_use) > 0L) {
            median(ref_use)
          } else {
            NA_real_
          },
          median_group = if (length(sg_use) > 0L) {
            median(sg_use)
          } else {
            NA_real_
          },
          p_value = p
        )
      }
    )
  ) %>%
    mutate(
      p_adj_BH = p.adjust(
        p_value,
        method = "BH"
      ),
      
      significance = vapply(
        p_adj_BH,
        p_to_stars,
        character(1)
      ),
      
      Group = factor(
        Group,
        levels = group_order
      )
    )
}


# Panels B/C: Fisher exact test, each foreground vs reference.
get_fisher_tests <- function(dat, feature_name) {
  
  ref <- dat %>%
    filter(
      as.character(Group) == reference_group
    )
  
  
  bind_rows(
    lapply(
      foreground_groups,
      function(grp) {
        
        sg <- dat %>%
          filter(
            as.character(Group) == grp
          )
        
        
        n1 <- if (nrow(ref) > 0L) {
          ref$n_total[1]
        } else {
          0L
        }
        
        p1 <- if (nrow(ref) > 0L) {
          ref$n_positive[1]
        } else {
          0L
        }
        
        n2 <- if (nrow(sg) > 0L) {
          sg$n_total[1]
        } else {
          0L
        }
        
        p2 <- if (nrow(sg) > 0L) {
          sg$n_positive[1]
        } else {
          0L
        }
        
        
        p <- NA_real_
        odds_ratio <- NA_real_
        ci_low <- NA_real_
        ci_high <- NA_real_
        
        
        if (n1 > 0L && n2 > 0L) {
          
          # Foreground in first row, reference group in second row.
          tab <- matrix(
            c(
              p2,
              n2 - p2,
              p1,
              n1 - p1
            ),
            nrow = 2,
            byrow = TRUE,
            dimnames = list(
              c(grp, reference_group),
              c("positive", "negative")
            )
          )
          
          
          ft <- tryCatch(
            stats::fisher.test(tab),
            error = function(e) NULL
          )
          
          
          if (!is.null(ft)) {
            
            p <- ft$p.value
            
            odds_ratio <- if (
              length(ft$estimate) == 1L
            ) {
              unname(ft$estimate)
            } else {
              NA_real_
            }
            
            
            if (length(ft$conf.int) == 2L) {
              ci_low <- unname(ft$conf.int[1])
              ci_high <- unname(ft$conf.int[2])
            }
          }
        }
        
        
        tibble(
          Feature = feature_name,
          Group = grp,
          reference_group = reference_group,
          n_reference = n1,
          reference_positive = p1,
          n_group = n2,
          group_positive = p2,
          odds_ratio_group_vs_reference = odds_ratio,
          CI_low = ci_low,
          CI_high = ci_high,
          p_value = p
        )
      }
    )
  ) %>%
    mutate(
      p_adj_BH = p.adjust(
        p_value,
        method = "BH"
      ),
      
      significance = vapply(
        p_adj_BH,
        p_to_stars,
        character(1)
      ),
      
      Group = factor(
        Group,
        levels = group_order
      )
    )
}


utr_tests <- get_wilcox_tests(
  utr_data,
  "max_utr3_length",
  "Maximum annotated 3'UTR length"
)


rep_utr_tests <- get_wilcox_tests(
  rep_utr_data,
  "representative_utr3_length",
  "Representative-transcript 3'UTR length"
)


intron_tests <- get_fisher_tests(
  intron_data,
  "3'UTR intron (representative transcript)"
)


intron_any_tests <- get_fisher_tests(
  intron_any_data,
  "3'UTR intron (any isoform)"
)


uorf_tests <- get_fisher_tests(
  uorf_data,
  "RiboTaper-supported uORF"
)


# Sensitivity tests for the other external uORF definitions.
uorf_robustness <- bind_rows(
  lapply(
    names(uorf_external_sets),
    function(definition) {
      
      flag_col <- make.names(definition)
      
      status <- uorf_df %>%
        transmute(
          gene,
          flag = as.logical(.data[[flag_col]])
        )
      
      
      dat <- calculate_percentage(
        status,
        "flag",
        definition
      )
      
      
      get_fisher_tests(
        dat,
        definition
      ) %>%
        mutate(
          definition = definition
        )
    }
  )
)


# ----------------------------- 12. Plot settings ------------------------
base_theme <- theme_classic(
  base_size = 16
) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      colour = "black",
      size = 27,
      face = "bold",
      angle = 22.5,
      hjust = 1,
      vjust = 1,
      lineheight = 0.95
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 31.5,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      size = 36,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 24,
      face = "bold",
      margin = margin(r = 6)
    ),
    
    plot.title = element_text(
      size = 40.5,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    
    plot.margin = margin(
      10,
      12,
      10,
      10
    )
  )


# ----------------------------- 13. Panel A: 3'UTR length ----------------

if (nrow(utr_data) == 0L) {
  stop(
    "No genes with annotated positive 3'UTR lengths are available for Panel A."
  )
}


utr_max_by_group <- utr_data %>%
  group_by(Group) %>%
  summarise(
    max_value = max(
      max_utr3_length,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


utr_label_positions <- utr_tests %>%
  left_join(
    utr_max_by_group,
    by = "Group"
  ) %>%
  mutate(
    label_y = max_value * 1.45
  )


finite_utr_y <- c(
  utr_data$max_utr3_length,
  utr_label_positions$label_y
)

finite_utr_y <- finite_utr_y[
  is.finite(finite_utr_y) & finite_utr_y > 0
]


if (length(finite_utr_y) == 0L) {
  stop("Cannot calculate finite y-axis limits for Panel A.")
}


len_y_upper <- max(finite_utr_y) * 1.15
len_y_upper <- max(len_y_upper, 10)


pA <- ggplot(
  utr_data,
  aes(
    x = Group,
    y = max_utr3_length,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.68,
    outlier.size = 1.0,
    outlier.alpha = 0.30,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_text(
    data = utr_label_positions,
    aes(
      x = Group,
      y = label_y,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 13.5,
    fontface = "bold",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = group_order,
    labels = group_labels,
    drop = FALSE
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  coord_cartesian(
    ylim = c(1, len_y_upper),
    clip = "off"
  ) +
  labs(
    title = "",
    x = NULL,
    y = "Maximum 3′UTR length (bp, log scale)"
  ) +
  base_theme


# ----------------------------- 14. Percentage-panel plotting function ---

make_percentage_plot <- function(
    dat,
    tests,
    y_title,
    panel_title
) {
  
  if (nrow(dat) == 0L) {
    stop("No data available for Panel: ", panel_title)
  }
  
  
  max_pct <- if (all(is.na(dat$percentage))) {
    0
  } else {
    max(dat$percentage, na.rm = TRUE)
  }
  
  
  star_offset <- max(
    8,
    max_pct * 0.08
  )
  
  
  pct_labels <- dat %>%
    mutate(
      label_y = percentage + 1.5,
      pct_label = ifelse(
        is.na(percentage),
        "",
        sprintf("%.1f%%", percentage)
      )
    )
  
  
  sig_labels <- tests %>%
    left_join(
      dat %>%
        dplyr::select(
          Group,
          percentage
        ),
      by = "Group"
    ) %>%
    mutate(
      label_y = percentage + star_offset + 2
    )
  
  
  plot_y_values <- c(
    dat$percentage,
    pct_labels$label_y,
    sig_labels$label_y
  )
  
  plot_y_values <- plot_y_values[
    is.finite(plot_y_values)
  ]
  
  
  y_upper <- if (length(plot_y_values) > 0L) {
    max(25, max(plot_y_values) + 8)
  } else {
    25
  }
  
  y_upper <- min(
    140,
    y_upper
  )
  
  
  ggplot(
    dat,
    aes(
      x = Group,
      y = percentage,
      fill = Group
    )
  ) +
    geom_col(
      width = 0.68,
      na.rm = TRUE
    ) +
    geom_text(
      data = pct_labels,
      aes(
        x = Group,
        y = label_y,
        label = pct_label
      ),
      inherit.aes = FALSE,
      size = 11.25,
      fontface = "bold",
      na.rm = TRUE
    ) +
    geom_text(
      data = sig_labels,
      aes(
        x = Group,
        y = label_y,
        label = significance
      ),
      inherit.aes = FALSE,
      size = 13.5,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = group_colors,
      drop = FALSE
    ) +
    scale_x_discrete(
      limits = group_order,
      labels = group_labels,
      drop = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(
        0,
        min(100, floor(y_upper / 10) * 10),
        by = 10
      ),
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +
    coord_cartesian(
      ylim = c(0, y_upper),
      clip = "off"
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = y_title
    ) +
    base_theme
}


# ----------------------------- 15. Panels B and C -----------------------

pB <- make_percentage_plot(
  intron_data,
  intron_tests,
  "Genes with a 3′UTR intron (%)",
  ""
)


pC <- make_percentage_plot(
  uorf_data,
  uorf_tests,
  "Genes with a RiboTaper-supported uORF (%)",
  ""
)


# ----------------------------- 16. Combine panels -----------------------

combined_plot <- pA + pB + pC +
  plot_layout(
    ncol = 3,
    widths = c(1.25, 1.10, 1.10)
  ) +
  plot_annotation(
    title = "Structural features of Col-0 HS responsive gene sets",
    caption = paste0(
      "Reference = other Col-0-expressed genes. ",
      "Induced: BH-adjusted p < ", alpha_cutoff,
      " and log2FC > ", lfc_cutoff, ". ",
      "Stars use BH-adjusted p-values; the four foreground-versus-reference ",
      "comparisons are corrected within each feature family. ",
      "Panel C has no 3′UTR annotation filter; negative uORF calls require ",
      "Ribo-seq TPM ≥ ", ribo_tpm_min, "."
    ),
    theme = theme(
      plot.title = element_text(
        size = 54,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      plot.caption = element_text(
        size = 25,
        hjust = 0
      )
    )
  )


print(combined_plot)


# ----------------------------- 17. Export tables ------------------------

write.csv(
  gene_groups,
  paste0(out_prefix, "_gene_group_assignments.csv"),
  row.names = FALSE
)


write.csv(
  group_sizes,
  paste0(out_prefix, "_gene_group_sizes.csv"),
  row.names = FALSE
)


write.csv(
  data.frame(gene = gene_universe),
  paste0(out_prefix, "_analysis_universe.csv"),
  row.names = FALSE
)


write.csv(
  gene_utr_max,
  paste0(out_prefix, "_utr_length_max_by_gene.csv"),
  row.names = FALSE
)


write.csv(
  gene_utr_length,
  paste0(out_prefix, "_utr_length_and_representative_transcript_by_gene.csv"),
  row.names = FALSE
)


write.csv(
  gene_features_3utr,
  paste0(out_prefix, "_representative_utr_intron_by_gene.csv"),
  row.names = FALSE
)


write.csv(
  uorf_df,
  paste0(out_prefix, "_uorf_status_by_gene.csv"),
  row.names = FALSE
)


write.csv(
  utr_data,
  paste0(out_prefix, "_utr_length_values.csv"),
  row.names = FALSE
)


write.csv(
  rep_utr_data,
  paste0(out_prefix, "_representative_utr_length_robustness_values.csv"),
  row.names = FALSE
)


write.csv(
  utr_tests,
  paste0(out_prefix, "_utr_length_tests.csv"),
  row.names = FALSE
)


write.csv(
  rep_utr_tests,
  paste0(out_prefix, "_representative_utr_length_robustness_tests.csv"),
  row.names = FALSE
)


write.csv(
  intron_data,
  paste0(out_prefix, "_utr_intron_percentages.csv"),
  row.names = FALSE
)


write.csv(
  intron_tests,
  paste0(out_prefix, "_utr_intron_tests.csv"),
  row.names = FALSE
)


write.csv(
  intron_any_data,
  paste0(out_prefix, "_utr_intron_any_isoform_robustness_percentages.csv"),
  row.names = FALSE
)


write.csv(
  intron_any_tests,
  paste0(out_prefix, "_utr_intron_any_isoform_robustness_tests.csv"),
  row.names = FALSE
)


write.csv(
  uorf_data,
  paste0(out_prefix, "_uorf_percentages.csv"),
  row.names = FALSE
)


write.csv(
  uorf_tests,
  paste0(out_prefix, "_uorf_tests.csv"),
  row.names = FALSE
)


write.csv(
  uorf_robustness,
  paste0(out_prefix, "_uorf_definition_robustness.csv"),
  row.names = FALSE
)


write.csv(
  coverage_table,
  paste0(out_prefix, "_panel_evaluable_denominators.csv"),
  row.names = FALSE
)


external_set_audit <- bind_rows(
  lapply(
    names(uorf_external_sets),
    function(nm) {
      
      s <- clean_gene_set(
        uorf_external_sets[[nm]]
      )
      
      tibble(
        definition = nm,
        n_unique_external_genes = length(s),
        n_in_analysis_universe = length(
          intersect(s, gene_universe)
        ),
        n_in_Ribo_detectable_set = length(
          intersect(s, ribo_detectable_genes)
        ),
        n_in_both_universe_and_detectable = length(
          intersect(
            intersect(s, gene_universe),
            ribo_detectable_genes
          )
        )
      )
    }
  )
)


write.csv(
  external_set_audit,
  paste0(out_prefix, "_external_uORF_set_audit.csv"),
  row.names = FALSE
)


# ----------------------------- 18. Save figures -------------------------
ggsave(
  paste0(out_prefix, ".png"),
  combined_plot,
  width = 30,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  paste0(out_prefix, ".pdf"),
  combined_plot,
  width = 18,
  height = 7.5,
  units = "in",
  bg = "white"
)

ggsave(
  paste0(out_prefix, "_panelA_UTR_length.png"),
  pA,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  paste0(out_prefix, "_panelA_UTR_length.pdf"),
  pA,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)


ggsave(
  paste0(out_prefix, "_panelB_UTR_intron.png"),
  pB,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  paste0(out_prefix, "_panelB_UTR_intron.pdf"),
  pB,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)


ggsave(
  paste0(out_prefix, "_panelC_uORF_RiboTaper.png"),
  pC,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  paste0(out_prefix, "_panelC_uORF_RiboTaper.pdf"),
  pC,
  width = 8,
  height = 6,
  units = "in",
  bg = "white"
)


# ----------------------------- 19. Final report ------------------------

message("")
message("Structural cis-feature comparison completed.")
message("HSR gene-set source: ", hsr_source)
message("GFF3 annotation: ", annotation_file)
message("Analysis universe: ", length(gene_universe), " genes")
message("Reference group: ", reference_group)
message("Significance labels: BH-adjusted p-values")
message("BH significance threshold: ", alpha_cutoff)
message("Ribo-seq evaluability threshold: TPM >= ", ribo_tpm_min)
message("")
message("Group sizes:")
print(group_sizes)

message("")
message("Genes with annotated 3′UTR lengths:")
print(
  utr_data %>%
    dplyr::count(Group, name = "n_genes_with_UTR")
)

message("")
message("Representative-transcript 3′UTR intron percentages:")
print(
  intron_data %>%
    dplyr::select(Group, n_total, n_positive, percentage)
)

message("")
message("RiboTaper uORF percentages:")
print(
  uorf_data %>%
    dplyr::select(Group, n_total, n_positive, percentage)
)

message("")
message("Main outputs:")
message(out_prefix, ".png")
message(out_prefix, ".pdf")
message(out_prefix, "_panelA_UTR_length.png / .pdf")
message(out_prefix, "_panelB_UTR_intron.png / .pdf")
message(out_prefix, "_panelC_uORF_RiboTaper.png / .pdf")
message("")
message(
  "Interpretation note: the broad HS-responsive gene set overlaps ",
  "the Acute, Sustained and Delayed sets by definition. ",
  "All foreground-vs-reference tests are gene-level comparisons, ",
  "not biological-replicate-level tests."
)
