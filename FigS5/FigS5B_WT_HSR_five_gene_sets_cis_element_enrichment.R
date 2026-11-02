# =============================================================================
# Cis-regulatory element enrichment in five Col-0 HSR gene sets

# Analysis:
#   - Extract full 1-kb promoters upstream of gene TSSs.
#   - Scan motifs on both DNA strands.
#   - One-sided Fisher's exact test for enrichment.
#   - BH correction across all available gene-set x motif tests.
#   - Significance stars are based on BH-adjusted p-values.
# =============================================================================


# ----------------------------- 0. Packages -----------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})


# ----------------------------- 1. Settings -----------------------------

alpha_cutoff <- 0.05
lfc_cutoff <- 1
promoter_length <- 1000L

annotation_candidates <- c(
  "Araport11_GFF3_genes_transposons.20241001.gff.gz",
  "Araport11_GFF3_genes_transposons.20230105.gff.gz"
)

output_prefix <- "Col_HSR_five_gene_sets_cis_element_enrichment"

# Odds ratios are visually capped at this value.
or_color_max <- 3

# Significance stars are based on BH-adjusted p-values.
# No raw-p-value significance mode is used.
significance_basis <- "BH-adjusted p-value"


# ----------------------------- 2. Motifs -------------------------------

motifs <- tibble::tribble(
  ~Element,          ~Pattern,      ~ElementClass,
  "W-box\n(WRKY)",   "TTGACY",      "Candidate element",
  "as-1\n(TGA)",     "TGACG",       "Candidate element",
  "HSE",             "NGAANNTTCN",  "Candidate element",
  "GCC-box",         "AGCCGCC",     "Candidate element",
  "ABRE",            "ACGTGKC",     "Candidate element",
  "DRE/CRT",         "RCCGAC",      "Candidate element",
  "TCA-element",     "CCATCTTTTT",  "Candidate element",
  "TATA-box",        "TATAWAWR",    "Negative control",
  "CAAT-box",        "CCAAT",       "Negative control",
  "GC-box",          "GGGCGG",      "Negative control"
) %>%
  mutate(
    ElementClass = factor(
      ElementClass,
      levels = c("Candidate element", "Negative control")
    ),
    Element = factor(
      Element,
      levels = Element
    )
  )


group_order <- c(
  "HS responsive genes",
  "Acute HS responsive genes",
  "Sustained HS responsive genes",
  "Delayed HS responsive genes",
  "Other expressed genes"
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
  
  unique(
    x[!is.na(x) & nzchar(x)]
  )
}


normalise_group_label <- function(x) {
  tolower(
    gsub("[^[:alnum:]]+", "", trimws(as.character(x)))
  )
}

map_group_label <- function(x) {
  
  z <- normalise_group_label(x)
  
  if (is.na(z) || !nzchar(z)) {
    return(NA_character_)
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
  
  if (grepl("otherexpress", z)) {
    return("Other expressed genes")
  }
  
  if (grepl("hsresponsive", z)) {
    return("HS responsive genes")
  }
  
  NA_character_
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


write_gene_list <- function(genes, filename) {
  write.csv(
    data.frame(gene = sort(unique(genes))),
    filename,
    row.names = FALSE
  )
}


# ----------------------------- 4. Recover current HSR gene sets --------

gene_sets_current <- NULL
gene_set_source <- NULL

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
      map_group_label,
      character(1)
    )
    
    recognised <- !is.na(mapped_groups)
    
    
    if (any(recognised)) {
      
      gene_sets_current <- tibble(
        gene = clean_gene_id(
          gs_candidate$gene[recognised]
        ),
        GeneGroup = mapped_groups[recognised]
      ) %>%
        filter(
          !is.na(gene),
          nzchar(gene)
        ) %>%
        distinct(
          gene,
          GeneGroup
        )
      
      gene_set_source <- "existing gene_sets data.frame"
    }
  }
}

if (
  is.null(gene_sets_current) &&
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
        "is_HSR_responsive",
        "is_acute",
        "is_sustained",
        "is_delayed"
      ) %in% names(gc)
    )
  ) {
    
    responsive_genes <- clean_gene_set(
      gc$gene[
        !is.na(gc$is_HSR_responsive) &
          gc$is_HSR_responsive
      ]
    )
    
    acute_genes <- clean_gene_set(
      gc$gene[
        !is.na(gc$is_acute) &
          gc$is_acute
      ]
    )
    
    sustained_genes <- clean_gene_set(
      gc$gene[
        !is.na(gc$is_sustained) &
          gc$is_sustained
      ]
    )
    
    delayed_genes <- clean_gene_set(
      gc$gene[
        !is.na(gc$is_delayed) &
          gc$is_delayed
      ]
    )
    
    
    if (exists("wt_expressed_genes", inherits = TRUE)) {
      
      expressed_genes <- clean_gene_set(
        get("wt_expressed_genes", inherits = TRUE)
      )
      
    } else if (exists("dds", inherits = TRUE)) {
      
      dds_obj <- get(
        "dds",
        inherits = TRUE
      )
      
      norm_counts <- DESeq2::counts(
        dds_obj,
        normalized = TRUE
      )
      
      rownames(norm_counts) <- clean_gene_id(
        rownames(norm_counts)
      )
      
      cd <- as.data.frame(
        SummarizedExperiment::colData(dds_obj)
      ) %>%
        tibble::rownames_to_column("Sample")
      
      
      if (!all(c("genotype", "treatment") %in% names(cd))) {
        stop(
          "colData(dds) must contain genotype and treatment columns ",
          "to identify the expressed Col-0 genes."
        )
      }
      
      
      genotype_norm <- tolower(
        gsub(
          "[^[:alnum:]]+",
          "",
          trimws(as.character(cd$genotype))
        )
      )
      
      
      wt_samples <- cd$Sample[
        genotype_norm %in% c("col", "col0", "wt", "wildtype")
      ]
      
      
      wt_samples <- intersect(
        wt_samples,
        colnames(norm_counts)
      )
      
      
      if (length(wt_samples) == 0L) {
        stop(
          "Could not identify Col-0 samples in the dds object."
        )
      }
      
      
      expressed_genes <- rownames(norm_counts)[
        rowSums(
          norm_counts[, wt_samples, drop = FALSE] > 0,
          na.rm = TRUE
        ) > 0
      ] %>%
        clean_gene_set()
      
    } else {
      
      stop(
        "Other expressed genes cannot be recovered: neither ",
        "wt_expressed_genes nor a dds object is available."
      )
    }
    
    
    other_genes <- setdiff(
      expressed_genes,
      responsive_genes
    )
    
    
    gene_sets_current <- bind_rows(
      
      tibble(
        gene = responsive_genes,
        GeneGroup = group_order[1]
      ),
      
      tibble(
        gene = acute_genes,
        GeneGroup = group_order[2]
      ),
      
      tibble(
        gene = sustained_genes,
        GeneGroup = group_order[3]
      ),
      
      tibble(
        gene = delayed_genes,
        GeneGroup = group_order[4]
      ),
      
      tibble(
        gene = other_genes,
        GeneGroup = group_order[5]
      )
      
    ) %>%
      filter(
        !is.na(gene),
        nzchar(gene)
      ) %>%
      distinct(
        gene,
        GeneGroup
      )
    
    
    gene_set_source <- "existing gene_classification plus expressed genes"
  }
}

if (
  is.null(gene_sets_current) &&
  exists("dds", inherits = TRUE)
) {
  
  dds_obj <- get(
    "dds",
    inherits = TRUE
  )
  
  
  if (!inherits(dds_obj, "DESeqDataSet")) {
    stop("'dds' must be a DESeqDataSet.")
  }
  
  
  cd <- as.data.frame(
    SummarizedExperiment::colData(dds_obj)
  )
  
  
  if (!all(c("genotype", "treatment") %in% names(cd))) {
    stop(
      "colData(dds) must contain genotype and treatment columns."
    )
  }
  
  
  if (
    !is.factor(cd$genotype) ||
    !is.factor(cd$treatment)
  ) {
    stop(
      "The genotype and treatment columns must be factors ",
      "to reconstruct the Col-0 HSR gene sets."
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
  
  
  genotype_ref <- normalise_genotype(
    levels(cd$genotype)[1]
  )
  
  
  treatment_ref <- toupper(
    trimws(as.character(levels(cd$treatment)[1]))
  )
  
  
  if (treatment_ref != "NT") {
    stop(
      "Cannot reconstruct the HSR classes: treatment reference must be NT."
    )
  }
  
  
  coef_names <- DESeq2::resultsNames(dds_obj)
  
  required_coefs <- c(
    "treatment_LAT_vs_NT",
    "treatment_REC_vs_NT"
  )
  
  
  missing_coefs <- setdiff(
    required_coefs,
    coef_names
  )
  
  
  if (length(missing_coefs) > 0L) {
    stop(
      "Missing expected DESeq2 coefficient(s): ",
      paste(missing_coefs, collapse = ", ")
    )
  }
  
  
  design_terms <- attr(
    stats::terms(DESeq2::design(dds_obj)),
    "term.labels"
  )
  
  
  has_interaction <- any(
    vapply(
      strsplit(design_terms, ":", fixed = TRUE),
      function(z) {
        all(c("genotype", "treatment") %in% z)
      },
      logical(1)
    )
  )
  
  
  observed_genotypes <- unique(
    normalise_genotype(cd$genotype)
  )
  
  observed_genotypes <- observed_genotypes[
    !is.na(observed_genotypes)
  ]
  
  
  if (
    length(observed_genotypes) > 1L &&
    !has_interaction
  ) {
    stop(
      "The dds object contains multiple genotypes without a ",
      "genotype:treatment interaction. The treatment main effects ",
      "cannot automatically be interpreted as Col-0-specific."
    )
  }
  
  
  if (
    length(observed_genotypes) > 1L &&
    (is.na(genotype_ref) || genotype_ref != "Col-0")
  ) {
    stop(
      "The reference genotype is not Col-0. ",
      "Cannot reconstruct Col-0-specific treatment effects."
    )
  }
  
  
  res_lat <- as.data.frame(
    DESeq2::results(
      dds_obj,
      name = "treatment_LAT_vs_NT",
      alpha = alpha_cutoff
    )
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene))
  
  
  res_rec <- as.data.frame(
    DESeq2::results(
      dds_obj,
      name = "treatment_REC_vs_NT",
      alpha = alpha_cutoff
    )
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = clean_gene_id(gene))
  
  
  lat_induced <- !is.na(res_lat$padj) &
    res_lat$padj < alpha_cutoff &
    !is.na(res_lat$log2FoldChange) &
    res_lat$log2FoldChange > lfc_cutoff
  
  
  rec_induced <- !is.na(res_rec$padj) &
    res_rec$padj < alpha_cutoff &
    !is.na(res_rec$log2FoldChange) &
    res_rec$log2FoldChange > lfc_cutoff
  
  
  lat_genes <- unique(res_lat$gene[lat_induced])
  rec_genes <- unique(res_rec$gene[rec_induced])
  
  
  responsive_genes <- clean_gene_set(
    union(lat_genes, rec_genes)
  )
  
  
  acute_genes <- clean_gene_set(
    setdiff(lat_genes, rec_genes)
  )
  
  
  sustained_genes <- clean_gene_set(
    intersect(lat_genes, rec_genes)
  )
  
  
  delayed_genes <- clean_gene_set(
    setdiff(rec_genes, lat_genes)
  )
  
  
  norm_counts <- DESeq2::counts(
    dds_obj,
    normalized = TRUE
  )
  
  
  rownames(norm_counts) <- clean_gene_id(
    rownames(norm_counts)
  )
  
  
  genotype_norm <- normalise_genotype(cd$genotype)
  
  
  wt_samples <- rownames(cd)[
    genotype_norm == "Col-0"
  ]
  
  
  wt_samples <- intersect(
    wt_samples,
    colnames(norm_counts)
  )
  
  
  if (length(wt_samples) == 0L) {
    stop("No Col-0 samples are available for expression filtering.")
  }
  
  
  expressed_genes <- rownames(norm_counts)[
    rowSums(
      norm_counts[, wt_samples, drop = FALSE] > 0,
      na.rm = TRUE
    ) > 0
  ] %>%
    clean_gene_set()
  
  
  other_genes <- setdiff(
    expressed_genes,
    responsive_genes
  )
  
  
  gene_sets_current <- bind_rows(
    
    tibble(
      gene = responsive_genes,
      GeneGroup = group_order[1]
    ),
    
    tibble(
      gene = acute_genes,
      GeneGroup = group_order[2]
    ),
    
    tibble(
      gene = sustained_genes,
      GeneGroup = group_order[3]
    ),
    
    tibble(
      gene = delayed_genes,
      GeneGroup = group_order[4]
    ),
    
    tibble(
      gene = other_genes,
      GeneGroup = group_order[5]
    )
    
  ) %>%
    distinct(
      gene,
      GeneGroup
    )
  
  
  gene_set_source <- "reconstructed from dds"
}


if (is.null(gene_sets_current)) {
  stop(
    "Could not recover the five Col-0 HSR gene sets. ",
    "Run the current HSR gene-set workflow first."
  )
}


message("Gene-set source: ", gene_set_source)

gene_sets_current <- gene_sets_current %>%
  mutate(
    gene = clean_gene_id(gene),
    GeneGroup = as.character(GeneGroup)
  ) %>%
  filter(
    !is.na(gene),
    nzchar(gene),
    GeneGroup %in% group_order
  ) %>%
  distinct(
    gene,
    GeneGroup
  )


gene_sets_list <- setNames(
  lapply(
    group_order,
    function(g) {
      clean_gene_set(
        gene_sets_current$gene[
          gene_sets_current$GeneGroup == g
        ]
      )
    }
  ),
  group_order
)


group_sizes <- tibble(
  GeneGroup = group_order,
  n_genes = vapply(
    gene_sets_list,
    length,
    integer(1)
  )
)


message("HSR gene-set sizes:")
print(group_sizes)


if (all(group_sizes$n_genes == 0L)) {
  stop("All five gene sets are empty.")
}


write.csv(
  gene_sets_current,
  paste0(output_prefix, "_gene_set_membership.csv"),
  row.names = FALSE
)


write.csv(
  group_sizes,
  paste0(output_prefix, "_gene_set_sizes.csv"),
  row.names = FALSE
)


# ----------------------------- 5. Load annotation -----------------------

found_annotations <- annotation_candidates[
  file.exists(annotation_candidates)
]


if (length(found_annotations) == 0L) {
  stop(
    "Could not find a GFF3 annotation file. Expected one of: ",
    paste(annotation_candidates, collapse = ", ")
  )
}


annotation_file <- found_annotations[1]


message("Using GFF3 annotation: ", annotation_file)


gff <- rtracklayer::import(annotation_file)
gff_df <- as.data.frame(gff)


if (!all(c("type", "ID", "seqnames", "start", "end", "strand") %in%
         names(gff_df))) {
  stop(
    "The imported GFF3 file is missing one or more required fields: ",
    "type, ID, seqnames, start, end, strand."
  )
}


genome <- BSgenome.Athaliana.TAIR.TAIR9

chr_lengths <- GenomeInfoDb::seqlengths(genome)
chr_available <- names(chr_lengths)


standardize_chr <- function(chr) {
  chr <- as.character(chr)
  
  map <- c(
    "1" = "Chr1",
    "2" = "Chr2",
    "3" = "Chr3",
    "4" = "Chr4",
    "5" = "Chr5",
    "Mt" = "ChrM",
    "Pt" = "ChrC",
    "M" = "ChrM",
    "C" = "ChrC"
  )
  
  out <- ifelse(
    chr %in% chr_available,
    chr,
    unname(map[chr])
  )
  
  out[!(out %in% chr_available)] <- NA_character_
  
  out
}


# ----------------------------- 6. Extract full 1-kb promoters ----------

gene_df <- gff_df %>%
  filter(
    tolower(type) == "gene"
  ) %>%
  transmute(
    gene = clean_gene_id(ID),
    chr = standardize_chr(seqnames),
    gene_start = as.integer(start),
    gene_end = as.integer(end),
    strand = as.character(strand)
  ) %>%
  filter(
    !is.na(gene),
    gene != "",
    !is.na(chr),
    !is.na(gene_start),
    !is.na(gene_end),
    gene_start >= 1L,
    gene_end >= gene_start,
    strand %in% c("+", "-")
  )


if (anyDuplicated(gene_df$gene)) {
  duplicated_ids <- unique(
    gene_df$gene[duplicated(gene_df$gene)]
  )
  
  warning(
    "Duplicate gene IDs in the GFF3 annotation; retaining one entry per ID. ",
    "Examples: ",
    paste(head(duplicated_ids, 10), collapse = ", ")
  )
  
  gene_df <- gene_df %>%
    distinct(
      gene,
      .keep_all = TRUE
    )
}


gene_df <- gene_df %>%
  mutate(
    chr_length = as.integer(
      chr_lengths[chr]
    ),
    
    tss = if_else(
      strand == "+",
      gene_start,
      gene_end
    ),
    
    promoter_start = if_else(
      strand == "+",
      tss - promoter_length,
      tss + 1L
    ),
    
    promoter_end = if_else(
      strand == "+",
      tss - 1L,
      tss + promoter_length
    )
  ) %>%
  filter(
    !is.na(chr_length),
    promoter_start >= 1L,
    promoter_end <= chr_length,
    promoter_end - promoter_start + 1L == promoter_length
  )


if (nrow(gene_df) == 0L) {
  stop(
    "No full-length 1-kb promoters could be extracted. ",
    "Check the GFF3 chromosome names and genome assembly."
  )
}


promoter_gr <- GenomicRanges::GRanges(
  seqnames = gene_df$chr,
  ranges = IRanges::IRanges(
    start = gene_df$promoter_start,
    end = gene_df$promoter_end
  ),
  strand = gene_df$strand
)


S4Vectors::mcols(promoter_gr)$gene <- gene_df$gene
names(promoter_gr) <- gene_df$gene


promoter_seqs <- Biostrings::getSeq(
  genome,
  promoter_gr
)


names(promoter_seqs) <- gene_df$gene


message("Genes with full-length 1-kb promoters: ",
        length(promoter_seqs))


write_gene_list(
  names(promoter_seqs),
  paste0(output_prefix, "_genes_with_full_promoters.csv")
)


# ----------------------------- 7. Scan promoter motifs -----------------

motif_presence <- tibble(
  gene = names(promoter_seqs)
)


for (i in seq_len(nrow(motifs))) {
  
  element_name <- as.character(
    motifs$Element[i]
  )
  
  pattern <- Biostrings::DNAString(
    as.character(motifs$Pattern[i])
  )
  
  pattern_rc <- Biostrings::reverseComplement(
    pattern
  )
  
  
  forward_hits <- Biostrings::vcountPattern(
    pattern,
    promoter_seqs,
    fixed = FALSE
  ) > 0L
  
  
  reverse_hits <- Biostrings::vcountPattern(
    pattern_rc,
    promoter_seqs,
    fixed = FALSE
  ) > 0L
  
  
  motif_presence[[element_name]] <- (
    forward_hits | reverse_hits
  )
}


presence_genes <- motif_presence$gene


write.csv(
  motif_presence,
  paste0(output_prefix, "_promoter_motif_presence.csv"),
  row.names = FALSE
)


# ----------------------------- 8. Enrichment tests ---------------------

results_list <- list()


for (group_name in group_order) {
  
  group_genes <- gene_sets_list[[group_name]]
  
  
  foreground <- intersect(
    group_genes,
    presence_genes
  )
  
  
  # Comparison set: every other gene with a valid full-length promoter.
  background <- setdiff(
    presence_genes,
    foreground
  )
  
  
  if (
    length(foreground) == 0L ||
    length(background) == 0L
  ) {
    
    message(
      "Group '", group_name,
      "' has no usable foreground or background genes; ",
      "motif tests will be marked NA."
    )
  }
  
  
  fg_rows <- match(
    foreground,
    motif_presence$gene
  )
  
  bg_rows <- match(
    background,
    motif_presence$gene
  )
  
  
  for (i in seq_len(nrow(motifs))) {
    
    element_name <- as.character(
      motifs$Element[i]
    )
    
    present_vec <- motif_presence[[element_name]]
    
    
    a <- if (length(fg_rows) > 0L) {
      sum(present_vec[fg_rows], na.rm = TRUE)
    } else {
      0L
    }
    
    
    b <- length(fg_rows) - a
    
    
    c <- if (length(bg_rows) > 0L) {
      sum(present_vec[bg_rows], na.rm = TRUE)
    } else {
      0L
    }
    
    
    d <- length(bg_rows) - c
    
    
    contingency <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("Gene_set", "Genome_background"),
        c("Motif_present", "Motif_absent")
      )
    )
    
    
    p_value <- if (
      length(fg_rows) > 0L &&
      length(bg_rows) > 0L
    ) {
      
      tryCatch(
        stats::fisher.test(
          contingency,
          alternative = "greater"
        )$p.value,
        
        error = function(e) NA_real_
      )
      
    } else {
      NA_real_
    }
    
    
    raw_or <- if (
      length(fg_rows) == 0L ||
      length(bg_rows) == 0L
    ) {
      
      NA_real_
      
    } else if ((b * c) == 0) {
      
      if ((a * d) > 0) {
        Inf
      } else {
        NA_real_
      }
      
    } else {
      
      (a * d) / (b * c)
    }
    
    
    # Haldane-Anscombe correction for displaying ORs when a cell is zero.
    or_display <- if (
      length(fg_rows) == 0L ||
      length(bg_rows) == 0L
    ) {
      
      NA_real_
      
    } else if (any(c(a, b, c, d) == 0)) {
      
      ((a + 0.5) * (d + 0.5)) /
        ((b + 0.5) * (c + 0.5))
      
    } else {
      
      raw_or
    }
    
    
    results_list[[length(results_list) + 1L]] <- tibble(
      GeneGroup = group_name,
      Element = element_name,
      Pattern = as.character(motifs$Pattern[i]),
      ElementClass = as.character(motifs$ElementClass[i]),
      
      n_foreground = length(fg_rows),
      foreground_with_motif = a,
      foreground_without_motif = b,
      
      n_background = length(bg_rows),
      background_with_motif = c,
      background_without_motif = d,
      
      odds_ratio_raw = raw_or,
      odds_ratio_display = or_display,
      p_value = p_value
    )
  }
}


enrichment <- bind_rows(results_list)


# ----------------------------- 9. BH correction and labels --------------

# Apply BH correction across all non-missing gene-set x motif p-values.
# The significance labels below are based only on p_adj_BH.
enrichment$p_adj_BH <- NA_real_

valid_p <- !is.na(enrichment$p_value)

enrichment$p_adj_BH[valid_p] <- p.adjust(
  enrichment$p_value[valid_p],
  method = "BH"
)


enrichment <- enrichment %>%
  mutate(
    significance = p_to_stars(p_adj_BH),
    
    GeneGroup = factor(
      GeneGroup,
      levels = rev(group_order)
    ),
    
    Element = factor(
      Element,
      levels = as.character(motifs$Element)
    ),
    
    ElementClass = factor(
      ElementClass,
      levels = c(
        "Candidate element",
        "Negative control"
      )
    ),
    
    odds_ratio_label = case_when(
      is.finite(odds_ratio_display) ~
        sprintf("%.2f", odds_ratio_display),
      
      is.infinite(odds_ratio_display) ~
        "Inf",
      
      TRUE ~
        "NA"
    ),
    
    tile_label = paste0(
      odds_ratio_label,
      "\n",
      significance
    )
  )


write.csv(
  enrichment,
  paste0(output_prefix, "_results.csv"),
  row.names = FALSE,
  quote = TRUE
)


message("BH-adjusted significant tests (adjp < ", alpha_cutoff, "): ",
        sum(enrichment$p_adj_BH < alpha_cutoff, na.rm = TRUE))


# ----------------------------- 10. Plot heatmap -------------------------

p_cis <- ggplot(
  enrichment,
  aes(
    x = Element,
    y = GeneGroup,
    fill = odds_ratio_display
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = tile_label),
    size = 5.2,
    lineheight = 0.95,
    fontface = "bold"
  ) +
  facet_grid(
    . ~ ElementClass,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#5878B7",
    mid = "#F7F7F7",
    high = "#E31A1C",
    midpoint = 1,
    limits = c(0, or_color_max),
    oob = scales::squish,
    na.value = "grey92",
    breaks = c(0, 0.5, 1, 1.5, 2, 3),
    name = "Odds ratio"
  ) +
  scale_x_discrete(
    drop = TRUE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    title = "Cis-regulatory element enrichment in Col-0 HSR gene sets",
    x = NULL,
    y = NULL,
    caption = paste0(
      "1-kb promoters; one-sided Fisher's exact test. ",
      "Stars indicate BH-adjusted p-values across the gene-set × motif tests: ",
      "* adjp < 0.05; ** adjp < 0.01; *** adjp < 0.001."
    )
  ) +
  theme_bw(
    base_size = 14
  ) +
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey35",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      colour = "black",
      size = 13
    ),
    
    axis.text.y = element_text(
      colour = "black",
      size = 14
    ),
    
    panel.grid = element_blank(),
    
    panel.spacing.x = grid::unit(
      0.25,
      "lines"
    ),
    
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    
    plot.caption = element_text(
      size = 11,
      hjust = 0
    ),
    
    plot.margin = margin(
      12,
      14,
      12,
      12
    )
  )


print(p_cis)


# ----------------------------- 11. Save figure -------------------------

png_file <- paste0(
  output_prefix,
  ".png"
)

pdf_file <- paste0(
  output_prefix,
  ".pdf"
)


ggsave(
  filename = png_file,
  plot = p_cis,
  width = 16,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  filename = pdf_file,
  plot = p_cis,
  width = 16,
  height = 8.5,
  units = "in",
  bg = "white"
)


# ----------------------------- 12. Final report ------------------------

message("")
message("Cis-regulatory element enrichment completed.")
message("Gene-set source: ", gene_set_source)
message("GFF3 annotation: ", annotation_file)
message("Promoter length: ", promoter_length, " bp")
message("BH significance threshold: ", alpha_cutoff)
message("Significance labels: BH-adjusted p-values")
message("")
message("Saved HSR gene-set membership: ",
        paste0(output_prefix, "_gene_set_membership.csv"))
message("Saved gene-set sizes: ",
        paste0(output_prefix, "_gene_set_sizes.csv"))
message("Saved genes with valid promoters: ",
        paste0(output_prefix, "_genes_with_full_promoters.csv"))
message("Saved promoter motif presence: ",
        paste0(output_prefix, "_promoter_motif_presence.csv"))
message("Saved enrichment results: ",
        paste0(output_prefix, "_results.csv"))
message("Saved PNG: ", png_file)
message("Saved PDF: ", pdf_file)
