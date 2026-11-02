############################################################################
# NMD cis-feature comparison: Cluster 5 subclusters vs background
#
# Panels:
#   A) Maximum annotated 3'UTR length across isoforms
#   B) 3'UTR intron status in the representative transcript per gene
#   C) RiboTaper-supported uORF prevalence; no 3'UTR annotation filter
#      is applied to this panel.
#
# Each Cluster 5 subcluster (5/1, 5/2, 5/3) is compared with Bg, defined as
# the selected analysis universe excluding all genes assigned to these three
# subclusters. BH-adjusted p-values are used for significance labels by default.
# RiboTaper-positive genes remain positive even below the TPM threshold;
# non-positive genes are evaluable negatives only when TPM >= ribo_tpm_min.
############################################################################

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(readxl)
})

# ----------------------------- Settings -------------------------------
parent_cluster <- 5
k_subclusters <- 3
annotation_file <- "Araport11_GFF3_genes_transposons.20241001.gff.gz"
ribotaper_file  <- "Dataset S1_RiboTaper.xlsx"
cips_file       <- "Dataset S5_CiPS.xlsx"
ribo_expr_file  <- "GSE183264_Ribo.genes.results"  # required for evaluable negatives
ribo_tpm_min    <- 1

# Background options:
#   "gff_annotated" = all GFF-annotated genes (genome-wide annotation background).
#   "dds_expressed" = GFF genes present in dds and with normalized count > 0
#                     in at least one sample (requires an existing dds object).
# Keep this choice explicit because it changes the biological comparison.
background_mode <- "dds_expressed"

# "BH" uses BH-adjusted p-values for significance labels; "raw" uses unadjusted p-values.
significance_basis <- "BH"

subcluster_file <- paste0(
  "NMD_cluster", parent_cluster, "_subclusters_k", k_subclusters,
  "_gene_assignments.csv"
)
out_prefix <- paste0("NMD_cluster", parent_cluster,
                     "_cis_feature_comparison_updated")
plot_group_order <- c("Bg", "5/1", "5/2", "5/3")
group_colors <- c("Bg" = "grey75", "5/1" = "#E60012",
                  "5/2" = "#164A91", "5/3" = "#26934F")
slide_width_cm <- 19.05
slide_height_cm <- 25.40
individual_panel_fraction <- 1 / 6
slide_area_in2 <- (slide_width_cm / 2.54) * (slide_height_cm / 2.54)
individual_panel_area_in2 <- slide_area_in2 * individual_panel_fraction
individual_plot_aspect_ratio <- 4.50 / 2.78
individual_plot_width_in <- sqrt(individual_panel_area_in2 * individual_plot_aspect_ratio)
individual_plot_height_in <- sqrt(individual_panel_area_in2 / individual_plot_aspect_ratio)

# ----------------------------- Helpers --------------------------------
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

clean_tx_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(transcript|mRNA|rna):", "", x, ignore.case = TRUE)
  trimws(x)
}

# Expand Parent attributes so each returned row has exactly one Parent.
expand_parent_rows <- function(df) {
  if (nrow(df) == 0L) return(df[0, , drop = FALSE])
  if (!"Parent" %in% names(df)) stop("The imported GFF has no 'Parent' column.")
  parent_values <- lapply(seq_len(nrow(df)), function(i) {
    value <- df$Parent[[i]]
    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
      return(character(0))
    }
    value <- paste(as.character(value), collapse = ",")
    vals <- trimws(unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE))
    vals[nzchar(vals) & !is.na(vals)]
  })
  keep <- lengths(parent_values) > 0L
  if (!any(keep)) return(df[0, , drop = FALSE])
  source_rows <- which(keep)
  out <- df[rep(source_rows, lengths(parent_values)), , drop = FALSE]
  out$Parent <- unlist(parent_values[keep], use.names = FALSE)
  rownames(out) <- NULL
  out
}

p_to_stars <- function(p) {
  if (length(p) == 0L || is.na(p) || !is.finite(p)) return("n/a")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

extract_agi <- function(df, col = NULL) {
  if (!is.null(col)) {
    if (!col %in% names(df)) {
      stop("Expected column '", col, "' not found. Available columns: ",
           paste(names(df), collapse = ", "))
    }
    df <- df[, col, drop = FALSE]
  }
  ids <- unlist(lapply(df, function(x) {
    if (is.character(x) || is.factor(x)) as.character(x) else character(0)
  }), use.names = FALSE)
  matches <- regmatches(ids, gregexpr("AT[1-5CM]G[0-9]{5}", ids, ignore.case = TRUE))
  unique(toupper(unlist(matches, use.names = FALSE)))
}

read_excel_table <- function(path, sheet, skip = 0) {
  as.data.frame(readxl::read_excel(path, sheet = sheet, skip = skip))
}

# ----------------------------- 1. Check inputs -------------------------
required_files <- c(annotation_file, subcluster_file, ribotaper_file,
                    cips_file, ribo_expr_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required input file(s) not found:\n - ",
       paste(missing_files, collapse = "\n - "),
       "\nCurrent working directory: ", getwd(),
       "\nThe Ribo-seq TPM table is required to define evaluable uORF-negative genes.")
}
if (length(parent_cluster) != 1L || is.na(parent_cluster) || parent_cluster < 1 ||
    parent_cluster != as.integer(parent_cluster)) {
  stop("parent_cluster must be a single positive integer.")
}
if (length(k_subclusters) != 1L || is.na(k_subclusters) || k_subclusters < 1 ||
    k_subclusters != as.integer(k_subclusters)) {
  stop("k_subclusters must be a single positive integer.")
}
if (!background_mode %in% c("dds_expressed", "gff_annotated")) {
  stop("background_mode must be 'dds_expressed' or 'gff_annotated'.")
}
if (length(significance_basis) != 1L || !significance_basis %in% c("raw", "BH")) {
  stop("significance_basis must be 'raw' or 'BH'.")
}
if (length(ribo_tpm_min) != 1L || is.na(ribo_tpm_min) ||
    !is.finite(ribo_tpm_min) || ribo_tpm_min < 0) {
  stop("ribo_tpm_min must be a single non-negative finite number.")
}
if (background_mode == "dds_expressed" && !exists("dds", inherits = TRUE)) {
  stop("background_mode='dds_expressed' requires the existing DESeq2 object 'dds'.")
}
if (background_mode == "dds_expressed" && !inherits(dds, "DESeqDataSet")) {
  stop("Object 'dds' must be a DESeqDataSet when background_mode='dds_expressed'.")
}

subclusters_raw <- read.csv(subcluster_file, stringsAsFactors = FALSE,
                            check.names = FALSE)
required_cols <- c("gene", "ParentCluster", "Subcluster")
if (!all(required_cols %in% names(subclusters_raw))) {
  stop("Subcluster file must contain: ", paste(required_cols, collapse = ", "))
}
subclusters_raw <- subclusters_raw %>%
  mutate(
    gene = clean_gene_id(gene),
    ParentCluster = suppressWarnings(as.integer(as.character(ParentCluster))),
    Subcluster = suppressWarnings(as.integer(as.character(Subcluster)))
  )
if (anyNA(subclusters_raw$gene) || any(!nzchar(subclusters_raw$gene))) {
  stop("The subcluster assignment file contains missing or empty gene IDs.")
}
if (anyNA(subclusters_raw$ParentCluster) || anyNA(subclusters_raw$Subcluster)) {
  stop("ParentCluster or Subcluster contains missing/non-integer values.")
}
subclusters_raw <- subclusters_raw %>%
  filter(ParentCluster == parent_cluster,
         Subcluster %in% seq_len(k_subclusters)) %>%
  distinct(gene, Subcluster)

if (!nrow(subclusters_raw)) {
  stop("No genes found for parent Cluster ", parent_cluster,
       ". Check the subcluster assignment file.")
}
conflicting_assignments <- subclusters_raw %>% count(gene, name = "n_assignments") %>%
  filter(n_assignments > 1L)
if (nrow(conflicting_assignments)) {
  stop("Some genes are assigned to more than one selected subcluster (examples: ",
       paste(head(conflicting_assignments$gene, 10), collapse = ", "), ").")
}
if (!all(seq_len(k_subclusters) %in% subclusters_raw$Subcluster)) {
  stop("The assignment file does not contain all subclusters 1, 2, and 3.")
}

# ----------------------------- 2. Import GFF and universe --------------
gff <- rtracklayer::import(annotation_file)
gff_df <- as.data.frame(gff)
if (!all(c("type", "ID", "start", "end") %in% names(gff_df))) {
  stop("GFF does not contain the expected type/ID/start/end fields.")
}
if (!"Parent" %in% names(gff_df)) {
  stop("Imported GFF does not contain a Parent column needed for transcript mapping.")
}

gff_gene_table <- gff_df %>%
  filter(tolower(type) == "gene") %>%
  transmute(gene = clean_gene_id(ID)) %>%
  filter(!is.na(gene), nzchar(gene))
if (anyDuplicated(gff_gene_table$gene)) {
  dup <- unique(gff_gene_table$gene[duplicated(gff_gene_table$gene)])
  stop("Duplicate gene IDs after cleaning in GFF gene features (examples: ",
       paste(head(dup, 10), collapse = ", "), ").")
}
gff_genes <- gff_gene_table$gene
if (!length(gff_genes)) stop("No gene features were found in the imported GFF.")

if (background_mode == "dds_expressed") {
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  count_gene_ids <- clean_gene_id(rownames(norm_counts))
  if (anyNA(count_gene_ids) || any(!nzchar(count_gene_ids))) {
    stop("The dds object contains missing/empty gene IDs after normalization.")
  }
  if (anyDuplicated(count_gene_ids)) {
    dup <- unique(count_gene_ids[duplicated(count_gene_ids)])
    stop("Gene IDs in dds are duplicated after normalization (examples: ",
         paste(head(dup, 10), collapse = ", "), ").")
  }
  rownames(norm_counts) <- count_gene_ids
  dds_expressed_genes <- rownames(norm_counts)[
    rowSums(norm_counts > 0, na.rm = TRUE) > 0
  ] %>% unique()
  gene_universe <- intersect(gff_genes, dds_expressed_genes)
  message("Background mode: GFF genes present in dds and expressed in at least one sample.")
} else {
  gene_universe <- gff_genes
  message("Background mode: all GFF-annotated genes.")
}
if (!length(gene_universe)) stop("The selected analysis universe contains no genes.")

subcluster_genes_original <- unique(subclusters_raw$gene)
not_in_universe <- setdiff(subcluster_genes_original, gene_universe)
if (length(not_in_universe)) {
  warning(length(not_in_universe),
          " subcluster genes are not in the selected analysis universe and will be omitted.")
}
subclusters <- subclusters_raw %>% filter(gene %in% gene_universe)
if (!nrow(subclusters)) stop("No Cluster 5 subcluster genes remain in the selected analysis universe.")
if (!all(seq_len(k_subclusters) %in% subclusters$Subcluster)) {
  stop("After applying the analysis universe, one or more subclusters are empty.")
}
subcluster_genes <- unique(subclusters$gene)
background_genes <- setdiff(gene_universe, subcluster_genes)
if (!length(background_genes)) stop("No background genes remain after excluding all three subclusters.")

genes_by_group <- c(
  list("Bg" = background_genes),
  list(
    "5/1" = unique(subclusters$gene[subclusters$Subcluster == 1]),
    "5/2" = unique(subclusters$gene[subclusters$Subcluster == 2]),
    "5/3" = unique(subclusters$gene[subclusters$Subcluster == 3])
  )
)

message("Genes assigned to each Cluster ", parent_cluster, " subcluster in the analysis universe:")
print(subclusters %>% count(Subcluster, name = "n_genes"))
message("Analysis-universe genes: ", length(gene_universe))
message("Background genes after excluding all three subclusters: ", length(background_genes))

# ----------------------------- 3. Map transcripts to genes -------------
mrna_raw <- gff_df %>%
  filter(tolower(type) %in% c("mrna", "transcript")) %>%
  dplyr::select(ID, Parent)
mrna_expanded <- expand_parent_rows(mrna_raw) %>%
  transmute(transcript = clean_tx_id(ID), gene = clean_gene_id(Parent)) %>%
  filter(!is.na(transcript), nzchar(transcript), !is.na(gene), nzchar(gene)) %>%
  distinct(transcript, gene)

map_tx_gene <- function(tx) {
  mapped <- mrna_expanded$gene[match(tx, mrna_expanded$transcript)]
  # Fallback only when the transcript ID itself clearly embeds an AGI gene ID.
  is_agi_tx <- grepl("^AT[1-5CM]G[0-9]{5}(\\.[0-9]+)?$", tx, ignore.case = TRUE)
  fallback <- ifelse(is_agi_tx, clean_gene_id(tx), NA_character_)
  ifelse(!is.na(mapped) & nzchar(mapped), mapped, fallback)
}

# ----------------------------- 4. 3'UTR length / representative isoform
# Panel A: maximum annotated 3'UTR length across isoforms.
# Panel B: longest-CDS transcript among those with an annotated 3'UTR;
# ties are resolved by longer 3'UTR and then transcript ID.

utr_raw <- gff_df %>%
  filter(tolower(type) == "three_prime_utr") %>%
  dplyr::select(Parent, start, end)
if (nrow(utr_raw) == 0L) stop("No three_prime_UTR features were found in the GFF.")

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

if (!nrow(utr_expanded)) stop("No valid 3'UTR segments could be mapped to genes.")

utr_by_tx <- utr_expanded %>%
  group_by(gene, transcript) %>%
  summarise(
    utr3_length = sum(width, na.rm = TRUE),
    n_utr_segments = n_distinct(paste(start, end, sep = ":")),
    .groups = "drop"
  )

gene_utr_max <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    max_utr3_length = max(utr3_length, na.rm = TRUE),
    n_annotated_utr_isoforms = n_distinct(transcript),
    .groups = "drop"
  ) %>%
  filter(is.finite(max_utr3_length), max_utr3_length > 0)

cds_raw <- gff_df %>%
  filter(tolower(type) == "cds") %>%
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

representative_tx <- utr_by_tx %>%
  left_join(cds_by_tx, by = "transcript") %>%
  mutate(cds_length = tidyr::replace_na(cds_length, 0)) %>%
  arrange(gene, desc(cds_length), desc(utr3_length), transcript) %>%
  group_by(gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  transmute(
    gene,
    representative_transcript = transcript,
    representative_utr3_length = utr3_length,
    representative_utr_segments = n_utr_segments,
    representative_cds_length = cds_length
  )

gene_utr_length <- gene_utr_max %>%
  left_join(representative_tx, by = "gene")

utr_intron_any <- utr_by_tx %>%
  group_by(gene) %>%
  summarise(
    has_3utr_intron_any = any(n_utr_segments > 1, na.rm = TRUE),
    .groups = "drop"
  )

gene_features_3utr <- representative_tx %>%
  transmute(
    gene,
    representative_transcript,
    has_3utr_intron = as.logical(representative_utr_segments > 1),
    representative_utr3_length,
    representative_utr_segments,
    representative_cds_length
  ) %>%
  left_join(utr_intron_any, by = "gene")

# ----------------------------- 5. uORF evidence + detectability --------
ribotaper_tab <- read_excel_table(ribotaper_file, sheet = "S1A_uORFs", skip = 0)
cips_tab <- read_excel_table(cips_file, sheet = "A.ALL", skip = 2)
ribotaper_genes <- extract_agi(ribotaper_tab, "gene_id")
cips_genes <- extract_agi(cips_tab, "gene_id")
if (!length(ribotaper_genes)) stop("No AGI gene IDs extracted from the RiboTaper table.")
if (!length(cips_genes)) warning("No AGI gene IDs extracted from the CiPS table; CiPS robustness analysis may be empty.")

if (!"AA" %in% names(cips_tab)) {
  warning("CiPS table has no 'AA' column; the >=10 aa sensitivity set will be empty.")
  cips_10aa_genes <- character(0)
} else {
  cips_10aa_tab <- cips_tab %>%
    mutate(.AA_numeric = suppressWarnings(as.numeric(as.character(AA)))) %>%
    filter(!is.na(.AA_numeric), .AA_numeric >= 10)
  cips_10aa_genes <- extract_agi(cips_10aa_tab, "gene_id")
}

uorf_external_sets <- list(
  "Ribo-seq uORF (RiboTaper)" = unique(clean_gene_id(ribotaper_genes)),
  "Ribo-seq uORF (CiPS)" = unique(clean_gene_id(cips_genes)),
  "Ribo-seq uORF (CiPS, >=10 aa)" = unique(clean_gene_id(cips_10aa_genes))
)
if (length(ribotaper_genes) && length(cips_genes)) {
  uorf_external_sets[["Ribo-seq uORF (RiboTaper or CiPS)"]] <-
    union(unique(clean_gene_id(ribotaper_genes)), unique(clean_gene_id(cips_genes)))
}
main_uorf_definition <- "Ribo-seq uORF (RiboTaper)"
main_uorf_col <- make.names(main_uorf_definition)

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
  pull(gene) %>% unique()
if (!length(ribo_detectable_genes)) {
  stop("No genes pass the Ribo-seq detectability threshold; check the TPM file and threshold.")
}

# For each uORF definition: positives remain TRUE even below the TPM threshold;
# only TPM-detectable non-positives are classified FALSE; all others are NA.
uorf_df <- tibble(GENEID = gene_universe)
for (nm in names(uorf_external_sets)) {
  target_set <- unique(uorf_external_sets[[nm]])
  flag <- rep(NA, nrow(uorf_df))
  if (length(target_set) > 0L) {
    positive <- uorf_df$GENEID %in% target_set
    detectable <- uorf_df$GENEID %in% ribo_detectable_genes
    flag[detectable & !positive] <- FALSE
    flag[positive] <- TRUE
  } else {
    message("No IDs were parsed for uORF definition '", nm,
            "'; its flags remain NA instead of assigning unsupported negatives.")
  }
  uorf_df[[make.names(nm)]] <- flag
}
uorf_df$has_uORF <- uorf_df[[main_uorf_col]]

ribo_positive_not_detectable <- intersect(
  ribotaper_genes, setdiff(gene_universe, ribo_detectable_genes)
)
message("RiboTaper unique genes: ", length(ribotaper_genes),
        "; TPM-detectable genes (TPM >= ", ribo_tpm_min, "): ", length(ribo_detectable_genes),
        "; RiboTaper positives below TPM threshold retained as positive: ",
        length(ribo_positive_not_detectable))
message("CiPS unique genes: ", length(cips_genes),
        "; CiPS >=10 aa unique genes: ", length(cips_10aa_genes))

# ----------------------------- 6. Define analysis groups --------------
# Bg excludes all genes from the three selected subclusters, so groups do not overlap.

# ----------------------------- 7. Assemble panel data -----------------
all_group_genes <- bind_rows(lapply(names(genes_by_group), function(grp) {
  tibble(gene = genes_by_group[[grp]], Group = grp)
})) %>%
  distinct(Group, gene) %>%
  mutate(Group = factor(Group, levels = plot_group_order))

# Panel A: maximum annotated 3'UTR per gene; representative length is a robustness output.
utr_data <- all_group_genes %>%
  inner_join(gene_utr_length, by = "gene") %>%
  filter(is.finite(max_utr3_length), max_utr3_length > 0)
rep_utr_data <- all_group_genes %>%
  inner_join(gene_utr_length, by = "gene") %>%
  filter(is.finite(representative_utr3_length), representative_utr3_length > 0)

# Panel B: representative-transcript intron call plus any-isoform status as a robustness check.
intron_gene_status <- gene_features_3utr %>%
  dplyr::select(gene, has_3utr_intron, has_3utr_intron_any) %>%
  filter(!is.na(has_3utr_intron))

calculate_percentage <- function(feature_status, flag_col, feature_name) {
  bind_rows(lapply(names(genes_by_group), function(grp) {
    genes <- unique(genes_by_group[[grp]])
    evaluable <- tibble(gene = genes) %>%
      left_join(feature_status, by = "gene") %>%
      filter(!is.na(.data[[flag_col]]))
    n_total <- nrow(evaluable)
    n_positive <- if (n_total > 0L) sum(evaluable[[flag_col]] %in% TRUE) else 0L
    tibble(
      Feature = feature_name,
      Group = grp,
      n_total = n_total,
      n_positive = n_positive,
      n_negative = n_total - n_positive,
      percentage = if (n_total > 0L) 100 * n_positive / n_total else NA_real_
    )
  })) %>% mutate(Group = factor(Group, levels = plot_group_order))
}

intron_data <- calculate_percentage(
  intron_gene_status, "has_3utr_intron", "3'UTR intron (representative transcript)"
)
intron_any_data <- calculate_percentage(
  intron_gene_status, "has_3utr_intron_any", "3'UTR intron (any isoform)"
)

# Panel C has no 3'UTR annotation filter. Evaluable genes have a positive call
# or pass the Ribo-seq TPM threshold, as defined above.
uorf_gene_status <- uorf_df %>% transmute(gene = GENEID, has_uORF = as.logical(has_uORF))
uorf_data <- calculate_percentage(uorf_gene_status, "has_uORF", "RiboTaper uORF")

# Report panel-specific denominators so differing coverage is transparent.
coverage_table <- bind_rows(lapply(names(genes_by_group), function(grp) {
  genes <- unique(genes_by_group[[grp]])
  tibble(
    Group = grp,
    n_genes_in_group = length(genes),
    n_with_annotated_3UTR = sum(genes %in% gene_utr_length$gene),
    n_with_representative_intron_call = sum(genes %in% gene_features_3utr$gene),
    n_RiboTaper_evaluable = sum(!is.na(uorf_df$has_uORF[match(genes, uorf_df$GENEID)])),
    n_RiboTaper_positive = sum(uorf_df$has_uORF[match(genes, uorf_df$GENEID)] %in% TRUE, na.rm = TRUE)
  )
}))

# ----------------------------- 8. Statistical tests -------------------
# Three planned comparisons per feature family: each subcluster vs Bg.
# BH adjustment is performed within each feature family across these three tests.
get_wilcox_tests <- function(dat, value_col, metric_label) {
  bg <- dat[[value_col]][dat$Group == "Bg"]
  out <- bind_rows(lapply(c("5/1", "5/2", "5/3"), function(grp) {
    sg <- dat[[value_col]][dat$Group == grp]
    bg_use <- bg[is.finite(bg) & bg > 0]
    sg_use <- sg[is.finite(sg) & sg > 0]
    p <- if (length(bg_use) >= 2L && length(sg_use) >= 2L) {
      tryCatch(stats::wilcox.test(bg_use, sg_use, exact = FALSE)$p.value,
               error = function(e) NA_real_)
    } else NA_real_
    tibble(
      Feature = metric_label,
      Group = factor(grp, levels = plot_group_order),
      n_background = length(bg_use),
      n_subcluster = length(sg_use),
      median_background = if (length(bg_use)) median(bg_use) else NA_real_,
      median_subcluster = if (length(sg_use)) median(sg_use) else NA_real_,
      p_value = p
    )
  }))
  out %>% mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    p_for_label = if (significance_basis == "BH") p_adj_BH else p_value,
    significance = vapply(p_for_label, p_to_stars, character(1))
  )
}

get_fisher_tests <- function(dat, feature_name) {
  bg <- dat %>% filter(Group == "Bg")
  bind_rows(lapply(c("5/1", "5/2", "5/3"), function(grp) {
    sg <- dat %>% filter(Group == grp)
    n1 <- if (nrow(bg)) bg$n_total[1] else 0L
    p1 <- if (nrow(bg)) bg$n_positive[1] else 0L
    n2 <- if (nrow(sg)) sg$n_total[1] else 0L
    p2 <- if (nrow(sg)) sg$n_positive[1] else 0L
    p <- or <- ci_low <- ci_high <- NA_real_
    if (n1 > 0L && n2 > 0L) {
      # Put the subcluster in row 1 so fisher.test's OR and CI are
      # oriented as subcluster odds / background odds.
      tab <- matrix(c(p2, n2 - p2, p1, n1 - p1), nrow = 2, byrow = TRUE,
                    dimnames = list(c(grp, "Bg"), c("positive", "negative")))
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
      Feature = feature_name,
      Group = factor(grp, levels = plot_group_order),
      n_background = n1,
      background_positive = p1,
      n_subcluster = n2,
      subcluster_positive = p2,
      odds_ratio_subcluster_vs_background = or,
      CI_low = ci_low,
      CI_high = ci_high,
      p_value = p
    )
  })) %>% mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    p_for_label = if (significance_basis == "BH") p_adj_BH else p_value,
    significance = vapply(p_for_label, p_to_stars, character(1))
  )
}

utr_tests <- get_wilcox_tests(utr_data, "max_utr3_length", "Maximum annotated 3'UTR length")
rep_utr_tests <- get_wilcox_tests(rep_utr_data, "representative_utr3_length", "Representative-transcript 3'UTR length")
intron_tests <- get_fisher_tests(intron_data, "3'UTR intron (representative transcript)")
intron_any_tests <- get_fisher_tests(intron_any_data, "3'UTR intron (any isoform)")
uorf_tests <- get_fisher_tests(uorf_data, "RiboTaper-supported uORF")

# Alternative uORF definitions are sensitivity analyses using the same evaluability rule.
uorf_robustness <- bind_rows(lapply(names(uorf_external_sets), function(definition) {
  flag_col <- make.names(definition)
  status <- uorf_df %>% transmute(gene = GENEID, flag = .data[[flag_col]])
  dat <- calculate_percentage(status, "flag", definition)
  get_fisher_tests(dat, definition) %>% mutate(definition = definition)
}))

# ----------------------------- 9. Plot settings ------------------------
base_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 4)),
    plot.title = element_blank(),
    plot.margin = margin(5, 6, 5, 5, unit = "pt")
  )

if (!nrow(utr_data)) stop("No genes with positive annotated 3'UTR lengths remain for Panel A.")

# Panel A: one significance label per subcluster comparison vs Bg.
utr_max_by_group <- utr_data %>%
  group_by(Group) %>%
  summarise(max_value = max(max_utr3_length, na.rm = TRUE), .groups = "drop")

utr_label_positions <- utr_tests %>%
  left_join(utr_max_by_group, by = "Group") %>%
  mutate(
    y_position = max_value * 1.18,
    label_y = y_position * 1.05
  )
len_y_upper_candidates <- c(utr_data$max_utr3_length, utr_label_positions$label_y)
len_y_upper <- max(len_y_upper_candidates[is.finite(len_y_upper_candidates)], na.rm = TRUE) * 1.08
len_y_upper <- max(len_y_upper, 100)

pA <- ggplot(utr_data, aes(x = Group, y = max_utr3_length, fill = Group)) +
  geom_boxplot(width = 0.68, outlier.size = 0.7, outlier.alpha = 0.25,
               linewidth = 0.55, na.rm = TRUE) +
  geom_text(data = utr_label_positions,
            aes(x = Group, y = label_y, label = significance),
            inherit.aes = FALSE, size = 6.0, fontface = "bold", na.rm = TRUE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_x_discrete(limits = plot_group_order, drop = FALSE) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(1, len_y_upper), clip = "off") +
  labs(title = "A. Maximum annotated 3'UTR length",
       x = NULL, y = "Maximum 3'UTR length (bp, log scale)") +
  base_theme

make_percentage_plot <- function(dat, y_title, panel_title) {
  if (!nrow(dat)) stop("No data available for: ", panel_title)
  sig_pos <- get_fisher_tests(dat, unique(dat$Feature)) %>%
    left_join(dat %>% dplyr::select(Group, percentage), by = "Group") %>%
    mutate(
      y_position = pmin(percentage + 8, 115),
      label_y = y_position + 1.5
    )
  candidates <- c(dat$percentage, sig_pos$label_y)
  candidates <- candidates[is.finite(candidates)]
  bar_y_max <- if (length(candidates)) min(125, max(25, max(candidates) + 4)) else 25

  ggplot(dat, aes(x = Group, y = percentage, fill = Group)) +
    geom_col(width = 0.68, na.rm = TRUE) +
    geom_text(aes(label = ifelse(is.na(percentage), "", sprintf("%.1f%%", percentage))),
              vjust = -0.30, size = 5.0, fontface = "bold", na.rm = TRUE) +
    geom_text(data = sig_pos,
              aes(x = Group, y = label_y, label = significance),
              inherit.aes = FALSE, size = 6.0, fontface = "bold", na.rm = TRUE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_x_discrete(limits = plot_group_order, drop = FALSE) +
    scale_y_continuous(limits = c(0, bar_y_max),
                       breaks = seq(0, min(100, floor(bar_y_max / 10) * 10), by = 10),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = panel_title, x = NULL, y = y_title) +
    base_theme
}

pB <- make_percentage_plot(
  intron_data, "Genes with a 3'UTR intron (%)",
  "B. 3'UTR intron status"
)
pC <- make_percentage_plot(
  uorf_data, "Genes with RiboTaper-supported uORF (%)",
  "C. RiboTaper-supported uORF"
)

combined_plot <- pA + pB + pC +
  plot_layout(ncol = 3, widths = c(1.2, 1, 1))
print(combined_plot)

# ----------------------------- 10. Export tables and figures ------------
write.csv(subclusters, paste0(out_prefix, "_subcluster_assignments_used.csv"), row.names = FALSE)
write.csv(data.frame(gene = gene_universe), paste0(out_prefix, "_analysis_universe.csv"), row.names = FALSE)
write.csv(data.frame(gene = background_genes), paste0(out_prefix, "_background_genes.csv"), row.names = FALSE)
write.csv(gene_utr_length, paste0(out_prefix, "_utr_length_by_gene.csv"), row.names = FALSE)
write.csv(gene_features_3utr, paste0(out_prefix, "_representative_utr_intron_by_gene.csv"), row.names = FALSE)
write.csv(uorf_df, paste0(out_prefix, "_uorf_status_by_gene.csv"), row.names = FALSE)
write.csv(utr_data, paste0(out_prefix, "_utr_length_values.csv"), row.names = FALSE)
write.csv(rep_utr_data, paste0(out_prefix, "_representative_utr_length_robustness_values.csv"), row.names = FALSE)
write.csv(utr_tests, paste0(out_prefix, "_utr_length_tests.csv"), row.names = FALSE)
write.csv(rep_utr_tests, paste0(out_prefix, "_representative_utr_length_robustness_tests.csv"), row.names = FALSE)
write.csv(intron_data, paste0(out_prefix, "_utr_intron_percentages.csv"), row.names = FALSE)
write.csv(intron_tests, paste0(out_prefix, "_utr_intron_tests.csv"), row.names = FALSE)
write.csv(intron_any_data, paste0(out_prefix, "_utr_intron_any_isoform_robustness_percentages.csv"), row.names = FALSE)
write.csv(intron_any_tests, paste0(out_prefix, "_utr_intron_any_isoform_robustness_tests.csv"), row.names = FALSE)
write.csv(uorf_data, paste0(out_prefix, "_uorf_percentages.csv"), row.names = FALSE)
write.csv(uorf_tests, paste0(out_prefix, "_uorf_tests.csv"), row.names = FALSE)
write.csv(uorf_robustness, paste0(out_prefix, "_uorf_definition_robustness.csv"), row.names = FALSE)
write.csv(coverage_table, paste0(out_prefix, "_panel_evaluable_denominators.csv"), row.names = FALSE)

external_set_audit <- bind_rows(lapply(names(uorf_external_sets), function(nm) {
  s <- unique(clean_gene_id(uorf_external_sets[[nm]]))
  tibble(
    definition = nm,
    n_unique_external_genes = length(s),
    n_in_analysis_universe = length(intersect(s, gene_universe)),
    n_in_Ribo_detectable_set = length(intersect(s, ribo_detectable_genes)),
    n_in_both_universe_and_detectable = length(intersect(intersect(s, gene_universe), ribo_detectable_genes))
  )
}))
write.csv(external_set_audit, paste0(out_prefix, "_external_uORF_set_audit.csv"), row.names = FALSE)

ggsave(paste0(out_prefix, ".png"), combined_plot,
       width = 13.5, height = 4.7, units = "in", dpi = 300, bg = "white")
ggsave(paste0(out_prefix, ".pdf"), combined_plot,
       width = 13.5, height = 4.7, units = "in", bg = "white")
ggsave(paste0(out_prefix, "_panelA_UTR_length.png"), pA,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_panelA_UTR_length.pdf"), pA,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", bg = "white")
ggsave(paste0(out_prefix, "_panelB_UTR_intron.png"), pB,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_panelB_UTR_intron.pdf"), pB,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", bg = "white")
ggsave(paste0(out_prefix, "_panelC_uORF_RiboTaper.png"), pC,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", dpi = 300, bg = "white")
ggsave(paste0(out_prefix, "_panelC_uORF_RiboTaper.pdf"), pC,
       width = individual_plot_width_in, height = individual_plot_height_in,
       units = "in", bg = "white")

message("Finished. Output prefix: ", out_prefix)
message("Background mode: ", background_mode)
message("Main uORF definition: ", main_uorf_definition)
message("Panel C is not filtered by 3'UTR annotation; evaluability is based on a positive uORF call or TPM >= ", ribo_tpm_min, ".")
message("Significance stars use ", significance_basis, " p-values; raw and BH-adjusted p-values are exported.")
message("Panel-specific evaluable denominator counts:")
print(coverage_table)

