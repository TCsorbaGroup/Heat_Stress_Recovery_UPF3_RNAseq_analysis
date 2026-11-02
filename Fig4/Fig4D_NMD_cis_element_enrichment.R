# Cis-regulatory element enrichment for the three subclusters of original Cluster 5
#
# Analysis outline:
#   * Extract 1-kb promoters immediately upstream of each gene's TSS.
#   * Search each promoter for the candidate cis-elements and reference promoter motifs below.
#   * For each Cluster 5 subcluster, test motif presence against all other annotated
#     Arabidopsis gene promoters using a one-sided Fisher exact test (enrichment).
#   * Adjust all available p-values across the motif-by-subcluster tests by
#     Benjamini-Hochberg. Plot stars default to raw p-values to preserve the
#     reference-figure behavior; raw and adjusted p-values are saved.
#
# IMPORTANT:
#   These are explicit consensus-sequence scans, not position-weight-matrix (PWM)
#   binding predictions. The as-1 test uses its TGACG core half-site, so it does not
#   require a paired/tandem architecture. Interpret these results as motif-occurrence
#   enrichment, not proof of TF binding or element activity.
#
# Required input files in the working directory:
#   NMD_cluster5_subclusters_k3_gene_assignments.csv
#   Araport11_GFF3_genes_transposons.20241001.gff.gz
#
# Required Bioconductor packages (install if missing):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("rtracklayer", "GenomicRanges", "Biostrings",
#                        "BSgenome.Athaliana.TAIR.TAIR9"))
# install.packages(c("dplyr", "tibble", "ggplot2", "scales"))

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ----------------------------- Settings -----------------------------
parent_cluster <- 5
k_subclusters <- 3
promoter_length <- 1000
annotation_file <- "Araport11_GFF3_genes_transposons.20241001.gff.gz"
subcluster_file <- paste0(
  "NMD_cluster", parent_cluster, "_subclusters_k", k_subclusters,
  "_gene_assignments.csv"
)
out_prefix <- paste0("NMD_cluster", parent_cluster,
                     "_subclusters_k", k_subclusters, "_cis_element_enrichment")

# "raw" preserves the supplied reference-figure behavior; "BH" uses BH-adjusted
# p-values for stars and is preferable for interpreting statistical significance
# across the motif-by-subcluster test family.
significance_basis <- "BH"
# Upper color limit for tile fill. Values above this are squished to the maximum color;
# the exact odds ratios are retained in the output CSV.
or_color_max <- 3

if (length(promoter_length) != 1L || is.na(promoter_length) ||
    !is.finite(promoter_length) || promoter_length < 1 ||
    promoter_length != as.integer(promoter_length)) {
  stop("promoter_length must be a single positive integer.")
}
if (length(significance_basis) != 1L ||
    !significance_basis %in% c("BH", "raw")) {
  stop("significance_basis must be either 'BH' or 'raw'.")
}
if (length(or_color_max) != 1L || is.na(or_color_max) ||
    !is.finite(or_color_max) || or_color_max <= 0) {
  stop("or_color_max must be a single positive finite number.")
}

# IUPAC consensus patterns. An element is present if at least one match is found
# on either DNA strand in the promoter of length promoter_length.
# Candidate elements are listed first. TATA/CAAT/GC motifs are shown as reference
# promoter motifs rather than true negative controls: these sequences can have
# biological promoter functions and should not be assumed to be inert.
motifs <- tibble::tribble(
  ~Element,                 ~Pattern,       ~ElementClass,
  "W-box\n(WRKY)",         "TTGACY",       "Candidate element",
  "as-1\n(TGA)",           "TGACG",        "Candidate element",
  "HSE",                    "NGAANNTTCN",   "Candidate element",
  "GCC-box",                "AGCCGCC",      "Candidate element",
  "ABRE",                   "ACGTGKC",      "Candidate element",
  "DRE/CRT",                "RCCGAC",       "Candidate element",
  "TCA-element",            "CCATCTTTTT",   "Candidate element",
  "TATA-box",               "TATAWAWR",     "Reference",
  "CAAT-box",               "CCAAT",        "Reference",
  "GC-box",                 "GGGCGG",       "Reference"
) %>%
  mutate(
    ElementClass = factor(ElementClass,
                          levels = c("Candidate element", "Reference")),
    Element = factor(Element, levels = Element)
  )

# ----------------------------- 1. Read subclusters ------------------
if (!file.exists(subcluster_file)) {
  stop("Cannot find ", subcluster_file,
       ". Run the Cluster 5 subclustering script first in this working directory.")
}
if (!file.exists(annotation_file)) {
  stop("Cannot find annotation file: ", annotation_file,
       ". Set annotation_file to the correct path.")
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(trimws(x))
}

subclusters <- read.csv(
  subcluster_file, stringsAsFactors = FALSE, check.names = FALSE
)
required_cols <- c("gene", "ParentCluster", "Subcluster")
if (!all(required_cols %in% names(subclusters))) {
  stop("The subcluster file must contain columns: ", paste(required_cols, collapse = ", "))
}

subclusters <- subclusters %>%
  mutate(
    gene = clean_gene_id(gene),
    ParentCluster = suppressWarnings(as.integer(as.character(ParentCluster))),
    Subcluster = suppressWarnings(as.integer(as.character(Subcluster)))
  )

if (anyNA(subclusters$gene) || any(!nzchar(subclusters$gene))) {
  stop("The subcluster assignment file contains missing or empty gene IDs.")
}
if (anyNA(subclusters$ParentCluster) || anyNA(subclusters$Subcluster)) {
  stop("ParentCluster or Subcluster contains missing/non-integer values.")
}
if (anyDuplicated(subclusters$gene)) {
  dup <- unique(subclusters$gene[duplicated(subclusters$gene)])
  stop("Duplicate gene IDs after cleaning in the subcluster assignment file (examples: ",
       paste(head(dup, 10), collapse = ", "), ").")
}

subclusters <- subclusters %>%
  filter(ParentCluster == parent_cluster,
         Subcluster %in% seq_len(k_subclusters))

if (nrow(subclusters) == 0) {
  stop("No genes found for Cluster ", parent_cluster,
       ". Check the subcluster file and parent_cluster setting.")
}

message("Cluster 5 subcluster gene counts (before promoter matching):")
print(subclusters %>% count(Subcluster, name = "n_genes"))

# ----------------------------- 2. Get gene TSSs ---------------------
gff <- rtracklayer::import(annotation_file)
gff_df <- as.data.frame(gff)

if (!all(c("type", "ID", "seqnames", "start", "end", "strand") %in% names(gff_df))) {
  stop("The GFF3 file does not contain the expected gene annotation columns.")
}

gene_df <- gff_df %>%
  filter(type == "gene") %>%
  transmute(
    gene = clean_gene_id(ID),
    chr = as.character(seqnames),
    start = as.integer(start),
    end = as.integer(end),
    strand = as.character(strand)
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(start), !is.na(end),
         start >= 1L, end >= start, strand %in% c("+", "-"))

if (anyDuplicated(gene_df$gene)) {
  dup <- unique(gene_df$gene[duplicated(gene_df$gene)])
  stop("Duplicate gene IDs after cleaning in the GFF3 gene features (examples: ",
       paste(head(dup, 10), collapse = ", "),
       "). Resolve annotation ID duplication rather than silently keeping one record.")
}

# Arabidopsis reference genome bundled with the BSgenome package (TAIR10).
genome <- BSgenome.Athaliana.TAIR.TAIR9
available_chromosomes <- names(GenomeInfoDb::seqlengths(genome))
n_gene_before_chr_filter <- nrow(gene_df)
gene_df <- gene_df %>% filter(chr %in% available_chromosomes)
message("Gene annotations excluded because chromosome names are absent from TAIR10: ",
        n_gene_before_chr_filter - nrow(gene_df))

if (nrow(gene_df) == 0) {
  stop("No gene annotations matched the chromosome names in the TAIR10 genome. ",
       "Check that the GFF3 annotation and genome build are compatible.")
}

# Build promoter coordinates numerically BEFORE constructing IRanges.
# GenomicRanges::promoters() can create a negative start for genes close to the
# beginning of a chromosome; IRanges rejects such ranges before they can be clipped.
# For a plus-strand gene, the promoter is [TSS - promoter_length, TSS - 1].
# For a minus-strand gene, it is [TSS + 1, TSS + promoter_length].
# The GFF3 gene boundary is used as a gene-level TSS proxy; alternative
# transcript-specific TSSs are not modeled separately.
chr_length <- as.integer(GenomeInfoDb::seqlengths(genome)[gene_df$chr])
tss <- ifelse(gene_df$strand == "+", gene_df$start, gene_df$end)
promoter_start <- ifelse(
  gene_df$strand == "+",
  tss - as.integer(promoter_length),
  tss + 1L
)
promoter_end <- ifelse(
  gene_df$strand == "+",
  tss - 1L,
  tss + as.integer(promoter_length)
)

# Keep only complete, in-bounds promoters. Do not truncate them: each sequence
# tested below must have exactly promoter_length bases for comparable detection.
keep_full_promoters <- !is.na(chr_length) &
  promoter_start >= 1L &
  promoter_end <= chr_length &
  (promoter_end - promoter_start + 1L) == promoter_length

n_promoters_excluded_boundary <- sum(!keep_full_promoters)
message("Genes excluded because a complete ", promoter_length,
        "-bp promoter extends beyond a chromosome boundary: ",
        n_promoters_excluded_boundary)

gene_promoter_df <- gene_df[keep_full_promoters, , drop = FALSE]
promoter_start <- promoter_start[keep_full_promoters]
promoter_end <- promoter_end[keep_full_promoters]

if (nrow(gene_promoter_df) == 0L) {
  stop("No full-length promoters could be extracted. Check annotation/genome chromosome names.")
}

promoter_gr <- GenomicRanges::GRanges(
  seqnames = gene_promoter_df$chr,
  ranges = IRanges::IRanges(start = promoter_start, end = promoter_end),
  strand = gene_promoter_df$strand
)
S4Vectors::mcols(promoter_gr)$gene <- gene_promoter_df$gene
promoter_ids <- gene_promoter_df$gene
names(promoter_gr) <- promoter_ids
promoter_seqs <- Biostrings::getSeq(genome, promoter_gr)
names(promoter_seqs) <- promoter_ids

message("Annotated genes with complete ", promoter_length,
        "-bp promoters: ", length(promoter_seqs))
message("Promoter matching denominator is restricted to genes with a complete promoter.")

# ----------------------------- 3. Scan promoter motifs --------------
# fixed = FALSE enables IUPAC degeneracy codes in the consensus pattern.
# Both orientations are checked because cis-elements may occur on either DNA strand.
motif_presence <- tibble(gene = names(promoter_seqs))

for (i in seq_len(nrow(motifs))) {
  element_name <- as.character(motifs$Element[i])
  pattern_string <- as.character(motifs$Pattern[i])
  pattern <- Biostrings::DNAString(pattern_string)
  pattern_rc <- Biostrings::reverseComplement(pattern)

  forward_hits <- Biostrings::vcountPattern(
    pattern, promoter_seqs, fixed = FALSE
  ) > 0L
  reverse_hits <- Biostrings::vcountPattern(
    pattern_rc, promoter_seqs, fixed = FALSE
  ) > 0L

  motif_presence[[element_name]] <- forward_hits | reverse_hits
}

# Save a gene-by-motif binary matrix so each promoter call is auditable.
write.csv(
  motif_presence,
  paste0(out_prefix, "_promoter_motif_presence.csv"),
  row.names = FALSE, quote = TRUE
)

# ----------------------------- 4. Fisher exact tests ----------------
# Foreground = genes in a given subcluster that have a valid 1-kb promoter.
# Background = all other annotated genes with a complete promoter in the same genome build.
# The Fisher test is one-sided (alternative = "greater") to test enrichment.
results_list <- list()
presence_genes <- motif_presence$gene

for (sc in seq_len(k_subclusters)) {
  foreground <- intersect(
    subclusters %>% filter(Subcluster == sc) %>% pull(gene) %>% unique(),
    presence_genes
  )
  background <- setdiff(presence_genes, foreground)

  message("Subcluster ", sc, ": ", length(foreground),
          " foreground genes with complete promoters; ", length(background),
          " background genes.")

  if (length(foreground) == 0 || length(background) == 0) {
    warning("Subcluster ", sc, " has no foreground or background genes after promoter matching.")
    next
  }

  fg_rows <- match(foreground, motif_presence$gene)
  bg_rows <- match(background, motif_presence$gene)

  for (i in seq_len(nrow(motifs))) {
    element_name <- as.character(motifs$Element[i])
    present_vec <- motif_presence[[element_name]]

    a <- sum(present_vec[fg_rows], na.rm = TRUE)  # foreground with motif
    b <- length(fg_rows) - a                     # foreground without motif
    c <- sum(present_vec[bg_rows], na.rm = TRUE)  # background with motif
    d <- length(bg_rows) - c                     # background without motif

    contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                          dimnames = list(c("Subcluster", "Genome_background"),
                                          c("Motif_present", "Motif_absent")))
    p_value <- tryCatch(
      stats::fisher.test(contingency, alternative = "greater")$p.value,
      error = function(e) NA_real_
    )

    raw_or <- if ((b * c) == 0) {
      if ((a * d) > 0) Inf else NA_real_
    } else {
      (a * d) / (b * c)
    }

    # Haldane-Anscombe 0.5 correction is used only for the displayed OR if
    # any cell is zero, avoiding infinite/undefined color values.
    or_display <- if (any(c(a, b, c, d) == 0)) {
      ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    } else {
      raw_or
    }

    results_list[[length(results_list) + 1L]] <- tibble(
      ParentCluster = parent_cluster,
      Subcluster = sc,
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

if (length(results_list) == 0) {
  stop("No enrichment tests could be calculated. Check Cluster 5 subcluster assignments.")
}

enrichment <- bind_rows(results_list) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    p_for_stars = if (significance_basis == "BH") p_adj_BH else p_value,
    significance = case_when(
      is.na(p_for_stars) ~ "NA",
      p_for_stars < 0.001 ~ "***",
      p_for_stars < 0.01 ~ "**",
      p_for_stars < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    SubclusterLabel = factor(
      paste0("Sc", Subcluster),
      levels = paste0("Sc", seq_len(k_subclusters))
    ),
    Element = factor(Element, levels = as.character(motifs$Element)),
    ElementClass = factor(ElementClass,
                          levels = c("Candidate element", "Reference")),
    tile_label = if_else(
      is.finite(odds_ratio_display),
      paste0(sprintf("%.2f", odds_ratio_display), "\n", significance),
      paste0("NA\n", significance)
    )
  )

write.csv(
  enrichment,
  paste0(out_prefix, "_results.csv"),
  row.names = FALSE, quote = TRUE
)

message("Significance stars are based on ", significance_basis,
        " p-values; BH-adjusted p-values are stored in the results CSV.")
message("Interpret motif hits as consensus-sequence occurrences, not evidence of TF binding.")

message("Genes included per subcluster after promoter matching:")
print(enrichment %>% distinct(Subcluster, n_foreground) %>% arrange(Subcluster))

# ----------------------------- 5. Plot -------------------------------
# Subclusters are displayed as columns; motif elements are displayed as rows.
# Only plot layout, text size/weight, and export dimensions are changed here.
# Wrap the caption so it stays within the exported figure boundaries.
plot_caption <- paste0(
  promoter_length, "-bp promoters; Fisher's exact test (enrichment). Stars use ",
  if (significance_basis == "BH") "BH-adjusted p-values" else "raw p-values",
  "; BH-adjusted p-values are also reported in the results table."
)
plot_caption <- paste(strwrap(plot_caption, width = 75), collapse = "\n")

p_cis <- ggplot(
  enrichment,
  aes(x = SubclusterLabel, y = Element, fill = odds_ratio_display)
) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(
    aes(label = tile_label),
    size = 4.2, lineheight = 0.95, fontface = "bold"
  ) +
  facet_grid(ElementClass ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_fill_gradient2(
    low = "#5878B7",
    mid = "#F7F7F7",
    high = "#E31A1C",
    midpoint = 1,
    limits = c(0, or_color_max),
    oob = scales::squish,
    breaks = c(0, 0.5, 1, 1.5, 2, 3),
    name = "Odds\nratio"
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  # Drop motif levels not used in each facet; keeping all levels duplicated labels
  # into the other facet and created large blank regions in the previous figure.
  scale_y_discrete(drop = TRUE) +
  labs(
    title = paste0("Cis-regulatory element enrichment: Cluster ",
                   parent_cluster, " subclusters"),
    x = NULL,
    y = NULL,
    caption = plot_caption
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey85", colour = "grey40", linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 15),
    axis.text.x = element_text(colour = "grey20", size = 20, face = "bold"),
    axis.text.x.top = element_text(colour = "grey20", size = 14, face = "bold"),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.y = element_text(colour = "grey15", size = 13, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing.y = grid::unit(0.25, "lines"),
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.caption = element_text(size = 15, face = "bold", hjust = 0, lineheight = 1.05),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_cis)

ggsave(
  paste0(out_prefix, ".png"), p_cis,
  width = 6, height = 6.8, units = "in", dpi = 300, bg = "white"
)
ggsave(
  paste0(out_prefix, ".pdf"), p_cis,
  width = 11, height = 6.8, units = "in", bg = "white"
)

# Save the motif definitions used in this analysis for reproducibility.
write.csv(
  motifs %>% mutate(Element = as.character(Element), ElementClass = as.character(ElementClass)),
  paste0(out_prefix, "_motif_definitions.csv"),
  row.names = FALSE, quote = TRUE
)

message("Done. Output prefix: ", out_prefix)
