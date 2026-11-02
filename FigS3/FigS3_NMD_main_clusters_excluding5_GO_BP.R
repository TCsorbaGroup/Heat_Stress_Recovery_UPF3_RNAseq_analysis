# =============================================================================
# GO Biological Process enrichment for the original k-means main clusters,
# excluding Cluster 5 (analysed separately).
#
# Analysis:
#   1. Read the original top-2000 k=6 cluster assignments.
#   2. Use all genes in the original top-2000 set as the background universe.
#   3. Run GO Biological Process over-representation analysis separately
#      for Clusters 1, 2, 3, 4 and 6.
#   4. Apply Benjamini-Hochberg correction within each cluster's GO analysis.
#   5. Reduce redundant significant GO terms using rrvgo.
#   6. Generate one faceted dotplot.
# =============================================================================


# ----------------------------- 0. Packages -----------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(rrvgo)
  library(dplyr)
  library(tibble)
  library(stringr)
})


# ----------------------------- 1. Settings -----------------------------

clusters_to_analyze <- c(1, 2, 3, 4, 6)

fdr_cutoff <- 0.05
rrvgo_threshold <- 0.7

# Maximum number of non-redundant significant GO terms shown per cluster.
# Set to Inf to plot all available terms.
max_terms_per_cluster <- 20

background_mode <- "top2000"

original_clusters_file <-
  "NMD_expression_trajectory_top2000_by_variance_k6_gene_clusters.csv"

out_prefix <- "NMD_main_clusters_excluding5_GO_BP"


# ----------------------------- 2. Helper functions --------------------
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}


# Return an empty tibble when a cluster has no results.
empty_enrichment_table <- function() {
  tibble::tibble()
}


# Format significant terms without semantic reduction.
# This is also used as a fallback if semantic-similarity reduction fails.
format_unreduced_terms <- function(terms) {
  if (nrow(terms) == 0L) {
    return(tibble::tibble())
  }

  terms %>%
    transmute(
      Cluster,
      GO_ID = ID,
      Description,
      p.adjust,
      Count,
      GeneRatio,
      n_terms_collapsed = 1L,
      log10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
    )
}


# ----------------------------- 3. Read and validate inputs ------------

if (!file.exists(original_clusters_file)) {
  stop(
    "Cannot find input file: ", original_clusters_file,
    "\nRun the original top-2000 k=6 clustering script first, ",
    "and use the same working directory."
  )
}


all_cluster_assignments <- read.csv(
  original_clusters_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (!all(c("gene", "Cluster") %in% names(all_cluster_assignments))) {
  stop(
    "The input file must contain the columns 'gene' and 'Cluster'."
  )
}


all_cluster_assignments <- all_cluster_assignments %>%
  mutate(
    gene = clean_gene_id(gene),
    Cluster = suppressWarnings(
      as.integer(as.character(Cluster))
    )
  )


if (anyNA(all_cluster_assignments$gene) ||
    any(!nzchar(all_cluster_assignments$gene))) {
  stop(
    "The original cluster file contains missing or empty gene IDs."
  )
}


if (anyNA(all_cluster_assignments$Cluster)) {
  stop(
    "Some original Cluster values are missing or cannot be converted ",
    "to integers."
  )
}


# Stop if canonicalisation produces duplicated gene identifiers.
# This catches both repeated rows and assignments that could otherwise
# silently merge different original identifiers.
if (anyDuplicated(all_cluster_assignments$gene)) {
  duplicated_genes <- unique(
    all_cluster_assignments$gene[
      duplicated(all_cluster_assignments$gene)
    ]
  )

  stop(
    "Duplicated gene IDs were found after identifier cleanup ",
    "in the original cluster file. Examples: ",
    paste(head(duplicated_genes, 10), collapse = ", "),
    ". Check the input before proceeding."
  )
}


missing_clusters <- setdiff(
  clusters_to_analyze,
  unique(all_cluster_assignments$Cluster)
)

if (length(missing_clusters) > 0L) {
  stop(
    "The input file does not contain expected cluster(s): ",
    paste(missing_clusters, collapse = ", ")
  )
}

if (!all(unique(all_cluster_assignments$Cluster) %in% 1:6)) {
  warning(
    "The input contains cluster labels outside 1:6. ",
    "Only the specified clusters will be tested."
  )
}


# ----------------------------- 4. Define the background ---------------

all_original_genes <- unique(all_cluster_assignments$gene)


universe_raw <- switch(
  background_mode,

  top2000 = all_original_genes,

  stop(
    "Unsupported background_mode: ", background_mode,
    ". This script is intended to use 'top2000'."
  )
)


# Restrict the analysis to recognised TAIR identifiers.
valid_tair_ids <- unique(
  toupper(
    AnnotationDbi::keys(
      org.At.tair.db,
      keytype = "TAIR"
    )
  )
)


universe <- intersect(
  unique(universe_raw),
  valid_tair_ids
)


if (length(universe) == 0L) {
  stop(
    "No valid TAIR identifiers remain in the selected background. ",
    "Check the input gene IDs and the organism annotation."
  )
}


# Restrict the tested sets to Clusters 1, 2, 3, 4 and 6.
cluster_assignments <- all_cluster_assignments %>%
  filter(
    Cluster %in% clusters_to_analyze,
    gene %in% universe
  )


if (nrow(cluster_assignments) == 0L) {
  stop(
    "No genes from the selected clusters remain after TAIR ID filtering."
  )
}


# Record the number of valid input genes for each selected cluster,
# including clusters with zero valid genes.
cluster_gene_counts <- cluster_assignments %>%
  count(Cluster, name = "n_valid_TAIR_genes") %>%
  right_join(
    tibble(Cluster = clusters_to_analyze),
    by = "Cluster"
  ) %>%
  mutate(
    n_valid_TAIR_genes = ifelse(
      is.na(n_valid_TAIR_genes),
      0L,
      n_valid_TAIR_genes
    )
  ) %>%
  arrange(match(Cluster, clusters_to_analyze))


counts_file <- paste0(
  out_prefix,
  "_cluster_gene_counts.csv"
)


write.csv(
  cluster_gene_counts,
  counts_file,
  row.names = FALSE
)


message("Clusters analysed: ",
        paste(clusters_to_analyze, collapse = ", "))

message("Background mode: ", background_mode)

message(
  "Genes in the background after TAIR ID filtering: ",
  length(universe)
)

message("Valid input genes by cluster:")
print(cluster_gene_counts)


# ----------------------------- 5. GO BP enrichment --------------------

run_enrichment <- function(cluster_id) {

  genes_cluster <- cluster_assignments %>%
    filter(Cluster == cluster_id) %>%
    pull(gene) %>%
    unique()


  if (length(genes_cluster) == 0L) {
    message(
      "Cluster ", cluster_id,
      ": no valid genes; skipping enrichment."
    )

    return(empty_enrichment_table())
  }


  message(
    "Running GO BP enrichment for Cluster ",
    cluster_id,
    " (", length(genes_cluster), " input genes)."
  )


  # BH correction is applied separately within each cluster's GO analysis.
  #
  # pvalueCutoff = 1 and qvalueCutoff = 1 preserve the returned result set
  # for subsequent explicit filtering at fdr_cutoff.
 
  ego <- clusterProfiler::enrichGO(
    gene = genes_cluster,
    universe = universe,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = FALSE
  )


  res <- as.data.frame(ego)


  if (nrow(res) == 0L) {
    message(
      "Cluster ", cluster_id,
      ": no GO BP terms were returned."
    )

    return(empty_enrichment_table())
  }


  res %>%
    mutate(
      Cluster = cluster_id,
      n_input_genes = length(genes_cluster),
      .before = 1
    )
}


all_go_results <- bind_rows(
  lapply(clusters_to_analyze, run_enrichment)
)


all_results_file <- paste0(
  out_prefix,
  "_all_enrichment.csv"
)


# Save the complete returned enrichment results, including non-significant
# terms, so the downstream filtering can be audited.
write.csv(
  all_go_results,
  all_results_file,
  row.names = FALSE
)


if (nrow(all_go_results) == 0L) {
  stop(
    "No GO BP terms were returned for any selected cluster. ",
    "The empty/full-result table was written to: ",
    all_results_file,
    ". Check gene IDs, annotation and background."
  )
}


# Explicit significance filter applied to the BH-adjusted p-values.
sig_go_results <- all_go_results %>%
  filter(
    !is.na(p.adjust),
    p.adjust < fdr_cutoff
  )


if (nrow(sig_go_results) == 0L) {
  stop(
    "No GO Biological Process terms passed BH-adjusted p < ",
    fdr_cutoff,
    ". The complete enrichment table was saved, but no reduced-term ",
    "table or dotplot was created."
  )
}


message(
  "Significant GO BP terms before semantic reduction: ",
  nrow(sig_go_results)
)


# ----------------------------- 6. Reduce redundant GO terms -----------

reduce_one_cluster <- function(cluster_id) {

  cluster_terms <- sig_go_results %>%
    filter(Cluster == cluster_id) %>%
    distinct(ID, .keep_all = TRUE)


  if (nrow(cluster_terms) == 0L) {
    message(
      "Cluster ", cluster_id,
      ": no significant GO BP terms to reduce."
    )

    return(tibble::tibble())
  }


  # A single term does not require semantic reduction.
  if (nrow(cluster_terms) == 1L) {
    return(format_unreduced_terms(cluster_terms))
  }


  # IMPORTANT:
  # keytype = "TAIR" matches the gene identifier type used by enrichGO.
  # The x argument itself contains GO term IDs, not gene IDs.
  sim_matrix <- tryCatch(
    rrvgo::calculateSimMatrix(
      unique(cluster_terms$ID),
      orgdb = "org.At.tair.db",
      keytype = "TAIR",
      ont = "BP",
      method = "Wang"
    ),

    error = function(e) {
      message(
        "Semantic-similarity calculation failed for Cluster ",
        cluster_id,
        ". Keeping significant terms without reduction. Error: ",
        conditionMessage(e)
      )

      NULL
    }
  )


  # If the similarity calculation failed, retain all significant terms.
  if (is.null(sim_matrix) || nrow(sim_matrix) == 0L) {
    return(format_unreduced_terms(cluster_terms))
  }


  # Separate terms that could and could not be represented in the
  # similarity matrix. Terms missing from the matrix must not disappear
  # silently from the final significant-term table.
  terms_with_similarity <- cluster_terms %>%
    filter(ID %in% rownames(sim_matrix))

  terms_without_similarity <- cluster_terms %>%
    filter(!ID %in% rownames(sim_matrix))


  if (nrow(terms_with_similarity) == 0L) {
    return(format_unreduced_terms(cluster_terms))
  }


  if (nrow(terms_with_similarity) == 1L) {
    return(
      bind_rows(
        format_unreduced_terms(terms_with_similarity),
        format_unreduced_terms(terms_without_similarity)
      )
    )
  }


  sim_matrix <- sim_matrix[
    terms_with_similarity$ID,
    terms_with_similarity$ID,
    drop = FALSE
  ]


  # Higher scores are better; use -log10(BH-adjusted p-values).
  scores <- setNames(
    -log10(
      pmax(
        terms_with_similarity$p.adjust,
        .Machine$double.xmin
      )
    ),
    terms_with_similarity$ID
  )


  reduced <- tryCatch(
    rrvgo::reduceSimMatrix(
      sim_matrix,
      scores = scores,
      threshold = rrvgo_threshold,
      orgdb = "org.At.tair.db"
    ) %>%
      as.data.frame(),

    error = function(e) {
      message(
        "Semantic-similarity reduction failed for Cluster ",
        cluster_id,
        ". Keeping significant terms without reduction. Error: ",
        conditionMessage(e)
      )

      NULL
    }
  )


  # Fallback if rrvgo fails or returns an unexpected structure.
  if (
    is.null(reduced) ||
    nrow(reduced) == 0L ||
    !all(c("go", "parent") %in% names(reduced))
  ) {
    return(format_unreduced_terms(cluster_terms))
  }


  # Count the number of GO terms represented by each cluster of similar
  # terms. The representative is the term selected by rrvgo.
  collapsed_counts <- reduced %>%
    group_by(parent) %>%
    summarise(
      n_terms_collapsed = n_distinct(go),
      .groups = "drop"
    )


  representatives <- reduced %>%
    filter(go == parent) %>%
    distinct(parent, .keep_all = TRUE) %>%
    dplyr::select(parent)


  representative_terms <- representatives %>%
    left_join(
      collapsed_counts,
      by = "parent"
    ) %>%
    left_join(
      terms_with_similarity,
      by = c("parent" = "ID")
    ) %>%
    transmute(
      Cluster = cluster_id,
      GO_ID = parent,
      Description,
      p.adjust,
      Count,
      GeneRatio,
      n_terms_collapsed = as.integer(n_terms_collapsed),
      log10FDR = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )


  # Preserve significant terms not present in the similarity matrix.
  bind_rows(
    representative_terms,
    format_unreduced_terms(terms_without_similarity)
  )
}


reduced_go <- bind_rows(
  lapply(clusters_to_analyze, reduce_one_cluster)
)


reduced_file <- paste0(
  out_prefix,
  "_reduced_terms.csv"
)


if (nrow(reduced_go) == 0L) {
  write.csv(
    reduced_go,
    reduced_file,
    row.names = FALSE
  )

  stop(
    "No representative GO terms remained after semantic-similarity ",
    "reduction. An empty reduced-term table was saved to: ",
    reduced_file
  )
}


write.csv(
  reduced_go,
  reduced_file,
  row.names = FALSE
)


message(
  "Representative GO BP terms after reduction: ",
  nrow(reduced_go)
)


# ----------------------------- 7. Prepare dotplot ----------------------

plot_data <- reduced_go %>%
  group_by(Cluster) %>%
  arrange(desc(log10FDR), .by_group = TRUE) %>%
  {
    if (is.finite(max_terms_per_cluster)) {
      slice_head(., n = max_terms_per_cluster)
    } else {
      .
    }
  } %>%
  ungroup() %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = clusters_to_analyze,
      labels = paste0("Sc", clusters_to_analyze)
    ),

    # Include GO ID and cluster in the internal label to ensure that
    # identical descriptions cannot create duplicate factor levels.
    TermKey = paste0(
      Description,
      " [GO:",
      GO_ID,
      ";Cluster:",
      as.character(Cluster),
      "]"
    )
  )


if (
  nrow(plot_data) == 0L ||
  all(!is.finite(plot_data$log10FDR))
) {
  stop(
    "No finite GO BP enrichment values are available for the dotplot."
  )
}


# Reverse order so stronger enrichment values appear toward the top
# of each panel after faceting.
plot_data$TermKey <- factor(
  plot_data$TermKey,
  levels = rev(unique(plot_data$TermKey))
)


p_go <- ggplot(
  plot_data,
  aes(
    x = log10FDR,
    y = TermKey
  )
) +
  geom_point(
    aes(
      size = Count,
      colour = log10FDR
    ),
    alpha = 0.95
  ) +
  facet_grid(
    rows = vars(Cluster),
    scales = "free_y",
    space = "free_y",
    drop = TRUE
  ) +
  scale_colour_gradient(
    low = "#3B0066",
    high = "#E60026",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(
    range = c(2.5, 8),
    breaks = pretty(plot_data$Count, n = 4),
    name = "Gene count"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_y_discrete(
    labels = function(x) {
      stringr::str_wrap(
        sub(" \\[GO:.*$", "", x),
        width = 40
      )
    }
  ) +
  labs(
    title = "GO BP enrichment: main clusters",
    subtitle = paste0(
      "GO Biological Process; BH-adjusted p < ",
      fdr_cutoff,
      "; background = original top-2000 gene set"
    ),
    x = expression(-log[10](FDR) ~ "of representative term"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 14,
      hjust = 0,
      face = "bold"
    ),
    plot.subtitle = element_blank()
    ,
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey45",
      linewidth = 0.5
    ),
    strip.text.y = element_text(
      angle = 270,
      size = 11,
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey88",
      linewidth = 0.3
    ),
    axis.text.y = element_text(
      size = 8.5,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 11
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.spacing.y = grid::unit(0.18, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )


print(p_go)


# ----------------------------- 8. Save outputs ------------------------

png_file <- paste0(
  out_prefix,
  "_dotplot.png"
)

pdf_file <- paste0(
  out_prefix,
  "_dotplot.pdf"
)

n_terms_plot <- nrow(plot_data)

plot_height <- max(
  8.0,
  1.8 + 0.24 * n_terms_plot
)


ggsave(
  filename = png_file,
  plot = p_go,
  width = 5.4,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)


ggsave(
  filename = pdf_file,
  plot = p_go,
  width = 8.2,
  height = plot_height,
  units = "in",
  bg = "white"
)


# ----------------------------- 9. Report ------------------------------

message("")
message("GO BP analysis completed.")
message("Background genes after TAIR ID filtering: ", length(universe))
message("Clusters analysed: ",
        paste(clusters_to_analyze, collapse = ", "))
message("BH-adjusted significance threshold: ", fdr_cutoff)
message("rrvgo semantic-similarity threshold: ", rrvgo_threshold)
message("")
message("Saved cluster gene counts: ", counts_file)
message("Saved full enrichment table: ", all_results_file)
message("Saved non-redundant GO terms: ", reduced_file)
message("Saved dotplot: ", png_file)
message("Saved dotplot: ", pdf_file)
message("")
message(
  "Interpretation note: the top2000 background asks which GO BP terms ",
  "are over-represented in each main cluster relative to the genes ",
  "selected for the original top-2000 clustering. It is not a test ",
  "against every gene in the Arabidopsis genome."
)
message(
  "Adjustment note: BH correction is performed separately for each ",
  "cluster's GO BP enrichment analysis."
)

