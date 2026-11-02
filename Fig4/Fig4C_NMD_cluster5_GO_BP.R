# GO Biological Process enrichment for the 3 subclusters of original Cluster 5
# Hypergeometric over-representation test with BH-adjusted p-values,
# followed by semantic-similarity reduction with rrvgo (threshold = 0.7).

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(rrvgo)
  library(dplyr)
  library(tidyr)
})

# ----------------------------- Settings -----------------------------
parent_cluster <- 5
k_subclusters <- 3
fdr_cutoff <- 0.05
rrvgo_threshold <- 0.7
max_terms_per_subcluster <- 20  # use Inf to show every non-redundant significant term

# Background interpretation:
# "top2000" compares each subcluster with all genes in the original k=6 analysis.
# "parent_cluster" compares each subcluster with all genes in original Cluster 5,
# which is more directly focused on terms over-represented within a subcluster
# relative to its parent cluster. The original setting is retained below.
background_mode <- "top2000"

subcluster_file <- paste0(
  "NMD_cluster", parent_cluster, "_subclusters_k", k_subclusters,
  "_gene_assignments.csv"
)
original_clusters_file <- paste0(
  "NMD_expression_trajectory_top2000_by_variance_k6_gene_clusters.csv"
)
out_prefix <- paste0(
  "NMD_cluster", parent_cluster, "_subclusters_k", k_subclusters, "_GO_BP"
)

# ----------------------------- 1. Read and validate inputs ----------
if (!file.exists(subcluster_file)) {
  stop("Cannot find ", subcluster_file,
       ". Run the Cluster 5 subclustering script first, and use the same working directory.")
}
if (!file.exists(original_clusters_file)) {
  stop("Cannot find ", original_clusters_file,
       ". Run the original top-2000 k=6 clustering script first.")
}

# Canonicalise identifiers consistently before matching them.
clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(gene|transcript):", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

subclusters_raw <- read.csv(
  subcluster_file, stringsAsFactors = FALSE, check.names = FALSE
)
all_cluster_assignments <- read.csv(
  original_clusters_file, stringsAsFactors = FALSE, check.names = FALSE
)

if (!all(c("gene", "ParentCluster", "Subcluster") %in% names(subclusters_raw))) {
  stop("The subcluster file must contain columns: gene, ParentCluster, Subcluster.")
}
if (!all(c("gene", "Cluster") %in% names(all_cluster_assignments))) {
  stop("The original cluster file must contain columns: gene, Cluster.")
}

subclusters <- subclusters_raw %>%
  mutate(
    gene = clean_gene_id(gene),
    ParentCluster = suppressWarnings(as.integer(as.character(ParentCluster))),
    Subcluster = suppressWarnings(as.integer(as.character(Subcluster)))
  )
all_cluster_assignments <- all_cluster_assignments %>%
  mutate(
    gene = clean_gene_id(gene),
    Cluster = suppressWarnings(as.integer(as.character(Cluster)))
  )

if (anyNA(subclusters$gene) || any(!nzchar(subclusters$gene))) {
  stop("The subcluster file contains missing or empty gene IDs.")
}
if (anyNA(all_cluster_assignments$gene) || any(!nzchar(all_cluster_assignments$gene))) {
  stop("The original cluster file contains missing or empty gene IDs.")
}
if (anyNA(subclusters$ParentCluster) || anyNA(subclusters$Subcluster)) {
  stop("Some ParentCluster or Subcluster values are missing or not integers.")
}
if (anyNA(all_cluster_assignments$Cluster)) {
  stop("Some original Cluster values are missing or not integers.")
}

if (anyDuplicated(subclusters$gene)) {
  examples <- unique(subclusters$gene[duplicated(subclusters$gene)])
  stop("Duplicated gene IDs in the subcluster file after ID cleanup (examples: ",
       paste(head(examples, 10), collapse = ", "), ").")
}
if (anyDuplicated(all_cluster_assignments$gene)) {
  examples <- unique(all_cluster_assignments$gene[duplicated(all_cluster_assignments$gene)])
  stop("Duplicated gene IDs in the original cluster file after ID cleanup (examples: ",
       paste(head(examples, 10), collapse = ", "), ").")
}

if (!all(subclusters$ParentCluster == parent_cluster)) {
  stop("The subcluster input contains rows not assigned to parent cluster ", parent_cluster, ".")
}
if (!all(subclusters$Subcluster %in% seq_len(k_subclusters))) {
  stop("The subcluster input contains Subcluster values outside 1..", k_subclusters, ".")
}

parent_genes <- all_cluster_assignments %>%
  filter(Cluster == parent_cluster) %>%
  pull(gene) %>%
  unique()

top2000_genes <- unique(all_cluster_assignments$gene)
if (length(parent_genes) == 0L) {
  stop("No genes found in original Cluster ", parent_cluster, ".")
}

not_in_parent <- setdiff(subclusters$gene, parent_genes)
if (length(not_in_parent) > 0L) {
  stop("Some subcluster genes are not assigned to original Cluster ", parent_cluster,
       " (examples: ", paste(head(not_in_parent, 10), collapse = ", "), ").")
}

universe_raw <- switch(
  background_mode,
  top2000 = top2000_genes,
  parent_cluster = parent_genes,
  stop("background_mode must be either 'top2000' or 'parent_cluster'.")
)

# Restrict both the background and tested gene lists to valid TAIR keys.
valid_tair_ids <- AnnotationDbi::keys(org.At.tair.db, keytype = "TAIR")
universe <- intersect(unique(universe_raw), valid_tair_ids)
subclusters <- subclusters %>% filter(gene %in% universe)

if (length(universe) == 0L) {
  stop("The selected background contains no valid TAIR IDs. Check gene-ID formatting and annotation keys.")
}
if (nrow(subclusters) == 0L) {
  stop("No subcluster genes remain after filtering to valid TAIR IDs.")
}

message("Parent cluster: ", parent_cluster)
message("Background mode: ", background_mode)
message("Genes in background after TAIR ID filtering: ", length(universe))
message("Valid input genes by subcluster:")
print(subclusters %>% count(Subcluster, name = "n_genes"))

# ----------------------------- 2. GO enrichment ---------------------
run_enrichment <- function(sc) {
  genes_sc <- subclusters %>%
    filter(Subcluster == sc) %>%
    pull(gene) %>%
    unique()

  if (length(genes_sc) == 0L) {
    message("Subcluster ", sc, ": no valid genes; skipping.")
    return(tibble())
  }

  message("Running GO BP enrichment for Subcluster ", sc,
          " (", length(genes_sc), " input genes).")

  # BH adjustment is performed by enrichGO separately for each subcluster.
  ego <- clusterProfiler::enrichGO(
    gene = genes_sc,
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
  if (nrow(res) == 0L) return(tibble())

  res %>% mutate(
    Subcluster = sc,
    n_input_genes = length(genes_sc),
    .before = 1
  )
}

all_go_results <- bind_rows(lapply(seq_len(k_subclusters), run_enrichment))
all_results_file <- paste0(out_prefix, "_all_enrichment.csv")

if (nrow(all_go_results) == 0L) {
  write.csv(all_go_results, all_results_file, row.names = FALSE)
  stop("No GO BP terms were returned. Check gene IDs and organism annotation.")
}

write.csv(all_go_results, all_results_file, row.names = FALSE)

sig_go_results <- all_go_results %>%
  filter(!is.na(p.adjust), p.adjust < fdr_cutoff)

if (nrow(sig_go_results) == 0L) {
  stop("No GO Biological Process terms passed BH-adjusted p < ", fdr_cutoff,
       ". The complete enrichment table was saved, but no reduced-term table or dotplot was created.")
}

# ----------------------------- 3. Reduce redundant GO terms ----------
format_unreduced_terms <- function(terms) {
  terms %>% transmute(
    Subcluster,
    GO_ID = ID,
    Description,
    p.adjust,
    Count,
    GeneRatio,
    n_terms_collapsed = 1L,
    log10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
  )
}

reduce_one_subcluster <- function(sc) {
  sc_terms <- sig_go_results %>%
    filter(Subcluster == sc) %>%
    distinct(ID, .keep_all = TRUE)

  if (nrow(sc_terms) == 0L) return(tibble())
  if (nrow(sc_terms) == 1L) return(format_unreduced_terms(sc_terms))

  sim_matrix <- tryCatch(
    rrvgo::calculateSimMatrix(
      unique(sc_terms$ID),
      orgdb = "org.At.tair.db",
      ont = "BP",
      method = "Wang"
    ),
    error = function(e) {
      message("Semantic-similarity calculation failed for Subcluster ", sc,
              "; retaining significant terms without reduction. Error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(sim_matrix) || nrow(sim_matrix) == 0L) {
    return(format_unreduced_terms(sc_terms))
  }

  # Retain terms for which a similarity score is unavailable instead of silently dropping them.
  terms_with_similarity <- sc_terms %>% filter(ID %in% rownames(sim_matrix))
  terms_without_similarity <- sc_terms %>% filter(!ID %in% rownames(sim_matrix))

  if (nrow(terms_with_similarity) == 0L) {
    return(format_unreduced_terms(sc_terms))
  }
  if (nrow(terms_with_similarity) == 1L) {
    return(bind_rows(
      format_unreduced_terms(terms_with_similarity),
      format_unreduced_terms(terms_without_similarity)
    ))
  }

  sim_matrix <- sim_matrix[terms_with_similarity$ID,
                           terms_with_similarity$ID, drop = FALSE]
  scores <- setNames(
    -log10(pmax(terms_with_similarity$p.adjust, .Machine$double.xmin)),
    terms_with_similarity$ID
  )

  reduced <- tryCatch(
    rrvgo::reduceSimMatrix(
      sim_matrix,
      scores = scores,
      threshold = rrvgo_threshold,
      orgdb = "org.At.tair.db"
    ) %>% as.data.frame(),
    error = function(e) {
      message("Semantic-similarity reduction failed for Subcluster ", sc,
              "; retaining significant terms without reduction. Error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(reduced) || nrow(reduced) == 0L ||
      !all(c("go", "parent") %in% names(reduced))) {
    return(format_unreduced_terms(sc_terms))
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
      Subcluster = sc,
      GO_ID = parent,
      Description,
      p.adjust,
      Count,
      GeneRatio,
      n_terms_collapsed = as.integer(n_terms_collapsed),
      log10FDR = -log10(pmax(p.adjust, .Machine$double.xmin))
    )

  bind_rows(representative_terms, format_unreduced_terms(terms_without_similarity))
}

reduced_go <- bind_rows(lapply(seq_len(k_subclusters), reduce_one_subcluster))
reduced_file <- paste0(out_prefix, "_reduced_terms.csv")

if (nrow(reduced_go) == 0L) {
  write.csv(reduced_go, reduced_file, row.names = FALSE)
  stop("No representative GO terms remained after semantic-similarity reduction.")
}

write.csv(reduced_go, reduced_file, row.names = FALSE)

# ----------------------------- 4. Prepare dotplot -------------------
plot_data <- reduced_go %>%
  group_by(Subcluster) %>%
  arrange(desc(log10FDR), .by_group = TRUE) %>%
  { if (is.finite(max_terms_per_subcluster)) slice_head(., n = max_terms_per_subcluster) else . } %>%
  ungroup() %>%
  mutate(
    Subcluster = factor(
      Subcluster,
      levels = seq_len(k_subclusters),
      labels = c("Sc1", "Subcluster 2", "Sc3")
    ),
    # Unique key prevents duplicate GO descriptions across facets from creating duplicate factor levels.
    TermKey = paste0(Description, " [GO:", GO_ID, ";SC:", as.integer(Subcluster), "]")
  )

plot_data$TermKey <- factor(plot_data$TermKey, levels = rev(unique(plot_data$TermKey)))

p_go <- ggplot(plot_data, aes(x = log10FDR, y = TermKey)) +
  geom_point(aes(size = Count, colour = log10FDR), alpha = 0.95) +
  facet_grid(
    rows = vars(Subcluster),
    scales = "free_y",
    space = "free_y",
    drop = FALSE
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
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(sub(" \\[GO:.*$", "", x), width = 28)
  ) +
  labs(
    title = paste0("GO BP: Cluster ", parent_cluster, " subclusters"),
    x = expression(-log[10](FDR) ~ "of representative term"),
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(size = 15, hjust = 0, face = "bold"),
    strip.background = element_rect(fill = "grey85", colour = "grey45", linewidth = 0.5),
    strip.text.y = element_text(angle = 270, size = 16, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
    axis.text.y = element_text(size = 14.5, colour = "black"),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 15),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    panel.spacing.y = grid::unit(0.18, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )

print(p_go)

# ----------------------------- 5. Save ------------------------------
n_terms_plot <- nrow(plot_data)
# Allocate vertical space for each displayed GO term to keep labels legible.
plot_height <- max(8.0, 1.8 + 0.24 * n_terms_plot)

png_file <- paste0(out_prefix, "_dotplot.png")
pdf_file <- paste0(out_prefix, "_dotplot.pdf")

ggsave(
  png_file, p_go,
  width = 6, height = plot_height, units = "in", dpi = 300, bg = "white"
)
ggsave(
  pdf_file, p_go,
  width = 3.5, height = plot_height, units = "in", bg = "white"
)

message("Saved GO enrichment table: ", all_results_file)
message("Saved non-redundant GO terms: ", reduced_file)
message("Saved dotplot: ", png_file, " and ", pdf_file)
message("Background note: 'top2000' assesses enrichment relative to all genes in the original k=6 clustering; 'parent_cluster' focuses more directly on enrichment relative to the parent cluster.")
