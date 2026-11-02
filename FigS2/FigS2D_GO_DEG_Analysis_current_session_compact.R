############################################################################
# GO analysis of DESeq2 DEGs: NT / LAT (P/T) / REC (P/T +1d), UPF3 vs Col-0
############################################################################

# Required packages must already be installed. Install separately if needed:
# BiocManager::install(c("clusterProfiler", "org.At.tair.db"))
# install.packages(c("ggplot2", "dplyr", "tidyr", "ggh4x"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggh4x)
})

# ----------------------------- Settings ------------------------------------
fdr_cutoff <- 0.05
lfc_cutoff <- 1
min_input_genes <- 5
min_gene_set_size <- 5
max_gene_set_size <- 500
max_terms_per_panel <- 10
simplify_cutoff <- 0.7
output_dir <- "."
figure_prefix <- "GO_dotplot_UPF3_vs_Col0"

# ----------------------------- Checks --------------------------------------
required_objects <- c("dds", "res_NT", "res_LAT", "res_REC")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects)) {
  stop("Missing required object(s) in the current R session: ",
       paste(missing_objects, collapse = ", "),
       ". Run the DESeq2 script first.")
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^gene:", "", x)
  x <- sub("\\.\\d+$", "", x)
  toupper(x)
}

# ------------------------ 1. Define DEG gene sets -------------------------
# DEG definition: adjusted p < 0.05 and |log2FC| >= 1.
prepare_go_genes <- function(res, treatment_name) {
  df <- as.data.frame(res)
  df$gene_id <- clean_gene_id(rownames(df))
  df <- df %>%
    filter(!is.na(padj), !is.na(log2FoldChange),
           padj < fdr_cutoff, abs(log2FoldChange) >= lfc_cutoff,
           !is.na(gene_id), gene_id != "")
  list(
    Up = unique(df$gene_id[df$log2FoldChange > 0]),
    Down = unique(df$gene_id[df$log2FoldChange < 0])
  )
}

go_genes_NT  <- prepare_go_genes(res_NT, "NT")
go_genes_LAT <- prepare_go_genes(res_LAT, "LAT")
go_genes_REC <- prepare_go_genes(res_REC, "REC")

gene_sets <- list(
  NT_Up = go_genes_NT$Up,
  NT_Down = go_genes_NT$Down,
  LAT_Up = go_genes_LAT$Up,
  LAT_Down = go_genes_LAT$Down,
  REC_Up = go_genes_REC$Up,
  REC_Down = go_genes_REC$Down
)

# Background: all genes that passed the original DESeq2 pre-filter.
universe_genes <- unique(clean_gene_id(rownames(dds)))
universe_genes <- universe_genes[!is.na(universe_genes) & universe_genes != ""]

# Save the exact input lists used for the GO tests.
gene_sets_table <- bind_rows(lapply(names(gene_sets), function(nm) {
  bits <- strsplit(nm, "_", fixed = TRUE)[[1]]
  ids <- unique(gene_sets[[nm]])
  if (!length(ids)) {
    return(data.frame(
      gene_set = character(), treatment = character(),
      direction = character(), gene_id = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    gene_set = rep(nm, length(ids)),
    treatment = rep(bits[1], length(ids)),
    direction = rep(bits[2], length(ids)),
    gene_id = ids,
    stringsAsFactors = FALSE
  )
}))
write.csv(gene_sets_table,
          file.path(output_dir, "GO_DEG_input_gene_lists.csv"),
          row.names = FALSE)

cat("DESeq2 background genes:", length(universe_genes), "\n")
cat("DEG counts by set:\n")
print(data.frame(gene_set = names(gene_sets), n_DEGs = lengths(gene_sets)), row.names = FALSE)

# ----------------------- 2. Run GO enrichment -----------------------------
run_go <- function(genes, universe, ontology) {
  genes <- intersect(unique(genes), universe)
  if (length(genes) < min_input_genes) return(NULL)

  clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    readable = FALSE
  )
}

ontologies <- c("BP", "CC", "MF")
go_results <- list()
for (set_name in names(gene_sets)) {
  for (ont in ontologies) {
    result_name <- paste(set_name, ont, sep = "_")
    go_results[[result_name]] <- run_go(gene_sets[[set_name]], universe_genes, ont)
    n_terms <- if (is.null(go_results[[result_name]])) {
      0L
    } else {
      nrow(as.data.frame(go_results[[result_name]]))
    }
    cat(result_name, ":", n_terms, "GO terms\n")
  }
}

# Convert results to a common table, retaining metadata even for empty analyses.
empty_go_table <- data.frame(
  ID = character(), Description = character(), GeneRatio = character(),
  BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
  qvalue = numeric(), geneID = character(), Count = integer(),
  stringsAsFactors = FALSE
)

results_to_table <- function(result_list, simplified = FALSE, significant_only = FALSE) {
  bind_rows(lapply(names(result_list), function(nm) {
    obj <- result_list[[nm]]
    if (is.null(obj)) return(NULL)
    df <- as.data.frame(obj)
    if (!nrow(df)) return(NULL)
    if (significant_only) {
      df <- df %>% filter(!is.na(p.adjust), p.adjust < fdr_cutoff)
      if (!nrow(df)) return(NULL)
    }
    bits <- strsplit(nm, "_", fixed = TRUE)[[1]]
    # First part is treatment; last part is ontology; middle is direction.
    ont <- bits[length(bits)]
    direction <- bits[length(bits) - 1L]
    treatment <- paste(bits[seq_len(length(bits) - 2L)], collapse = "_")
    df %>% mutate(
      analysis = nm,
      treatment_code = treatment,
      treatment_label = recode(treatment, NT = "NT", LAT = "P/T", REC = "P/T +1d"),
      direction = direction,
      ontology = ont,
      logFDR = -log10(pmax(p.adjust, 1e-300)),
      .before = 1
    )
  }))
}

all_go_table <- results_to_table(go_results)
if (is.null(all_go_table) || !nrow(all_go_table)) {
  all_go_table <- empty_go_table %>% mutate(
    analysis = character(), treatment_code = character(), treatment_label = character(),
    direction = character(), ontology = character(), logFDR = numeric(), .before = 1
  )
}
write.csv(all_go_table, file.path(output_dir, "GO_results_all_tests.csv"), row.names = FALSE)

# Save summary for all 18 set/ontology combinations.
go_summary <- bind_rows(lapply(names(go_results), function(nm) {
  obj <- go_results[[nm]]
  df <- if (is.null(obj)) empty_go_table else as.data.frame(obj)
  bits <- strsplit(nm, "_", fixed = TRUE)[[1]]
  data.frame(
    analysis = nm,
    treatment_code = bits[1],
    treatment_label = recode(bits[1], NT = "NT", LAT = "P/T", REC = "P/T +1d"),
    direction = bits[2], ontology = bits[3],
    input_DEGs = length(intersect(gene_sets[[paste(bits[1], bits[2], sep = "_")]], universe_genes)),
    GO_terms_returned = nrow(df),
    significant_GO_terms_FDR_lt_0_05 = if (nrow(df)) sum(!is.na(df$p.adjust) & df$p.adjust < fdr_cutoff) else 0L,
    stringsAsFactors = FALSE
  )
}))
write.csv(go_summary, file.path(output_dir, "GO_analysis_summary.csv"), row.names = FALSE)

# ----------------- 3. Simplify redundant GO terms -------------------------
go_results_simplified <- lapply(go_results, function(obj) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (!nrow(df)) return(obj)
  tryCatch(
    clusterProfiler::simplify(
      obj, cutoff = simplify_cutoff, by = "p.adjust",
      select_fun = min, measure = "Wang"
    ),
    error = function(e) {
      warning("GO simplify failed for one analysis; retaining original terms: ", conditionMessage(e))
      obj
    }
  )
})

term_comparison <- bind_rows(lapply(names(go_results), function(nm) {
  original <- go_results[[nm]]
  reduced <- go_results_simplified[[nm]]
  n_orig <- if (is.null(original)) 0L else nrow(as.data.frame(original))
  n_red <- if (is.null(reduced)) 0L else nrow(as.data.frame(reduced))
  data.frame(analysis = nm, Original = n_orig, Simplified = n_red,
             Removed = n_orig - n_red, stringsAsFactors = FALSE)
}))
write.csv(term_comparison,
          file.path(output_dir, "GO_term_simplification_comparison.csv"),
          row.names = FALSE)

sig_simplified_table <- results_to_table(go_results_simplified, significant_only = TRUE)
if (is.null(sig_simplified_table) || !nrow(sig_simplified_table)) {
  sig_simplified_table <- all_go_table[0, , drop = FALSE]
}
write.csv(sig_simplified_table,
          file.path(output_dir, "GO_results_significant_simplified.csv"),
          row.names = FALSE)

# -------------------------- 4. Prepare plot -------------------------------
plot_data <- bind_rows(lapply(names(go_results_simplified), function(nm) {
  obj <- go_results_simplified[[nm]]
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj) %>%
    filter(!is.na(p.adjust), p.adjust < fdr_cutoff)
  if (!nrow(df)) return(NULL)

  bits <- strsplit(nm, "_", fixed = TRUE)[[1]]
  treatment <- bits[1]
  direction <- bits[2]
  ontology <- bits[3]

  df %>%
    mutate(
      treatment_code = treatment,
      Treatment = recode(treatment, NT = "NT", LAT = "P/T", REC = "P/T +1d"),
      Direction = direction,
      Ontology = ontology,
      GeneCount = Count,
      logFDR = -log10(pmax(p.adjust, 1e-300))
    ) %>%
    group_by(Treatment, Direction, Ontology) %>%
    slice_min(order_by = p.adjust, n = max_terms_per_panel, with_ties = FALSE) %>%
    ungroup()
}))

if (is.null(plot_data) || !nrow(plot_data)) {
  stop("No GO terms passed adjusted p < ", fdr_cutoff,
       ". CSV results were saved, but no dotplot can be drawn.")
}

plot_data <- plot_data %>%
  mutate(
    Treatment = factor(Treatment, levels = c("NT", "P/T", "P/T +1d")),
    Direction = factor(Direction, levels = c("Down", "Up")),
    Ontology = factor(Ontology, levels = c("BP", "CC", "MF"))
  )

# Stable term ordering; each facet can use only terms present in its panel.
ordered_terms <- plot_data %>%
  group_by(Ontology, Description) %>%
  summarise(best_FDR = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  arrange(Ontology, best_FDR, Description) %>%
  pull(Description) %>% as.character() %>% unique()
plot_data <- plot_data %>%
  mutate(Description = factor(Description, levels = rev(ordered_terms)))

write.csv(
  plot_data %>% mutate(Description = as.character(Description)),
  file.path(output_dir, "GO_results_plotted_terms.csv"),
  row.names = FALSE
)

# ---------------------------- 5. Plot --------------------------------------
p_go_compact <- ggplot(
  plot_data,
  aes(x = logFDR, y = Description, size = GeneCount, color = logFDR)
) +
  geom_point(alpha = 0.85, stroke = 0) +
  ggh4x::facet_grid2(
    rows = vars(Ontology),
    cols = vars(Treatment, Direction),
    scales = "free_y",
    space = "free_y",
    independent = "none",
    strip = ggh4x::strip_vanilla()
  ) +
  scale_x_continuous(
    name = expression(-log[10](FDR)),
    breaks = scales::pretty_breaks(n = 3),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_color_gradient(
    low = "purple4", high = "red",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(range = c(0.8, 2.7), name = "Number of genes") +
  labs(
    title = expression("GO terms " * italic(upf3) * " vs Col-0"),
    x = expression(-log[10](FDR)),
    y = NULL
  ) +
  theme_bw(base_size = 5.4) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8.4, face = "plain",
                              margin = margin(b = 2)),
    strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.25),
    strip.text.x = element_text(size = 5.2, face = "bold", margin = margin(1, 1, 1, 1)),
    strip.text.y = element_text(size = 5.2, face = "bold", margin = margin(1, 1, 1, 1)),
    axis.text.x = element_text(size = 4.5, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "grey25"),
    axis.title.x = element_text(size = 5.2),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.spacing.x = grid::unit(0.04, "cm"),
    panel.spacing.y = grid::unit(0.08, "cm"),
    legend.position = "right",
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key.height = grid::unit(0.22, "cm"),
    legend.key.width = grid::unit(0.22, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

print(p_go_compact)

plot_width_in <- 5
plot_height_in <- 6

ggsave(
  filename = file.path(output_dir, paste0(figure_prefix, ".png")),
  plot = p_go_compact,
  width = plot_width_in,
  height = plot_height_in,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, paste0(figure_prefix, ".pdf")),
  plot = p_go_compact,
  width = plot_width_in,
  height = plot_height_in,
  units = "in",
  bg = "white"
)

cat("\nSaved files in:", normalizePath(output_dir), "\n")
cat(" - ", paste0(figure_prefix, ".png"), "\n", sep = "")
cat(" - ", paste0(figure_prefix, ".pdf"), "\n", sep = "")
cat(" - GO_DEG_input_gene_lists.csv\n")
cat(" - GO_results_all_tests.csv\n")
cat(" - GO_analysis_summary.csv\n")
cat(" - GO_term_simplification_comparison.csv\n")
cat(" - GO_results_significant_simplified.csv\n")
cat(" - GO_results_plotted_terms.csv\n")
