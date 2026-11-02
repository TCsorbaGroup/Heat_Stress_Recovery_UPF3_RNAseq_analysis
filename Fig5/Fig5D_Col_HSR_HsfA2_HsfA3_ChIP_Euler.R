###############################################################################
# Euler diagram: mutually exclusive Col-0 HSR classes vs HsfA2/HsfA3 ChIP targets
#
# Uses the existing `gene_sets` object from the Col_HSR_five_gene_sets workflow
# whenever possible. The overlapping union class "HS responsive genes" and the
# non-responsive "Other expressed genes" class are intentionally excluded.
#
# HSR classes:
#   Acute    = induced at LAT, not induced at REC
#   Sustained= induced at LAT and REC
#   Delayed  = induced at REC, not induced at LAT
#
# ChIP sets:
#   Broad = clusters c7/c11/c12/c14
#   Core  = clusters c12/c14
#
# This is a set-overlap visualization, not a formal enrichment test. The HSR
# classes are required to be mutually exclusive; each can overlap a ChIP set.
###############################################################################

if (!requireNamespace("eulerr", quietly = TRUE)) {
  stop("The R package 'eulerr' is required. Install it before running this script.")
}

suppressPackageStartupMessages({
  library(eulerr)
})

# ------------------------------ Settings -------------------------------------
chip_file  <- "13059_2023_2970_MOESM3_ESM.xlsx"
chip_sheet <- "gene-clusters-with-short-descri"
output_prefix <- "Col_HSR_three_classes_vs_HsfA23_ChIP"
message("Output directory: ", normalizePath(getwd(), mustWork = TRUE))

acute_label     <- "Acute HS responsive genes"
sustained_label <- "Sustained HS responsive genes"
delayed_label   <- "Delayed HS responsive genes"

# --------------------------- Helpers and checks ------------------------------
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

extract_agi_ids <- function(x) {
  x <- as.character(x)
  matches <- regmatches(
    x,
    gregexpr("AT[1-5CM]G[0-9]{5}(?:\\.[0-9]+)?", x,
             ignore.case = TRUE, perl = TRUE)
  )
  clean_gene_set(unlist(matches, use.names = FALSE))
}

normalise_label <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", trimws(as.character(x))))
}

write_gene_set <- function(genes, path) {
  write.csv(data.frame(gene = sort(unique(genes))), path, row.names = FALSE)
}

# --------------------- 1. Extract the three HSR classes ----------------------
hsr_classes <- NULL
hsr_source <- NULL

gs_ok <- exists("gene_sets", inherits = TRUE)
if (gs_ok) {
  gs_candidate <- get("gene_sets", inherits = TRUE)
  gs_ok <- is.data.frame(gs_candidate) && all(c("gene", "GeneGroup") %in% names(gs_candidate))
} else {
  gs_candidate <- NULL
}

if (gs_ok) {
  gs <- gs_candidate
  groups_chr <- trimws(as.character(gs$GeneGroup))
  labels_present <- unique(groups_chr[!is.na(groups_chr) & nzchar(groups_chr)])
  labels_norm <- normalise_label(labels_present)

  aliases <- list(
    Acute = normalise_label(c(acute_label, "Acute", "Transient", "Transient acute", "LAT only")),
    Sustained = normalise_label(c(sustained_label, "Sustained", "Sustained HSR", "LAT and REC")),
    Delayed = normalise_label(c(delayed_label, "Delayed", "Delayed HSR", "REC only"))
  )
  matched <- lapply(aliases, function(a) labels_present[labels_norm %in% a])
  if (!length(matched$Acute)) matched$Acute <- labels_present[grepl("acute|transient", labels_norm)]
  if (!length(matched$Sustained)) matched$Sustained <- labels_present[grepl("sustain", labels_norm)]
  if (!length(matched$Delayed)) matched$Delayed <- labels_present[grepl("delayed", labels_norm)]

  if (all(lengths(matched) == 1L)) {
    hsr_classes <- list(
      Acute = clean_gene_set(gs$gene[groups_chr == matched$Acute]),
      Sustained = clean_gene_set(gs$gene[groups_chr == matched$Sustained]),
      Delayed = clean_gene_set(gs$gene[groups_chr == matched$Delayed])
    )
    hsr_source <- "existing gene_sets table"
  } else {
    message("`gene_sets` exists, but its three HSR class labels were not identified uniquely; trying alternative sources.")
  }
} else {
  message("`gene_sets` is unavailable or lacks gene/GeneGroup columns; trying alternative sources.")
}

if (is.null(hsr_classes) && exists("gene_classification", inherits = TRUE)) {
  gc <- get("gene_classification", inherits = TRUE)
  required_gc <- c("gene", "is_acute", "is_sustained", "is_delayed")
  if (is.data.frame(gc) && all(required_gc %in% names(gc))) {
    hsr_classes <- list(
      Acute = clean_gene_set(gc$gene[!is.na(gc$is_acute) & gc$is_acute]),
      Sustained = clean_gene_set(gc$gene[!is.na(gc$is_sustained) & gc$is_sustained]),
      Delayed = clean_gene_set(gc$gene[!is.na(gc$is_delayed) & gc$is_delayed])
    )
    hsr_source <- "existing gene_classification flags"
  }
}

if (is.null(hsr_classes) && exists("dds", inherits = TRUE)) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("`gene_sets` and `gene_classification` are unavailable; package DESeq2 is required to rebuild the classes from `dds`.")
  }
  dds_obj <- get("dds", inherits = TRUE)
  if (!inherits(dds_obj, "DESeqDataSet")) stop("Object `dds` must be a DESeqDataSet.")
  cd <- as.data.frame(SummarizedExperiment::colData(dds_obj))
  if (!all(c("genotype", "treatment") %in% names(cd))) {
    stop("colData(dds) must contain 'genotype' and 'treatment' to reconstruct the Col-0 HSR classes.")
  }
  if (!is.factor(cd$genotype) || !is.factor(cd$treatment)) {
    stop("To reconstruct HSR classes from `dds`, genotype and treatment must be factors with verified reference levels.")
  }
  genotype_ref <- tolower(as.character(levels(cd$genotype)[1]))
  treatment_ref <- as.character(levels(cd$treatment)[1])
  if (!genotype_ref %in% c("col", "col-0", "wt", "wildtype", "wild_type")) {
    stop("Cannot reconstruct Col-0 treatment response: the reference genotype is '",
         levels(cd$genotype)[1], "', not Col-0/WT.")
  }
  if (treatment_ref != "NT") {
    stop("Cannot reconstruct treatment induction classes: the reference treatment is '",
         treatment_ref, "', not NT.")
  }

  alpha_cutoff <- 0.05
  lfc_cutoff <- 1
  required_coef_names <- c("treatment_LAT_vs_NT", "treatment_REC_vs_NT")
  missing_coef_names <- setdiff(required_coef_names, DESeq2::resultsNames(dds_obj))
  if (length(missing_coef_names)) {
    stop("Cannot rebuild classes from dds; missing DESeq2 coefficient(s): ",
         paste(missing_coef_names, collapse = ", "))
  }
  rlat <- as.data.frame(DESeq2::results(dds_obj, name = "treatment_LAT_vs_NT"))
  rrec <- as.data.frame(DESeq2::results(dds_obj, name = "treatment_REC_vs_NT"))
  rlat$gene <- clean_gene_id(rownames(rlat))
  rrec$gene <- clean_gene_id(rownames(rrec))
  if (anyDuplicated(rlat$gene) || anyDuplicated(rrec$gene)) {
    stop("Duplicate gene IDs appeared while reconstructing HSR classes from DESeq2 results.")
  }
  lat_induced <- !is.na(rlat$padj) & rlat$padj < alpha_cutoff &
    !is.na(rlat$log2FoldChange) & rlat$log2FoldChange > lfc_cutoff
  rec_induced <- !is.na(rrec$padj) & rrec$padj < alpha_cutoff &
    !is.na(rrec$log2FoldChange) & rrec$log2FoldChange > lfc_cutoff
  lat_genes <- unique(rlat$gene[lat_induced])
  rec_genes <- unique(rrec$gene[rec_induced])
  hsr_classes <- list(
    Acute = clean_gene_set(setdiff(lat_genes, rec_genes)),
    Sustained = clean_gene_set(intersect(lat_genes, rec_genes)),
    Delayed = clean_gene_set(setdiff(rec_genes, lat_genes))
  )
  hsr_source <- "reconstructed from dds: padj < 0.05 and log2FC > 1, Col-0 vs NT reference"
}

if (is.null(hsr_classes)) {
  stop(
    "Could not recover gene membership for the three HSR classes. `gene_set_sizes` only contains counts, not gene IDs. ",
    "Keep `gene_sets` or `gene_classification` in the global environment, or make `dds` available."
  )
}

hsr_classes <- lapply(hsr_classes, clean_gene_set)
if (!identical(names(hsr_classes), c("Acute", "Sustained", "Delayed"))) {
  hsr_classes <- hsr_classes[c("Acute", "Sustained", "Delayed")]
}
if (any(lengths(hsr_classes) == 0L)) {
  print(lengths(hsr_classes))
  stop("At least one HSR class is empty; check the source objects and induction criteria.")
}

cat("HSR class source: ", hsr_source, "\n", sep = "")
cat("HSR class sizes: ", paste(names(hsr_classes), lengths(hsr_classes), sep = " = ", collapse = "; "), "\n", sep = "")

class_pairs <- combn(names(hsr_classes), 2, simplify = FALSE)
class_pair_overlaps <- vapply(class_pairs, function(pair) {
  length(intersect(hsr_classes[[pair[1]]], hsr_classes[[pair[2]]]))
}, integer(1))
names(class_pair_overlaps) <- vapply(class_pairs, function(pair) paste(pair, collapse = " & "), character(1))
if (any(class_pair_overlaps != 0L)) {
  print(class_pair_overlaps)
  stop("At least two HSR classes overlap after gene-ID cleaning.")
}

# Export the exact HSR membership used by the figures.
hsr_membership <- bind_rows(lapply(names(hsr_classes), function(nm) {
  data.frame(HSR_class = nm, gene = hsr_classes[[nm]], stringsAsFactors = FALSE)
}))
write.csv(hsr_membership, paste0(output_prefix, "_HSR_class_membership.csv"), row.names = FALSE)

# -------------------------- 2. Load ChIP gene sets ---------------------------
if (exists("hsfa23_all_relevant", inherits = TRUE) && exists("hsfa23_core", inherits = TRUE)) {
  chip_all <- clean_gene_set(get("hsfa23_all_relevant", inherits = TRUE))
  chip_core <- clean_gene_set(get("hsfa23_core", inherits = TRUE))
  chip_source <- "Existing objects: hsfa23_all_relevant / hsfa23_core"
} else {
  if (!file.exists(chip_file)) {
    stop("Neither `hsfa23_all_relevant`/`hsfa23_core` nor the ChIP Excel file were found.\nExpected file: ",
         normalizePath(chip_file, mustWork = FALSE),
         "\nSet `chip_file` or load the ChIP gene sets first.")
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("The package 'readxl' is needed to read the ChIP Excel file. Install it with install.packages('readxl').")
  }
  chip_clusters <- as.data.frame(readxl::read_excel(chip_file, sheet = chip_sheet), stringsAsFactors = FALSE)
  if (ncol(chip_clusters) < 2L) stop("The ChIP worksheet must have gene ID and cluster columns.")
  chip_gene_raw <- chip_clusters[[1]]
  chip_cluster_raw <- trimws(as.character(chip_clusters[[2]]))
  cluster_text <- ifelse(grepl("[0-9]+", chip_cluster_raw),
                         sub("^[^0-9]*([0-9]+).*", "\\1", chip_cluster_raw), NA_character_)
  chip_cluster_num <- suppressWarnings(as.integer(cluster_text))
  chip_all <- extract_agi_ids(chip_gene_raw[chip_cluster_num %in% c(7, 11, 12, 14)])
  chip_core <- extract_agi_ids(chip_gene_raw[chip_cluster_num %in% c(12, 14)])
  chip_source <- paste0("Excel: ", chip_file, " [", chip_sheet, "]")
}

if (length(chip_all) == 0L) stop("The broad ChIP target set is empty; check gene IDs and cluster labels.")
if (length(chip_core) == 0L) warning("The core ChIP target set is empty; its diagram may not be informative.")
core_not_broad <- setdiff(chip_core, chip_all)
if (length(core_not_broad)) {
  stop("The core ChIP set contains genes absent from the broad c7/c11/c12/c14 set; check the ChIP sets. Examples: ",
       paste(head(core_not_broad, 10), collapse = ", "))
}

cat("ChIP source: ", chip_source, "\n", sep = "")
cat("Broad ChIP genes: ", length(chip_all), "; core ChIP genes: ", length(chip_core), "\n", sep = "")
write_gene_set(chip_all, paste0(output_prefix, "_ChIP_broad_gene_set.csv"))
write_gene_set(chip_core, paste0(output_prefix, "_ChIP_core_gene_set.csv"))

# ------------------------- 3. Define disjoint regions ------------------------
make_euler_regions <- function(chip_set) {
  acute <- hsr_classes$Acute
  sustained <- hsr_classes$Sustained
  delayed <- hsr_classes$Delayed
  hsr_union <- unique(c(acute, sustained, delayed))

  # With input = "disjoint", values are counts in exact/non-overlapping regions.
  c(
    Acute = length(setdiff(acute, chip_set)),
    Sustained = length(setdiff(sustained, chip_set)),
    Delayed = length(setdiff(delayed, chip_set)),
    ChIP = length(setdiff(chip_set, hsr_union)),
    "Acute&ChIP" = length(intersect(acute, chip_set)),
    "Sustained&ChIP" = length(intersect(sustained, chip_set)),
    "Delayed&ChIP" = length(intersect(delayed, chip_set)),
    "Acute&Sustained" = 0,
    "Acute&Delayed" = 0,
    "Sustained&Delayed" = 0,
    "Acute&Sustained&Delayed" = 0,
    "Acute&Sustained&ChIP" = 0,
    "Acute&Delayed&ChIP" = 0,
    "Sustained&Delayed&ChIP" = 0,
    "Acute&Sustained&Delayed&ChIP" = 0
  )
}

make_overlap_summary <- function(chip_set, chip_label) {
  data.frame(
    ChIP_set = chip_label,
    HSR_class = names(hsr_classes),
    n_HSR_genes = as.integer(lengths(hsr_classes)),
    n_in_ChIP = vapply(hsr_classes, function(gs) length(intersect(gs, chip_set)), integer(1)),
    n_not_in_ChIP = vapply(hsr_classes, function(gs) length(setdiff(gs, chip_set)), integer(1)),
    pct_in_ChIP = vapply(hsr_classes, function(gs) {
      if (!length(gs)) NA_real_ else 100 * length(intersect(gs, chip_set)) / length(gs)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

# ----------------------------- 4. Draw diagrams ------------------------------
plot_euler_diagram <- function(chip_set, suffix, title) {
  regions <- make_euler_regions(chip_set)
  fit <- eulerr::euler(regions, input = "disjoint", shape = "ellipse")
  message("\n", title, " | diagError = ", signif(fit$diagError, 3),
          "; stress = ", signif(fit$stress, 3))

  png_file <- file.path(getwd(), paste0(output_prefix, "_", suffix, ".png"))
  pdf_file <- file.path(getwd(), paste0(output_prefix, "_", suffix, ".pdf"))

  labels <- c(
    "Acute HS responsive",
    "Sustained HS responsive",
    "Delayed HS responsive",
    "HsfA2/HsfA3 ChIP"
  )
  fills <- list(fill = c("#E69F00", "#D55E00", "#0072B2", "#7B61A8"), alpha = 0.42)

  draw_diagram <- function() {
    p <- plot(
      fit,
      fills = fills,
      edges = list(col = "grey25", lwd = 1.3),
      labels = labels,
      quantities = list(font = 2, cex = 0.9),
      main = title,
      legend = FALSE,
      complement = FALSE
    )
    print(p)
  }

  grDevices::png(png_file, width = 2400, height = 1900, res = 300, bg = "white")
  tryCatch(
    draw_diagram(),
    error = function(e) {
      grDevices::dev.off()
      stop("Failed while drawing PNG diagram '", png_file, "': ", conditionMessage(e), call. = FALSE)
    }
  )
  grDevices::dev.off()

  grDevices::pdf(pdf_file, width = 8, height = 6.4, useDingbats = FALSE, onefile = TRUE)
  tryCatch(
    draw_diagram(),
    error = function(e) {
      grDevices::dev.off()
      stop("Failed while drawing PDF diagram '", pdf_file, "': ", conditionMessage(e), call. = FALSE)
    }
  )
  grDevices::dev.off()

  for (f in c(png_file, pdf_file)) {
    if (!file.exists(f)) {
      stop("The diagram file was not created: ", f, call. = FALSE)
    }
    file_size <- file.info(f)$size
    if (is.na(file_size) || file_size < 1000) {
      stop(
        "The diagram file exists but appears empty or incomplete: ", f,
        " (size = ", ifelse(is.na(file_size), "NA", paste(file_size, "bytes")),
        "). Ensure the output drive has free space and permission to write files.",
        call. = FALSE
      )
    }
    message("Saved diagram: ", normalizePath(f, mustWork = TRUE),
            " (", round(file_size / 1024, 1), " KB)")
  }

  write.csv(
    data.frame(ChIP_set = suffix, Region = names(regions), n_genes = as.integer(regions)),
    file.path(getwd(), paste0(output_prefix, "_", suffix, "_region_counts.csv")),
    row.names = FALSE
  )
  invisible(fit)
}

# ----------------------------- 5. Run and export -----------------------------
cat("\nGene-set sizes used in diagrams:\n")
print(data.frame(HSR_class = names(hsr_classes), n_genes = as.integer(lengths(hsr_classes))))
cat("\nPairwise HSR intersections (must all be zero):\n")
print(class_pair_overlaps)

chip_all_summary <- make_overlap_summary(chip_all, "all relevant: c7/c11/c12/c14")
chip_core_summary <- make_overlap_summary(chip_core, "core: c12/c14")
overlap_summary <- rbind(chip_all_summary, chip_core_summary)
print(overlap_summary, row.names = FALSE)
write.csv(overlap_summary, paste0(output_prefix, "_overlap_summary.csv"), row.names = FALSE)

fit_all <- plot_euler_diagram(
  chip_all, "all_relevant",
  "Three Col-0 HSR classes and all relevant HsfA2/HsfA3 ChIP targets"
)
if (length(chip_core) > 0L) {
  fit_core <- plot_euler_diagram(
    chip_core, "core",
    "Three Col-0 HSR classes and core HsfA2/HsfA3 ChIP targets"
  )
}

message("Done. Output prefix: ", output_prefix)
