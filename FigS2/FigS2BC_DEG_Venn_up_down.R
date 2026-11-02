###############################################################################
# DEG Venn diagrams (NT, P/T, P/T+1d) and overlap with empirical NMD targets
# DEG definition agreed for this analysis:
#   adjusted p-value (padj) < 0.05 AND abs(log2FoldChange) >= 1
# Positive log2FC = higher in UPF3; negative log2FC = higher in Col-0.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------- Settings ------------------------------------
alpha_cutoff <- 0.05
lfc_cutoff <- 1
output_dir <- getwd()
figure_width_in <- 7.5
figure_height_in <- 3.33
png_dpi <- 300
tiff_dpi <- 600

out_path <- function(filename) file.path(output_dir, filename)

# ------------------------- Check / obtain DESeq2 results ----------------------
# Prefer the res_* objects already created by the DESeq2 script, to preserve the
# exact results used in the current session. If any are missing, reconstruct
# the missing result objects from dds using the same genotype contrasts.
if (!exists("dds", inherits = TRUE) &&
    !all(vapply(c("res_NT", "res_LAT", "res_REC"), exists,
                logical(1), inherits = TRUE))) {
  stop("Cannot find the DESeq2 object `dds` or all of `res_NT`, `res_LAT`, `res_REC`. Run the first DESeq2 script first.")
}

if (!exists("res_NT", inherits = TRUE)) {
  if (!exists("dds", inherits = TRUE)) stop("`res_NT` is missing and `dds` is not available.")
  res_NT <- DESeq2::results(dds, name = "genotype_UPF3_vs_Col", alpha = alpha_cutoff)
}
if (!exists("res_LAT", inherits = TRUE)) {
  if (!exists("dds", inherits = TRUE)) stop("`res_LAT` is missing and `dds` is not available.")
  res_LAT <- DESeq2::results(
    dds,
    contrast = list(c("genotype_UPF3_vs_Col", "genotypeUPF3.treatmentLAT")),
    alpha = alpha_cutoff
  )
}
if (!exists("res_REC", inherits = TRUE)) {
  if (!exists("dds", inherits = TRUE)) stop("`res_REC` is missing and `dds` is not available.")
  res_REC <- DESeq2::results(
    dds,
    contrast = list(c("genotype_UPF3_vs_Col", "genotypeUPF3.treatmentREC")),
    alpha = alpha_cutoff
  )
}

clean_gene_id <- function(x) {
  x <- sub("^gene:", "", as.character(x), ignore.case = TRUE)
  x <- sub("\\.\\d+$", "", x)
  toupper(trimws(x))
}

# ---------------------- Read empirical NMD target list -----------------------
# Expected input: NMD_target_genes.txt in the working/output directory.
# Supports one gene ID per line or IDs separated by whitespace, commas,
# semicolons, or tabs. Only Arabidopsis TAIR locus IDs are retained; optional
# gene: prefixes and numeric version suffixes are normalized.
nmd_target_file <- file.path(output_dir, "NMD_target_genes.txt")
nmd_target_file_loaded <- FALSE
external_target_ids <- character(0)

if (!file.exists(nmd_target_file)) {
  message(
    "NMD target file not found: ", nmd_target_file,
    "\nThe DEG Venn diagrams/counts will still be created, but the NMD overlap analysis will be skipped."
  )
} else {
  nmd_lines <- readLines(nmd_target_file, warn = FALSE, encoding = "UTF-8")
  nmd_lines <- sub("^\ufeff", "", nmd_lines)
  nmd_lines <- sub("#.*$", "", nmd_lines)  # remove comments
  nmd_text <- paste(nmd_lines, collapse = " ")
  nmd_tokens <- unlist(strsplit(trimws(nmd_text), "[[:space:],;]+"), use.names = FALSE)
  nmd_tokens <- nmd_tokens[nzchar(nmd_tokens)]
  nmd_tokens <- sub("^>", "", nmd_tokens)  # tolerate FASTA-style prefixes
  nmd_tokens <- sub("^gene:", "", nmd_tokens, ignore.case = TRUE)

  # Nuclear, chloroplast, and mitochondrial Arabidopsis locus IDs; optional version suffix.
  is_tair_id <- grepl("^AT([1-5]G|CG|MG)[0-9]{5}(\\.\\d+)?$", nmd_tokens, ignore.case = TRUE)
  external_target_ids <- unique(clean_gene_id(nmd_tokens[is_tair_id]))
  external_target_ids <- external_target_ids[!is.na(external_target_ids) & nzchar(external_target_ids)]

  if (length(external_target_ids) == 0L) {
    warning(
      "No valid Arabidopsis TAIR gene IDs were recognized in ", nmd_target_file,
      ". Check that the file contains gene IDs such as AT1G01010. NMD overlap analysis will be skipped."
    )
  } else {
    nmd_target_file_loaded <- TRUE
    cat("Loaded NMD target IDs from ", normalizePath(nmd_target_file), ": ",
        length(external_target_ids), " unique genes.\n", sep = "")
    if (any(!is_tair_id)) {
      message("Ignored ", sum(!is_tair_id), " non-TAIR token(s) in the NMD target file.")
    }
  }
}

# Get unique gene IDs meeting the agreed DEG definition and split by direction.
get_deg_directions <- function(res) {
  df <- as.data.frame(res)
  df$gene_id <- clean_gene_id(rownames(df))
  eligible <- !is.na(df$padj) & df$padj < alpha_cutoff &
    !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= lfc_cutoff
  list(
    up = unique(df$gene_id[eligible & df$log2FoldChange >= lfc_cutoff]),
    down = unique(df$gene_id[eligible & df$log2FoldChange <= -lfc_cutoff])
  )
}

deg_nt  <- get_deg_directions(res_NT)
deg_lat <- get_deg_directions(res_LAT)
deg_rec <- get_deg_directions(res_REC)

venn_sets <- list(
  Up = list("NT" = deg_nt$up, "P/T" = deg_lat$up, "P/T +1d" = deg_rec$up),
  Down = list("NT" = deg_nt$down, "P/T" = deg_lat$down, "P/T +1d" = deg_rec$down)
)

cat("DEG counts using padj <", alpha_cutoff, " and |log2FC| >=", lfc_cutoff, ":\n", sep = "")
print(data.frame(
  treatment = c("NT", "P/T", "P/T +1d"),
  up_in_UPF3 = c(length(deg_nt$up), length(deg_lat$up), length(deg_rec$up)),
  down_in_UPF3 = c(length(deg_nt$down), length(deg_lat$down), length(deg_rec$down))
), row.names = FALSE)

# ------------------------------- Venn plotting -------------------------------
venn_cols <- c("NT" = "#4DAF4A", "P/T" = "#E3120B", "P/T +1d" = "#377EB8")
venn_regions <- c("100", "010", "001", "110", "101", "011", "111")

venn_circle <- function(cx, cy, r = 1, n = 240) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}

venn_ctr <- data.frame(
  set = names(venn_cols),
  cx = c(-0.55, 0.55, 0),
  cy = c(0.30, 0.30, -0.62),
  stringsAsFactors = FALSE
)

venn_lab_pos <- data.frame(
  region = venn_regions,
  x = c(-1.05, 1.05, 0, 0, -0.62, 0.62, 0),
  y = c(0.62, 0.62, -1.25, 0.82, -0.30, -0.30, 0.05),
  stringsAsFactors = FALSE
)

venn_one <- function(sets, title_text) {
  # Normalize and de-duplicate IDs before calculating membership.
  sets <- lapply(sets, function(x) unique(clean_gene_id(x[!is.na(x) & nzchar(x)])))
  ids <- unique(unlist(sets, use.names = FALSE))

  if (length(ids)) {
    memb <- vapply(sets, function(s) ids %in% s, logical(length(ids)))
    colnames(memb) <- names(sets)
    code <- apply(memb, 1, function(r) paste0(as.integer(r), collapse = ""))
  } else {
    code <- character(0)
  }

  cnt <- table(factor(code, levels = venn_regions))
  lab <- venn_lab_pos
  lab$n <- as.integer(cnt[lab$region])

  circ <- do.call(rbind, lapply(seq_len(nrow(venn_ctr)), function(i) {
    cdat <- venn_circle(venn_ctr$cx[i], venn_ctr$cy[i])
    cdat$set <- venn_ctr$set[i]
    cdat
  }))
  circ$set <- factor(circ$set, levels = names(venn_cols))

  setlab <- data.frame(
    set = names(sets),
    n = vapply(sets, length, integer(1)),
    x = c(-1.45, 1.45, 0),
    y = c(1.45, 1.45, -1.85),
    stringsAsFactors = FALSE
  )
  setlab$text <- paste0(setlab$set, "\n(n = ", setlab$n, ")")
  setlab$set <- factor(setlab$set, levels = names(venn_cols))

  title_display <- bquote(
    .(title_text) ~ italic(upf3) ~
      "(" * "total unique: " * .(length(ids)) * ")"
  )
  p <- ggplot() +
    geom_polygon(
      data = circ,
      aes(x = x, y = y, fill = set, colour = set, group = set),
      alpha = 0.25, linewidth = 0.65
    ) +
    geom_text(data = lab, aes(x = x, y = y, label = n), size = 3.0) +
    geom_text(
      data = setlab,
      aes(x = x, y = y, label = text, colour = set),
      fontface = "bold", size = 2.55, lineheight = 0.9
    ) +
    scale_fill_manual(values = venn_cols, guide = "none", drop = FALSE) +
    scale_colour_manual(values = venn_cols, guide = "none", drop = FALSE) +
    coord_equal(xlim = c(-2.1, 2.1), ylim = c(-2.2, 1.9), clip = "off") +
    theme_void(base_size = 8) +
    labs(title = title_display) +
    theme(
      plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5,
                                margin = margin(b = 1.5)),
      plot.margin = margin(1.5, 2, 1.5, 2, unit = "mm")
    )

  counts <- data.frame(
    direction = title_text,
    region_code = lab$region,
    n_genes = lab$n,
    total_unique_genes = length(ids),
    stringsAsFactors = FALSE
  )
  list(plot = p, counts = counts, sets = sets)
}

venn_up <- venn_one(venn_sets$Up, "DEGs-UP in")
venn_dn <- venn_one(venn_sets$Down, "DEGs-DOWN in")

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package `patchwork` is required to arrange the Up and Down Venn diagrams. Install it with install.packages('patchwork').")
}
p_venn <- venn_up$plot + venn_dn$plot + patchwork::plot_layout(ncol = 2)
print(p_venn)

ggsave(out_path("Figure_Venn_DEGs_up_down.png"), plot = p_venn,
       width = figure_width_in, height = figure_height_in, units = "in",
       dpi = png_dpi, bg = "white")
ggsave(out_path("Figure_Venn_DEGs_up_down.pdf"), plot = p_venn,
       width = figure_width_in, height = figure_height_in, units = "in",
       bg = "white")
ggsave(out_path("Figure_Venn_DEGs_up_down.tiff"), plot = p_venn,
       width = figure_width_in, height = figure_height_in, units = "in",
       dpi = tiff_dpi, compression = "lzw", bg = "white")

venn_cnt <- rbind(venn_up$counts, venn_dn$counts)
region_labels <- c(
  "100" = "NT only",
  "010" = "P/T only",
  "001" = "P/T +1d only",
  "110" = "NT + P/T",
  "101" = "NT + P/T+1d",
  "011" = "P/T + P/T+1d",
  "111" = "all three"
)
venn_cnt$region <- unname(region_labels[venn_cnt$region_code])
write.csv(venn_cnt, out_path("DEG_overlap_counts_Venn.csv"), row.names = FALSE)
cat("\nVenn region counts:\n")
print(venn_cnt, row.names = FALSE)
cat("\nSaved figure files to: ", normalizePath(output_dir), "\n", sep = "")

# --------------------- Optional overlap with empirical NMD targets -----------
# The empirical NMD target list is read from NMD_target_genes.txt above.
# If the file is absent or contains no recognized TAIR IDs, Venn plots/counts
# are still made and the NMD-target section is skipped with an explicit message.
if (nmd_target_file_loaded) {
  if (exists("dds", inherits = TRUE)) {
    venn_bg <- unique(clean_gene_id(rownames(dds)))
  } else {
    venn_bg <- unique(clean_gene_id(c(
      rownames(as.data.frame(res_NT)),
      rownames(as.data.frame(res_LAT)),
      rownames(as.data.frame(res_REC))
    )))
  }
  venn_bg <- venn_bg[!is.na(venn_bg) & nzchar(venn_bg)]
  nmd_set <- intersect(unique(clean_gene_id(external_target_ids)), venn_bg)

  fisher_overlap_row <- function(gene_set, background, nmd_targets, direction, state) {
    gene_set <- intersect(unique(clean_gene_id(gene_set)), background)
    in_set <- background %in% gene_set
    in_nmd <- background %in% nmd_targets
    tab <- table(
      factor(in_set, levels = c(TRUE, FALSE)),
      factor(in_nmd, levels = c(TRUE, FALSE))
    )
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    n_in_set <- sum(in_set)
    n_overlap <- sum(in_set & in_nmd)
    data.frame(
      direction = direction,
      state = state,
      n_DEGs = n_in_set,
      n_on_NMD_list = n_overlap,
      pct_on_NMD_list = if (n_in_set > 0) round(100 * n_overlap / n_in_set, 1) else NA_real_,
      n_background = length(background),
      n_NMD_targets_in_background = sum(in_nmd),
      odds_ratio_vs_background = if (is.null(ft)) NA_real_ else unname(ft$estimate),
      CI_low = if (is.null(ft)) NA_real_ else unname(ft$conf.int[1]),
      CI_high = if (is.null(ft)) NA_real_ else unname(ft$conf.int[2]),
      p_fisher = if (is.null(ft)) NA_real_ else ft$p.value,
      stringsAsFactors = FALSE
    )
  }

  venn_nmd <- do.call(rbind, lapply(c("Up", "Down"), function(direction) {
    s <- venn_sets[[direction]]
    state_names <- c(names(s), "union of the three states")
    do.call(rbind, lapply(state_names, function(state) {
      g <- if (state == "union of the three states") unique(unlist(s, use.names = FALSE)) else s[[state]]
      fisher_overlap_row(g, venn_bg, nmd_set, direction, state)
    }))
  }))
  venn_nmd$p_adj_BH <- p.adjust(venn_nmd$p_fisher, method = "BH")
  write.csv(venn_nmd, out_path("DEG_NMD_target_overlap.csv"), row.names = FALSE)
  cat("\nEmpirical NMD-target overlap results:\n")
  print(venn_nmd, row.names = FALSE)
  cat("NMD targets in the tested background: ", length(nmd_set), " of ", length(venn_bg), " background genes.\n", sep = "")
} else {
  message("NMD target list was not loaded from NMD_target_genes.txt; DEG Venn figure/counts were saved, but NMD-target overlap was skipped.")
}
