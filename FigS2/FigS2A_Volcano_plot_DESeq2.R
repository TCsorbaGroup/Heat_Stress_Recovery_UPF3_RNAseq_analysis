# -----------------------------------------------------------------------------
# DESeq2 volcano plots: NT, LAT and REC
# Requires res_NT, res_LAT and res_REC from the DESeq2 analysis script.
# Output: Volcano_NT_LAT_REC.png
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- Settings ----------------------------------------------------------------
output_file <- "Volcano_NT_LAT_REC.png"
plot_width_in <- 7.2
plot_height_in <- 3.1  # upper third of a 10-inch-tall slide, with margins
plot_dpi <- 300

# ---- Check inputs -------------------------------------------------------------
required_results <- c("res_NT", "res_LAT", "res_REC")
missing_results <- required_results[!vapply(required_results, exists, logical(1), inherits = TRUE)]
if (length(missing_results) > 0) {
  stop(
    "Missing DESeq2 result object(s): ", paste(missing_results, collapse = ", "),
    ". Run the DESeq2 analysis script first."
  )
}

# ---- Prepare data -------------------------------------------------------------
prepare_volcano <- function(res, treatment_name) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)

  # A volcano plot needs a finite fold-change coordinate.
  df <- df[!is.na(df$log2FoldChange) & is.finite(df$log2FoldChange), , drop = FALSE]

  # Keep genes whose adjusted p-value is NA; display them along the bottom in a
  # separate colour. Their y-coordinate is only a visual placeholder, not a p-value.
  df$padj_NA <- is.na(df$padj)
  df$minus_log10_padj <- 0

  valid_padj <- !is.na(df$padj) & is.finite(df$padj) & df$padj > 0
  df$minus_log10_padj[valid_padj] <- -log10(df$padj[valid_padj])

  # A padj of exactly zero cannot be plotted as -log10(0). Use a finite value
  # slightly above the largest plotted significance. If all non-NA padj values
  # are zero, use a conservative finite fallback rather than an extreme value.
  zero_padj <- !is.na(df$padj) & is.finite(df$padj) & df$padj == 0
  if (any(zero_padj)) {
    if (any(valid_padj)) {
      pseudo_padj <- max(min(df$padj[valid_padj], na.rm = TRUE) / 10,
                         .Machine$double.xmin)
    } else {
      pseudo_padj <- .Machine$double.eps
    }
    df$minus_log10_padj[zero_padj] <- -log10(pseudo_padj)
  }

  # Mutually exclusive display categories; significance refers to adjusted p-value.
  df$category <- "NS"
  df$category[!df$padj_NA & abs(df$log2FoldChange) >= 1] <- "Log2 FC only"
  df$category[!df$padj_NA & !is.na(df$padj) & df$padj < 0.05] <- "Adjusted p-value only"
  df$category[!df$padj_NA & !is.na(df$padj) & df$padj < 0.05 &
                abs(df$log2FoldChange) >= 1] <- "Adjusted p-value and log2 FC"
  df$category[df$padj_NA] <- "padj NA"

  df$treatment <- treatment_name
  df
}

volcano_NT  <- prepare_volcano(res_NT,  "NT")
volcano_LAT <- prepare_volcano(res_LAT, "P/T")
volcano_REC <- prepare_volcano(res_REC, "P/T +1d")

volcano_all <- rbind(volcano_NT, volcano_LAT, volcano_REC)
volcano_all$treatment <- factor(
  volcano_all$treatment,
  levels = c("NT", "P/T", "P/T +1d")
)
volcano_all$category <- factor(
  volcano_all$category,
  levels = c("NS", "Log2 FC only", "Adjusted p-value only",
             "Adjusted p-value and log2 FC", "padj NA")
)

# ---- Count significantly up/downregulated genes -------------------------------
# Consistent with the accepted DEG definition: BH-adjusted p < 0.05 and |log2FC| >= 1.
count_one_treatment <- function(df) {
  data.frame(
    treatment = as.character(df$treatment[1]),
    down = sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange <= -1),
    up   = sum(!is.na(df$padj) & df$padj < 0.05 & df$log2FoldChange >= 1)
  )
}

gene_counts <- do.call(
  rbind,
  lapply(split(volcano_all, volcano_all$treatment, drop = TRUE), count_one_treatment)
)
gene_counts$treatment <- factor(gene_counts$treatment, levels = c("NT", "P/T", "P/T +1d"))

count_labels_down <- data.frame(
  treatment = gene_counts$treatment,
  label = paste0("Down: n = ", gene_counts$down)
)
count_labels_up <- data.frame(
  treatment = gene_counts$treatment,
  label = paste0("Up: n = ", gene_counts$up)
)

# Print the counts for verification in the R console.
print(gene_counts, row.names = FALSE)

# ---- Plot ---------------------------------------------------------------------
p <- ggplot(
  volcano_all,
  aes(x = log2FoldChange, y = minus_log10_padj, colour = category)
) +
  geom_point(size = 0.48, alpha = 0.72) +
  geom_vline(
    xintercept = c(-1, 1), linetype = "dashed",
    colour = "grey45", linewidth = 0.35
  ) +
  geom_hline(
    yintercept = -log10(0.05), linetype = "dashed",
    colour = "grey45", linewidth = 0.35
  ) +
  facet_wrap(
    ~ treatment,
    nrow = 1,
    labeller = labeller(treatment = as_labeller(
      c(
        "NT" = "'NT' ~ '|' ~ italic(upf3) ~ 'vs' ~ 'Col-0'",
        "P/T" = "'P/T' ~ '|' ~ italic(upf3) ~ 'vs' ~ 'Col-0'",
        "P/T +1d" = "'P/T +1d' ~ '|' ~ italic(upf3) ~ 'vs' ~ 'Col-0'"
      ),
      default = label_parsed
    ))
  ) +
  # Downregulated count: top-left of each panel.
  geom_text(
    data = count_labels_down,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.06, vjust = 1.15,
    size = 3.1, fontface = "bold"
  ) +
  # Upregulated count: top-right of each panel.
  geom_text(
    data = count_labels_up,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = 1.06, vjust = 1.15,
    size = 3.1, fontface = "bold"
  ) +
  scale_colour_manual(
    values = c(
      "NS" = "grey30",
      "Log2 FC only" = "#4DAF4A",
      "Adjusted p-value only" = "#4169E1",
      "Adjusted p-value and log2 FC" = "#F03B3B",
      "padj NA" = "grey75"
    ),
    breaks = c(
      "NS", "Log2 FC only", "Adjusted p-value only",
      "Adjusted p-value and log2 FC", "padj NA"
    ),
    labels = c(
      "NS", "|log2FC| >= 1 only", "padj < 0.05 only",
      "padj < 0.05 + FC", "padj NA"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Differential expression by treatment",
    x = expression(Log[2] ~ "fold change"),
    y = expression(-Log[10] ~ "(adjusted p-value)"),
    colour = NULL,
    caption = "Dashed lines: |log2FC| = 1 and adjusted p = 0.05. padj = NA genes are retained at y = 0."
  ) +
  guides(
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(size = 1.8, alpha = 1)
    )
  ) +
  theme_bw(base_size = 9.5) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "plain",
          size = 10.5,
          margin = margin(b = 2)
        ),
        plot.caption = element_text(
          size = 7.2,
          hjust = 0,
          margin = margin(t = 2)
        ),
        strip.background = element_blank(),
        strip.text = element_text(
          size = 8.2,
          face = "bold"
        ),
        legend.position = "top",
        legend.text = element_text(size = 6.8),
        legend.key.size = grid::unit(0.23, "cm"),
        legend.spacing.x = grid::unit(0.10, "cm"),
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        axis.title = element_text(size = 9.2),
        axis.text = element_text(size = 8.0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(
          colour = "grey90",
          linewidth = 0.25
        ),
        
        panel.spacing = grid::unit(0.16, "cm"),
        plot.margin = margin(3, 5, 3, 5)
      )

print(p)

volcano_gene_counts <- c(
  NT = sum(!is.na(res_NT$log2FoldChange)),
  LAT = sum(!is.na(res_LAT$log2FoldChange)),
  REC = sum(!is.na(res_REC$log2FoldChange))
)

print(volcano_gene_counts)

ggsave(
  filename = output_file,
  plot = p,
  width = plot_width_in,
  height = plot_height_in,
  units = "in",
  dpi = plot_dpi,
  bg = "white",
  limitsize = FALSE
)

message(
  "Saved volcano plot: ", normalizePath(output_file, mustWork = FALSE),
  " (", plot_width_in, " x ", plot_height_in, " inches; ", plot_dpi, " dpi)"
)
