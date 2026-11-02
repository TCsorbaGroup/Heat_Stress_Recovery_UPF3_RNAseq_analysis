library(DESeq2)

counts <- read.delim(
  "Arabidopsis_gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

count_matrix <- counts[, 7:24]
rownames(count_matrix) <- counts$Geneid

count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

sample_info <- data.frame(
  genotype = c(
    rep("Col", 3),
    rep("UPF3", 3),
    rep("Col", 3),
    rep("UPF3", 3),
    rep("Col", 3),
    rep("UPF3", 3)
  ),
  treatment = c(
    rep("NT", 6),
    rep("LAT", 6),
    rep("REC", 6)
  )
)

rownames(sample_info) <- colnames(count_matrix)

sample_info$genotype <- factor(
  sample_info$genotype,
  levels = c("Col", "UPF3")
)

sample_info$treatment <- factor(
  sample_info$treatment,
  levels = c("NT", "LAT", "REC")
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ genotype * treatment
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

res_NT <- results(
  dds,
  name = "genotype_UPF3_vs_Col",
  alpha = 0.05
)

res_LAT <- results(
  dds,
  contrast = list(
    c("genotype_UPF3_vs_Col",
      "genotypeUPF3.treatmentLAT")
  ),
  alpha = 0.05
  
)

res_REC <- results(
  dds,
  contrast = list(
    c("genotype_UPF3_vs_Col",
      "genotypeUPF3.treatmentREC")
  ),
  alpha = 0.05
)

res_NT_full <- as.data.frame(res_NT)
res_NT_full$gene_id <- rownames(res_NT_full)

res_NT_full <- res_NT_full[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.csv(
  res_NT_full,
  "DESeq2_UPF3_vs_Col_NT_full.csv",
  row.names = FALSE
)

res_LAT_full <- as.data.frame(res_LAT)
res_LAT_full$gene_id <- rownames(res_LAT_full)

res_LAT_full <- res_LAT_full[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.csv(
  res_LAT_full,
  "DESeq2_UPF3_vs_Col_LAT_full.csv",
  row.names = FALSE
)

res_REC_full <- as.data.frame(res_REC)
res_REC_full$gene_id <- rownames(res_REC_full)

res_REC_full <- res_REC_full[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.csv(
  res_REC_full,
  "DESeq2_UPF3_vs_Col_REC_full.csv",
  row.names = FALSE
)

deg_NT <- as.data.frame(res_NT)
deg_NT$gene_id <- rownames(deg_NT)

deg_NT <- deg_NT[
  !is.na(deg_NT$padj) &
    deg_NT$padj < 0.05 &
    abs(deg_NT$log2FoldChange) >= 1,
]

deg_NT$direction <- ifelse(
  deg_NT$log2FoldChange > 0,
  "UPF3_up",
  "Col_up"
)

deg_NT <- deg_NT[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "direction"
)]

deg_LAT <- as.data.frame(res_LAT)
deg_LAT$gene_id <- rownames(deg_LAT)

deg_LAT <- deg_LAT[
  !is.na(deg_LAT$padj) &
    deg_LAT$padj < 0.05 &
    abs(deg_LAT$log2FoldChange) >= 1,
]

deg_LAT$direction <- ifelse(
  deg_LAT$log2FoldChange > 0,
  "UPF3_up",
  "Col_up"
)

deg_LAT <- deg_LAT[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "direction"
)]

deg_REC <- as.data.frame(res_REC)
deg_REC$gene_id <- rownames(deg_REC)

deg_REC <- deg_REC[
  !is.na(deg_REC$padj) &
    deg_REC$padj < 0.05 &
    abs(deg_REC$log2FoldChange) >= 1,
]

deg_REC$direction <- ifelse(
  deg_REC$log2FoldChange > 0,
  "UPF3_up",
  "Col_up"
)

deg_REC <- deg_REC[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "direction"
)]

write.csv(deg_NT, "DEG_UPF3_vs_Col_NT.csv", row.names = FALSE)
write.csv(deg_LAT, "DEG_UPF3_vs_Col_LAT.csv", row.names = FALSE)
write.csv(deg_REC, "DEG_UPF3_vs_Col_REC.csv", row.names = FALSE)



