#!/usr/bin/env Rscript
#
# Adult Differential Expression Analysis with DESeq2
#
# Purpose: Identify differentially expressed genes between diapause and non-diapause
#          conditions in adult female mosquitoes
#
# Priority: HIGH - needed for GWAS candidate validation
#
# Input:
#   - results/count_matrices/adults_counts.tsv
#   - results/count_matrices/adults_metadata.csv
#   - data/metadata/gwas_candidates_final_34.txt (if exists)
#
# Output:
#   - results/differential_expression/adults/deseq2_results_all.tsv
#   - results/differential_expression/adults/significant_genes_fdr05.tsv
#   - results/differential_expression/adults/gwas_candidates_expression.tsv
#   - results/differential_expression/adults/statistics_summary.txt
#   - results/differential_expression/adults/volcano_plot.pdf
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
})

cat("===============================================================\n")
cat("Adult Differential Expression Analysis with DESeq2\n")
cat("===============================================================\n\n")

# Set paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
count_dir <- file.path(project_base, "results/count_matrices")
results_dir <- file.path(project_base, "results/differential_expression/adults")
metadata_dir <- file.path(project_base, "data/metadata")

# Create output directory
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("Paths:\n")
cat("  Count directory:", count_dir, "\n")
cat("  Results directory:", results_dir, "\n\n")

# ============================================================================
# 1. Load count matrix
# ============================================================================

cat("1. Loading count matrix...\n")
counts_file <- file.path(count_dir, "adults_counts.tsv")
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)

cat("  Shape:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("  Sample names:", paste(head(colnames(counts)), collapse = ", "), "...\n\n")

# ============================================================================
# 2. Load metadata
# ============================================================================

cat("2. Loading metadata...\n")
metadata_file <- file.path(count_dir, "adults_metadata.csv")
metadata <- read.csv(metadata_file, row.names = 1)

cat("  Samples:", nrow(metadata), "\n")
cat("  Columns:", paste(colnames(metadata), collapse = ", "), "\n")

# Verify condition column exists
if (!"condition" %in% colnames(metadata)) {
  stop("Error: 'condition' column not found in metadata!")
}

# Show condition summary
cat("\n  Condition summary:\n")
print(table(metadata$condition))
cat("\n")

# Verify samples match between counts and metadata
if (!all(rownames(metadata) == colnames(counts))) {
  cat("  WARNING: Sample order mismatch, reordering counts to match metadata\n")
  counts <- counts[, rownames(metadata)]
}

# ============================================================================
# 3. Create DESeq2 dataset
# ============================================================================

cat("3. Creating DESeq2 dataset...\n")

# Convert condition to factor (DESeq2 requirement)
metadata$condition <- factor(metadata$condition)

# Set reference level (non_diapause as reference)
if ("non_diapause" %in% levels(metadata$condition)) {
  metadata$condition <- relevel(metadata$condition, ref = "non_diapause")
  cat("  Reference level: non_diapause\n")
  cat("  Comparison: diapause vs non_diapause\n")
} else {
  cat("  WARNING: 'non_diapause' not found in conditions\n")
  cat("  Using default reference:", levels(metadata$condition)[1], "\n")
}

# Create DESeq2 dataset
# Design: ~ condition (no batch correction - stage-specific analysis)
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),  # DESeq2 requires integer counts
  colData = metadata,
  design = ~ condition
)

cat("  DESeq2 dataset created\n")
cat("  Design formula: ~ condition\n\n")

# ============================================================================
# 4. Pre-filtering
# ============================================================================

cat("4. Pre-filtering low-count genes...\n")

# Keep genes with at least 10 reads total across all samples
keep <- rowSums(counts(dds)) >= 10
dds_filtered <- dds[keep, ]

cat("  Genes before filtering:", nrow(dds), "\n")
cat("  Genes after filtering:", nrow(dds_filtered), "\n")
cat("  Genes removed:", nrow(dds) - nrow(dds_filtered), "\n\n")

# ============================================================================
# 5. Run DESeq2
# ============================================================================

cat("5. Running DESeq2 differential expression...\n")
cat("  This may take several minutes...\n")

system.time({
  dds_deseq <- DESeq(dds_filtered)
})

cat("  DESeq2 analysis complete\n\n")

# ============================================================================
# 6. Extract results
# ============================================================================

cat("6. Extracting results...\n")

# Get results (default: last coefficient in design, which is diapause vs non_diapause)
res <- results(dds_deseq, alpha = 0.05)

cat("  Total genes tested:", nrow(res), "\n")
cat("  Genes with adjusted p-value:", sum(!is.na(res$padj)), "\n\n")

# Summary statistics
cat("  Summary at FDR < 0.05:\n")
print(summary(res))
cat("\n")

# ============================================================================
# 7. Count DEGs at different thresholds
# ============================================================================

cat("7. Counting DEGs at different thresholds...\n\n")

thresholds <- c(0.01, 0.05, 0.1)
deg_counts <- data.frame(
  threshold = thresholds,
  up = NA,
  down = NA,
  total = NA
)

for (i in seq_along(thresholds)) {
  fdr <- thresholds[i]
  sig <- res[which(res$padj < fdr & !is.na(res$padj)), ]
  up <- sum(sig$log2FoldChange > 0)
  down <- sum(sig$log2FoldChange < 0)

  deg_counts$up[i] <- up
  deg_counts$down[i] <- down
  deg_counts$total[i] <- nrow(sig)

  cat(sprintf("  FDR < %.2f: %d DEGs (%d up, %d down)\n",
    fdr, nrow(sig), up, down
  ))
}
cat("\n")

# ============================================================================
# 8. Save results
# ============================================================================

cat("8. Saving results...\n")

# Convert results to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c(
  "gene_id", "baseMean", "log2FoldChange",
  "lfcSE", "stat", "pvalue", "padj"
)]

# Save all results
all_results_file <- file.path(results_dir, "deseq2_results_all.tsv")
write.table(res_df, all_results_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  All results:", all_results_file, "\n")

# Save significant genes (FDR < 0.05)
sig_genes <- res_df[which(res_df$padj < 0.05 & !is.na(res_df$padj)), ]
sig_genes <- sig_genes[order(sig_genes$padj), ]
sig_file <- file.path(results_dir, "significant_genes_fdr05.tsv")
write.table(sig_genes, sig_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Significant genes (FDR < 0.05):", sig_file, "\n")
cat("    Total:", nrow(sig_genes), "\n\n")

# ============================================================================
# 9. Extract GWAS candidate expression
# ============================================================================

cat("9. Extracting GWAS candidate gene expression...\n")

# Check if GWAS candidate file exists
gwas_file <- file.path(metadata_dir, "gwas_candidates_final_34.txt")
if (file.exists(gwas_file)) {
  cat("  Loading GWAS candidates from:", gwas_file, "\n")

  gwas_candidates <- read.table(gwas_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(gwas_candidates) <- "gene_id"

  cat("  Total GWAS candidates:", nrow(gwas_candidates), "\n")

  # Extract expression for candidates
  candidate_expr <- res_df[res_df$gene_id %in% gwas_candidates$gene_id, ]

  if (nrow(candidate_expr) > 0) {
    cat("  Candidates found in expression data:", nrow(candidate_expr), "\n")

    # Count how many are significant
    sig_candidates <- candidate_expr[which(candidate_expr$padj < 0.05 & !is.na(candidate_expr$padj)), ]
    cat("  Candidates with FDR < 0.05:", nrow(sig_candidates), "\n")

    # Save candidate expression
    cand_file <- file.path(results_dir, "gwas_candidates_expression.tsv")
    write.table(candidate_expr, cand_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Saved:", cand_file, "\n")
  } else {
    cat("  WARNING: No GWAS candidates found in expression data\n")
    cat("  This might indicate a gene ID mismatch\n")
  }
} else {
  cat("  GWAS candidate file not found:", gwas_file, "\n")
  cat("  Skipping candidate extraction\n")
}

cat("\n")

# ============================================================================
# 10. Save summary statistics
# ============================================================================

cat("10. Saving summary statistics...\n")

summary_file <- file.path(results_dir, "statistics_summary.txt")
sink(summary_file)

cat("Adult Differential Expression Analysis Summary\n")
cat("===============================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Sample Information:\n")
cat("  Total samples:", ncol(counts), "\n")
cat("  Diapause:", sum(metadata$condition == "diapause"), "\n")
cat("  Non-diapause:", sum(metadata$condition == "non_diapause"), "\n\n")

cat("Gene Filtering:\n")
cat("  Total genes in count matrix:", nrow(dds), "\n")
cat("  Genes passing filter (≥10 reads):", nrow(dds_filtered), "\n")
cat("  Genes removed:", nrow(dds) - nrow(dds_filtered), "\n\n")

cat("Differential Expression Results:\n")
cat("  Total genes tested:", nrow(res), "\n")
cat("  Genes with valid adjusted p-value:", sum(!is.na(res$padj)), "\n\n")

cat("DEGs at Different Thresholds:\n")
print(deg_counts)
cat("\n")

if (exists("candidate_expr") && nrow(candidate_expr) > 0) {
  cat("GWAS Candidate Genes:\n")
  cat("  Total candidates:", nrow(gwas_candidates), "\n")
  cat("  Candidates in expression data:", nrow(candidate_expr), "\n")
  cat("  Candidates with FDR < 0.05:", nrow(sig_candidates), "\n")
  cat("  Enrichment proportion:", round(nrow(sig_candidates) / nrow(candidate_expr), 3), "\n")
}

sink()

cat("  Summary saved:", summary_file, "\n\n")

# ============================================================================
# 11. Create volcano plot
# ============================================================================

cat("11. Creating volcano plot...\n")

# Prepare data for plotting
plot_data <- res_df
plot_data$sig <- ifelse(plot_data$padj < 0.05 & !is.na(plot_data$padj), "FDR < 0.05", "Not significant")
plot_data$log10p <- -log10(plot_data$pvalue)

# Volcano plot
volcano <- ggplot(plot_data, aes(x = log2FoldChange, y = log10p, color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("FDR < 0.05" = "red", "Not significant" = "gray")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(
    title = "Adult: Diapause vs Non-Diapause",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

volcano_file <- file.path(results_dir, "volcano_plot.pdf")
ggsave(volcano_file, volcano, width = 8, height = 6)
cat("  Volcano plot saved:", volcano_file, "\n\n")

# ============================================================================
# Summary
# ============================================================================

cat("===============================================================\n")
cat("Analysis Complete!\n")
cat("===============================================================\n\n")

cat("Output files:\n")
cat("  1. All results:", basename(all_results_file), "\n")
cat("  2. Significant genes:", basename(sig_file), "\n")
cat("  3. Statistics summary:", basename(summary_file), "\n")
cat("  4. Volcano plot:", basename(volcano_file), "\n")
if (exists("cand_file")) {
  cat("  5. GWAS candidates:", basename(cand_file), "\n")
}

cat("\nKey Statistics:\n")
cat("  Genes tested:", nrow(res), "\n")
cat("  DEGs (FDR < 0.05):", nrow(sig_genes), "(", round(100 * nrow(sig_genes) / nrow(res), 2), "%)\n")

if (exists("sig_candidates")) {
  cat("  GWAS candidates DE:", nrow(sig_candidates), "/", nrow(candidate_expr), "\n")
}

cat("\n")
cat("Next steps:\n")
cat("  - Review volcano plot and significant genes\n")
cat("  - Run enrichment test (Fisher's exact) for GWAS candidates\n")
cat("  - Repeat analysis for embryos and larvae\n\n")

cat("===============================================================\n")