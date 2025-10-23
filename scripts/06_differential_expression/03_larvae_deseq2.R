#!/usr/bin/env Rscript
#
# Larvae Differential Expression Analysis with DESeq2
#
# Purpose: Identify differentially expressed genes between diapause and non-diapause
#          conditions in larvae
#
# Input:
#   - results/count_matrices/larvae_counts.tsv
#   - results/count_matrices/larvae_metadata.csv
#
# Output:
#   - results/differential_expression/larvae/deseq2_results_all.tsv
#   - results/differential_expression/larvae/significant_genes_fdr05.tsv
#   - results/differential_expression/larvae/statistics_summary.txt
#
# Author: LCosme
# Date: October 21, 2025

# This is a simplified version of the adult script
# See 01_adults_deseq2.R for detailed comments

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

cat("===============================================================\n")
cat("Larvae Differential Expression Analysis with DESeq2\n")
cat("===============================================================\n\n")

# Paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
count_dir <- file.path(project_base, "results/count_matrices")
results_dir <- file.path(project_base, "results/differential_expression/larvae")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
cat("Loading count matrix and metadata...\n")
counts <- read.delim(file.path(count_dir, "larvae_counts.tsv"), row.names = 1, check.names = FALSE)
metadata <- read.csv(file.path(count_dir, "larvae_metadata.csv"), row.names = 1)

cat("  Samples:", ncol(counts), "\n")
cat("  Genes:", nrow(counts), "\n")
cat("  Conditions:", paste(table(metadata$condition), collapse = " vs "), "\n\n")

# Create DESeq2 dataset
metadata$condition <- factor(metadata$condition)
if ("non_diapause" %in% levels(metadata$condition)) {
  metadata$condition <- relevel(metadata$condition, ref = "non_diapause")
}

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = metadata,
  design = ~ condition
)

# Filter and run
cat("Running DESeq2...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("  Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

cat("\nResults:\n")
print(summary(res))

# Count DEGs
deg_fdr05 <- sum(res$padj < 0.05 & !is.na(res$padj))
cat("\nDEGs at FDR < 0.05:", deg_fdr05, "/", nrow(res), "(", round(100 * deg_fdr05 / nrow(res), 2), "%)\n\n")

# Save
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.table(res_df, file.path(results_dir, "deseq2_results_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sig_genes <- res_df[which(res_df$padj < 0.05 & !is.na(res_df$padj)), ]
sig_genes <- sig_genes[order(sig_genes$padj), ]
write.table(sig_genes, file.path(results_dir, "significant_genes_fdr05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Results saved to:", results_dir, "\n")
cat("===============================================================\n")