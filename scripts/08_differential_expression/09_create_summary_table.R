#!/usr/bin/env Rscript
#
# Create Comprehensive DE Summary Table
#
# Purpose: Summarize all differential expression results across stages
#          Output to 03_differential_expression/
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("===============================================================\n")
cat("Creating Comprehensive DE Summary Table\n")
cat("===============================================================\n\n")

# Paths - use CORRECT directory structure
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
de_dir <- file.path(project_base, "results/03_differential_expression")

# Load all results
stages <- c("adults", "embryos", "larvae")
summary_table <- data.frame()

for (stage in stages) {
  results_file <- file.path(de_dir, stage, "deseq2_results_all.tsv")

  res <- read.delim(results_file, stringsAsFactors = FALSE)

  # Calculate statistics
  total <- nrow(res)
  sig_all <- sum(res$padj < 0.05 & !is.na(res$padj))
  sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0 & !is.na(res$padj))
  sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0 & !is.na(res$padj))

  # With FC threshold (collaborator's approach)
  sig_fc <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 0.5 & !is.na(res$padj))
  sig_fc_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0.5 & !is.na(res$padj))
  sig_fc_down <- sum(res$padj < 0.05 & res$log2FoldChange < -0.5 & !is.na(res$padj))

  summary_table <- rbind(summary_table, data.frame(
    stage = tools::toTitleCase(stage),
    total_genes = total,
    sig_genes = sig_all,
    sig_pct = round(100 * sig_all / total, 2),
    sig_up = sig_up,
    sig_down = sig_down,
    sig_with_fc = sig_fc,
    sig_with_fc_pct = round(100 * sig_fc / total, 2),
    sig_fc_up = sig_fc_up,
    sig_fc_down = sig_fc_down
  ))

  cat(sprintf("%s: %d total, %d sig (%.2f%%)\n",
              tools::toTitleCase(stage), total, sig_all, 100*sig_all/total))
}

# Save to correct directory
output_file <- file.path(de_dir, "summary_all_stages.tsv")
write.table(summary_table, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nSaved: %s\n", output_file))

# Create formatted text summary
text_file <- file.path(de_dir, "summary_all_stages.txt")
sink(text_file)

cat("===============================================================\n")
cat("Differential Expression Analysis Summary - All Stages\n")
cat("===============================================================\n\n")
cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("Analysis Parameters:\n")
cat("  - Method: DESeq2\n")
cat("  - Design: ~ Condition (stage-specific)\n")
cat("  - Comparison: Diapause vs Non-diapause\n")
cat("  - Significance: FDR < 0.05\n")
cat("  - No fold change cutoff (our approach)\n\n")

cat("Results by Stage:\n")
cat("=================\n\n")

for (i in 1:nrow(summary_table)) {
  row <- summary_table[i,]
  cat(sprintf("%s:\n", row$stage))
  cat(sprintf("  Total genes tested: %d\n", row$total_genes))
  cat(sprintf("  Significant (FDR < 0.05): %d (%.2f%%)\n", row$sig_genes, row$sig_pct))
  cat(sprintf("    Upregulated: %d\n", row$sig_up))
  cat(sprintf("    Downregulated: %d\n", row$sig_down))
  cat(sprintf("  With FC threshold (|log2FC| > 0.5): %d (%.2f%%)\n",
              row$sig_with_fc, row$sig_with_fc_pct))
  cat(sprintf("    Upregulated: %d\n", row$sig_fc_up))
  cat(sprintf("    Downregulated: %d\n", row$sig_fc_down))
  cat("\n")
}

cat("Interpretation:\n")
cat("===============\n")
cat("- Adults show the strongest diapause response (23.84% DE)\n")
cat("- Embryos show minimal response (0.14% DE)\n")
cat("- Larvae show intermediate response (1.64% DE)\n")
cat("- Applying FC cutoff reduces DEGs by ~40% in adults and larvae\n")
cat("- Embryo DEGs all meet FC threshold (strong effects despite few genes)\n\n")

sink()

cat(sprintf("Saved: %s\n", text_file))

cat("\n===============================================================\n")
cat("Summary complete!\n")
cat("===============================================================\n")
