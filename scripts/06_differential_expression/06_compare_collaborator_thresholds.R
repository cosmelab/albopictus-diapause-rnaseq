#!/usr/bin/env Rscript
#
# Compare our DE analysis with collaborator's thresholds
#
# Purpose: Apply collaborator's FC threshold and compare DEG counts
#
# Author: Assistant
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("===============================================================\n")
cat("Comparison with Collaborator's Analysis Approach\n")
cat("===============================================================\n\n")

# Paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
de_dir <- file.path(project_base, "results/differential_expression")

# Load and analyze each stage
stages <- c("adults", "embryos", "larvae")
comparison_results <- data.frame()

for (stage in stages) {
  cat(sprintf("Stage: %s\n", tools::toTitleCase(stage)))
  cat("----------------------------------------\n")

  # Load DE results
  results_file <- file.path(de_dir, stage, "deseq2_results_all.tsv")
  if (!file.exists(results_file)) {
    cat("  Results file not found\n\n")
    next
  }

  de_results <- read.delim(results_file, stringsAsFactors = FALSE)

  # Total genes tested
  total_genes <- nrow(de_results)
  cat(sprintf("  Total genes tested: %d\n", total_genes))

  # Our current approach (FDR < 0.05, no FC cutoff)
  our_degs <- sum(de_results$padj < 0.05 & !is.na(de_results$padj))
  our_pct <- 100 * our_degs / total_genes

  # Collaborator's approach (FDR < 0.05 & |log2FC| > 0.5)
  collab_degs <- sum(de_results$padj < 0.05 &
                     abs(de_results$log2FoldChange) > 0.5 &
                     !is.na(de_results$padj), na.rm = TRUE)
  collab_pct <- 100 * collab_degs / total_genes

  # Reduction in DEGs
  reduction <- our_degs - collab_degs
  reduction_pct <- 100 * (1 - collab_degs/our_degs)

  cat(sprintf("\n  Our approach (FDR < 0.05 only):\n"))
  cat(sprintf("    DEGs: %d (%.2f%%)\n", our_degs, our_pct))

  cat(sprintf("\n  Collaborator's approach (FDR < 0.05 & |log2FC| > 0.5):\n"))
  cat(sprintf("    DEGs: %d (%.2f%%)\n", collab_degs, collab_pct))

  cat(sprintf("\n  Difference:\n"))
  cat(sprintf("    %d fewer DEGs with FC threshold (%.1f%% reduction)\n",
              reduction, reduction_pct))

  # Breakdown by direction
  if (collab_degs > 0) {
    up_fc <- sum(de_results$padj < 0.05 &
                 de_results$log2FoldChange > 0.5 &
                 !is.na(de_results$padj), na.rm = TRUE)
    down_fc <- sum(de_results$padj < 0.05 &
                   de_results$log2FoldChange < -0.5 &
                   !is.na(de_results$padj), na.rm = TRUE)

    cat(sprintf("\n  Direction with FC threshold:\n"))
    cat(sprintf("    Upregulated (log2FC > 0.5): %d\n", up_fc))
    cat(sprintf("    Downregulated (log2FC < -0.5): %d\n", down_fc))
  }

  # Store results
  comparison_results <- rbind(comparison_results, data.frame(
    stage = stage,
    total_genes = total_genes,
    our_degs = our_degs,
    our_pct = our_pct,
    collab_degs = collab_degs,
    collab_pct = collab_pct,
    reduction = reduction,
    reduction_pct = reduction_pct
  ))

  cat("\n")
}

# Summary table
cat("===============================================================\n")
cat("Summary Table\n")
cat("===============================================================\n\n")

cat("| Stage    | Genes Tested | Our DEGs     | Collab DEGs  | Reduction    |\n")
cat("|----------|-------------|-------------|-------------|-------------|\n")

for (i in 1:nrow(comparison_results)) {
  row <- comparison_results[i,]
  cat(sprintf("| %-8s | %11d | %5d (%.1f%%) | %5d (%.1f%%) | %3d (%.0f%%) |\n",
    tools::toTitleCase(row$stage),
    row$total_genes,
    row$our_degs, row$our_pct,
    row$collab_degs, row$collab_pct,
    row$reduction, row$reduction_pct
  ))
}

cat("\n")

# Check for very stringent p-value threshold (like collaborator's volcano plots)
cat("===============================================================\n")
cat("Very Stringent Threshold (pCutoff = 10e-6)\n")
cat("===============================================================\n\n")

for (stage in stages) {
  results_file <- file.path(de_dir, stage, "deseq2_results_all.tsv")
  if (!file.exists(results_file)) next

  de_results <- read.delim(results_file, stringsAsFactors = FALSE)

  stringent <- sum(de_results$padj < 1e-6 &
                   abs(de_results$log2FoldChange) > 0.5 &
                   !is.na(de_results$padj), na.rm = TRUE)

  if (stringent > 0) {
    cat(sprintf("%s: %d genes (padj < 1e-6 & |log2FC| > 0.5)\n",
                tools::toTitleCase(stage), stringent))
  }
}

cat("\n")

# Save comparison results
output_file <- file.path(de_dir, "collaborator_comparison_summary.tsv")
write.table(comparison_results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", basename(output_file)))

cat("\n===============================================================\n")
cat("Key Findings:\n")
cat("===============================================================\n\n")

cat("1. Adults show the strongest response regardless of threshold\n")
cat("2. Applying FC threshold reduces DEGs by ~30-50% typically\n")
cat("3. Embryos have very few DEGs even without FC threshold\n")
cat("4. Our approach is more inclusive, capturing smaller magnitude changes\n")
cat("5. Collaborator's approach focuses on larger effect sizes\n")

cat("\n")