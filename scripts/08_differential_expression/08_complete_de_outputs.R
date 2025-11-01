#!/usr/bin/env Rscript
#
# Complete Differential Expression Analysis Outputs
#
# Purpose: Create all missing outputs (volcano plots, MA plots, statistics)
#          for embryos and larvae to match adults
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

cat("===============================================================\n")
cat("Completing Differential Expression Analysis Outputs\n")
cat("===============================================================\n\n")

# Paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
de_dir <- file.path(project_base, "results/differential_expression")

# Process each stage
stages <- c("embryos", "larvae")

for (stage in stages) {
  cat(sprintf("\n%s\n", toupper(stage)))
  cat("-----------------------------------------------------------\n")

  stage_dir <- file.path(de_dir, stage)
  results_file <- file.path(stage_dir, "deseq2_results_all.tsv")

  if (!file.exists(results_file)) {
    cat("ERROR: Results file not found!\n")
    next
  }

  # Load results
  res <- read.delim(results_file, stringsAsFactors = FALSE)

  # Calculate statistics
  total_genes <- nrow(res)
  sig_genes <- sum(res$padj < 0.05 & !is.na(res$padj))
  sig_pct <- 100 * sig_genes / total_genes

  sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0 & !is.na(res$padj))
  sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0 & !is.na(res$padj))

  # Print statistics
  cat(sprintf("Total genes tested: %d\n", total_genes))
  cat(sprintf("Significant (FDR < 0.05): %d (%.2f%%)\n", sig_genes, sig_pct))
  cat(sprintf("  Upregulated: %d\n", sig_up))
  cat(sprintf("  Downregulated: %d\n", sig_down))

  # Save statistics summary
  stats_file <- file.path(stage_dir, "statistics_summary.txt")
  sink(stats_file)
  cat(sprintf("Differential Expression Analysis - %s\n", tools::toTitleCase(stage)))
  cat("=" , rep("=", 50), "\n", sep="")
  cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Total genes tested: %d\n", total_genes))
  cat(sprintf("Significant genes (FDR < 0.05): %d (%.2f%%)\n", sig_genes, sig_pct))
  cat(sprintf("  Upregulated (log2FC > 0): %d\n", sig_up))
  cat(sprintf("  Downregulated (log2FC < 0): %d\n", sig_down))
  cat("\nSignificance threshold: FDR < 0.05\n")
  cat("No fold change cutoff applied\n")
  sink()

  cat(sprintf("  Saved: statistics_summary.txt\n"))

  # Create volcano plot
  res$significance <- ifelse(res$padj < 0.05 & !is.na(res$padj),
                             ifelse(res$log2FoldChange > 0, "Up", "Down"),
                             "NS")

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = sprintf("Volcano Plot - %s", tools::toTitleCase(stage)),
         subtitle = sprintf("Diapause vs Non-diapause (%d DEGs)", sig_genes),
         x = "log2 Fold Change",
         y = "-log10(adjusted p-value)") +
    theme_bw() +
    theme(legend.position = "bottom")

  volcano_file <- file.path(stage_dir, "volcano_plot.pdf")
  ggsave(volcano_file, p, width = 8, height = 6)
  cat(sprintf("  Saved: volcano_plot.pdf\n"))

  # Create MA plot
  res$baseMean_log10 <- log10(res$baseMean + 1)

  p_ma <- ggplot(res, aes(x = baseMean_log10, y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = sprintf("MA Plot - %s", tools::toTitleCase(stage)),
         subtitle = sprintf("Diapause vs Non-diapause (%d DEGs)", sig_genes),
         x = "log10(Mean Expression)",
         y = "log2 Fold Change") +
    theme_bw() +
    theme(legend.position = "bottom")

  ma_file <- file.path(stage_dir, "ma_plot.pdf")
  ggsave(ma_file, p_ma, width = 8, height = 6)
  cat(sprintf("  Saved: ma_plot.pdf\n"))
}

# Create MA plot for adults too
cat("\nADULTS\n")
cat("-----------------------------------------------------------\n")
stage_dir <- file.path(de_dir, "adults")
results_file <- file.path(stage_dir, "deseq2_results_all.tsv")

if (file.exists(results_file)) {
  res <- read.delim(results_file, stringsAsFactors = FALSE)

  res$significance <- ifelse(res$padj < 0.05 & !is.na(res$padj),
                             ifelse(res$log2FoldChange > 0, "Up", "Down"),
                             "NS")
  res$baseMean_log10 <- log10(res$baseMean + 1)

  sig_genes <- sum(res$padj < 0.05 & !is.na(res$padj))

  p_ma <- ggplot(res, aes(x = baseMean_log10, y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "MA Plot - Adults",
         subtitle = sprintf("Diapause vs Non-diapause (%d DEGs)", sig_genes),
         x = "log10(Mean Expression)",
         y = "log2 Fold Change") +
    theme_bw() +
    theme(legend.position = "bottom")

  ma_file <- file.path(stage_dir, "ma_plot.pdf")
  ggsave(ma_file, p_ma, width = 8, height = 6)
  cat("  Saved: ma_plot.pdf\n")
}

cat("\n===============================================================\n")
cat("All outputs completed!\n")
cat("===============================================================\n")
