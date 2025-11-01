#!/usr/bin/env Rscript
#
# Create Stage-Specific PCA Plots
#
# Purpose: Generate PCA plots for each developmental stage to assess
#          sample clustering and validate our stage-specific analysis approach
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(tidyverse)
})

cat("===============================================================\n")
cat("Creating Stage-Specific PCA Plots\n")
cat("===============================================================\n\n")

# Paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
count_dir <- file.path(project_base, "results/count_matrices")
results_dir <- file.path(project_base, "results/qc_analysis")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Process each stage
stages <- c("adults", "embryos", "larvae")

for (stage in stages) {
  cat(sprintf("\nProcessing %s...\n", stage))

  # Load data
  counts_file <- file.path(count_dir, paste0(stage, "_counts.tsv"))
  metadata_file <- file.path(count_dir, paste0(stage, "_metadata.csv"))

  if (!file.exists(counts_file) || !file.exists(metadata_file)) {
    cat("  Files not found, skipping\n")
    next
  }

  counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  metadata <- read.csv(metadata_file, row.names = 1)

  cat(sprintf("  Samples: %d\n", ncol(counts)))
  cat(sprintf("  Genes: %d\n", nrow(counts)))

  # Create DESeq2 dataset for normalization
  metadata$condition <- factor(metadata$condition)
  if ("non_diapause" %in% levels(metadata$condition)) {
    metadata$condition <- relevel(metadata$condition, ref = "non_diapause")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~ condition
  )

  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  # Variance stabilizing transformation
  vst <- vst(dds, blind = TRUE)

  # PCA
  pca_data <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  # Create plot
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    labs(
      title = paste(tools::toTitleCase(stage), "- PCA Plot"),
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = "Condition"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )

  # Save
  out_file <- file.path(results_dir, paste0(stage, "_pca.pdf"))
  ggsave(out_file, p, width = 7, height = 5)
  cat(sprintf("  Saved: %s\n", basename(out_file)))
}

cat("\n===============================================================\n")
cat("PCA plots complete!\n")
cat("Output directory:", results_dir, "\n")
cat("===============================================================\n")