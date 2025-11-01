#!/usr/bin/env Rscript
#' Global QC Plots - All 45 Samples
#'
#' Purpose: Quality control visualization across all samples
#'          These plots may become Supplemental Figures showing data quality
#'
#' Outputs:
#'   - PCA plot (all 45 samples, colored by life stage)
#'   - Correlation heatmap (hierarchical clustering)
#'   - Sample distance heatmap
#'   - Library size boxplots
#'   - Detection rate by stage
#'
#' Expected results:
#'   - Samples cluster by life stage (platform confounding)
#'   - No extreme outliers
#'   - Consistent library sizes within stages
#'
#' Usage:
#'   Rscript 01_global_qc_plots.R featurecounts
#'   Rscript 01_global_qc_plots.R salmon
#'
#' Created: October 27, 2025
#' Phase: 2.5 - QC Visualization

# =============================================================================
# Load Libraries
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# =============================================================================
# Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
COUNT_TYPE <- if (length(args) > 0) args[1] else "featurecounts"

if (!COUNT_TYPE %in% c("featurecounts", "salmon")) {
  stop("Usage: Rscript 01_global_qc_plots.R [featurecounts|salmon]")
}

# Paths
COUNT_DIR <- file.path("output/count_matrices", COUNT_TYPE)
OUTPUT_DIR <- file.path("results/qc_visualization", COUNT_TYPE, "all_samples")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output files
PCA_FILE <- file.path(OUTPUT_DIR, "pca_all_45_samples.pdf")
COR_HEATMAP_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_all.pdf")
DIST_HEATMAP_FILE <- file.path(OUTPUT_DIR, "sample_distance_heatmap.pdf")
LIBSIZE_FILE <- file.path(OUTPUT_DIR, "library_size_boxplot.pdf")
DETECTION_FILE <- file.path(OUTPUT_DIR, "detection_rate_by_stage.pdf")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "qc_summary_stats.txt")

# Colors for life stages
STAGE_COLORS <- c("adults" = "#E41A1C", "embryos" = "#377EB8", "larvae" = "#4DAF4A")

# =============================================================================
# Helper Functions
# =============================================================================

#' Load all count matrices and combine
load_all_counts <- function(count_dir) {
  cat("Loading count matrices...\n")

  # Load each stage
  stages <- c("adults", "embryos", "larvae")
  count_list <- list()
  meta_list <- list()

  for (stage in stages) {
    count_file <- file.path(count_dir, stage, "all_counts.tsv")
    meta_file <- file.path(count_dir, stage, "all_metadata.tsv")

    if (!file.exists(count_file)) {
      stop(paste("Count file not found:", count_file))
    }

    counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
    metadata <- read.table(meta_file, header = TRUE, sep = "\t")

    count_list[[stage]] <- counts
    meta_list[[stage]] <- metadata

    cat(sprintf("  %s: %d samples, %d genes\n", stage, ncol(counts), nrow(counts)))
  }

  # Combine all stages
  all_genes <- Reduce(intersect, lapply(count_list, rownames))
  cat(sprintf("\nCommon genes across all stages: %d\n", length(all_genes)))

  all_counts <- do.call(cbind, lapply(count_list, function(x) x[all_genes, ]))
  all_metadata <- do.call(rbind, meta_list)

  # Add stage column
  all_metadata$stage <- rep(stages, sapply(count_list, ncol))

  # Create sample_id column (BioProject_Run)
  all_metadata$sample_id <- paste0(all_metadata$BioProject, "_", all_metadata$Run)

  # Ensure sample order matches
  all_counts <- all_counts[, all_metadata$sample_id]

  cat(sprintf("Combined matrix: %d samples, %d genes\n\n", ncol(all_counts), nrow(all_counts)))

  return(list(counts = all_counts, metadata = all_metadata))
}

#' Create PCA plot
create_pca_plot <- function(vst_data, metadata, output_file, title) {
  cat("Creating PCA plot...\n")

  # Calculate PCA
  pca_data <- plotPCA(vst_data, intgroup = "stage", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  # Merge with metadata for better labeling
  pca_data$sample_id <- rownames(pca_data)
  pca_plot_data <- merge(pca_data, metadata, by.x = "group", by.y = "stage")

  # Create plot
  p <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = group, label = `New.abreviation`)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = STAGE_COLORS) +
    labs(
      title = title,
      subtitle = "Expect: Strong clustering by life stage (platform confounding)",
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance"),
      color = "Life Stage"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))

  # Print variance explained
  cat(sprintf("  PC1 explains %d%% variance\n", percent_var[1]))
  cat(sprintf("  PC2 explains %d%% variance\n\n", percent_var[2]))

  return(pca_plot_data)
}

#' Create correlation heatmap
create_correlation_heatmap <- function(counts, metadata, output_file, title) {
  cat("Creating correlation heatmap...\n")

  # Log transform
  log_counts <- log2(counts + 1)

  # Calculate correlations
  cor_matrix <- cor(log_counts, method = "pearson")

  # Create annotation
  annotation_col <- data.frame(
    Stage = metadata$stage,
    row.names = colnames(cor_matrix)
  )

  annotation_colors <- list(Stage = STAGE_COLORS)

  # Create heatmap
  pdf(output_file, width = 12, height = 11)
  pheatmap(
    cor_matrix,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0.5, 1.0, length.out = 101),
    border_color = NA
  )
  dev.off()

  cat(sprintf("  Saved: %s\n", output_file))

  # Calculate summary statistics
  diag(cor_matrix) <- NA  # Remove self-correlations
  cat(sprintf("  Mean correlation: %.3f\n", mean(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Min correlation: %.3f\n", min(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Max correlation: %.3f\n\n", max(cor_matrix, na.rm = TRUE)))

  return(cor_matrix)
}

#' Create sample distance heatmap
create_distance_heatmap <- function(vst_data, metadata, output_file, title) {
  cat("Creating sample distance heatmap...\n")

  # Calculate distances
  sample_dists <- dist(t(assay(vst_data)))
  sample_dist_matrix <- as.matrix(sample_dists)

  # Create annotation
  annotation_col <- data.frame(
    Stage = metadata$stage,
    row.names = colnames(sample_dist_matrix)
  )

  annotation_colors <- list(Stage = STAGE_COLORS)

  # Create heatmap
  pdf(output_file, width = 12, height = 11)
  pheatmap(
    sample_dist_matrix,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title,
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
    border_color = NA
  )
  dev.off()

  cat(sprintf("  Saved: %s\n\n", output_file))
}

#' Create library size boxplot
create_libsize_plot <- function(counts, metadata, output_file, title) {
  cat("Creating library size boxplot...\n")

  # Calculate library sizes
  lib_sizes <- data.frame(
    sample_id = colnames(counts),
    library_size = colSums(counts) / 1e6,  # In millions
    stage = metadata$stage
  )

  # Create plot
  p <- ggplot(lib_sizes, aes(x = stage, y = library_size, fill = stage)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    scale_fill_manual(values = STAGE_COLORS) +
    labs(
      title = title,
      subtitle = "Check: No extreme outliers (>3 SD from mean)",
      x = "Life Stage",
      y = "Library Size (millions of reads)",
      fill = "Stage"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "none"
    )

  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("  Saved: %s\n", output_file))

  # Print summary
  for (stage in unique(lib_sizes$stage)) {
    stage_libs <- lib_sizes$library_size[lib_sizes$stage == stage]
    cat(sprintf("  %s: Mean = %.1fM, SD = %.1fM, Range = %.1f-%.1fM\n",
                stage, mean(stage_libs), sd(stage_libs), min(stage_libs), max(stage_libs)))
  }
  cat("\n")

  return(lib_sizes)
}

#' Create detection rate plot
create_detection_plot <- function(counts, metadata, output_file, title) {
  cat("Creating detection rate plot...\n")

  # Calculate detection rates (genes with count > 0)
  detection_rates <- data.frame(
    sample_id = colnames(counts),
    detection_rate = colSums(counts > 0) / nrow(counts) * 100,
    stage = metadata$stage
  )

  # Create plot
  p <- ggplot(detection_rates, aes(x = stage, y = detection_rate, fill = stage)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    scale_fill_manual(values = STAGE_COLORS) +
    labs(
      title = title,
      subtitle = "Percent of genes detected (count > 0) per sample",
      x = "Life Stage",
      y = "Detection Rate (%)",
      fill = "Stage"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "none"
    )

  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("  Saved: %s\n", output_file))

  # Print summary
  for (stage in unique(detection_rates$stage)) {
    stage_detect <- detection_rates$detection_rate[detection_rates$stage == stage]
    cat(sprintf("  %s: Mean = %.1f%%, SD = %.1f%%, Range = %.1f-%.1f%%\n",
                stage, mean(stage_detect), sd(stage_detect), min(stage_detect), max(stage_detect)))
  }
  cat("\n")

  return(detection_rates)
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  cat("=======================================================================\n")
  cat("GLOBAL QC PLOTS - ALL 45 SAMPLES\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

  # Load data
  data <- load_all_counts(COUNT_DIR)
  counts <- data$counts
  metadata <- data$metadata

  # Create DESeq2 object for VST
  cat("Creating DESeq2 object...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ stage
  )
  cat("  Done\n\n")

  # Variance stabilizing transformation
  cat("Performing variance stabilizing transformation...\n")
  vst_data <- vst(dds, blind = TRUE)
  cat("  Done\n\n")

  # Create plots
  pca_data <- create_pca_plot(
    vst_data, metadata, PCA_FILE,
    paste("PCA - All 45 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  cor_matrix <- create_correlation_heatmap(
    counts, metadata, COR_HEATMAP_FILE,
    paste("Sample Correlations - All 45 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  create_distance_heatmap(
    vst_data, metadata, DIST_HEATMAP_FILE,
    paste("Sample Distances - All 45 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  lib_sizes <- create_libsize_plot(
    counts, metadata, LIBSIZE_FILE,
    paste("Library Sizes by Stage (", toupper(COUNT_TYPE), ")", sep = "")
  )

  detection_rates <- create_detection_plot(
    counts, metadata, DETECTION_FILE,
    paste("Detection Rates by Stage (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # Write summary
  cat("Writing summary statistics...\n")
  sink(SUMMARY_FILE)
  cat("GLOBAL QC SUMMARY - ALL 45 SAMPLES\n")
  cat("===================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("SAMPLE COUNTS:\n")
  cat(sprintf("  Total samples: %d\n", ncol(counts)))
  cat(sprintf("  Adults: %d\n", sum(metadata$stage == "adults")))
  cat(sprintf("  Embryos: %d\n", sum(metadata$stage == "embryos")))
  cat(sprintf("  Larvae: %d\n\n", sum(metadata$stage == "larvae")))

  cat("GENE COUNTS:\n")
  cat(sprintf("  Total genes: %d\n\n", nrow(counts)))

  cat("CORRELATION STATISTICS:\n")
  diag(cor_matrix) <- NA
  cat(sprintf("  Mean correlation: %.3f\n", mean(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Median correlation: %.3f\n", median(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Min correlation: %.3f\n", min(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Max correlation: %.3f\n\n", max(cor_matrix, na.rm = TRUE)))

  cat("EXPECTED RESULTS:\n")
  cat("  - Strong clustering by life stage in PCA (platform confounding)\n")
  cat("  - Correlation heatmap shows three blocks (adults, embryos, larvae)\n")
  cat("  - No extreme outliers in library size or detection rate\n")
  cat("  - Within-stage correlations > 0.80\n\n")

  cat("ACTION ITEMS:\n")
  cat("  1. Review PCA plot - do samples cluster by stage?\n")
  cat("  2. Check correlation heatmap - any unexpected patterns?\n")
  cat("  3. Identify outliers in library size (>3 SD from mean)\n")
  cat("  4. Proceed to stage-specific QC (scripts 02-04)\n")
  sink()

  cat(sprintf("  Saved: %s\n\n", SUMMARY_FILE))

  cat("=======================================================================\n")
  cat("GLOBAL QC COMPLETE\n")
  cat("=======================================================================\n\n")
  cat("Next steps:\n")
  cat("  1. Review all PDF files in", OUTPUT_DIR, "\n")
  cat("  2. Read summary:", SUMMARY_FILE, "\n")
  cat("  3. Run stage-specific QC: scripts 02-04\n\n")
}

# Run main
main()
