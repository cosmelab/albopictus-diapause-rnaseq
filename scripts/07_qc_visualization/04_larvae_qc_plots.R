#!/usr/bin/env Rscript
#' Larvae QC Plots - With Degradation Warning
#'
#' Purpose: Stage-specific QC for larvae (17 samples)
#'          CRITICAL: LD samples at 21d/40d are degraded (dead larvae)
#'          Publication-quality figures with appropriate warnings
#'
#' Experimental Design:
#'   - 2×3 factorial: Photoperiod (SD/LD) × Time (11d/21d/40d)
#'   - 2-3 replicates per group (unbalanced)
#'   - VALID comparison: 11d only (SD vs LD, before LD death)
#'   - SD time series: 11d → 21d → 40d (diapause deepening)
#'   - LD 21d/40d: Biological artifacts (degraded tissue)
#'
#' Outputs:
#'   - Valid comparison PCA (11d only)
#'   - SD time series PCA (diapause progression)
#'   - LD degradation warning plot
#'   - Global larvae PCA (showing degradation)
#'   - Correlation heatmaps
#'
#' Expected Results:
#'   - 11d: Clear SD vs LD separation
#'   - SD time series: Gradual progression (not abrupt)
#'   - LD 21d/40d: Visible outliers (expected degradation)
#'
#' Usage:
#'   Rscript 04_larvae_qc_plots.R featurecounts
#'   Rscript 04_larvae_qc_plots.R salmon
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
  library(gridExtra)
})

# =============================================================================
# Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
COUNT_TYPE <- if (length(args) > 0) args[1] else "featurecounts"

if (!COUNT_TYPE %in% c("featurecounts", "salmon")) {
  stop("Usage: Rscript 04_larvae_qc_plots.R [featurecounts|salmon]")
}

# Paths
COUNT_DIR <- file.path("output/count_matrices", COUNT_TYPE, "larvae")
OUTPUT_DIR <- file.path("results/qc_visualization", COUNT_TYPE, "larvae")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output files
PCA_11D_FILE <- file.path(OUTPUT_DIR, "pca_11d_valid_comparison.pdf")
PCA_SD_TIMESERIES_FILE <- file.path(OUTPUT_DIR, "pca_SD_timeseries.pdf")
PCA_LD_WARNING_FILE <- file.path(OUTPUT_DIR, "pca_LD_degradation_warning.pdf")
PCA_ALL_FILE <- file.path(OUTPUT_DIR, "pca_all_larvae.pdf")
COR_11D_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_11d.pdf")
COR_SD_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_SD_timeseries.pdf")
COR_ALL_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_larvae_with_warning.pdf")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "qc_summary_larvae.txt")

# Colors
PHOTOPERIOD_COLORS <- c("SD" = "#D95F02", "LD" = "#1B9E77")
TIMEPOINT_COLORS <- c("11d" = "#66C2A5", "21d" = "#FC8D62", "40d" = "#8DA0CB")
DEGRADATION_COLORS <- c("Valid" = "#1B9E77", "Degraded" = "#E41A1C")

# =============================================================================
# Helper Functions
# =============================================================================

#' Load larvae count data
load_larvae_data <- function(count_dir) {
  cat("Loading larvae count matrices...\n")

  # Load all larvae
  all_counts <- read.table(
    file.path(count_dir, "all_counts.tsv"),
    header = TRUE, row.names = 1, sep = "\t"
  )
  all_meta <- read.table(
    file.path(count_dir, "all_metadata.tsv"),
    header = TRUE, sep = "\t"
  )

  # Extract photoperiod and timepoint from sample names
  all_meta$photoperiod <- ifelse(grepl("^SD", all_meta$New.abreviation), "SD", "LD")
  all_meta$timepoint <- ifelse(grepl("11d", all_meta$New.abreviation), "11d",
                                ifelse(grepl("21d", all_meta$New.abreviation), "21d", "40d"))
  all_meta$sample_id <- paste0(all_meta$BioProject, "_", all_meta$Run)

  # Mark degraded samples
  all_meta$quality <- ifelse(all_meta$photoperiod == "LD" & all_meta$timepoint %in% c("21d", "40d"),
                              "Degraded", "Valid")

  # Ensure sample order matches
  all_counts <- all_counts[, all_meta$sample_id]

  cat(sprintf("  All larvae: %d samples, %d genes\n", ncol(all_counts), nrow(all_counts)))
  cat(sprintf("    SD: %d, LD: %d\n", sum(all_meta$photoperiod == "SD"), sum(all_meta$photoperiod == "LD")))
  cat(sprintf("    11d: %d, 21d: %d, 40d: %d\n",
              sum(all_meta$timepoint == "11d"),
              sum(all_meta$timepoint == "21d"),
              sum(all_meta$timepoint == "40d")))
  cat(sprintf("    Valid: %d, Degraded: %d\n\n",
              sum(all_meta$quality == "Valid"),
              sum(all_meta$quality == "Degraded")))

  return(list(counts = all_counts, metadata = all_meta))
}

#' Create 11d comparison PCA (valid comparison)
create_11d_pca <- function(counts, metadata, output_file, title) {
  cat("Creating 11d PCA (valid SD vs LD comparison)...\n")

  # Subset to 11d only
  samples_11d <- metadata$timepoint == "11d"
  counts_11d <- counts[, samples_11d]
  meta_11d <- metadata[samples_11d, ]

  cat(sprintf("  Samples: %d (SD: %d, LD: %d)\n",
              nrow(meta_11d),
              sum(meta_11d$photoperiod == "SD"),
              sum(meta_11d$photoperiod == "LD")))

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_11d,
    colData = meta_11d,
    design = ~ photoperiod
  )

  # VST
  vst_data <- vst(dds, blind = TRUE)

  # PCA
  pcaData <- plotPCA(vst_data, intgroup = "photoperiod", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, meta_11d, by = "sample_id")

  # Plot
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = photoperiod,
                              label = New.abreviation)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -1, size = 3) +
    scale_color_manual(values = PHOTOPERIOD_COLORS) +
    labs(
      title = title,
      subtitle = "VALID comparison - LD larvae not yet degraded at 11d",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Photoperiod"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "darkgreen"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  PC1: %d%%, PC2: %d%%\n\n", percentVar[1], percentVar[2]))

  return(list(vst = vst_data, pca = plot_data))
}

#' Create SD time series PCA
create_sd_timeseries_pca <- function(counts, metadata, output_file, title) {
  cat("Creating SD time series PCA (diapause progression)...\n")

  # Subset to SD only
  samples_sd <- metadata$photoperiod == "SD"
  counts_sd <- counts[, samples_sd]
  meta_sd <- metadata[samples_sd, ]

  cat(sprintf("  Samples: %d (11d: %d, 21d: %d, 40d: %d)\n",
              nrow(meta_sd),
              sum(meta_sd$timepoint == "11d"),
              sum(meta_sd$timepoint == "21d"),
              sum(meta_sd$timepoint == "40d")))

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sd,
    colData = meta_sd,
    design = ~ timepoint
  )

  # VST
  vst_data <- vst(dds, blind = TRUE)

  # PCA
  pcaData <- plotPCA(vst_data, intgroup = "timepoint", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, meta_sd, by = "sample_id")

  # Plot with trajectory
  centroids <- plot_data %>%
    group_by(timepoint) %>%
    summarize(PC1_mean = mean(PC1), PC2_mean = mean(PC2), .groups = "drop")

  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = timepoint)) +
    geom_point(size = 4, alpha = 0.6) +
    geom_point(data = centroids, aes(x = PC1_mean, y = PC2_mean), size = 6, alpha = 0.9) +
    geom_path(data = centroids, aes(x = PC1_mean, y = PC2_mean, group = 1),
              arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
              size = 1.2, alpha = 0.7) +
    scale_color_manual(values = TIMEPOINT_COLORS) +
    labs(
      title = title,
      subtitle = "Diapause deepening trajectory: 11d → 21d → 40d (SD only)",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Timepoint"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  PC1: %d%%, PC2: %d%%\n", percentVar[1], percentVar[2]))
  cat("  Expect: Gradual progression, not abrupt jumps\n\n")

  return(list(vst = vst_data, pca = plot_data))
}

#' Create LD degradation warning plot
create_ld_warning_pca <- function(counts, metadata, output_file, title) {
  cat("Creating LD degradation warning plot...\n")

  # Subset to LD only
  samples_ld <- metadata$photoperiod == "LD"
  counts_ld <- counts[, samples_ld]
  meta_ld <- metadata[samples_ld, ]

  cat(sprintf("  LD samples: %d (11d: %d, 21d: %d, 40d: %d)\n",
              nrow(meta_ld),
              sum(meta_ld$timepoint == "11d"),
              sum(meta_ld$timepoint == "21d"),
              sum(meta_ld$timepoint == "40d")))

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_ld,
    colData = meta_ld,
    design = ~ timepoint
  )

  # VST
  vst_data <- vst(dds, blind = TRUE)

  # PCA
  pcaData <- plotPCA(vst_data, intgroup = c("timepoint", "quality"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, meta_ld, by = "sample_id")

  # Plot with degradation color coding
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = quality, shape = timepoint)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_color_manual(values = DEGRADATION_COLORS) +
    scale_shape_manual(values = c("11d" = 16, "21d" = 17, "40d" = 15)) +
    annotate("text", x = min(plot_data$PC1), y = max(plot_data$PC2),
             label = "WARNING: Red points are degraded samples\n(dead larvae, not valid for analysis)",
             hjust = 0, vjust = 1, size = 4, color = "#E41A1C", fontface = "bold") +
    labs(
      title = title,
      subtitle = "LD larvae hatch ~3-4 days, then die without food/water",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Sample Quality",
      shape = "Timepoint"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "#E41A1C"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  PC1: %d%%, PC2: %d%%\n", percentVar[1], percentVar[2]))
  cat("  WARNING: LD 21d/40d samples are biological artifacts\n\n")
}

#' Create global larvae PCA showing degradation
create_global_pca <- function(vst_data, metadata, output_file, title) {
  cat("Creating global larvae PCA...\n")

  # PCA
  pcaData <- plotPCA(vst_data, intgroup = c("photoperiod", "quality"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, metadata, by = "sample_id")

  # Plot
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = quality, shape = photoperiod)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = DEGRADATION_COLORS) +
    scale_shape_manual(values = c("SD" = 16, "LD" = 17)) +
    labs(
      title = title,
      subtitle = "Degraded LD samples (21d/40d) separate from valid samples",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Sample Quality",
      shape = "Photoperiod"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  PC1: %d%%, PC2: %d%%\n\n", percentVar[1], percentVar[2]))
}

#' Create correlation heatmap
create_cor_heatmap <- function(counts, metadata, ann_cols, output_file, title, add_warning = FALSE) {
  cat(sprintf("Creating correlation heatmap: %s...\n", basename(output_file)))

  # Log transform
  log_counts <- log2(counts + 1)

  # Calculate correlations
  cor_matrix <- cor(log_counts, method = "pearson")

  # Create annotation
  annotation_df <- as.data.frame(metadata[, ann_cols, drop = FALSE])
  rownames(annotation_df) <- colnames(cor_matrix)

  # Colors
  ann_colors <- list()
  if ("photoperiod" %in% ann_cols) ann_colors$photoperiod <- PHOTOPERIOD_COLORS
  if ("timepoint" %in% ann_cols) ann_colors$timepoint <- TIMEPOINT_COLORS
  if ("quality" %in% ann_cols) ann_colors$quality <- DEGRADATION_COLORS

  # Create heatmap
  pdf(output_file, width = 11, height = 10)
  pheatmap(
    cor_matrix,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0.5, 1.0, length.out = 101),
    border_color = NA
  )
  dev.off()

  cat(sprintf("  Saved: %s\n", output_file))

  # Summary stats
  diag(cor_matrix) <- NA
  cat(sprintf("  Mean correlation: %.3f\n", mean(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Min correlation: %.3f\n", min(cor_matrix, na.rm = TRUE)))

  if ("quality" %in% metadata) {
    # Check degraded sample correlations
    degraded_samples <- metadata$sample_id[metadata$quality == "Degraded"]
    if (length(degraded_samples) > 0) {
      valid_samples <- metadata$sample_id[metadata$quality == "Valid"]
      cross_cor <- cor_matrix[valid_samples, degraded_samples]
      cat(sprintf("  Valid-Degraded correlation: %.3f (expect lower)\n", mean(cross_cor, na.rm = TRUE)))
    }
  }
  cat("\n")

  return(cor_matrix)
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  cat("=======================================================================\n")
  cat("LARVAE QC PLOTS - WITH DEGRADATION WARNING\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

  # Load data
  data <- load_larvae_data(COUNT_DIR)
  counts <- data$counts
  metadata <- data$metadata

  # === VALID 11D COMPARISON ===
  cat("--- 11d Comparison (VALID SD vs LD) ---\n")

  data_11d <- create_11d_pca(
    counts, metadata,
    PCA_11D_FILE,
    paste("Larvae 11d - SD vs LD [VALID] (", toupper(COUNT_TYPE), ")", sep = "")
  )

  samples_11d <- metadata$timepoint == "11d"
  cor_11d <- create_cor_heatmap(
    counts[, samples_11d], metadata[samples_11d, ], c("photoperiod"),
    COR_11D_FILE,
    paste("Larvae 11d Correlation Heatmap [VALID] (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === SD TIME SERIES ===
  cat("--- SD Time Series (Diapause Progression) ---\n")

  data_sd <- create_sd_timeseries_pca(
    counts, metadata,
    PCA_SD_TIMESERIES_FILE,
    paste("Larvae SD Time Series - Diapause Deepening (", toupper(COUNT_TYPE), ")", sep = "")
  )

  samples_sd <- metadata$photoperiod == "SD"
  cor_sd <- create_cor_heatmap(
    counts[, samples_sd], metadata[samples_sd, ], c("timepoint"),
    COR_SD_FILE,
    paste("Larvae SD Time Series Correlation Heatmap (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === LD DEGRADATION WARNING ===
  cat("--- LD Samples (Degradation Warning) ---\n")

  create_ld_warning_pca(
    counts, metadata,
    PCA_LD_WARNING_FILE,
    paste("Larvae LD Samples - DEGRADATION WARNING (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === GLOBAL LARVAE ===
  cat("--- Global Larvae (All Samples) ---\n")

  # Create DESeq2 for all
  dds_all <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ quality
  )
  vst_all <- vst(dds_all, blind = TRUE)

  create_global_pca(
    vst_all, metadata,
    PCA_ALL_FILE,
    paste("Larvae PCA - All 17 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  cor_all <- create_cor_heatmap(
    counts, metadata, c("photoperiod", "quality"),
    COR_ALL_FILE,
    paste("Larvae Correlation Heatmap - All Samples (", toupper(COUNT_TYPE), ")", sep = ""),
    add_warning = TRUE
  )

  # === SUMMARY ===
  cat("Writing summary statistics...\n")
  sink(SUMMARY_FILE)
  cat("LARVAE QC SUMMARY\n")
  cat("=================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("SAMPLE DISTRIBUTION:\n")
  cat(sprintf("  Total: %d\n", ncol(counts)))
  cat(sprintf("  SD_11d: %d\n", sum(metadata$photoperiod == "SD" & metadata$timepoint == "11d")))
  cat(sprintf("  LD_11d: %d\n", sum(metadata$photoperiod == "LD" & metadata$timepoint == "11d")))
  cat(sprintf("  SD_21d: %d\n", sum(metadata$photoperiod == "SD" & metadata$timepoint == "21d")))
  cat(sprintf("  LD_21d: %d (DEGRADED)\n", sum(metadata$photoperiod == "LD" & metadata$timepoint == "21d")))
  cat(sprintf("  SD_40d: %d\n", sum(metadata$photoperiod == "SD" & metadata$timepoint == "40d")))
  cat(sprintf("  LD_40d: %d (DEGRADED)\n\n", sum(metadata$photoperiod == "LD" & metadata$timepoint == "40d")))

  cat("SAMPLE QUALITY:\n")
  cat(sprintf("  Valid: %d\n", sum(metadata$quality == "Valid")))
  cat(sprintf("  Degraded: %d (LD 21d/40d)\n\n", sum(metadata$quality == "Degraded")))

  cat("CORRELATION SUMMARY:\n")
  cat("  11d comparison (valid):\n")
  diag(cor_11d) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_11d, na.rm = TRUE), min(cor_11d, na.rm = TRUE)))

  cat("  SD time series:\n")
  diag(cor_sd) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_sd, na.rm = TRUE), min(cor_sd, na.rm = TRUE)))

  cat("  All larvae:\n")
  diag(cor_all) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n\n", mean(cor_all, na.rm = TRUE), min(cor_all, na.rm = TRUE)))

  cat("EXPECTED RESULTS:\n")
  cat("  - 11d: Clear SD vs LD separation\n")
  cat("  - SD time series: Gradual progression (not abrupt)\n")
  cat("  - LD 21d/40d: Visible outliers (expected degradation)\n")
  cat("  - Degraded samples have lower correlation with valid samples\n\n")

  cat("BIOLOGICAL INTERPRETATION:\n")
  cat("  - 11d: Early diapause maintenance (VALID comparison)\n")
  cat("  - SD time series: Diapause deepening trajectory\n")
  cat("  - LD 21d/40d: Dead larvae, degraded tissue (DO NOT analyze as normal)\n\n")

  cat("CRITICAL WARNINGS:\n")
  cat("  - DO NOT analyze LD 21d/40d vs SD as 'diapause vs non-diapause'\n")
  cat("  - LD larvae hatch ~3-4 days, then die without food/water\n")
  cat("  - Expression at 21d/40d reflects decomposition, not development\n")
  cat("  - ONLY VALID cross-photoperiod comparison: 11d\n\n")

  cat("QC CHECKLIST:\n")
  cat("  [ ] 11d shows SD vs LD separation\n")
  cat("  [ ] SD time series shows gradual changes\n")
  cat("  [ ] LD 21d/40d are visible outliers (expected)\n")
  cat("  [ ] Within-group correlations > 0.9 (for valid samples)\n")
  cat("  [ ] Proceed to MultiQC integration (script 05)\n")

  sink()
  cat(sprintf("  Saved: %s\n\n", SUMMARY_FILE))

  cat("=======================================================================\n")
  cat("LARVAE QC COMPLETE\n")
  cat("=======================================================================\n\n")
  cat("Review:\n")
  cat("  - Valid comparison (11d): ", PCA_11D_FILE, "\n", sep = "")
  cat("  - Degradation warning: ", PCA_LD_WARNING_FILE, "\n", sep = "")
  cat("  - All plots in:", OUTPUT_DIR, "\n")
  cat("  - Summary:", SUMMARY_FILE, "\n\n")
  cat("Next: Run script 05 (MultiQC integration)\n\n")
}

# Run main
main()
