#!/usr/bin/env Rscript
#' Embryos QC Plots - Time Series Analysis
#'
#' Purpose: Stage-specific QC for embryos (12 samples) with time series
#'          Publication-quality figures showing diapause programming dynamics
#'
#' Experimental Design:
#'   - 2×2 factorial: Photoperiod (SD/LD) × Time (72h/135h)
#'   - 3 replicates per group
#'   - 72h: Early development (before diapause entry)
#'   - 135h: Diapause entry (~120h in SD)
#'   - Question: Does photoperiod effect INCREASE over time?
#'
#' Outputs:
#'   - Global embryos PCA (all 12)
#'   - 72h-specific PCA
#'   - 135h-specific PCA
#'   - Interaction plot (photoperiod × time)
#'   - Time series trajectory plot
#'   - Correlation heatmaps
#'
#' Expected Results:
#'   - Photoperiod separation at both timepoints
#'   - Possible time effect (developmental progression)
#'   - Smooth trajectory (not abrupt changes)
#'   - Within-group correlations > 0.9
#'
#' Usage:
#'   Rscript 03_embryos_qc_plots.R featurecounts
#'   Rscript 03_embryos_qc_plots.R salmon
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
  stop("Usage: Rscript 03_embryos_qc_plots.R [featurecounts|salmon]")
}

# Paths
COUNT_DIR <- file.path("output/count_matrices", COUNT_TYPE, "embryos")
OUTPUT_DIR <- file.path("results/qc_visualization", COUNT_TYPE, "embryos")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output files
PCA_ALL_FILE <- file.path(OUTPUT_DIR, "pca_all_embryos.pdf")
PCA_72H_FILE <- file.path(OUTPUT_DIR, "pca_72h_only.pdf")
PCA_135H_FILE <- file.path(OUTPUT_DIR, "pca_135h_only.pdf")
PCA_PHOTOPERIOD_FILE <- file.path(OUTPUT_DIR, "pca_by_photoperiod.pdf")
PCA_TIMEPOINT_FILE <- file.path(OUTPUT_DIR, "pca_by_timepoint.pdf")
PCA_INTERACTION_FILE <- file.path(OUTPUT_DIR, "pca_interaction.pdf")
TRAJECTORY_FILE <- file.path(OUTPUT_DIR, "time_series_trajectory.pdf")
COR_ALL_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_embryos.pdf")
COR_72H_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_72h.pdf")
COR_135H_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_135h.pdf")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "qc_summary_embryos.txt")

# Colors
PHOTOPERIOD_COLORS <- c("SD" = "#D95F02", "LD" = "#1B9E77")
TIMEPOINT_COLORS <- c("72h" = "#8DA0CB", "135h" = "#FC8D62")

# =============================================================================
# Helper Functions
# =============================================================================

#' Load embryos count data
load_embryos_data <- function(count_dir) {
  cat("Loading embryos count matrices...\n")

  # Load all embryos
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
  all_meta$timepoint <- ifelse(grepl("72h", all_meta$New.abreviation), "72h", "135h")
  all_meta$sample_id <- paste0(all_meta$BioProject, "_", all_meta$Run)

  # Ensure sample order matches
  all_counts <- all_counts[, all_meta$sample_id]

  cat(sprintf("  All embryos: %d samples, %d genes\n", ncol(all_counts), nrow(all_counts)))
  cat(sprintf("    SD: %d, LD: %d\n", sum(all_meta$photoperiod == "SD"), sum(all_meta$photoperiod == "LD")))
  cat(sprintf("    72h: %d, 135h: %d\n\n", sum(all_meta$timepoint == "72h"), sum(all_meta$timepoint == "135h")))

  return(list(counts = all_counts, metadata = all_meta))
}

#' Create PCA with custom coloring
create_pca_custom <- function(vst_data, metadata, color_by, shape_by = NULL,
                               output_file, title, color_palette, shape_palette = NULL) {
  cat(sprintf("Creating PCA: %s...\n", basename(output_file)))

  # Get PCA data
  intgroup <- if (is.null(shape_by)) c(color_by) else c(color_by, shape_by)
  pcaData <- plotPCA(vst_data, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, metadata, by = "sample_id")

  # Create plot
  if (is.null(shape_by)) {
    p <- ggplot(plot_data, aes_string(x = "PC1", y = "PC2", color = color_by)) +
      geom_point(size = 4, alpha = 0.8)
  } else {
    p <- ggplot(plot_data, aes_string(x = "PC1", y = "PC2", color = color_by, shape = shape_by)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_shape_manual(values = c(16, 17))
  }

  p <- p +
    scale_color_manual(values = color_palette) +
    labs(
      title = title,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = color_by,
      shape = shape_by
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  PC1: %d%%, PC2: %d%%\n\n", percentVar[1], percentVar[2]))

  return(plot_data)
}

#' Create timepoint-specific PCA
create_timepoint_pca <- function(counts, metadata, timepoint_val, output_file, title) {
  cat(sprintf("Creating %s PCA...\n", timepoint_val))

  # Subset data
  subset_samples <- metadata$timepoint == timepoint_val
  subset_counts <- counts[, subset_samples]
  subset_meta <- metadata[subset_samples, ]

  cat(sprintf("  Samples: %d\n", ncol(subset_counts)))

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = subset_counts,
    colData = subset_meta,
    design = ~ photoperiod
  )

  # VST
  vst_data <- vst(dds, blind = TRUE)

  # Create PCA
  pcaData <- plotPCA(vst_data, intgroup = "photoperiod", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, subset_meta, by = "sample_id")

  # Plot
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = photoperiod,
                              label = New.abreviation)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -1, size = 3) +
    scale_color_manual(values = PHOTOPERIOD_COLORS) +
    labs(
      title = title,
      subtitle = sprintf("n = %d per group", sum(subset_meta$photoperiod == "SD")),
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Photoperiod"
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

  return(list(vst = vst_data, pca = plot_data))
}

#' Create time series trajectory plot
create_trajectory_plot <- function(vst_data, metadata, output_file, title) {
  cat("Creating time series trajectory plot...\n")

  # Get PCA data
  pcaData <- plotPCA(vst_data, intgroup = c("photoperiod", "timepoint"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, metadata, by = "sample_id")

  # Calculate group centroids
  centroids <- plot_data %>%
    group_by(photoperiod, timepoint) %>%
    summarize(
      PC1_mean = mean(PC1),
      PC2_mean = mean(PC2),
      .groups = "drop"
    )

  # Create trajectory lines (72h → 135h for each photoperiod)
  trajectory_lines <- centroids %>%
    arrange(photoperiod, timepoint) %>%
    group_by(photoperiod) %>%
    mutate(
      PC1_start = lag(PC1_mean),
      PC2_start = lag(PC2_mean)
    ) %>%
    filter(!is.na(PC1_start))

  # Plot
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = photoperiod, shape = timepoint)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_point(data = centroids, aes(x = PC1_mean, y = PC2_mean), size = 6, alpha = 0.9) +
    geom_segment(data = trajectory_lines,
                 aes(x = PC1_start, y = PC2_start, xend = PC1_mean, yend = PC2_mean),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 size = 1.5, alpha = 0.8) +
    scale_color_manual(values = PHOTOPERIOD_COLORS) +
    scale_shape_manual(values = c("72h" = 16, "135h" = 17)) +
    labs(
      title = title,
      subtitle = "Arrows show 72h → 135h developmental trajectory",
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Photoperiod",
      shape = "Timepoint"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  Trajectories show developmental progression from 72h to 135h\n\n"))
}

#' Create correlation heatmap
create_cor_heatmap <- function(counts, metadata, group_cols, output_file, title) {
  cat(sprintf("Creating correlation heatmap: %s...\n", basename(output_file)))

  # Log transform
  log_counts <- log2(counts + 1)

  # Calculate correlations
  cor_matrix <- cor(log_counts, method = "pearson")

  # Create annotation
  annotation_df <- as.data.frame(metadata[, group_cols, drop = FALSE])
  rownames(annotation_df) <- colnames(cor_matrix)

  # Colors
  ann_colors <- list(
    photoperiod = PHOTOPERIOD_COLORS,
    timepoint = TIMEPOINT_COLORS
  )

  # Create heatmap
  pdf(output_file, width = 10, height = 9)
  pheatmap(
    cor_matrix,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0.7, 1.0, length.out = 101),
    border_color = NA
  )
  dev.off()

  cat(sprintf("  Saved: %s\n", output_file))

  # Summary stats
  diag(cor_matrix) <- NA
  cat(sprintf("  Mean correlation: %.3f\n", mean(cor_matrix, na.rm = TRUE)))
  cat(sprintf("  Min correlation: %.3f\n\n", min(cor_matrix, na.rm = TRUE)))

  return(cor_matrix)
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  cat("=======================================================================\n")
  cat("EMBRYOS QC PLOTS - TIME SERIES ANALYSIS\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

  # Load data
  data <- load_embryos_data(COUNT_DIR)
  counts <- data$counts
  metadata <- data$metadata

  # Create DESeq2 object for all embryos
  cat("Creating DESeq2 object (all embryos)...\n")
  dds_all <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ photoperiod + timepoint + photoperiod:timepoint
  )
  vst_all <- vst(dds_all, blind = TRUE)
  cat("  Done\n\n")

  # === GLOBAL EMBRYOS PCA ===
  cat("--- Global Embryos Analysis ---\n")

  pca_all <- create_pca_custom(
    vst_all, metadata, "photoperiod", "timepoint",
    PCA_ALL_FILE,
    paste("Embryos PCA - All 12 Samples (", toupper(COUNT_TYPE), ")", sep = ""),
    PHOTOPERIOD_COLORS
  )

  pca_photo <- create_pca_custom(
    vst_all, metadata, "photoperiod", NULL,
    PCA_PHOTOPERIOD_FILE,
    paste("Embryos PCA - Colored by Photoperiod (", toupper(COUNT_TYPE), ")", sep = ""),
    PHOTOPERIOD_COLORS
  )

  pca_time <- create_pca_custom(
    vst_all, metadata, "timepoint", NULL,
    PCA_TIMEPOINT_FILE,
    paste("Embryos PCA - Colored by Timepoint (", toupper(COUNT_TYPE), ")", sep = ""),
    TIMEPOINT_COLORS
  )

  # Create interaction plot (same as all but emphasized)
  pca_interaction <- create_pca_custom(
    vst_all, metadata, "photoperiod", "timepoint",
    PCA_INTERACTION_FILE,
    paste("Embryos Interaction - Photoperiod × Time (", toupper(COUNT_TYPE), ")", sep = ""),
    PHOTOPERIOD_COLORS
  )

  # Create trajectory plot
  create_trajectory_plot(
    vst_all, metadata,
    TRAJECTORY_FILE,
    paste("Embryos Developmental Trajectory (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # Correlation heatmap (all embryos)
  cor_all <- create_cor_heatmap(
    counts, metadata, c("photoperiod", "timepoint"),
    COR_ALL_FILE,
    paste("Embryos Correlation Heatmap - All 12 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === 72H ANALYSIS ===
  cat("--- 72h Timepoint (Early Development) ---\n")

  data_72h <- create_timepoint_pca(
    counts, metadata, "72h",
    PCA_72H_FILE,
    paste("Embryos 72h - SD vs LD (", toupper(COUNT_TYPE), ")", sep = "")
  )

  samples_72h <- metadata$timepoint == "72h"
  cor_72h <- create_cor_heatmap(
    counts[, samples_72h], metadata[samples_72h, ], c("photoperiod"),
    COR_72H_FILE,
    paste("Embryos 72h Correlation Heatmap (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === 135H ANALYSIS ===
  cat("--- 135h Timepoint (Diapause Entry) ---\n")

  data_135h <- create_timepoint_pca(
    counts, metadata, "135h",
    PCA_135H_FILE,
    paste("Embryos 135h - SD vs LD (", toupper(COUNT_TYPE), ")", sep = "")
  )

  samples_135h <- metadata$timepoint == "135h"
  cor_135h <- create_cor_heatmap(
    counts[, samples_135h], metadata[samples_135h, ], c("photoperiod"),
    COR_135H_FILE,
    paste("Embryos 135h Correlation Heatmap (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === SUMMARY ===
  cat("Writing summary statistics...\n")
  sink(SUMMARY_FILE)
  cat("EMBRYOS QC SUMMARY\n")
  cat("==================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("SAMPLE DISTRIBUTION:\n")
  cat(sprintf("  Total: %d\n", ncol(counts)))
  cat(sprintf("  SD_72h: %d\n", sum(metadata$photoperiod == "SD" & metadata$timepoint == "72h")))
  cat(sprintf("  LD_72h: %d\n", sum(metadata$photoperiod == "LD" & metadata$timepoint == "72h")))
  cat(sprintf("  SD_135h: %d\n", sum(metadata$photoperiod == "SD" & metadata$timepoint == "135h")))
  cat(sprintf("  LD_135h: %d\n\n", sum(metadata$photoperiod == "LD" & metadata$timepoint == "135h")))

  cat("CORRELATION SUMMARY:\n")
  cat("  All embryos:\n")
  diag(cor_all) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_all, na.rm = TRUE), min(cor_all, na.rm = TRUE)))

  cat("  72h timepoint:\n")
  diag(cor_72h) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_72h, na.rm = TRUE), min(cor_72h, na.rm = TRUE)))

  cat("  135h timepoint:\n")
  diag(cor_135h) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n\n", mean(cor_135h, na.rm = TRUE), min(cor_135h, na.rm = TRUE)))

  cat("EXPECTED RESULTS:\n")
  cat("  - Photoperiod separation at both timepoints\n")
  cat("  - Possible time effect (72h vs 135h)\n")
  cat("  - Interaction: Photoperiod effect may increase over time\n")
  cat("  - Smooth developmental trajectory (not abrupt jumps)\n")
  cat("  - Within-group correlations > 0.9\n\n")

  cat("BIOLOGICAL INTERPRETATION:\n")
  cat("  - 72h: Early diapause programming\n")
  cat("  - 135h: Diapause entry (~120h in SD embryos)\n")
  cat("  - Trajectory shows diapause preparation vs normal development\n\n")

  cat("QC CHECKLIST:\n")
  cat("  [ ] Photoperiod separation visible at both timepoints\n")
  cat("  [ ] Within-timepoint correlations > 0.9\n")
  cat("  [ ] Smooth trajectory (no abrupt changes)\n")
  cat("  [ ] No extreme outliers\n")
  cat("  [ ] Proceed to larvae QC (script 04)\n")

  sink()
  cat(sprintf("  Saved: %s\n\n", SUMMARY_FILE))

  cat("=======================================================================\n")
  cat("EMBRYOS QC COMPLETE\n")
  cat("=======================================================================\n\n")
  cat("Review:\n")
  cat("  - Time series trajectory: ", TRAJECTORY_FILE, "\n", sep = "")
  cat("  - All plots in:", OUTPUT_DIR, "\n")
  cat("  - Summary:", SUMMARY_FILE, "\n\n")
  cat("Next: Run script 04 (larvae QC)\n\n")
}

# Run main
main()
