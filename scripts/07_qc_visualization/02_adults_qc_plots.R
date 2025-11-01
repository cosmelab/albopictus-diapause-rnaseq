#!/usr/bin/env Rscript
#' Adults QC Plots - Split by Blood Meal Status
#'
#' Purpose: Stage-specific QC for adults (16 samples)
#'          Split analysis by blood meal status (confounding factor)
#'          Publication-quality figures for supplemental materials
#'
#' Experimental Design:
#'   - 2×2 factorial: Photoperiod (SD/LD) × Blood Meal (Unfed/Fed)
#'   - 4 replicates per group
#'   - Primary comparison: SD_Fed vs LD_Fed (diapause vs non-diapause oocytes)
#'   - Secondary: SD_Unfed vs LD_Unfed (photoperiod without reproduction)
#'
#' Outputs:
#'   - Global adults PCA (all 16)
#'   - Unfed-specific PCA (8 samples)
#'   - Fed-specific PCA (8 samples) - PRIMARY
#'   - Correlation heatmaps for each comparison
#'   - Boxplots of normalized counts
#'
#' Expected Results:
#'   - Replicates cluster together (within-group r > 0.9)
#'   - Clear SD vs LD separation
#'   - Blood meal effect visible
#'   - No extreme outliers
#'
#' Usage:
#'   Rscript 02_adults_qc_plots.R featurecounts
#'   Rscript 02_adults_qc_plots.R salmon
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
  stop("Usage: Rscript 02_adults_qc_plots.R [featurecounts|salmon]")
}

# Paths
COUNT_DIR <- file.path("output/count_matrices", COUNT_TYPE, "adults")
OUTPUT_DIR <- file.path("results/qc_visualization", COUNT_TYPE, "adults")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Output files
PCA_ALL_FILE <- file.path(OUTPUT_DIR, "pca_all_adults.pdf")
PCA_UNFED_FILE <- file.path(OUTPUT_DIR, "pca_unfed_only.pdf")
PCA_FED_FILE <- file.path(OUTPUT_DIR, "pca_fed_only.pdf")
PCA_PHOTOPERIOD_FILE <- file.path(OUTPUT_DIR, "pca_by_photoperiod.pdf")
PCA_BLOODMEAL_FILE <- file.path(OUTPUT_DIR, "pca_by_blood_meal.pdf")
COR_ALL_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_adults.pdf")
COR_UNFED_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_unfed.pdf")
COR_FED_FILE <- file.path(OUTPUT_DIR, "correlation_heatmap_fed.pdf")
BOXPLOT_FILE <- file.path(OUTPUT_DIR, "boxplot_by_group.pdf")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "qc_summary_adults.txt")

# Colors
PHOTOPERIOD_COLORS <- c("SD" = "#D95F02", "LD" = "#1B9E77")
BLOODMEAL_COLORS <- c("Unfed" = "#7570B3", "Fed" = "#E7298A")

# =============================================================================
# Helper Functions
# =============================================================================

#' Load adults count data
load_adults_data <- function(count_dir) {
  cat("Loading adults count matrices...\n")

  # Load all adults
  all_counts <- read.table(
    file.path(count_dir, "all_counts.tsv"),
    header = TRUE, row.names = 1, sep = "\t"
  )
  all_meta <- read.table(
    file.path(count_dir, "all_metadata.tsv"),
    header = TRUE, sep = "\t"
  )

  # Extract photoperiod and blood meal from sample names
  all_meta$photoperiod <- ifelse(grepl("^SD", all_meta$New.abreviation), "SD", "LD")
  all_meta$blood_meal <- ifelse(grepl("Unfed", all_meta$New.abreviation), "Unfed", "Fed")
  all_meta$sample_id <- paste0(all_meta$BioProject, "_", all_meta$Run)

  # Ensure sample order matches
  all_counts <- all_counts[, all_meta$sample_id]

  cat(sprintf("  All adults: %d samples, %d genes\n", ncol(all_counts), nrow(all_counts)))
  cat(sprintf("    SD: %d, LD: %d\n", sum(all_meta$photoperiod == "SD"), sum(all_meta$photoperiod == "LD")))
  cat(sprintf("    Unfed: %d, Fed: %d\n\n", sum(all_meta$blood_meal == "Unfed"), sum(all_meta$blood_meal == "Fed")))

  return(list(counts = all_counts, metadata = all_meta))
}

#' Create PCA plot with custom coloring
create_pca_custom <- function(vst_data, metadata, color_by, shape_by = NULL,
                               output_file, title, color_palette) {
  cat(sprintf("Creating PCA plot: %s...\n", basename(output_file)))

  # Get PCA data
  pcaData <- plotPCA(vst_data, intgroup = c(color_by), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Merge with full metadata
  pcaData$sample_id <- rownames(pcaData)
  plot_data <- merge(pcaData, metadata, by = "sample_id")

  # Create base plot
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

#' Create subset PCA (unfed or fed only)
create_subset_pca <- function(counts, metadata, subset_col, subset_val,
                               output_file, title) {
  cat(sprintf("Creating %s PCA...\n", subset_val))

  # Subset data
  subset_samples <- metadata[[subset_col]] == subset_val
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

#' Create correlation heatmap
create_cor_heatmap <- function(counts, metadata, group_col, output_file, title) {
  cat(sprintf("Creating correlation heatmap: %s...\n", basename(output_file)))

  # Log transform
  log_counts <- log2(counts + 1)

  # Calculate correlations
  cor_matrix <- cor(log_counts, method = "pearson")

  # Create annotation
  if (!is.null(group_col)) {
    annotation_df <- as.data.frame(metadata[, group_col, drop = FALSE])
    rownames(annotation_df) <- colnames(cor_matrix)

    # Colors
    if (group_col == "photoperiod") {
      ann_colors <- list(photoperiod = PHOTOPERIOD_COLORS)
    } else if (group_col == "blood_meal") {
      ann_colors <- list(blood_meal = BLOODMEAL_COLORS)
    } else {
      ann_colors <- NULL
    }
  } else {
    annotation_df <- NULL
    ann_colors <- NULL
  }

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
  cat(sprintf("  Min correlation: %.3f\n", min(cor_matrix, na.rm = TRUE)))

  # Within-group correlations
  if (!is.null(group_col)) {
    for (group in unique(metadata[[group_col]])) {
      group_samples <- metadata$sample_id[metadata[[group_col]] == group]
      group_cor <- cor_matrix[group_samples, group_samples]
      diag(group_cor) <- NA
      cat(sprintf("  Within-%s correlation: %.3f (min: %.3f)\n",
                  group, mean(group_cor, na.rm = TRUE), min(group_cor, na.rm = TRUE)))
    }
  }
  cat("\n")

  return(cor_matrix)
}

#' Create boxplot of normalized counts
create_boxplot <- function(vst_data, metadata, output_file, title) {
  cat("Creating boxplot of normalized counts...\n")

  # Get VST counts
  vst_counts <- assay(vst_data)

  # Reshape for plotting
  vst_long <- as.data.frame(vst_counts) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "vst_count")

  # Merge with metadata
  vst_long <- merge(vst_long, metadata[, c("sample_id", "photoperiod", "blood_meal")], by = "sample_id")

  # Create combined group
  vst_long$group <- paste(vst_long$photoperiod, vst_long$blood_meal, sep = "_")

  # Plot
  p <- ggplot(vst_long, aes(x = group, y = vst_count, fill = photoperiod)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = PHOTOPERIOD_COLORS) +
    labs(
      title = title,
      subtitle = "VST-normalized expression distribution by treatment group",
      x = "Treatment Group",
      y = "VST-normalized Expression",
      fill = "Photoperiod"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 7)
  cat(sprintf("  Saved: %s\n\n", output_file))
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  cat("=======================================================================\n")
  cat("ADULTS QC PLOTS - SPLIT BY BLOOD MEAL STATUS\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

  # Load data
  data <- load_adults_data(COUNT_DIR)
  counts <- data$counts
  metadata <- data$metadata

  # Create DESeq2 object for all adults
  cat("Creating DESeq2 object (all adults)...\n")
  dds_all <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ photoperiod + blood_meal
  )
  vst_all <- vst(dds_all, blind = TRUE)
  cat("  Done\n\n")

  # === GLOBAL ADULTS PCA ===
  cat("--- Global Adults Analysis ---\n")

  pca_all <- create_pca_custom(
    vst_all, metadata, "photoperiod", "blood_meal",
    PCA_ALL_FILE,
    paste("Adults PCA - All 16 Samples (", toupper(COUNT_TYPE), ")", sep = ""),
    PHOTOPERIOD_COLORS
  )

  pca_photo <- create_pca_custom(
    vst_all, metadata, "photoperiod", NULL,
    PCA_PHOTOPERIOD_FILE,
    paste("Adults PCA - Colored by Photoperiod (", toupper(COUNT_TYPE), ")", sep = ""),
    PHOTOPERIOD_COLORS
  )

  pca_blood <- create_pca_custom(
    vst_all, metadata, "blood_meal", NULL,
    PCA_BLOODMEAL_FILE,
    paste("Adults PCA - Colored by Blood Meal Status (", toupper(COUNT_TYPE), ")", sep = ""),
    BLOODMEAL_COLORS
  )

  cor_all <- create_cor_heatmap(
    counts, metadata, "photoperiod",
    COR_ALL_FILE,
    paste("Adults Correlation Heatmap - All 16 Samples (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === UNFED ANALYSIS ===
  cat("--- Unfed Analysis (Secondary) ---\n")

  unfed_data <- create_subset_pca(
    counts, metadata, "blood_meal", "Unfed",
    PCA_UNFED_FILE,
    paste("Adults Unfed - SD vs LD (", toupper(COUNT_TYPE), ")", sep = "")
  )

  unfed_samples <- metadata$blood_meal == "Unfed"
  cor_unfed <- create_cor_heatmap(
    counts[, unfed_samples], metadata[unfed_samples, ], "photoperiod",
    COR_UNFED_FILE,
    paste("Adults Unfed Correlation Heatmap (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === FED ANALYSIS (PRIMARY) ===
  cat("--- Fed Analysis (PRIMARY COMPARISON) ---\n")

  fed_data <- create_subset_pca(
    counts, metadata, "blood_meal", "Fed",
    PCA_FED_FILE,
    paste("Adults Fed - SD vs LD [PRIMARY] (", toupper(COUNT_TYPE), ")", sep = "")
  )

  fed_samples <- metadata$blood_meal == "Fed"
  cor_fed <- create_cor_heatmap(
    counts[, fed_samples], metadata[fed_samples, ], "photoperiod",
    COR_FED_FILE,
    paste("Adults Fed Correlation Heatmap [PRIMARY] (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === BOXPLOT ===
  create_boxplot(
    vst_all, metadata,
    BOXPLOT_FILE,
    paste("Adults Normalized Counts by Group (", toupper(COUNT_TYPE), ")", sep = "")
  )

  # === SUMMARY ===
  cat("Writing summary statistics...\n")
  sink(SUMMARY_FILE)
  cat("ADULTS QC SUMMARY\n")
  cat("=================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("SAMPLE DISTRIBUTION:\n")
  cat(sprintf("  Total: %d\n", ncol(counts)))
  cat(sprintf("  SD_Unfed: %d\n", sum(metadata$photoperiod == "SD" & metadata$blood_meal == "Unfed")))
  cat(sprintf("  LD_Unfed: %d\n", sum(metadata$photoperiod == "LD" & metadata$blood_meal == "Unfed")))
  cat(sprintf("  SD_Fed: %d\n", sum(metadata$photoperiod == "SD" & metadata$blood_meal == "Fed")))
  cat(sprintf("  LD_Fed: %d\n\n", sum(metadata$photoperiod == "LD" & metadata$blood_meal == "Fed")))

  cat("CORRELATION SUMMARY:\n")
  cat("  All adults:\n")
  diag(cor_all) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_all, na.rm = TRUE), min(cor_all, na.rm = TRUE)))

  cat("  Unfed only:\n")
  diag(cor_unfed) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n", mean(cor_unfed, na.rm = TRUE), min(cor_unfed, na.rm = TRUE)))

  cat("  Fed only (PRIMARY):\n")
  diag(cor_fed) <- NA
  cat(sprintf("    Mean: %.3f, Min: %.3f\n\n", mean(cor_fed, na.rm = TRUE), min(cor_fed, na.rm = TRUE)))

  cat("EXPECTED RESULTS:\n")
  cat("  - Replicates cluster together (within-group r > 0.9)\n")
  cat("  - Clear SD vs LD separation in PCA\n")
  cat("  - Blood meal effect visible\n")
  cat("  - Fed comparison shows strongest separation (diapause signal)\n\n")

  cat("QC CHECKLIST:\n")
  cat("  [ ] All replicates cluster together in PCA\n")
  cat("  [ ] Within-group correlations > 0.9\n")
  cat("  [ ] SD vs LD separation visible in Fed samples\n")
  cat("  [ ] No extreme outliers\n")
  cat("  [ ] Proceed to embryos QC (script 03)\n")

  sink()
  cat(sprintf("  Saved: %s\n\n", SUMMARY_FILE))

  cat("=======================================================================\n")
  cat("ADULTS QC COMPLETE\n")
  cat("=======================================================================\n\n")
  cat("Review:\n")
  cat("  - Primary comparison (Fed): ", PCA_FED_FILE, "\n", sep = "")
  cat("  - All plots in:", OUTPUT_DIR, "\n")
  cat("  - Summary:", SUMMARY_FILE, "\n\n")
  cat("Next: Run script 03 (embryos QC)\n\n")
}

# Run main
main()
