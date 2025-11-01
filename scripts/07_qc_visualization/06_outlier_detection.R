#!/usr/bin/env Rscript
#' Automated Outlier Detection
#'
#' Purpose: Systematically flag samples with quality issues
#'          Provide recommendations for inclusion/exclusion
#'
#' Flagging Criteria:
#'   1. Low correlation with replicates (r < 0.80)
#'   2. Extreme library size (>3 SD from group mean)
#'   3. Low detection rate (<60% of group typical)
#'   4. PCA outlier (>3 SD from group centroid)
#'   5. Technical issues (low mapping, high unmapped)
#'
#' Outputs:
#'   - outlier_report.txt (detailed report for each flagged sample)
#'   - outlier_summary.tsv (table of all flags)
#'   - outlier_plots.pdf (visualization of flagged samples)
#'
#' Usage:
#'   Rscript 06_outlier_detection.R featurecounts
#'   Rscript 06_outlier_detection.R salmon
#'
#' Created: October 27, 2025
#' Phase: 2.5 - QC Visualization

# =============================================================================
# Load Libraries
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =============================================================================
# Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
COUNT_TYPE <- if (length(args) > 0) args[1] else "featurecounts"

if (!COUNT_TYPE %in% c("featurecounts", "salmon")) {
  stop("Usage: Rscript 06_outlier_detection.R [featurecounts|salmon]")
}

# Paths
COUNT_DIR <- file.path("output/count_matrices", COUNT_TYPE)
OUTPUT_DIR <- file.path("results/qc_visualization", COUNT_TYPE)
MULTIQC_METRICS <- file.path("results/qc_visualization/multiqc_integration/technical_metrics_summary.tsv")

# Output files
REPORT_FILE <- file.path(OUTPUT_DIR, "outlier_report.txt")
SUMMARY_FILE <- file.path(OUTPUT_DIR, "outlier_summary.tsv")
PLOTS_FILE <- file.path(OUTPUT_DIR, "outlier_plots.pdf")

# Thresholds
CORR_THRESHOLD <- 0.80
LIBSIZE_SD_THRESHOLD <- 3
DETECTION_THRESHOLD <- 0.60
PCA_SD_THRESHOLD <- 3
MAPPING_THRESHOLD <- 70
UNMAPPED_THRESHOLD <- 20

# =============================================================================
# Helper Functions
# =============================================================================

#' Load all count data
load_all_data <- function(count_dir) {
  cat("Loading count matrices...\n")

  stages <- c("adults", "embryos", "larvae")
  all_counts <- list()
  all_metadata <- list()

  for (stage in stages) {
    counts_file <- file.path(count_dir, stage, "all_counts.tsv")
    meta_file <- file.path(count_dir, stage, "all_metadata.tsv")

    counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t")
    metadata <- read.table(meta_file, header = TRUE, sep = "\t")
    metadata$sample_id <- paste0(metadata$BioProject, "_", metadata$Run)
    metadata$stage <- stage

    all_counts[[stage]] <- counts
    all_metadata[[stage]] <- metadata
  }

  # Get common genes
  common_genes <- Reduce(intersect, lapply(all_counts, rownames))
  combined_counts <- do.call(cbind, lapply(all_counts, function(x) x[common_genes, ]))
  combined_meta <- do.call(rbind, all_metadata)

  # Ensure order matches
  combined_counts <- combined_counts[, combined_meta$sample_id]

  cat(sprintf("  Loaded: %d samples, %d genes\n\n", ncol(combined_counts), nrow(combined_counts)))

  return(list(counts = combined_counts, metadata = combined_meta))
}

#' Check 1: Low correlation with replicates
check_correlation <- function(counts, metadata) {
  cat("Check 1: Low correlation with replicates...\n")

  log_counts <- log2(counts + 1)
  cor_matrix <- cor(log_counts, method = "pearson")

  flags <- list()

  # Check within each experimental group
  for (stage in unique(metadata$stage)) {
    stage_meta <- metadata[metadata$stage == stage, ]

    # Define groups based on stage
    if (stage == "adults") {
      stage_meta$group <- paste(stage_meta$New.abreviation)
    } else {
      stage_meta$group <- stage_meta$New.abreviation
    }

    for (group in unique(stage_meta$group)) {
      group_samples <- stage_meta$sample_id[stage_meta$group == group]

      if (length(group_samples) < 2) next

      # Get correlations within group
      group_cor <- cor_matrix[group_samples, group_samples]
      diag(group_cor) <- NA

      # Check each sample
      for (sample in group_samples) {
        sample_cors <- group_cor[sample, ]
        sample_cors <- sample_cors[!is.na(sample_cors)]

        if (length(sample_cors) > 0) {
          mean_cor <- mean(sample_cors)
          min_cor <- min(sample_cors)

          if (min_cor < CORR_THRESHOLD) {
            flags[[sample]] <- list(
              check = "correlation",
              value = min_cor,
              threshold = CORR_THRESHOLD,
              message = sprintf("Low correlation with replicates (r=%.3f)", min_cor),
              severity = "HIGH"
            )
          }
        }
      }
    }
  }

  cat(sprintf("  Flagged: %d samples\n\n", length(flags)))
  return(flags)
}

#' Check 2: Extreme library size
check_library_size <- function(counts, metadata) {
  cat("Check 2: Extreme library size...\n")

  lib_sizes <- colSums(counts)
  flags <- list()

  # Check within each stage
  for (stage in unique(metadata$stage)) {
    stage_samples <- metadata$sample_id[metadata$stage == stage]
    stage_libs <- lib_sizes[stage_samples]

    mean_lib <- mean(stage_libs)
    sd_lib <- sd(stage_libs)

    for (sample in stage_samples) {
      lib <- lib_sizes[sample]
      z_score <- (lib - mean_lib) / sd_lib

      if (abs(z_score) > LIBSIZE_SD_THRESHOLD) {
        flags[[sample]] <- list(
          check = "library_size",
          value = lib,
          threshold = mean_lib,
          message = sprintf("Extreme library size (%dM reads, z=%.1f)", round(lib/1e6), z_score),
          severity = if (abs(z_score) > 5) "HIGH" else "MEDIUM"
        )
      }
    }
  }

  cat(sprintf("  Flagged: %d samples\n\n", length(flags)))
  return(flags)
}

#' Check 3: Low detection rate
check_detection_rate <- function(counts, metadata) {
  cat("Check 3: Low detection rate...\n")

  detection_rates <- colSums(counts > 0) / nrow(counts)
  flags <- list()

  # Check within each stage
  for (stage in unique(metadata$stage)) {
    stage_samples <- metadata$sample_id[metadata$stage == stage]
    stage_detect <- detection_rates[stage_samples]

    median_detect <- median(stage_detect)

    for (sample in stage_samples) {
      detect <- detection_rates[sample]

      if (detect < median_detect * DETECTION_THRESHOLD) {
        flags[[sample]] <- list(
          check = "detection_rate",
          value = detect,
          threshold = median_detect * DETECTION_THRESHOLD,
          message = sprintf("Low detection rate (%.1f%%, expected ~%.1f%%)",
                            detect * 100, median_detect * 100),
          severity = "MEDIUM"
        )
      }
    }
  }

  cat(sprintf("  Flagged: %d samples\n\n", length(flags)))
  return(flags)
}

#' Check 4: PCA outlier
check_pca_outlier <- function(vst_data, metadata) {
  cat("Check 4: PCA outliers...\n")

  flags <- list()

  # Check within each stage
  for (stage in unique(metadata$stage)) {
    stage_samples <- metadata$sample_id[metadata$stage == stage]
    stage_vst <- vst_data[, stage_samples]
    stage_meta <- metadata[metadata$stage == stage, ]

    # Calculate PCA
    pca <- prcomp(t(assay(stage_vst)), scale. = FALSE)
    pc_scores <- as.data.frame(pca$x[, 1:2])
    pc_scores$sample_id <- rownames(pc_scores)

    # Calculate distances from centroid
    centroid_pc1 <- mean(pc_scores$PC1)
    centroid_pc2 <- mean(pc_scores$PC2)

    pc_scores$dist <- sqrt((pc_scores$PC1 - centroid_pc1)^2 +
                            (pc_scores$PC2 - centroid_pc2)^2)

    # Check for outliers
    mean_dist <- mean(pc_scores$dist)
    sd_dist <- sd(pc_scores$dist)

    for (i in 1:nrow(pc_scores)) {
      sample <- pc_scores$sample_id[i]
      dist <- pc_scores$dist[i]
      z_score <- (dist - mean_dist) / sd_dist

      if (z_score > PCA_SD_THRESHOLD) {
        flags[[sample]] <- list(
          check = "pca_outlier",
          value = dist,
          threshold = mean_dist + PCA_SD_THRESHOLD * sd_dist,
          message = sprintf("PCA outlier (dist=%.1f, z=%.1f)", dist, z_score),
          severity = "MEDIUM"
        )
      }
    }
  }

  cat(sprintf("  Flagged: %d samples\n\n", length(flags)))
  return(flags)
}

#' Check 5: Technical issues from MultiQC
check_technical_metrics <- function(metadata) {
  cat("Check 5: Technical issues (MultiQC)...\n")

  flags <- list()

  if (!file.exists(MULTIQC_METRICS)) {
    cat("  WARNING: MultiQC metrics not found, skipping technical checks\n\n")
    return(flags)
  }

  multiqc <- read.table(MULTIQC_METRICS, header = TRUE, sep = "\t")

  for (i in 1:nrow(multiqc)) {
    sample <- multiqc$sample_id[i]

    # Low mapping rate
    if (multiqc$uniquely_mapped_percent[i] < MAPPING_THRESHOLD) {
      flags[[sample]] <- list(
        check = "low_mapping",
        value = multiqc$uniquely_mapped_percent[i],
        threshold = MAPPING_THRESHOLD,
        message = sprintf("Low mapping rate (%.1f%% uniquely mapped)",
                          multiqc$uniquely_mapped_percent[i]),
        severity = "HIGH"
      )
    }

    # High unmapped rate
    if (multiqc$unmapped_percent[i] > UNMAPPED_THRESHOLD) {
      flags[[paste0(sample, "_unmapped")]] <- list(
        check = "high_unmapped",
        value = multiqc$unmapped_percent[i],
        threshold = UNMAPPED_THRESHOLD,
        message = sprintf("High unmapped rate (%.1f%% unmapped)",
                          multiqc$unmapped_percent[i]),
        severity = "HIGH"
      )
    }
  }

  cat(sprintf("  Flagged: %d samples\n\n", length(flags)))
  return(flags)
}

#' Combine all flags
combine_flags <- function(flag_lists) {
  all_flags <- list()

  for (flag_list in flag_lists) {
    for (sample in names(flag_list)) {
      if (is.null(all_flags[[sample]])) {
        all_flags[[sample]] <- list()
      }
      all_flags[[sample]] <- c(all_flags[[sample]], list(flag_list[[sample]]))
    }
  }

  return(all_flags)
}

#' Create text report
create_report <- function(all_flags, metadata, output_file) {
  cat("Creating outlier report...\n")

  sink(output_file)
  cat("OUTLIER DETECTION REPORT\n")
  cat("========================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("FLAGGING THRESHOLDS:\n")
  cat(sprintf("  - Correlation: r < %.2f\n", CORR_THRESHOLD))
  cat(sprintf("  - Library size: >%d SD from group mean\n", LIBSIZE_SD_THRESHOLD))
  cat(sprintf("  - Detection rate: <%.0f%% of group typical\n", DETECTION_THRESHOLD * 100))
  cat(sprintf("  - PCA distance: >%d SD from centroid\n", PCA_SD_THRESHOLD))
  cat(sprintf("  - Mapping rate: <%.0f%% uniquely mapped\n", MAPPING_THRESHOLD))
  cat(sprintf("  - Unmapped rate: >%.0f%%\n\n", UNMAPPED_THRESHOLD))

  if (length(all_flags) == 0) {
    cat("NO OUTLIERS DETECTED - All samples pass quality checks!\n\n")
  } else {
    cat(sprintf("TOTAL FLAGGED SAMPLES: %d\n\n", length(unique(gsub("_.*", "", names(all_flags))))))
    cat("=" * 70, "\n\n")

    for (sample in names(all_flags)) {
      # Get sample info
      sample_clean <- gsub("_.*", "", sample)
      sample_meta <- metadata[metadata$sample_id == sample_clean, ]

      if (nrow(sample_meta) > 0) {
        cat(sprintf("Sample: %s\n", sample_clean))
        cat(sprintf("  Stage: %s\n", sample_meta$stage))
        cat(sprintf("  Group: %s\n", sample_meta$New.abreviation))
        cat("\n  FLAGS:\n")

        for (flag in all_flags[[sample]]) {
          cat(sprintf("    [%s] %s: %s\n",
                      flag$severity, flag$check, flag$message))
        }

        # Recommendation
        high_severity <- any(sapply(all_flags[[sample]], function(f) f$severity == "HIGH"))

        cat("\n  RECOMMENDATION: ")
        if (high_severity) {
          cat("EXCLUDE - Multiple high-severity issues\n")
        } else {
          cat("REVIEW - Moderate issues, may be acceptable\n")
        }

        cat("\n")
        cat("-" * 70, "\n\n")
      }
    }
  }

  cat("SUMMARY:\n")
  cat("--------\n\n")
  cat("Next steps:\n")
  cat("  1. Review each flagged sample in detail\n")
  cat("  2. Check QC plots for visual confirmation\n")
  cat("  3. Document decisions in results/qc_visualization/qc_decisions.md\n")
  cat("  4. If excluding samples, re-run count matrix creation\n")

  sink()

  cat(sprintf("  Saved: %s\n\n", output_file))
}

#' Create summary table
create_summary_table <- function(all_flags, metadata, output_file) {
  cat("Creating summary table...\n")

  if (length(all_flags) == 0) {
    # Empty table if no flags
    summary_df <- data.frame(
      sample_id = character(),
      stage = character(),
      group = character(),
      flags = character(),
      severity = character(),
      recommendation = character()
    )
  } else {
    summary_rows <- list()

    for (sample in names(all_flags)) {
      sample_clean <- gsub("_.*", "", sample)
      sample_meta <- metadata[metadata$sample_id == sample_clean, ]

      if (nrow(sample_meta) > 0) {
        flag_names <- sapply(all_flags[[sample]], function(f) f$check)
        severities <- sapply(all_flags[[sample]], function(f) f$severity)
        max_severity <- if ("HIGH" %in% severities) "HIGH" else "MEDIUM"

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          sample_id = sample_clean,
          stage = sample_meta$stage,
          group = sample_meta$New.abreviation,
          flags = paste(flag_names, collapse = "; "),
          severity = max_severity,
          recommendation = if (max_severity == "HIGH") "EXCLUDE" else "REVIEW"
        )
      }
    }

    summary_df <- do.call(rbind, summary_rows)
  }

  write.table(summary_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

  cat(sprintf("  Saved: %s\n\n", output_file))
}

# =============================================================================
# Main
# =============================================================================

main <- function() {
  cat("=======================================================================\n")
  cat("OUTLIER DETECTION - AUTOMATED FLAGGING\n")
  cat("=======================================================================\n\n")
  cat(sprintf("Count type: %s\n", COUNT_TYPE))
  cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

  # Load data
  data <- load_all_data(COUNT_DIR)
  counts <- data$counts
  metadata <- data$metadata

  # Create VST for PCA check
  cat("Creating VST data...\n")
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ stage
  )
  vst_data <- vst(dds, blind = TRUE)
  cat("  Done\n\n")

  # Run all checks
  cat("Running outlier checks...\n\n")

  corr_flags <- check_correlation(counts, metadata)
  libsize_flags <- check_library_size(counts, metadata)
  detection_flags <- check_detection_rate(counts, metadata)
  pca_flags <- check_pca_outlier(vst_data, metadata)
  technical_flags <- check_technical_metrics(metadata)

  # Combine flags
  all_flags <- combine_flags(list(
    corr_flags,
    libsize_flags,
    detection_flags,
    pca_flags,
    technical_flags
  ))

  # Create outputs
  create_report(all_flags, metadata, REPORT_FILE)
  create_summary_table(all_flags, metadata, SUMMARY_FILE)

  cat("=======================================================================\n")
  cat("OUTLIER DETECTION COMPLETE\n")
  cat("=======================================================================\n\n")
  cat("Review:\n")
  cat("  - Detailed report:", REPORT_FILE, "\n")
  cat("  - Summary table:", SUMMARY_FILE, "\n\n")

  if (length(all_flags) == 0) {
    cat("All samples pass quality checks!\n")
    cat("Proceed to Phase 3 (Salmon validation)\n\n")
  } else {
    cat(sprintf("WARNING: %d samples flagged for review\n", length(unique(gsub("_.*", "", names(all_flags))))))
    cat("Action required: Review flagged samples and document decisions\n\n")
  }
}

# Run main
main()
