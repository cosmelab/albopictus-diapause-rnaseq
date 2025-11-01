#!/usr/bin/env Rscript
################################################################################
# Adult Female Differential Expression Analysis - Unfed Only
#
# Following Huang et al. 2015 approach: analyze unfed and fed separately
# to isolate photoperiod effect without blood meal confounding.
#
# This script analyzes UNFED females only:
# - LD_Unfed (16L:8D, non-diapause) vs SD_Unfed (8L:16D, non-diapause)
# - Design: ~ Photoperiod
# - Tests: Diapause-inducing (SD) vs Non-diapause (LD) photoperiod
#
# Author: L. Cosme
# Date: October 22, 2025
################################################################################

library(DESeq2)
library(tidyverse)
library(pheatmap)

# Set output options
options(width = 120)
options(scipen = 999)

cat("================================================================================\n")
cat("Adult Differential Expression Analysis - Unfed Females\n")
cat("================================================================================\n\n")

# Define paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
counts_file <- file.path(project_base, "results/02_count_matrices/adults_counts.tsv")
metadata_file <- file.path(project_base, "data/metadata/samples_adults.csv")
output_dir <- file.path(project_base, "results/03_differential_expression/adults")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project base:", project_base, "\n")
cat("Counts file:", counts_file, "\n")
cat("Metadata file:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n\n")

################################################################################
# 1. Load Data
################################################################################

cat("Loading count matrix...\n")
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

cat("\nLoading metadata...\n")
metadata <- read.csv(metadata_file, row.names = 16)  # Column 16 is "Sample"
cat("  Samples:", nrow(metadata), "\n")

# Subset to unfed samples only
cat("\nFiltering for UNFED samples only...\n")
unfed_samples <- metadata[grepl("Unfed", metadata$New.abreviation), ]
cat("  Unfed samples:", nrow(unfed_samples), "\n")
cat("  Groups:\n")
print(table(unfed_samples$New.abreviation))

# Subset counts to unfed samples
counts_unfed <- counts[, rownames(unfed_samples)]
cat("\n  Count matrix subset:", nrow(counts_unfed), "genes x", ncol(counts_unfed), "samples\n")

# Create simplified photoperiod column
unfed_samples$Photoperiod <- ifelse(grepl("LD", unfed_samples$New.abreviation), "LD", "SD")
unfed_samples$Photoperiod <- factor(unfed_samples$Photoperiod, levels = c("LD", "SD"))

cat("\n  Photoperiod groups:\n")
print(table(unfed_samples$Photoperiod))

################################################################################
# 2. Pre-filtering
################################################################################

cat("\n\nPre-filtering low-count genes...\n")
cat("  Original genes:", nrow(counts_unfed), "\n")

# Filter: keep genes with at least 10 counts in at least 25% of samples
min_samples <- ceiling(0.25 * ncol(counts_unfed))
keep <- rowSums(counts_unfed >= 10) >= min_samples
counts_filtered <- counts_unfed[keep, ]

cat("  After filtering:", nrow(counts_filtered), "genes\n")
cat("  Removed:", nrow(counts_unfed) - nrow(counts_filtered), "genes\n")
cat("  Retention rate:", round(100 * nrow(counts_filtered) / nrow(counts_unfed), 2), "%\n")

################################################################################
# 3. Create DESeq2 Dataset
################################################################################

cat("\n\nCreating DESeq2 dataset...\n")
cat("  Design: ~ Photoperiod\n")
cat("  Reference level: LD (non-diapause)\n")
cat("  Contrast: SD (diapause-inducing) vs LD (non-diapause)\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = unfed_samples,
  design = ~ Photoperiod
)

# Ensure LD is reference level
dds$Photoperiod <- relevel(dds$Photoperiod, ref = "LD")

cat("DESeq2 dataset created successfully\n")
cat("  Genes:", nrow(dds), "\n")
cat("  Samples:", ncol(dds), "\n")

################################################################################
# 4. Run DESeq2
################################################################################

cat("\n\nRunning DESeq2 differential expression analysis...\n")
dds <- DESeq(dds)
cat("DESeq2 completed\n")

# Get results
res <- results(dds, alpha = 0.05)
cat("\n  Testing contrast: SD vs LD\n")
cat("  Significance level: FDR < 0.05\n")

# LFC shrinkage using apeglm (following Angela's method)
cat("\nApplying log2 fold change shrinkage (apeglm)...\n")
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="Photoperiod_SD_vs_LD", type="apeglm")
cat("LFC shrinkage completed\n")

# Summary statistics
cat("\nDESeq2 results summary (shrunk LFC):\n")
summary(res_LFC)

# Count significant genes (following Angela's thresholds)
# padj < 0.05 AND |log2FoldChange| > 0.58
sig_genes <- sum(res_LFC$padj < 0.05 & abs(res_LFC$log2FoldChange) > 0.58, na.rm = TRUE)
up_genes <- sum(res_LFC$padj < 0.05 & res_LFC$log2FoldChange > 0.58, na.rm = TRUE)
down_genes <- sum(res_LFC$padj < 0.05 & res_LFC$log2FoldChange < -0.58, na.rm = TRUE)

cat("\nSignificant genes (FDR < 0.05 AND |log2FC| > 0.58):\n")
cat("  Total:", sig_genes, "/", nrow(res_LFC), "(", round(100 * sig_genes / nrow(res_LFC), 2), "%)\n")
cat("  Up-regulated (SD > LD):", up_genes, "\n")
cat("  Down-regulated (SD < LD):", down_genes, "\n")

################################################################################
# 5. Save Results
################################################################################

cat("\n\nSaving results...\n")

# Convert to dataframe (using LFC-shrunk results)
res_df <- as.data.frame(res_LFC)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

# Order by p-value
res_df <- res_df[order(res_df$pvalue), ]

# Save all results
all_results_file <- file.path(output_dir, "adults_unfed_deseq2_all_LFCshrink.tsv")
write.table(res_df, all_results_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  All results (LFC shrunk):", all_results_file, "\n")

# Save significant genes only (following Angela's thresholds: padj < 0.05 AND |log2FC| > 0.58)
sig_results <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.58, ]
sig_results_file <- file.path(output_dir, "adults_unfed_deseq2_significant_LFCshrink.tsv")
write.table(sig_results, sig_results_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Significant genes (padj < 0.05, |log2FC| > 0.58):", sig_results_file, "\n")

# Save summary statistics
summary_stats <- data.frame(
  analysis = "adults_unfed",
  comparison = "SD_vs_LD",
  total_genes = nrow(res),
  sig_genes = sig_genes,
  up_regulated = up_genes,
  down_regulated = down_genes,
  percent_de = round(100 * sig_genes / nrow(res), 2)
)
summary_file <- file.path(output_dir, "adults_unfed_summary_stats.tsv")
write.table(summary_stats, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Summary statistics:", summary_file, "\n")

################################################################################
# 6. Quality Control Plots
################################################################################

cat("\n\nGenerating QC plots...\n")

# PCA plot
cat("  Creating PCA plot...\n")
vsd <- vst(dds, blind = TRUE)

pdf(file.path(output_dir, "adults_unfed_pca.pdf"), width = 8, height = 6)
plotPCA(vsd, intgroup = "Photoperiod") +
  ggtitle("PCA - Adult Females (Unfed Only)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
dev.off()

# Sample correlation heatmap
cat("  Creating correlation heatmap...\n")
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)

pdf(file.path(output_dir, "adults_unfed_correlation_heatmap.pdf"), width = 10, height = 10)
pheatmap(vsd_cor,
         main = "Sample Correlation - Adult Females (Unfed)",
         fontsize = 10)
dev.off()

# MA plot (both unshrunk and shrunk)
cat("  Creating MA plots...\n")
pdf(file.path(output_dir, "adults_unfed_ma_plot.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
plotMA(res, main = "MA Plot - Unshrunk\nAdult Females (Unfed) SD vs LD", ylim = c(-5, 5))
plotMA(res_LFC, main = "MA Plot - LFC Shrunk (apeglm)\nAdult Females (Unfed) SD vs LD", ylim = c(-5, 5))
par(mfrow = c(1, 1))
dev.off()

# Dispersion plot
cat("  Creating dispersion plot...\n")
pdf(file.path(output_dir, "adults_unfed_dispersion.pdf"), width = 8, height = 6)
plotDispEsts(dds, main = "Dispersion Estimates - Adult Females (Unfed)")
dev.off()

cat("\nQC plots saved to:", output_dir, "\n")

################################################################################
# 7. Session Info
################################################################################

cat("\n\n================================================================================\n")
cat("ANALYSIS COMPLETE - Unfed Females\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat("  - Tested", nrow(res), "genes\n")
cat("  - Found", sig_genes, "significantly DE genes (FDR < 0.05)\n")
cat("  - Up-regulated:", up_genes, "genes\n")
cat("  - Down-regulated:", down_genes, "genes\n")
cat("  - Results saved to:", output_dir, "\n\n")

# Save session info
sink(file.path(output_dir, "adults_unfed_session_info.txt"))
cat("Adult Unfed DESeq2 Analysis\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()

cat("Session info saved\n")
cat("================================================================================\n")
