#!/usr/bin/env Rscript

################################################################################
# DESeq2 Differential Expression Analysis
#
# Purpose: Perform differential expression analysis on batch-corrected counts
# Input: ComBat-seq corrected count matrix
# Output: DE results, normalized counts, visualizations
#
# Author: RNA-seq Analysis Pipeline
# Date: October 20, 2025
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(ggrepel)
})

cat("====================================================================== \n")
cat("DESeq2 Differential Expression Analysis\n")
cat("====================================================================== \n\n")

# Set paths
base_dir <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
setwd(base_dir)

# Input files
count_file <- "results/count_matrices/salmon_gene_counts_combat_corrected.tsv"
metadata_file <- "data/metadata/samples.csv"

# Output directories
results_dir <- "results/differential_expression"
figures_dir <- "results/de_figures"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading batch-corrected count data...\n")
counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat(sprintf("  Loaded: %d genes × %d samples\n", nrow(counts), ncol(counts)))

# Load metadata
cat("\nLoading sample metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1)
metadata$Condition <- factor(metadata$Condition,
                             levels = c("non-diapause-inducing", "diapause-inducing"))
cat(sprintf("  Loaded metadata for %d samples\n", nrow(metadata)))
cat(sprintf("  Conditions: %s\n", paste(levels(metadata$Condition), collapse = ", ")))

# Ensure sample order matches
metadata <- metadata[colnames(counts), ]
stopifnot(all(rownames(metadata) == colnames(counts)))

# Create DESeq2 dataset
cat("\n====================================================================== \n")
cat("Running DESeq2 Analysis\n")
cat("====================================================================== \n\n")

# Since counts are batch-corrected, use simple design
cat("Creating DESeq2 dataset with design: ~ Condition\n")
cat("(Platform batch effects already removed by ComBat-seq)\n\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),  # Ensure integer counts
  colData = metadata,
  design = ~ Condition
)

# Filter low-count genes (optional - may already be filtered)
cat("Filtering genes...\n")
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
cat(sprintf("  Genes after filtering: %d\n", nrow(dds)))

# Run DESeq2
cat("\nRunning DESeq2 normalization and testing...\n")
dds <- DESeq(dds)

# Extract results
cat("\nExtracting differential expression results...\n")
res <- results(dds,
               contrast = c("Condition", "diapause-inducing", "non-diapause-inducing"),
               alpha = 0.05)

# Add gene symbols if available (placeholder for now)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Order by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Summary statistics
cat("\n====================================================================== \n")
cat("Differential Expression Summary\n")
cat("====================================================================== \n\n")

sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("Total genes tested: %d\n", sum(!is.na(res$padj))))
cat(sprintf("Significantly DE genes (padj < 0.05): %d\n", sig_genes))
cat(sprintf("  - Upregulated in diapause: %d\n", sig_up))
cat(sprintf("  - Downregulated in diapause: %d\n", sig_down))
cat("\n")

# More stringent cutoffs
sig_genes_fc <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)
cat(sprintf("DE genes with |log2FC| > 1: %d\n", sig_genes_fc))

sig_genes_strict <- sum(res$padj < 0.01, na.rm = TRUE)
cat(sprintf("DE genes with padj < 0.01: %d\n", sig_genes_strict))

# Save results
cat("\n====================================================================== \n")
cat("Saving Results\n")
cat("====================================================================== \n\n")

# Full results
output_file <- file.path(results_dir, "deseq2_results_all_genes.tsv")
write.table(res_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved all results to: %s\n", output_file))

# Significant genes only
sig_df <- res_df[which(res_df$padj < 0.05), ]
output_file <- file.path(results_dir, "deseq2_results_significant.tsv")
write.table(sig_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved significant genes to: %s\n", output_file))

# Top 100 genes
top_df <- head(res_df, 100)
output_file <- file.path(results_dir, "deseq2_results_top100.tsv")
write.table(top_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved top 100 genes to: %s\n", output_file))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
output_file <- file.path(results_dir, "normalized_counts.tsv")
write.table(norm_counts, output_file, sep = "\t", quote = FALSE)
cat(sprintf("Saved normalized counts to: %s\n", output_file))

# VST transformed counts (for visualization)
vst_counts <- assay(vst(dds, blind = FALSE))
output_file <- file.path(results_dir, "vst_transformed_counts.tsv")
write.table(vst_counts, output_file, sep = "\t", quote = FALSE)
cat(sprintf("Saved VST counts to: %s\n", output_file))

# Save DESeq2 object
output_file <- file.path(results_dir, "deseq2_dds_object.rds")
saveRDS(dds, output_file)
cat(sprintf("Saved DESeq2 object to: %s\n", output_file))

cat("\n====================================================================== \n")
cat("Creating Visualizations\n")
cat("====================================================================== \n\n")

# 1. MA Plot
cat("Creating MA plot...\n")
pdf(file.path(figures_dir, "ma_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), alpha = 0.05,
       main = "MA Plot: Diapause vs Non-Diapause",
       cex = 0.5)
dev.off()

# 2. Volcano Plot
cat("Creating volcano plot...\n")
pdf(file.path(figures_dir, "volcano_plot.pdf"), width = 10, height = 8)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Diapause vs Non-Diapause',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colAlpha = 0.8,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)
dev.off()

# 3. PCA of VST counts
cat("Creating PCA plot...\n")
vst_data <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vst_data, intgroup = "Condition", returnData = TRUE)

pdf(file.path(figures_dir, "pca_vst_counts.pdf"), width = 8, height = 6)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples (After Batch Correction)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 4. Heatmap of top DE genes
cat("Creating heatmap of top 50 DE genes...\n")
top_genes <- head(res_df$gene_id[!is.na(res_df$padj)], 50)
top_vst <- vst_counts[top_genes, ]

# Scale by row for visualization
top_vst_scaled <- t(scale(t(top_vst)))

# Create annotation for samples
annotation_col <- data.frame(
  Condition = metadata$Condition,
  row.names = rownames(metadata)
)

# Color palette
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)

pdf(file.path(figures_dir, "heatmap_top50_genes.pdf"), width = 10, height = 12)
pheatmap(top_vst_scaled,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         annotation_col = annotation_col,
         color = colors,
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "Top 50 Differentially Expressed Genes",
         fontsize_col = 8)
dev.off()

# 5. Sample-to-sample distance heatmap
cat("Creating sample distance heatmap...\n")
sample_dists <- dist(t(vst_counts))
sample_dist_matrix <- as.matrix(sample_dists)

pdf(file.path(figures_dir, "sample_distance_heatmap.pdf"), width = 10, height = 8)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Sample-to-Sample Distances")
dev.off()

cat("\nAll visualizations saved to:", figures_dir, "\n")

# Final report
cat("\n====================================================================== \n")
cat("Analysis Complete!\n")
cat("====================================================================== \n\n")
cat("Key outputs:\n")
cat("  - DE results: results/differential_expression/deseq2_results_*.tsv\n")
cat("  - Normalized counts: results/differential_expression/normalized_counts.tsv\n")
cat("  - Visualizations: results/de_figures/*.pdf\n")
cat("\nNext steps:\n")
cat("  1. Review significant genes\n")
cat("  2. Check GWAS candidate genes in results\n")
cat("  3. Perform GO enrichment analysis\n")
cat("  4. Generate custom visualizations for manuscript\n")

# Session info
cat("\n====================================================================== \n")
cat("Session Information\n")
cat("====================================================================== \n\n")
sessionInfo()