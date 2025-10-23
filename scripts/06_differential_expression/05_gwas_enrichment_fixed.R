#!/usr/bin/env Rscript
#
# GWAS Candidate Enrichment Test (Fixed for ID matching)
#
# Purpose: Test if GWAS candidates are enriched among differentially expressed genes
#          Handles different ID formats (LOC vs gene-LOC)
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("===============================================================\n")
cat("GWAS Candidate Enrichment Test (Fixed)\n")
cat("===============================================================\n\n")

# Paths
project_base <- "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
de_dir <- file.path(project_base, "results/differential_expression")
metadata_dir <- file.path(project_base, "data/metadata")
results_dir <- file.path(project_base, "results/gwas_validation")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load GWAS candidates
# ============================================================================

cat("1. Loading GWAS candidate list...\n")
gwas_file <- file.path(metadata_dir, "gwas_candidates_final_34.txt")

if (!file.exists(gwas_file)) {
  cat("ERROR: GWAS candidate file not found!\n")
  quit(status = 1)
}

candidates <- read.table(gwas_file, header = FALSE, stringsAsFactors = FALSE)
colnames(candidates) <- "gene_id"

# Add both formats for matching
candidates$gene_id_with_prefix <- paste0("gene-", candidates$gene_id)
candidates$gene_id_original <- candidates$gene_id

cat("  Total GWAS candidates:", nrow(candidates), "\n")
cat("  First few IDs:", paste(head(candidates$gene_id, 3), collapse=", "), "\n\n")

# ============================================================================
# 2. Load DE results and match candidates
# ============================================================================

cat("2. Loading differential expression results...\n")

stages <- c("adults", "embryos", "larvae")
de_results <- list()
candidate_expr <- list()

for (stage in stages) {
  results_file <- file.path(de_dir, stage, "deseq2_results_all.tsv")

  if (!file.exists(results_file)) {
    cat(sprintf("  WARNING: Results not found for %s\n", stage))
    next
  }

  # Load DE results
  de_results[[stage]] <- read.delim(results_file, stringsAsFactors = FALSE)

  # Try to match candidates with both formats
  candidate_expr[[stage]] <- de_results[[stage]] %>%
    filter(gene_id %in% candidates$gene_id_with_prefix |
           gene_id %in% candidates$gene_id_original |
           gsub("gene-", "", gene_id) %in% candidates$gene_id_original) %>%
    mutate(stage = stage)

  cat(sprintf("  %s: %d genes tested, %d candidates found\n",
    tools::toTitleCase(stage),
    nrow(de_results[[stage]]),
    nrow(candidate_expr[[stage]])
  ))

  # Show which candidates were found
  if (nrow(candidate_expr[[stage]]) > 0) {
    cat("    Found:", paste(head(candidate_expr[[stage]]$gene_id, 3), collapse=", "), "...\n")
  }
}

cat("\n")

# ============================================================================
# 3. Fisher's Exact Test
# ============================================================================

cat("3. Performing Fisher's exact test for enrichment...\n\n")

enrichment_results <- data.frame(
  stage = character(),
  genes_tested = integer(),
  genome_de = integer(),
  genome_de_pct = numeric(),
  candidates_total = integer(),
  candidates_found = integer(),
  candidates_de = integer(),
  candidates_de_pct = numeric(),
  odds_ratio = numeric(),
  p_value = numeric(),
  significant = character(),
  stringsAsFactors = FALSE
)

for (stage in names(de_results)) {
  cat(sprintf("Testing %s...\n", tools::toTitleCase(stage)))

  # Get candidate expression for this stage
  cand_stage <- candidate_expr[[stage]]

  # Count how many candidates we actually found in the data
  candidates_found <- nrow(cand_stage)

  if (candidates_found == 0) {
    cat("  No candidates found in expression data\n")
    cat("  Check gene ID format mismatch\n\n")

    # Add NA row for this stage
    enrichment_results <- rbind(enrichment_results, data.frame(
      stage = stage,
      genes_tested = nrow(de_results[[stage]]),
      genome_de = sum(de_results[[stage]]$padj < 0.05 & !is.na(de_results[[stage]]$padj)),
      genome_de_pct = 100 * sum(de_results[[stage]]$padj < 0.05 & !is.na(de_results[[stage]]$padj)) / nrow(de_results[[stage]]),
      candidates_total = nrow(candidates),
      candidates_found = 0,
      candidates_de = NA,
      candidates_de_pct = NA,
      odds_ratio = NA,
      p_value = NA,
      significant = "No candidates found"
    ))
    next
  }

  # Count DEGs genome-wide
  genome_tested <- nrow(de_results[[stage]])
  genome_de <- sum(de_results[[stage]]$padj < 0.05 & !is.na(de_results[[stage]]$padj))
  genome_de_pct <- 100 * genome_de / genome_tested

  # Count DEGs among candidates
  candidates_de <- sum(cand_stage$padj < 0.05 & !is.na(cand_stage$padj))
  candidates_de_pct <- 100 * candidates_de / candidates_found

  cat(sprintf("  Genome-wide: %d / %d DEGs (%.2f%%)\n",
    genome_de, genome_tested, genome_de_pct
  ))
  cat(sprintf("  Candidates: %d / %d DEGs (%.2f%%)\n",
    candidates_de, candidates_found, candidates_de_pct
  ))

  # Fisher's exact test
  if (candidates_found > 0) {
    a <- candidates_de
    b <- candidates_found - candidates_de
    c <- genome_de - candidates_de
    d <- (genome_tested - genome_de) - (candidates_found - candidates_de)

    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(contingency_table, alternative = "greater")

    cat(sprintf("  Odds ratio: %.2f\n", fisher_result$estimate))
    cat(sprintf("  P-value: %.2e\n", fisher_result$p.value))
    cat(sprintf("  Significant: %s\n",
      ifelse(fisher_result$p.value < 0.05, "YES", "NO")
    ))

    # Add to results
    enrichment_results <- rbind(enrichment_results, data.frame(
      stage = stage,
      genes_tested = genome_tested,
      genome_de = genome_de,
      genome_de_pct = genome_de_pct,
      candidates_total = nrow(candidates),
      candidates_found = candidates_found,
      candidates_de = candidates_de,
      candidates_de_pct = candidates_de_pct,
      odds_ratio = as.numeric(fisher_result$estimate),
      p_value = fisher_result$p.value,
      significant = ifelse(fisher_result$p.value < 0.05, "YES", "NO")
    ))
  }

  cat("\n")
}

# ============================================================================
# 4. Apply collaborator's threshold (|log2FC| > 0.5)
# ============================================================================

cat("4. Testing with collaborator's threshold (padj < 0.05 & |log2FC| > 0.5)...\n\n")

for (stage in names(de_results)) {
  cat(sprintf("%s:\n", tools::toTitleCase(stage)))

  # Count with FC threshold
  genome_de_fc <- sum(de_results[[stage]]$padj < 0.05 &
                      abs(de_results[[stage]]$log2FoldChange) > 0.5 &
                      !is.na(de_results[[stage]]$padj), na.rm = TRUE)

  if (nrow(candidate_expr[[stage]]) > 0) {
    candidates_de_fc <- sum(candidate_expr[[stage]]$padj < 0.05 &
                            abs(candidate_expr[[stage]]$log2FoldChange) > 0.5 &
                            !is.na(candidate_expr[[stage]]$padj), na.rm = TRUE)

    cat(sprintf("  Genome-wide with FC cutoff: %d / %d (%.2f%%)\n",
      genome_de_fc, nrow(de_results[[stage]]),
      100 * genome_de_fc / nrow(de_results[[stage]])
    ))
    cat(sprintf("  Candidates with FC cutoff: %d / %d (%.2f%%)\n",
      candidates_de_fc, nrow(candidate_expr[[stage]]),
      100 * candidates_de_fc / nrow(candidate_expr[[stage]])
    ))
  } else {
    cat("  No candidates found\n")
  }
  cat("\n")
}

# ============================================================================
# 5. Save results
# ============================================================================

cat("5. Saving results...\n")

# Save enrichment results
enrichment_file <- file.path(results_dir, "enrichment_test_results_fixed.txt")
sink(enrichment_file)

cat("GWAS Candidate Enrichment Test Results\n")
cat("=====================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Gene ID Matching:\n")
cat("  Original candidate format: LOC109...\n")
cat("  DE results format: gene-LOC109...\n")
cat("  Matching attempted with both formats\n\n")

cat("Results by Stage:\n")
cat("=================\n\n")

for (i in 1:nrow(enrichment_results)) {
  row <- enrichment_results[i, ]

  cat(sprintf("%s:\n", tools::toTitleCase(row$stage)))
  cat(sprintf("  Genes tested: %d\n", row$genes_tested))
  cat(sprintf("  Genome-wide DEGs: %d (%.2f%%)\n", row$genome_de, row$genome_de_pct))
  cat(sprintf("  GWAS candidates: %d total, %d found in data\n",
    row$candidates_total, row$candidates_found))

  if (!is.na(row$candidates_de)) {
    cat(sprintf("  Candidates DE: %d (%.2f%%)\n", row$candidates_de, row$candidates_de_pct))
    cat(sprintf("  Odds Ratio: %.2f\n", row$odds_ratio))
    cat(sprintf("  P-value: %.2e\n", row$p_value))
    cat(sprintf("  Significant at α=0.05: %s\n", row$significant))
  } else {
    cat("  No enrichment test performed (no candidates found)\n")
  }
  cat("\n")
}

sink()

cat("  Saved:", basename(enrichment_file), "\n")

# Save tables
write.table(enrichment_results,
  file.path(results_dir, "enrichment_summary_fixed.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

# Save candidate expression
if (length(candidate_expr) > 0) {
  combined_candidates <- bind_rows(candidate_expr)
  if (nrow(combined_candidates) > 0) {
    write.table(combined_candidates,
      file.path(results_dir, "candidate_expression_fixed.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

cat("\n===============================================================\n")
cat("Analysis Complete!\n")
cat("===============================================================\n\n")

# Summary
if (sum(sapply(candidate_expr, nrow)) == 0) {
  cat("⚠ WARNING: No GWAS candidates found in expression data!\n")
  cat("This indicates a gene ID format mismatch.\n")
  cat("Check that gene IDs match between:\n")
  cat("  - GWAS candidates file\n")
  cat("  - GTF annotation used for quantification\n")
  cat("  - Count matrices\n")
} else {
  cat("Key Findings:\n")
  for (i in 1:nrow(enrichment_results)) {
    if (enrichment_results$candidates_found[i] > 0) {
      cat(sprintf("  %s: %d/%d candidates DE (OR=%.2f, p=%.2e) - %s\n",
        tools::toTitleCase(enrichment_results$stage[i]),
        enrichment_results$candidates_de[i],
        enrichment_results$candidates_found[i],
        enrichment_results$odds_ratio[i],
        enrichment_results$p_value[i],
        enrichment_results$significant[i]
      ))
    }
  }
}

cat("\n")