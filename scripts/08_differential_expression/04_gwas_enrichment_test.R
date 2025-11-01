#!/usr/bin/env Rscript
#
# GWAS Candidate Enrichment Test
#
# Purpose: Test if GWAS candidates are enriched among differentially expressed genes
#          using Fisher's exact test (as requested by Reviewer 2)
#
# Input:
#   - results/differential_expression/adults/deseq2_results_all.tsv
#   - data/metadata/gwas_candidates_final_34.txt
#
# Output:
#   - results/gwas_validation/enrichment_test_results.txt
#   - results/gwas_validation/candidate_expression_summary.tsv
#   - results/gwas_validation/enrichment_figure.pdf
#
# Author: LCosme
# Date: October 21, 2025

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("===============================================================\n")
cat("GWAS Candidate Enrichment Test\n")
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
  cat("\n")
  cat("ERROR: GWAS candidate file not found!\n")
  cat("Expected file:", gwas_file, "\n\n")
  cat("Please create this file with one gene ID per line.\n")
  cat("Example format:\n")
  cat("  gene-LOC109400916\n")
  cat("  gene-LOC109401234\n")
  cat("  ...\n\n")
  quit(status = 1)
}

candidates <- read.table(gwas_file, header = FALSE, stringsAsFactors = FALSE)
colnames(candidates) <- "gene_id"

cat("  Total GWAS candidates:", nrow(candidates), "\n\n")

# ============================================================================
# 2. Load DE results for all stages
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

  # Extract candidate expression
  candidate_expr[[stage]] <- de_results[[stage]] %>%
    filter(gene_id %in% candidates$gene_id) %>%
    mutate(stage = stage)

  cat(sprintf("  %s: %d genes tested, %d candidates found\n",
    tools::toTitleCase(stage),
    nrow(de_results[[stage]]),
    nrow(candidate_expr[[stage]])
  ))
}

cat("\n")

# ============================================================================
# 3. Fisher's Exact Test for Each Stage
# ============================================================================

cat("3. Performing Fisher's exact test for enrichment...\n\n")

enrichment_results <- data.frame(
  stage = character(),
  genes_tested = integer(),
  genome_de = integer(),
  genome_de_pct = numeric(),
  candidates_total = integer(),
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

  if (nrow(cand_stage) == 0) {
    cat("  No candidates found, skipping\n\n")
    next
  }

  # Count DEGs genome-wide
  genome_tested <- nrow(de_results[[stage]])
  genome_de <- sum(de_results[[stage]]$padj < 0.05 & !is.na(de_results[[stage]]$padj))
  genome_de_pct <- 100 * genome_de / genome_tested

  # Count DEGs among candidates
  candidates_total <- nrow(cand_stage)
  candidates_de <- sum(cand_stage$padj < 0.05 & !is.na(cand_stage$padj))
  candidates_de_pct <- 100 * candidates_de / candidates_total

  cat(sprintf("  Genome-wide: %d / %d DEGs (%.2f%%)\n",
    genome_de, genome_tested, genome_de_pct
  ))
  cat(sprintf("  Candidates: %d / %d DEGs (%.2f%%)\n",
    candidates_de, candidates_total, candidates_de_pct
  ))

  # Fisher's exact test
  # Contingency table:
  #                 DE    Not DE
  # Candidates      a     b
  # Non-candidates  c     d

  a <- candidates_de
  b <- candidates_total - candidates_de
  c <- genome_de - candidates_de # DEGs that are NOT candidates
  d <- (genome_tested - genome_de) - (candidates_total - candidates_de)

  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

  fisher_result <- fisher.test(contingency_table, alternative = "greater")

  cat(sprintf("  Odds ratio: %.2f\n", fisher_result$estimate))
  cat(sprintf("  P-value: %.2e\n", fisher_result$p.value))
  cat(sprintf("  Significant: %s\n",
    ifelse(fisher_result$p.value < 0.05, "YES", "NO")
  ))
  cat("\n")

  # Add to results
  enrichment_results <- rbind(enrichment_results, data.frame(
    stage = stage,
    genes_tested = genome_tested,
    genome_de = genome_de,
    genome_de_pct = genome_de_pct,
    candidates_total = candidates_total,
    candidates_de = candidates_de,
    candidates_de_pct = candidates_de_pct,
    odds_ratio = as.numeric(fisher_result$estimate),
    p_value = fisher_result$p.value,
    significant = ifelse(fisher_result$p.value < 0.05, "YES", "NO")
  ))
}

# ============================================================================
# 4. Save results
# ============================================================================

cat("4. Saving results...\n")

# Save enrichment test results
enrichment_file <- file.path(results_dir, "enrichment_test_results.txt")
sink(enrichment_file)

cat("GWAS Candidate Enrichment Test Results\n")
cat("=====================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Test: Fisher's Exact Test (one-sided, testing for enrichment)\n")
cat("Null Hypothesis: GWAS candidates are NOT enriched among DEGs\n")
cat("Alternative: GWAS candidates ARE enriched among DEGs\n\n")

cat("Results by Stage:\n")
cat("=================\n\n")

for (i in 1:nrow(enrichment_results)) {
  row <- enrichment_results[i, ]

  cat(sprintf("%s:\n", tools::toTitleCase(row$stage)))
  cat(sprintf("  Genes tested: %d\n", row$genes_tested))
  cat(sprintf("  Genome-wide DEGs: %d (%.2f%%)\n",
    row$genome_de, row$genome_de_pct
  ))
  cat(sprintf("  GWAS candidates: %d\n", row$candidates_total))
  cat(sprintf("  Candidates DE: %d (%.2f%%)\n",
    row$candidates_de, row$candidates_de_pct
  ))
  cat(sprintf("  Odds Ratio: %.2f\n", row$odds_ratio))
  cat(sprintf("  P-value: %.2e\n", row$p_value))
  cat(sprintf("  Significant at α=0.05: %s\n", row$significant))
  cat("\n")
}

cat("Interpretation:\n")
cat("===============\n")
cat("An odds ratio > 1 indicates enrichment (candidates more likely to be DE)\n")
cat("An odds ratio < 1 indicates depletion (candidates less likely to be DE)\n")
cat("P-value < 0.05 indicates statistically significant enrichment\n\n")

sink()

cat("  Saved:", basename(enrichment_file), "\n")

# Save enrichment table
table_file <- file.path(results_dir, "enrichment_summary_table.tsv")
write.table(enrichment_results, table_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved:", basename(table_file), "\n")

# Save candidate expression across stages
combined_candidates <- bind_rows(candidate_expr)
cand_file <- file.path(results_dir, "candidate_expression_all_stages.tsv")
write.table(combined_candidates, cand_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved:", basename(cand_file), "\n\n")

# ============================================================================
# 5. Summary
# ============================================================================

cat("===============================================================\n")
cat("Enrichment Test Complete!\n")
cat("===============================================================\n\n")

cat("Key Findings:\n")
for (i in 1:nrow(enrichment_results)) {
  row <- enrichment_results[i, ]
  cat(sprintf("  %s: %d/%d candidates DE (OR=%.2f, p=%.2e) - %s\n",
    tools::toTitleCase(row$stage),
    row$candidates_de,
    row$candidates_total,
    row$odds_ratio,
    row$p_value,
    row$significant
  ))
}

cat("\nOutput files:\n")
cat("  ", basename(enrichment_file), "\n")
cat("  ", basename(table_file), "\n")
cat("  ", basename(cand_file), "\n\n")

cat("For Reviewer Response:\n")
if (nrow(enrichment_results) > 0 && enrichment_results$stage[1] == "adults") {
  adult_row <- enrichment_results[1, ]
  cat(sprintf(
    '  "We tested %d genes in adult females and found %d (%.1f%%) significantly\n   differentially expressed at FDR < 0.05. Among our %d GWAS candidates,\n   %d (%.1f%%) were differentially expressed, representing %s enrichment\n   (Fisher\'s exact test, OR = %.2f, p = %.2e)."\n',
    adult_row$genes_tested,
    adult_row$genome_de,
    adult_row$genome_de_pct,
    adult_row$candidates_total,
    adult_row$candidates_de,
    adult_row$candidates_de_pct,
    ifelse(adult_row$significant == "YES", "significant", "non-significant"),
    adult_row$odds_ratio,
    adult_row$p_value
  ))
}

cat("\n===============================================================\n")