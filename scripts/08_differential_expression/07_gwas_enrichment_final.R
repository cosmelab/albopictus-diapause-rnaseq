#!/usr/bin/env Rscript
#
# GWAS Candidate Enrichment Test - WORKING VERSION
#
# Purpose: Test if GWAS candidates are enriched among differentially expressed genes
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
  cat("ERROR: GWAS candidate file not found!\n")
  quit(status = 1)
}

# Read and clean the file
candidates_raw <- readLines(gwas_file)
candidates_clean <- candidates_raw[candidates_raw != "" & !is.na(candidates_raw)]

candidates <- data.frame(
  gene_id_original = candidates_clean,
  gene_id_with_prefix = paste0("gene-", candidates_clean),
  stringsAsFactors = FALSE
)

cat(sprintf("  Total GWAS candidates: %d\n", nrow(candidates)))
cat(sprintf("  First few IDs: %s\n\n", paste(head(candidates$gene_id_original, 3), collapse=", ")))

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

  # Match candidates - try with gene- prefix
  candidate_expr[[stage]] <- de_results[[stage]] %>%
    filter(gene_id %in% candidates$gene_id_with_prefix) %>%
    mutate(stage = stage)

  cat(sprintf("  %s: %d genes tested, %d candidates found\n",
    tools::toTitleCase(stage),
    nrow(de_results[[stage]]),
    nrow(candidate_expr[[stage]])
  ))

  # Show which candidates were found
  if (nrow(candidate_expr[[stage]]) > 0) {
    found_ids <- gsub("gene-", "", candidate_expr[[stage]]$gene_id[1:min(3, nrow(candidate_expr[[stage]]))])
    cat(sprintf("    Found: %s...\n", paste(found_ids, collapse=", ")))
  }
}

cat("\n")

# ============================================================================
# 3. Fisher's Exact Test
# ============================================================================

cat("3. Performing Fisher's exact test for enrichment...\n\n")

enrichment_results <- data.frame()

for (stage in names(de_results)) {
  cat(sprintf("Testing %s...\n", tools::toTitleCase(stage)))

  cand_stage <- candidate_expr[[stage]]
  candidates_found <- nrow(cand_stage)

  if (candidates_found == 0) {
    cat("  No candidates found in expression data\n\n")
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

  # Store results
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

  cat("\n")
}

# ============================================================================
# 4. Save results
# ============================================================================

cat("4. Saving results...\n")

# Save enrichment summary
write.table(enrichment_results,
  file.path(results_dir, "enrichment_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

# Save candidate expression
if (length(candidate_expr) > 0) {
  combined_candidates <- bind_rows(candidate_expr)
  if (nrow(combined_candidates) > 0) {
    write.table(combined_candidates,
      file.path(results_dir, "candidate_expression.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

cat("  Saved enrichment_summary.tsv\n")
cat("  Saved candidate_expression.tsv\n")

cat("\n===============================================================\n")
cat("Analysis Complete!\n")
cat("===============================================================\n\n")

if (nrow(enrichment_results) > 0) {
  cat("Summary:\n")
  for (i in 1:nrow(enrichment_results)) {
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

cat("\n")
