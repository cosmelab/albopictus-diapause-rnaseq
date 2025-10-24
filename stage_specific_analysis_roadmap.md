# Stage-Specific Analysis Roadmap

**Created:** October 21, 2025
**Purpose:** Step-by-step guide for analyzing each developmental stage independently

---

## Overview

Due to complete confounding (Platform = Stage = Study), we must:
1. Analyze each developmental stage completely independently
2. Use meta-analysis to combine results
3. Focus on adults for primary GWAS validation

---

## Phase 1: Data Preparation ✅ COMPLETED

### Files Already Split by Stage
```
results/count_matrices/
├── adults_counts.tsv    # 16 samples (PRJNA268379)
├── embryos_counts.tsv   # 12 samples (PRJNA158021)
└── larvae_counts.tsv     # 17 samples (PRJNA187045)
```

### Next: Create Stage-Specific Metadata

#### Step 1.1: Split Metadata File
```python
# Script: scripts/06_differential_expression/01_split_metadata.py
import pandas as pd

# Read full metadata
metadata = pd.read_csv('data/metadata/samples.csv')

# Split by developmental stage
adults_meta = metadata[metadata['Stage'] == 'Adult']
embryos_meta = metadata[metadata['Stage'] == 'Embryo']
larvae_meta = metadata[metadata['Stage'] == 'Pharate_Larvae']

# Save separately
adults_meta.to_csv('data/metadata/samples_adults.csv', index=False)
embryos_meta.to_csv('data/metadata/samples_embryos.csv', index=False)
larvae_meta.to_csv('data/metadata/samples_larvae.csv', index=False)
```

---

## Phase 2: Adults Analysis (PRIORITY - Most relevant for GWAS)

### Step 2.1: QC for Adults Only
```python
# Script: scripts/06_differential_expression/02_adults_qc.py
# - PCA plot (16 samples, colored by condition)
# - Sample correlation heatmap
# - Gene filtering (low counts)
# - Check for outliers
```

### Step 2.2: DESeq2 for Adults
```r
# Script: scripts/06_differential_expression/03_adults_deseq2.R
library(DESeq2)

# Load data
counts <- read.table('results/count_matrices/adults_counts.tsv', header=TRUE, row.names=1)
metadata <- read.csv('data/metadata/samples_adults.csv')

# Create DESeq object - SIMPLE DESIGN, NO BATCH
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Condition  # Just condition, no platform covariate
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 4  # At least 10 counts in 25% of samples
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results
results <- results(dds, contrast=c("Condition", "diapause", "non_diapause"))

# Save all results
write.csv(results, 'results/differential_expression/adults_all_genes.csv')
```

### Step 2.3: Extract GWAS Candidates
```r
# Script: scripts/06_differential_expression/04_adults_gwas_candidates.R

# List of 34 GWAS candidate gene IDs
gwas_candidates <- c("LOC109431832", "LOC109432456", ...)  # Full list

# Extract results for candidates
candidate_results <- results[rownames(results) %in% gwas_candidates, ]

# Create summary table
summary_table <- data.frame(
  Gene = rownames(candidate_results),
  log2FC = candidate_results$log2FoldChange,
  pvalue = candidate_results$pvalue,
  padj = candidate_results$padj,
  Significant = candidate_results$padj < 0.05
)

write.csv(summary_table, 'results/differential_expression/adults_gwas_candidates.csv')
```

### Step 2.4: Visualizations for Adults
```r
# Script: scripts/06_differential_expression/05_adults_visualization.R

# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(results,
  lab = rownames(results),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Adults: Diapause vs Non-Diapause',
  selectLab = gwas_candidates  # Highlight GWAS candidates
)

# MA plot
plotMA(results, main="Adults MA Plot")

# Heatmap of top DEGs
top_genes <- head(order(results$padj), 50)
pheatmap(assay(vst(dds))[top_genes,])
```

---

## Phase 3: Embryos Analysis

### Step 3.1: QC for Embryos
```python
# Script: scripts/06_differential_expression/06_embryos_qc.py
# Same structure as adults QC
```

### Step 3.2: DESeq2 for Embryos
```r
# Script: scripts/06_differential_expression/07_embryos_deseq2.R
# Identical structure to adults, but with embryo data
# design = ~ Condition  # No batch correction
```

### Step 3.3: Extract GWAS Candidates (Embryos)
```r
# Script: scripts/06_differential_expression/08_embryos_gwas_candidates.R
# Check which GWAS candidates replicate in embryos
```

---

## Phase 4: Larvae Analysis

### Step 4.1: QC for Larvae
```python
# Script: scripts/06_differential_expression/09_larvae_qc.py
```

### Step 4.2: DESeq2 for Larvae
```r
# Script: scripts/06_differential_expression/10_larvae_deseq2.R
# design = ~ Condition  # No batch correction
```

### Step 4.3: Extract GWAS Candidates (Larvae)
```r
# Script: scripts/06_differential_expression/11_larvae_gwas_candidates.R
```

---

## Phase 5: Meta-Analysis

### Step 5.1: Combine P-values
```r
# Script: scripts/06_differential_expression/12_meta_analysis.R
library(metap)

# Load results from all three stages
adults_res <- read.csv('results/differential_expression/adults_all_genes.csv')
embryos_res <- read.csv('results/differential_expression/embryos_all_genes.csv')
larvae_res <- read.csv('results/differential_expression/larvae_all_genes.csv')

# Get common genes
common_genes <- Reduce(intersect, list(
  adults_res$Gene,
  embryos_res$Gene,
  larvae_res$Gene
))

# For each common gene
meta_results <- data.frame()
for (gene in common_genes) {
  p_vals <- c(
    adults_res[adults_res$Gene == gene, 'pvalue'],
    embryos_res[embryos_res$Gene == gene, 'pvalue'],
    larvae_res[larvae_res$Gene == gene, 'pvalue']
  )

  # Fisher's method
  fisher_p <- sumlog(p_vals)$p

  # Stouffer's weighted method (by sample size)
  weights <- c(16, 12, 17)  # Sample sizes
  stouffer_p <- sumz(p_vals, weights=weights)$p

  meta_results <- rbind(meta_results, data.frame(
    Gene = gene,
    fisher_p = fisher_p,
    stouffer_p = stouffer_p
  ))
}
```

### Step 5.2: Cross-Stage Summary
```r
# Script: scripts/06_differential_expression/13_cross_stage_summary.R

# Create summary of GWAS candidates across stages
gwas_summary <- data.frame(
  Gene = gwas_candidates,
  Adults_log2FC = adults_gwas$log2FC,
  Adults_padj = adults_gwas$padj,
  Embryos_log2FC = embryos_gwas$log2FC,
  Embryos_padj = embryos_gwas$padj,
  Larvae_log2FC = larvae_gwas$log2FC,
  Larvae_padj = larvae_gwas$padj,
  Consistent_Direction = sign(Adults_log2FC) == sign(Embryos_log2FC) &
                         sign(Adults_log2FC) == sign(Larvae_log2FC),
  Stages_Significant = sum(c(Adults_padj < 0.05,
                             Embryos_padj < 0.05,
                             Larvae_padj < 0.05))
)

write.csv(gwas_summary, 'results/differential_expression/gwas_cross_stage_summary.csv')
```

---

## Phase 6: Final Visualizations

### Step 6.1: Cross-Stage Heatmap
```r
# Script: scripts/06_differential_expression/14_cross_stage_heatmap.R
# Heatmap showing log2FC for GWAS candidates across all 3 stages
```

### Step 6.2: Summary Figure
```r
# Script: scripts/06_differential_expression/15_summary_figure.R
# Multi-panel figure:
# A. PCA plots (3 panels, one per stage)
# B. Volcano plots (3 panels)
# C. GWAS candidate heatmap
# D. Venn diagram of DEGs across stages
```

---

## Expected Outputs

### Per Stage
```
results/differential_expression/
├── adults_all_genes.csv           # All gene results
├── adults_gwas_candidates.csv     # 34 GWAS candidates
├── adults_volcano.pdf
├── adults_ma_plot.pdf
├── adults_pca.pdf
├── adults_top50_heatmap.pdf
├── embryos_all_genes.csv
├── embryos_gwas_candidates.csv
├── [embryo plots...]
├── larvae_all_genes.csv
├── larvae_gwas_candidates.csv
└── [larvae plots...]
```

### Cross-Stage
```
results/differential_expression/
├── meta_analysis_fisher.csv
├── meta_analysis_stouffer.csv
├── gwas_cross_stage_summary.csv
├── cross_stage_heatmap.pdf
└── summary_figure.pdf
```

---

## Timeline

### Day 1 (Today)
- [ ] Split metadata files
- [ ] Run adults QC
- [ ] Run adults DESeq2
- [ ] Extract adults GWAS results

### Day 2
- [ ] Adults visualizations
- [ ] Run embryos analysis
- [ ] Run larvae analysis

### Day 3
- [ ] Meta-analysis
- [ ] Cross-stage summaries
- [ ] Final figures

### Day 4
- [ ] Write methods section
- [ ] Prepare reviewer response
- [ ] Documentation cleanup

---

## Key Commands

### Running Scripts in Container
```bash
# For R scripts
singularity exec albopictus-diapause-rnaseq.sif \
  Rscript scripts/06_differential_expression/03_adults_deseq2.R

# For Python scripts
singularity exec albopictus-diapause-rnaseq.sif \
  python scripts/06_differential_expression/02_adults_qc.py
```

### Submitting to SLURM
```bash
# Create submission script
sbatch --wrap="singularity exec albopictus-diapause-rnaseq.sif \
  Rscript scripts/06_differential_expression/03_adults_deseq2.R" \
  --job-name=adults_de \
  --output=logs/adults_de_%j.log \
  --cpus-per-task=8 \
  --mem=32G \
  --time=2:00:00
```

---

## Critical Reminders

1. **NO BATCH CORRECTION** - Each stage is analyzed independently
2. **Use Container** - All analyses must use the .sif file
3. **Adults First** - Most important for GWAS validation
4. **Simple Design** - Just `~ Condition`, no covariates
5. **Document Everything** - This unusual situation needs clear explanation

---

## Success Criteria

✅ Analysis is successful if:
1. Adults show clear differential expression for some GWAS candidates
2. No batch correction was attempted
3. Each stage analyzed independently
4. Meta-analysis combines results appropriately
5. Methods clearly explain the confounding issue

❌ Analysis fails if:
1. Attempting any batch correction
2. Combining stages before analysis
3. Using platform as covariate
4. Mixing samples across stages

---

**This roadmap bypasses the impossible batch correction and provides a statistically valid approach to the confounded data.**