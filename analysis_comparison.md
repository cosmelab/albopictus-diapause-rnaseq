# Analysis Comparison: Our Approach vs Collaborators

## Overview
The collaborators remapped all three diapause studies to a newer reference genome. We're doing the same but with modern nf-core pipeline.

## Key Differences

### 1. Pipeline Framework
- **Collaborators**: Manual scripts (Trimmomatic → STAR → HTSeq)
- **Our approach**: nf-core/rnaseq v3.19.0 (automated, containerized)

### 2. Quantification Method
- **Collaborators**: STAR alignment + HTSeq counting (gene-level)
- **Our approach**: Salmon pseudo-alignment (faster, includes TPM)

### 3. Quality Control
- **Collaborators**: FastQC + manual inspection
- **Our approach**: MultiQC + Qualimap + RSeQC (comprehensive)

### 4. Trimming Strategy
- **Collaborators**: Trimmomatic with HEADCROP:15 (remove first 15 bases)
- **Our approach**: nf-core default trimming (adaptive)

### 5. Reference Genome
- **Collaborators**: New PI group assembly (not public?)
- **Our approach**: Public reference (data/references/albo.fasta.gz)

## Project Assignments
Based on their GitHub structure:
- **Angela** → PRJNA268379 (Adult females)
- **Mackenzie** → PRJNA158021 (Embryos)  
- **Sarah** → PRJNA187045 (Pharate larvae)

## Statistical Models
They used project-specific DESeq2 models:
- Adults: Factorial design (photoperiod × blood meal)
- Embryos: Time course with condition
- Pharate larvae: Time course comparison

## Key Scripts to Review

### QC and Counting
- `scripts/Sarah/CountReadsAssignedToGenes.R` - How they summarized counts
- `misc/TestingStrandedness.md` - Important for count accuracy

### Differential Expression
- `scripts/Angela/Angela_DESeq.R` - Adult analysis
- `scripts/Mackenzie/Mackenzie_DESeq.R` - Embryo analysis
- `scripts/Sarah/DESeq.R` - Pharate larvae analysis

### Downstream Analysis
- GO enrichment scripts
- KEGG pathway analysis

## Advantages of Our Approach

1. **Reproducibility**: Containerized pipeline with locked versions
2. **Speed**: Salmon is much faster than STAR+HTSeq
3. **QC**: More comprehensive metrics from nf-core
4. **Automation**: Less manual intervention, fewer errors
5. **Standards**: Following current best practices

## How to Address Reviewers

When reviewers ask about differences from original studies:
1. We used the same samples and experimental designs
2. Modern pipeline improves accuracy and reproducibility
3. Salmon has been validated to produce comparable results to STAR+HTSeq
4. We can show correlation between our counts and theirs
5. Any differences in DEGs likely due to:
   - Better QC filtering
   - Improved quantification algorithms
   - Different reference annotation version

## Next Steps

1. Extract key parameters from their scripts (e.g., HTSeq strandedness)
2. Run our QC analysis
3. Compare gene counts between methods
4. Document all differences for transparency