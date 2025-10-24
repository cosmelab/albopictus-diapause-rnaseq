#!/bin/bash
# Create comprehensive directory structure for RNA-seq analysis results
# Author: L. Cosme
# Date: October 2025

# Set the base directory
BASE_DIR="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
RESULTS_DIR="${BASE_DIR}/results"

echo "Creating comprehensive results directory structure..."
echo "Base directory: ${RESULTS_DIR}"

# Create main results directories
mkdir -p ${RESULTS_DIR}/{01_quality_control,02_count_matrices,03_differential_expression,04_gwas_validation,05_meta_analysis,06_functional_enrichment,07_method_comparison,08_manuscript_figures,09_supplementary_tables,10_session_info}

# 01_quality_control subdirectories
mkdir -p ${RESULTS_DIR}/01_quality_control/{stage_specific,combined}
mkdir -p ${RESULTS_DIR}/01_quality_control/stage_specific/{adults,embryos,larvae}

# 02_count_matrices subdirectories
mkdir -p ${RESULTS_DIR}/02_count_matrices/{raw,filtered,normalized}

# 03_differential_expression subdirectories
mkdir -p ${RESULTS_DIR}/03_differential_expression/{adults,embryos,larvae}

# 04_gwas_validation directory already created

# 05_meta_analysis directory already created

# 06_functional_enrichment subdirectories
mkdir -p ${RESULTS_DIR}/06_functional_enrichment/{adults,embryos,larvae,cross_stage}

# 07_method_comparison subdirectories
mkdir -p ${RESULTS_DIR}/07_method_comparison/{salmon_vs_featurecounts,salmon_vs_htseq}

# 08_manuscript_figures subdirectories
mkdir -p ${RESULTS_DIR}/08_manuscript_figures/{main_figures,supplementary_figures}

# 09_supplementary_tables directory already created

# 10_session_info directory already created

# Create logs and checkpoint directories
mkdir -p ${BASE_DIR}/logs/sessions
mkdir -p ${BASE_DIR}/checkpoints

# Create scripts directories
mkdir -p ${BASE_DIR}/scripts/{06_differential_expression,07_visualization,08_functional_analysis,09_method_validation,10_manuscript_figures,utils}

echo "Directory structure created successfully!"

# Create directory tree visualization
echo ""
echo "Directory structure:"
tree ${RESULTS_DIR} -d -L 3

# Create initial README files
echo "Creating README templates..."

# Main results README
cat > ${RESULTS_DIR}/README.md << 'EOF'
# Results Directory Structure

## Overview
This directory contains all analysis results for the Aedes albopictus diapause RNA-seq project.

## Directory Organization

### 01_quality_control/
Quality control metrics and visualizations for all samples

### 02_count_matrices/
Raw, filtered, and normalized expression matrices

### 03_differential_expression/
Stage-specific differential expression results

### 04_gwas_validation/
GWAS candidate gene validation results

### 05_meta_analysis/
Cross-stage meta-analysis results

### 06_functional_enrichment/
GO and KEGG enrichment analyses

### 07_method_comparison/
Validation of quantification methods

### 08_manuscript_figures/
Publication-ready figures

### 09_supplementary_tables/
Supplementary data tables

### 10_session_info/
Software versions and session information

## Generation
Created: $(date)
Script: scripts/utils/create_directory_structure.sh
EOF

# QC README
cat > ${RESULTS_DIR}/01_quality_control/README.md << 'EOF'
# Quality Control Results

## Contents
- `stage_specific/`: QC metrics separated by developmental stage
- `combined/`: QC metrics for all samples together

## Key Metrics
- Mapping rates
- Library sizes
- rRNA contamination
- Duplication rates
- Sample correlations

## Generation
Scripts: `scripts/04_qc_analysis/`
EOF

# Count matrices README
cat > ${RESULTS_DIR}/02_count_matrices/README.md << 'EOF'
# Expression Count Matrices

## Directories

### raw/
Unfiltered count matrices directly from quantification tools

### filtered/
Count matrices after low-count gene filtering

### normalized/
Normalized expression values (rlog, VST)

## File Naming Convention
`[stage]_[tool]_[type].tsv`
- stage: adults, embryos, larvae, or all
- tool: salmon, featurecounts, htseq
- type: counts, tpm, fpkm

## Generation
Scripts: `scripts/04_qc_analysis/01_gather_count_matrices.py`
EOF

# DE README
cat > ${RESULTS_DIR}/03_differential_expression/README.md << 'EOF'
# Differential Expression Results

## Organization
Separate directories for each developmental stage:
- `adults/`: Adult female analysis (priority for GWAS)
- `embryos/`: Embryonic stage analysis
- `larvae/`: Pharate larvae analysis

## Files per Stage
- `deseq2_results_all.tsv`: Complete DE results
- `significant_genes_fdr05.tsv`: Filtered for FDR < 0.05
- `gwas_candidates_results.tsv`: 34 GWAS genes only
- Visualization plots (volcano, MA, heatmaps)

## Statistical Approach
- Method: DESeq2
- Design: ~ Condition (no batch correction)
- Significance: FDR < 0.05

## Generation
Scripts: `scripts/06_differential_expression/`
EOF

echo "README templates created!"

# Create progress tracker
cat > ${BASE_DIR}/progress_tracker.md << 'EOF'
# Analysis Progress Tracker

## Data Processing
- [x] Raw data download
- [x] Reference genome preparation
- [x] GTF annotation preprocessing
- [x] RNA-seq pipeline (nf-core)
- [x] Count matrix generation

## Quality Control
- [x] Initial QC metrics extraction
- [ ] Stage-specific QC analysis
- [ ] QC visualization generation

## Differential Expression
- [ ] Adults analysis
- [ ] Embryos analysis
- [ ] Larvae analysis

## GWAS Validation
- [ ] Extract candidate results
- [ ] Cross-stage comparison
- [ ] Validation summary

## Meta-Analysis
- [ ] Fisher's method
- [ ] Stouffer's method
- [ ] Consistency assessment

## Functional Analysis
- [ ] GO enrichment
- [ ] KEGG pathways
- [ ] Network analysis

## Manuscript
- [ ] Main figures
- [ ] Supplementary figures
- [ ] Tables
- [ ] Methods section
- [ ] Results section

## Last Updated
$(date)
EOF

echo ""
echo "✓ Directory structure created"
echo "✓ README templates generated"
echo "✓ Progress tracker initialized"
echo ""
echo "Next steps:"
echo "1. Run metadata splitting: python scripts/06_differential_expression/01_split_metadata.py"
echo "2. Run adult DE analysis: Rscript scripts/06_differential_expression/03_adults_deseq2.R"