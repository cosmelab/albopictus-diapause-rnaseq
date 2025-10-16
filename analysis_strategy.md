# Comprehensive Analysis Strategy for Albopictus Diapause RNA-seq

## Overview
Build a fully reproducible analysis addressing molecular ecology reviewers' concerns while validating GWAS candidate genes.

## Key Goals
1. **Transparency**: Every step documented and reproducible
2. **Validation**: Compare with collaborators' published results
3. **Robustness**: Address reviewer concerns about batch effects and statistical models
4. **Accessibility**: Container includes everything needed to reproduce

## Phase 1: Setup and Reference Collection

### 1.1 Clone Collaborator Repositories
```bash
mkdir -p references/collaborator_repos
cd references/collaborator_repos

# Clone relevant repos (need to identify URLs):
# - Huang et al. 2015 (PRJNA268379 - Adults)
# - Poelchau et al. 2013a (PRJNA158021 - Embryos)
# - Poelchau et al. 2013b (PRJNA187045 - Pharate larvae)
```

### 1.2 Extract Key Information
- Original DE analysis scripts
- Gene lists and annotations used
- Statistical models applied
- QC thresholds and filtering criteria
- Any custom code or parameters

### 1.3 Document Differences
Create comparison table:
- Our approach (nf-core) vs their approach
- Modern best practices vs 2013-2015 methods
- How we address reviewer concerns

## Phase 2: Data Organization and QC

### 2.1 Directory Restructuring
```
albopictus-diapause-rnaseq/
├── analysis/                    # All analysis code
│   ├── 01_qc_extraction/       # Extract from nf-core outputs
│   ├── 02_count_matrices/      # Prepare for DE
│   ├── 03_de_analysis/         # Per-project DE
│   ├── 04_meta_analysis/       # Cross-project
│   └── 05_publication/         # Final outputs
├── data/                       # Raw and processed data
│   ├── nf_core_outputs/        # Symlinks to output/
│   ├── count_matrices/         # Combined counts
│   └── metadata/               # Sample sheets, gene lists
├── references/                 # External resources
│   ├── collaborator_repos/     # Their code
│   ├── manuscripts/            # PDFs of papers
│   └── reviewer_comments/      # What to address
├── results/                    # Our analysis outputs
│   ├── 01_qc_reports/
│   ├── 02_de_results/
│   └── 03_figures/
└── container/                  # Docker/Singularity files
    ├── Dockerfile
    ├── requirements.txt
    └── data_snapshot/          # Processed data for container
```

### 2.2 QC Strategy
1. Extract all QC metrics from nf-core outputs
2. Create summary table comparing to published QC standards
3. Flag any samples below thresholds
4. Document decisions for reviewer transparency

## Phase 3: Differential Expression Analysis

### 3.1 Statistical Models (Addressing Reviewer Concerns)
```r
# Document each model choice with justification

# Model 1: Batch effect correction
# "Reviewers concerned about platform differences"
# Include platform as covariate or use ComBat-seq

# Model 2: Appropriate design matrices
# Show we're using correct factorial/time-series models

# Model 3: Multiple testing correction
# Use FDR, document threshold choices
```

### 3.2 Validation Strategy
1. Run our pipeline on same samples
2. Compare DEG lists with published results
3. Calculate overlap statistics
4. If differences, explain why (better methods, QC, etc.)

## Phase 4: Container Development

### 4.1 Enhanced Dockerfile
```dockerfile
# Base image with all tools
FROM mambaorg/micromamba:1.5.1

# Copy analysis code
COPY analysis/ /opt/analysis/

# Copy processed data (optional - for full reproducibility)
COPY data/count_matrices/ /opt/data/count_matrices/
COPY results/01_qc_reports/ /opt/data/qc_reports/

# Include manuscript generation
COPY manuscript/ /opt/manuscript/
```

### 4.2 Container Versions
1. **Light version**: Tools only, user provides data
2. **Full version**: Includes processed data, can reproduce from counts
3. **Complete version**: Everything including nf-core outputs (large!)

### 4.3 Full Reproducibility Package
Include everything needed to replicate from scratch:
```dockerfile
# Include SRA download scripts
COPY scripts/01_download_data/ /opt/scripts/data_acquisition/

# Include nf-core pipeline scripts  
COPY scripts/02_run_rnaseq/ /opt/scripts/nf_core_pipeline/

# Include sample metadata
COPY data/metadata/samples_ncbi.txt /opt/data/metadata/
```

This allows users to:
1. Download raw data from SRA using our exact scripts
2. Run nf-core pipeline with our exact parameters
3. Reproduce all downstream analysis
4. Verify every step of our process

## Phase 5: Documentation and Reviewer Response

### 5.1 Methods Documentation
- Detailed methods.md with every parameter
- Comparison table: our methods vs original papers
- Justification for each analytical choice

### 5.2 Reviewer Response Document
```markdown
# Response to Molecular Ecology Reviewers

## Concern 1: Batch Effects
- We implemented [specific approach]
- Results show [evidence of correction]

## Concern 2: Statistical Models
- Models chosen based on [experimental design]
- Validated against published results
```

## Implementation Timeline

### Week 1: Foundation
- [ ] Clone collaborator repos
- [ ] Run QC extraction
- [ ] Reorganize directory structure

### Week 2: Analysis
- [ ] Generate count matrices
- [ ] Run DE analysis per project
- [ ] Compare with published results

### Week 3: Integration
- [ ] Meta-analysis across projects
- [ ] Focus on candidate genes
- [ ] Create publication figures

### Week 4: Package
- [ ] Update Dockerfile
- [ ] Build new container
- [ ] Write methods and response

## Success Criteria
1. Can reproduce published DEG lists (or explain differences)
2. Container allows full reproduction from counts
3. All reviewer concerns addressed with evidence
4. Clear documentation for next researcher

## Notes
- Keep molecular ecology standards in mind
- Every decision should be justified
- Transparency > perfection
- Make it easy for reviewers to verify our work