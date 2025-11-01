# Differential Expression Analysis Scripts

**Purpose:** Stage-specific differential expression analysis using DESeq2

**Author:** L. Cosme
**Date:** October 22, 2025

---

## Overview

This directory contains all scripts for differential expression (DE) analysis of the albopictus diapause RNA-seq data. Following the discovery that Platform = Stage = Study (complete confounding), we analyze each developmental stage independently.

---

## Script Organization

### Adults (PRJNA268379) - 16 samples

Following Huang et al. 2015 methodology, adult females are analyzed separately by blood meal status to isolate photoperiod effects:

**01. `01_adults_deseq2_unfed.R`**
- Analysis: Unfed females only (8 samples)
- Comparison: SD_Unfed (diapause-inducing) vs LD_Unfed (non-diapause)
- Design: `~ Photoperiod`
- Batch script: `run_adults_unfed_deseq2.sh`
- Job: 20616099 (submitted Oct 22, 2025)

**02. `02_adults_deseq2_fed.R`**
- Analysis: Blood-fed females only (8 samples)
- Comparison: SD_Fed (diapause) vs LD_Fed (non-diapause)
- Design: `~ Photoperiod`
- Batch script: `run_adults_fed_deseq2.sh` (to be created)

### Embryos (PRJNA158021) - 12 samples

**03. `03_embryos_deseq2.R`** (to be created)
- Analysis: All embryo samples
- Comparison: Diapause-destined vs Non-diapause
- Design: `~ Condition`
- No blood meal confounding

### Larvae (PRJNA187045) - 17 samples

**04. `04_larvae_deseq2.R`** (to be created)
- Analysis: All pharate larval samples
- Comparison: Diapause vs Non-diapause
- Design: `~ Condition`
- No blood meal confounding

### GWAS Validation

**05. `05_extract_gwas_candidates.R`** (to be created)
- Extract results for 34 GWAS candidate genes
- Combine across all stages
- Create summary table

**06. `06_fishers_exact_enrichment.R`** (to be created)
- Test if GWAS candidates are enriched among DEGs
- Fisher's exact test
- Calculate odds ratios and confidence intervals

### Visualization

**07. `07_create_volcano_plots.R`** (to be created)
- EnhancedVolcano plots for each analysis
- Highlight GWAS candidates
- Publication-quality figures

**08. `08_create_heatmaps.R`** (to be created)
- Cross-stage expression heatmaps
- GWAS candidate heatmaps
- Clustering analyses

### Tables

**09. `09_generate_supplementary_tables.R`** (to be created)
- Table S1: Sample information
- Table S2: GWAS candidate expression
- Table S3: Genome-wide DE statistics
- Table S4: Enrichment test results

---

## Analysis Parameters

### Pre-filtering
- Minimum counts: ≥10 in ≥25% of samples
- Applied per-stage (not across all samples)

### DESeq2 Settings
- Significance threshold: FDR < 0.05
- Test: Wald test (default)
- Multiple testing correction: Benjamini-Hochberg
- Shrinkage: None yet (TODO: add lfcShrink with apeglm)

### Design Formulas

| Stage | Samples | Design | Notes |
|-------|---------|--------|-------|
| Adults (unfed) | 8 | `~ Photoperiod` | LD vs SD, no blood meal |
| Adults (fed) | 8 | `~ Photoperiod` | LD vs SD, with blood meal |
| Embryos | 12 | `~ Condition` | Diapause vs non-diapause |
| Larvae | 17 | `~ Condition` | Diapause vs non-diapause |

---

## Running the Scripts

### Individual Script Execution

Using Singularity container (required):
```bash
# Load module
module load singularity-ce/3.9.3

# Run R script
singularity exec --cleanenv --no-home \
  albopictus-diapause-rnaseq.sif \
  Rscript scripts/06_differential_expression/01_adults_deseq2_unfed.R
```

### Batch Submission (Recommended)

```bash
# Adults unfed
sbatch scripts/06_differential_expression/run_adults_unfed_deseq2.sh

# Adults fed
sbatch scripts/06_differential_expression/run_adults_fed_deseq2.sh

# Check job status
squeue -u $USER
```

---

## Output Organization

Results are saved to `results/03_differential_expression/`:

```
results/03_differential_expression/
├── adults/
│   ├── adults_unfed_deseq2_all.tsv           # All genes
│   ├── adults_unfed_deseq2_significant.tsv   # FDR < 0.05 only
│   ├── adults_unfed_summary_stats.tsv        # Summary statistics
│   ├── adults_unfed_pca.pdf                  # QC plot
│   ├── adults_unfed_correlation_heatmap.pdf  # QC plot
│   ├── adults_unfed_ma_plot.pdf              # MA plot
│   ├── adults_unfed_dispersion.pdf           # Dispersion estimates
│   ├── adults_unfed_session_info.txt         # R session info
│   ├── adults_fed_deseq2_all.tsv
│   ├── adults_fed_deseq2_significant.tsv
│   └── ... (similar files for fed)
├── embryos/
│   └── ... (similar structure)
├── larvae/
│   └── ... (similar structure)
└── combined/
    ├── gwas_candidates_all_stages.tsv
    ├── enrichment_test_results.tsv
    └── summary_statistics.tsv
```

---

## Job Tracking

| Job ID | Script | Status | Submitted | Completed | Notes |
|--------|--------|--------|-----------|-----------|-------|
| 20616099 | 01_adults_deseq2_unfed.R | Running | Oct 22, 2025 | - | Unfed females |
| - | 02_adults_deseq2_fed.R | Pending | - | - | To be submitted |
| - | 03_embryos_deseq2.R | Not created | - | - | - |
| - | 04_larvae_deseq2.R | Not created | - | - | - |

---

## Expected Results

Based on preliminary analysis (archived docs):

| Stage | Total Genes | DEGs (FDR<0.05) | % DE | Notes |
|-------|-------------|-----------------|------|-------|
| Adults (unfed) | ~28,000 | TBD | TBD | Most relevant for GWAS |
| Adults (fed) | ~28,000 | TBD | TBD | Blood meal effect |
| Embryos | ~25,000 | ~35 | 0.14% | Very few DEGs expected |
| Larvae | ~25,000 | ~400 | 1.6% | Moderate DEGs |

---

## Next Steps

1. ✅ Adults unfed analysis (Job 20616099 running)
2. ⏳ Adults fed analysis (script created, needs submission)
3. ⏳ Embryos analysis (script to be created)
4. ⏳ Larvae analysis (script to be created)
5. ⏳ GWAS candidate extraction
6. ⏳ Enrichment testing
7. ⏳ Visualization and tables

---

## Notes

- **Container-only policy:** All analyses MUST run in Singularity container
- **No batch correction:** Each stage analyzed independently (platform confounding)
- **Blood meal splitting:** Adults require special handling due to experimental design
- **GWAS priority:** Adults are most important for validating GWAS candidates

---

## Troubleshooting

### Common Issues

**Issue:** "cannot find function 'DESeq'"
**Solution:** Ensure running inside container with DESeq2 installed

**Issue:** "No such file or directory"
**Solution:** Check paths are absolute, not relative

**Issue:** matplotlib/fontconfig warnings
**Solution:** These are harmless; plots are still generated

### Getting Help

- Check session logs in `logs/`
- Review session info files in output directories
- Consult `archived_docs/` for historical context

---

**Last updated:** October 22, 2025 12:30
