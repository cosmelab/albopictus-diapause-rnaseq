# QC Analysis Plan and Directory Structure

## Agreed Directory Structure (from project_rules.md)

```
results/
├── organized/          # All organized outputs go here
│   ├── qc_reports/     # QC analysis outputs
│   │   ├── metrics/    # Extracted QC metrics
│   │   ├── figures/    # QC plots
│   │   └── tables/     # Summary tables
│   ├── gene_counts/    # For future: combined count matrices
│   └── gene_tpm/       # For future: combined TPM matrices
└── analysis/           # Downstream analysis (DE, etc.)
```

## Scripts Organization

All scripts in `scripts/03_qc_analysis/` will be numbered:

1. **01_extract_qc_metrics.py**
   - Input: `output/PRJNA*/SRR*/` (nf-core outputs)
   - Output: `results/organized/qc_reports/metrics/`
   - Creates: `all_samples_qc_metrics.csv` and per-sample JSONs

2. **02_generate_qc_figures.py**
   - Input: `results/organized/qc_reports/metrics/all_samples_qc_metrics.csv`
   - Output: `results/organized/qc_reports/figures/`
   - Creates: figure_1_qc_overview.pdf, figure_2_detailed_qc.pdf, etc.

3. **03_create_summary_tables.py** (if needed)
   - Input: `results/organized/qc_reports/metrics/all_samples_qc_metrics.csv`
   - Output: `results/organized/qc_reports/tables/`
   - Creates: Publication-ready tables

## Current Status
- Created extraction and figure scripts but used wrong output directory
- Need to update scripts to use `results/organized/qc_reports/`
- Need to rename scripts to follow numbering convention

## Next Steps
1. Update scripts to use correct directory structure
2. Rename scripts with 01_, 02_ prefix
3. Update run script to reflect correct paths
4. Test the pipeline