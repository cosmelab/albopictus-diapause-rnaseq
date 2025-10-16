# Scripts Reorganization Plan

## Current Structure Issues:
1. Numbered directories (01_, 02_) are rigid - hard to insert new steps
2. "03_qc_analysis" has too many old scripts mixed with new ones
3. "05_visualization" overlaps with QC and DE analysis
4. Directory names don't match our analysis phases

## Proposed New Structure:
```
scripts/
├── data_acquisition/          # Previously 01_download_data
│   └── sra_download/
├── nf_core_pipeline/         # Previously 02_run_rnaseq
│   └── array_job_scripts/
├── qc_extraction/            # Extract QC from nf-core outputs
│   ├── extract_metrics.py
│   ├── generate_figures.py
│   └── run_qc_analysis.sh
├── count_matrices/           # Prepare count data for DE
│   ├── combine_salmon_counts.py
│   └── prepare_tpm_matrices.py
├── differential_expression/  # DE analysis by project
│   ├── PRJNA268379_adults/
│   ├── PRJNA158021_embryos/
│   └── PRJNA187045_pharate/
├── meta_analysis/           # Cross-project integration
│   ├── candidate_gene_analysis.py
│   └── batch_correction.py
├── publication/             # Final figures and tables
│   ├── figures/
│   └── tables/
└── utils/                   # Helper scripts
    └── container_tools/
```

## Benefits:
1. Follows our analysis workflow phases
2. Clear separation of concerns
3. Easy to find relevant scripts
4. No rigid numbering system
5. Matches the results/ directory structure

## Migration Plan:
1. Create new directory structure
2. Move only the scripts we need (ignore old/outdated ones)
3. Update paths in scripts as needed
4. Archive old structure for reference