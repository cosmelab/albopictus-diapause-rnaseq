# QC Analysis Workflow

## Step 1: Extract QC Metrics (Per Sample)
Create a script that:
1. Reads MultiQC general_stats.txt for each sample
2. Extracts key metrics from individual QC files (qualimap, rseqc)
3. Creates a per-sample JSON/CSV with all metrics

## Step 2: Aggregate by Project
1. Combine all samples within each project (PRJNA*)
2. Calculate mean, SD, min, max for each metric
3. Identify potential outliers

## Step 3: Create Master Summary
1. Combine all projects into one master table
2. Add metadata (photoperiod, time point, etc.)
3. Format for manuscript supplementary table

## Step 4: Generate Figures
1. Basic QC plots (mapping, read distribution)
2. Advanced analysis (PCA, correlations)
3. Project-specific comparisons

## Key Files to Extract From:
- `multiqc_general_stats.txt` - Main metrics
- `rnaseq_qc_results.txt` - Qualimap RNA-seq metrics
- `*.infer_experiment.txt` - Strand specificity
- `*.read_distribution.txt` - Genomic feature distribution
- `salmon.merged.gene_counts.tsv` - For gene detection stats

## Output Structure:
```
qc_results/
├── per_sample/
│   ├── PRJNA268379_SRR1663685_qc.json
│   ├── PRJNA268379_SRR1663687_qc.json
│   └── ...
├── per_project/
│   ├── PRJNA268379_summary.csv
│   ├── PRJNA158021_summary.csv
│   └── PRJNA187045_summary.csv
├── combined/
│   ├── all_samples_qc_metrics.csv
│   └── manuscript_table_s1_qc.csv
└── figures/
    ├── figure_1_qc_overview.pdf
    └── figure_s1_detailed_qc.pdf
```