# QC Analysis Scripts

This directory contains scripts for extracting and visualizing quality control metrics from nf-core RNA-seq outputs.

## Scripts Overview

### 1. `extract_qc_metrics.py`
Extracts QC metrics from nf-core/rnaseq pipeline outputs.

**Input Requirements:**
- Expects nf-core output structure in `output/` directory:
  ```
  output/
  ├── PRJNA268379/
  │   ├── SRR1663685/
  │   │   ├── multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt
  │   │   ├── star_salmon/qualimap/*/rnaseq_qc_results.txt
  │   │   ├── star_salmon/salmon.merged.gene_counts.tsv
  │   │   └── star_salmon/rseqc/infer_experiment/*.infer_experiment.txt
  │   └── ...
  └── ...
  ```

**Metrics Extracted:**
- **Sequencing metrics**: Total reads, read length, GC%, duplication rate
- **Mapping metrics**: % mapped, % uniquely mapped, % multi-mapped
- **RNA-seq metrics**: % exonic/intronic/intergenic reads, rRNA%, 5'/3' bias
- **Quality metrics**: Error rate, % properly paired, insert size
- **Gene detection**: Number of genes with >0, >10, >100, >1000 counts
- **Strand specificity**: Forward/reverse strand percentages

**Outputs:**
```
qc_results/
├── all_samples_qc_metrics.csv          # Master table with all samples
└── per_sample/
    ├── PRJNA268379_SRR1663685_qc.json  # Individual sample JSON
    ├── PRJNA268379_SRR1663687_qc.json
    └── ...
```

### 2. `generate_qc_figures.py`
Creates publication-ready QC figures from extracted metrics.

**Input Requirements:**
- `qc_results/all_samples_qc_metrics.csv` (created by extract_qc_metrics.py)

**Figures Generated:**

1. **Figure 1: QC Overview** (4 panels)
   - Panel A: Sequencing depth bar chart
   - Panel B: Mapping rates stacked bar chart
   - Panel C: Read distribution (exonic/intronic/intergenic)
   - Panel D: Gene detection boxplots by project

2. **Figure 2: Detailed QC** (4 panels)
   - Panel A: 5'/3' bias scatter plot
   - Panel B: rRNA contamination boxplots
   - Panel C: Insert size violin plots
   - Panel D: Duplication rate vs library size

3. **Supplementary Figure S1**: Correlation heatmap of all QC metrics

**Tables Generated:**
- `table_1_qc_summary.csv`: Summary statistics by life stage (mean ± SD)
- `table_1_qc_summary.tex`: LaTeX formatted table for manuscripts

**Outputs:**
```
qc_results/
├── figures/
│   ├── figure_1_qc_overview.pdf
│   ├── figure_1_qc_overview.png
│   ├── figure_2_detailed_qc.pdf
│   ├── figure_2_detailed_qc.png
│   ├── figure_s1_correlation_heatmap.pdf
│   └── figure_s1_correlation_heatmap.png
└── tables/
    ├── table_1_qc_summary.csv
    └── table_1_qc_summary.tex
```

### 3. `run_qc_analysis.sh`
Automated script to run both Python scripts using Singularity container.

**Usage:**
```bash
./scripts/03_qc_analysis/run_qc_analysis.sh
```

**What it does:**
1. Loads singularity module
2. Checks if container exists
3. Runs extract_qc_metrics.py
4. Checks if extraction succeeded
5. Runs generate_qc_figures.py
6. Reports results location

## Manual Usage

### Running in Singularity Container

```bash
# Load module
module load singularity/3.9.3

# Interactive mode
singularity shell --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif
cd /proj
python scripts/03_qc_analysis/extract_qc_metrics.py
python scripts/03_qc_analysis/generate_qc_figures.py
exit

# Or execute directly
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
    bash -c "cd /proj && python scripts/03_qc_analysis/extract_qc_metrics.py"
```

## Output Interpretation

### all_samples_qc_metrics.csv
CSV file with one row per sample containing ~30 QC metrics:
- `project_id`: PRJNA number
- `sample_id`: SRR accession
- `total_reads_millions`: Library size
- `uniquely_mapped_percent`: Should be >70%
- `exonic_percent`: Should be >50% for good RNA-seq
- `rrna_percent`: Should be <5% ideally
- `genes_detected_gt10`: Typical range 10,000-15,000

### Quality Thresholds
Good quality RNA-seq should have:
- Total reads: >20 million
- Uniquely mapped: >70%
- Exonic reads: >50%
- rRNA contamination: <5%
- 5'/3' bias: Close to 1.0
- Genes detected (>10 counts): >10,000

### Troubleshooting

**Error: "QC metrics file not found"**
- Run extract_qc_metrics.py first

**Error: "No such file or directory"**
- Check that nf-core pipeline completed successfully
- Verify output directory structure matches expected format

**Missing metrics**
- Some metrics may be NA if specific QC tools were skipped in nf-core
- Check multiqc_general_stats.txt for available columns

## Customization

To add new metrics:
1. Edit `extract_qc_metrics.py` and add extraction logic
2. Update `generate_qc_figures.py` if you want to plot new metrics
3. Add metric description to this README

## Dependencies
- Python 3.10+
- pandas
- matplotlib
- seaborn
- numpy
- scikit-learn (for PCA in correlation analysis)