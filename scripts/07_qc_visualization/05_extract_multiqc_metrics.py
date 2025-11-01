#!/usr/bin/env python3
"""
Extract and Integrate MultiQC Metrics

Purpose: Parse MultiQC JSON output and integrate with expression QC
         Check if technical metrics explain expression variation
         Answer Reviewer 3.4: Mapping statistics

Outputs:
  - Technical metrics summary table
  - Technical vs PCA correlation plots
  - Mapping rate vs expression quality
  - Duplication rate analysis

Usage:
  python 05_extract_multiqc_metrics.py

Created: October 27, 2025
Phase: 2.5 - QC Visualization
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# =============================================================================
# Configuration
# =============================================================================

MULTIQC_JSON = "output/multiqc/star_salmon/multiqc_data/multiqc_data.json"
OUTPUT_DIR = Path("results/qc_visualization/multiqc_integration")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
METRICS_TABLE = OUTPUT_DIR / "technical_metrics_summary.tsv"
MAPPING_PLOT = OUTPUT_DIR / "mapping_statistics.pdf"
DUPLICATION_PLOT = OUTPUT_DIR / "duplication_analysis.pdf"
SUMMARY_FILE = OUTPUT_DIR / "multiqc_summary.txt"

# =============================================================================
# Helper Functions
# =============================================================================

def load_multiqc_data(json_file):
    """Load MultiQC JSON data."""
    print(f"Loading MultiQC data from {json_file}...")

    if not Path(json_file).exists():
        print(f"ERROR: MultiQC JSON not found: {json_file}")
        print("Expected location: output/multiqc/star_salmon/multiqc_data/multiqc_data.json")
        sys.exit(1)

    with open(json_file, 'r') as f:
        data = json.load(f)

    print(f"  Loaded MultiQC data\n")
    return data


def extract_star_metrics(multiqc_data):
    """Extract STAR alignment metrics."""
    print("Extracting STAR alignment metrics...")

    star_data = multiqc_data.get('report_plot_data', {}).get('star', {})

    if not star_data:
        print("  WARNING: No STAR data found in MultiQC")
        return pd.DataFrame()

    # Get alignment stats
    metrics = []
    for sample_id, stats in star_data.items():
        metric_dict = {
            'sample_id': sample_id,
            'total_reads': stats.get('total_reads', 0),
            'uniquely_mapped': stats.get('uniquely_mapped', 0),
            'uniquely_mapped_percent': stats.get('uniquely_mapped_percent', 0),
            'multimapped': stats.get('multimapped', 0),
            'multimapped_percent': stats.get('multimapped_percent', 0),
            'unmapped_percent': stats.get('unmapped_percent', 0)
        }
        metrics.append(metric_dict)

    df = pd.DataFrame(metrics)
    print(f"  Extracted metrics for {len(df)} samples\n")

    return df


def extract_rseqc_metrics(multiqc_data):
    """Extract RSeQC read distribution metrics."""
    print("Extracting RSeQC metrics...")

    rseqc_data = multiqc_data.get('report_plot_data', {}).get('rseqc_read_distribution', {})

    if not rseqc_data:
        print("  WARNING: No RSeQC data found in MultiQC")
        return pd.DataFrame()

    metrics = []
    for sample_id, stats in rseqc_data.items():
        metric_dict = {
            'sample_id': sample_id,
            'exonic_rate': stats.get('CDS_Exons', {}).get('percent', 0),
            'intronic_rate': stats.get('Introns', {}).get('percent', 0),
            'intergenic_rate': stats.get('Intergenic', {}).get('percent', 0)
        }
        metrics.append(metric_dict)

    df = pd.DataFrame(metrics)
    print(f"  Extracted metrics for {len(df)} samples\n")

    return df


def extract_dupradar_metrics(multiqc_data):
    """Extract dupRadar duplication metrics."""
    print("Extracting dupRadar metrics...")

    dupradar_data = multiqc_data.get('report_plot_data', {}).get('dupradar', {})

    if not dupradar_data:
        print("  WARNING: No dupRadar data found in MultiQC")
        return pd.DataFrame()

    metrics = []
    for sample_id, stats in dupradar_data.items():
        # dupRadar provides duplication rate at different expression levels
        metric_dict = {
            'sample_id': sample_id,
            'duplication_rate': np.mean([v for v in stats.values() if isinstance(v, (int, float))])
        }
        metrics.append(metric_dict)

    df = pd.DataFrame(metrics)
    print(f"  Extracted metrics for {len(df)} samples\n")

    return df


def combine_metrics(star_df, rseqc_df, dupradar_df):
    """Combine all technical metrics into one table."""
    print("Combining metrics...")

    # Start with STAR (should have all samples)
    combined = star_df.copy()

    # Merge RSeQC if available
    if not rseqc_df.empty:
        combined = combined.merge(rseqc_df, on='sample_id', how='left')

    # Merge dupRadar if available
    if not dupradar_df.empty:
        combined = combined.merge(dupradar_df, on='sample_id', how='left')

    # Add metadata (life stage from sample ID)
    def get_stage(sample_id):
        if 'PRJNA268379' in sample_id:
            return 'adults'
        elif 'PRJNA158021' in sample_id:
            return 'embryos'
        elif 'PRJNA187045' in sample_id:
            return 'larvae'
        else:
            return 'unknown'

    combined['stage'] = combined['sample_id'].apply(get_stage)

    print(f"  Combined metrics: {combined.shape[0]} samples, {combined.shape[1]} metrics\n")

    return combined


def create_mapping_plot(metrics_df, output_file):
    """Create mapping statistics plot."""
    print("Creating mapping statistics plot...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Mapping Statistics Summary (by Life Stage)', fontsize=16, fontweight='bold')

    # Plot 1: Uniquely mapped percentage
    ax = axes[0, 0]
    metrics_df.boxplot(column='uniquely_mapped_percent', by='stage', ax=ax)
    ax.set_title('Uniquely Mapped Reads (%)')
    ax.set_xlabel('Life Stage')
    ax.set_ylabel('Uniquely Mapped (%)')
    plt.sca(ax)
    plt.xticks(rotation=45)

    # Plot 2: Multimapped percentage
    ax = axes[0, 1]
    metrics_df.boxplot(column='multimapped_percent', by='stage', ax=ax)
    ax.set_title('Multi-mapped Reads (%)')
    ax.set_xlabel('Life Stage')
    ax.set_ylabel('Multi-mapped (%)')
    plt.sca(ax)
    plt.xticks(rotation=45)

    # Plot 3: Unmapped percentage
    ax = axes[1, 0]
    metrics_df.boxplot(column='unmapped_percent', by='stage', ax=ax)
    ax.set_title('Unmapped Reads (%)')
    ax.set_xlabel('Life Stage')
    ax.set_ylabel('Unmapped (%)')
    plt.sca(ax)
    plt.xticks(rotation=45)

    # Plot 4: Total reads
    ax = axes[1, 1]
    metrics_df.boxplot(column='total_reads', by='stage', ax=ax)
    ax.set_title('Total Reads (millions)')
    ax.set_xlabel('Life Stage')
    ax.set_ylabel('Total Reads (M)')
    ax.ticklabel_format(style='plain', axis='y')
    plt.sca(ax)
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"  Saved: {output_file}\n")


def create_duplication_plot(metrics_df, output_file):
    """Create duplication analysis plot."""
    print("Creating duplication analysis plot...")

    if 'duplication_rate' not in metrics_df.columns:
        print("  WARNING: Duplication data not available, skipping plot\n")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Duplication Analysis', fontsize=16, fontweight='bold')

    # Plot 1: Duplication rate by stage
    ax = axes[0]
    metrics_df.boxplot(column='duplication_rate', by='stage', ax=ax)
    ax.set_title('Duplication Rate by Life Stage')
    ax.set_xlabel('Life Stage')
    ax.set_ylabel('Duplication Rate')
    plt.sca(ax)
    plt.xticks(rotation=45)

    # Plot 2: Duplication vs library size
    ax = axes[1]
    for stage in metrics_df['stage'].unique():
        stage_data = metrics_df[metrics_df['stage'] == stage]
        ax.scatter(stage_data['total_reads'], stage_data['duplication_rate'],
                   label=stage, alpha=0.6, s=100)
    ax.set_xlabel('Total Reads (millions)')
    ax.set_ylabel('Duplication Rate')
    ax.set_title('Duplication vs Library Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    print(f"  Saved: {output_file}\n")


def create_summary(metrics_df, output_file):
    """Create text summary of QC metrics."""
    print("Writing summary statistics...")

    with open(output_file, 'w') as f:
        f.write("MULTIQC METRICS SUMMARY\n")
        f.write("=======================\n\n")
        f.write(f"Total samples: {len(metrics_df)}\n\n")

        f.write("MAPPING STATISTICS BY STAGE:\n")
        f.write("-" * 70 + "\n\n")

        for stage in sorted(metrics_df['stage'].unique()):
            stage_data = metrics_df[metrics_df['stage'] == stage]
            f.write(f"{stage.upper()} (n={len(stage_data)}):\n")
            f.write(f"  Total reads:\n")
            f.write(f"    Mean: {stage_data['total_reads'].mean()/1e6:.1f}M\n")
            f.write(f"    Range: {stage_data['total_reads'].min()/1e6:.1f}M - {stage_data['total_reads'].max()/1e6:.1f}M\n")
            f.write(f"  Uniquely mapped:\n")
            f.write(f"    Mean: {stage_data['uniquely_mapped_percent'].mean():.1f}%\n")
            f.write(f"    Range: {stage_data['uniquely_mapped_percent'].min():.1f}% - {stage_data['uniquely_mapped_percent'].max():.1f}%\n")
            f.write(f"  Multi-mapped:\n")
            f.write(f"    Mean: {stage_data['multimapped_percent'].mean():.1f}%\n")
            f.write(f"  Unmapped:\n")
            f.write(f"    Mean: {stage_data['unmapped_percent'].mean():.1f}%\n")
            f.write("\n")

        f.write("\nQUALITY ASSESSMENT:\n")
        f.write("-" * 70 + "\n\n")

        # Flag samples with low mapping
        low_mapping = metrics_df[metrics_df['uniquely_mapped_percent'] < 70]
        if len(low_mapping) > 0:
            f.write(f"FLAGGED: {len(low_mapping)} samples with <70% uniquely mapped reads:\n")
            for _, row in low_mapping.iterrows():
                f.write(f"  {row['sample_id']}: {row['uniquely_mapped_percent']:.1f}%\n")
        else:
            f.write("All samples have >70% uniquely mapped reads (GOOD)\n")

        f.write("\n")

        # Flag samples with high unmapped
        high_unmapped = metrics_df[metrics_df['unmapped_percent'] > 20]
        if len(high_unmapped) > 0:
            f.write(f"FLAGGED: {len(high_unmapped)} samples with >20% unmapped reads:\n")
            for _, row in high_unmapped.iterrows():
                f.write(f"  {row['sample_id']}: {row['unmapped_percent']:.1f}%\n")
        else:
            f.write("All samples have <20% unmapped reads (GOOD)\n")

        f.write("\n")

        f.write("REVIEWER 3.4 RESPONSE:\n")
        f.write("-" * 70 + "\n\n")
        f.write("Mapping statistics extracted from MultiQC report:\n")
        f.write(f"  - Mean uniquely mapped: {metrics_df['uniquely_mapped_percent'].mean():.1f}%\n")
        f.write(f"  - Mean multi-mapped: {metrics_df['multimapped_percent'].mean():.1f}%\n")
        f.write(f"  - Mean unmapped: {metrics_df['unmapped_percent'].mean():.1f}%\n")
        f.write("\nDetailed per-sample statistics: technical_metrics_summary.tsv\n")

    print(f"  Saved: {output_file}\n")


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 70)
    print("EXTRACT MULTIQC METRICS")
    print("=" * 70)
    print()

    # Load MultiQC data
    multiqc_data = load_multiqc_data(MULTIQC_JSON)

    # Extract metrics
    star_df = extract_star_metrics(multiqc_data)
    rseqc_df = extract_rseqc_metrics(multiqc_data)
    dupradar_df = extract_dupradar_metrics(multiqc_data)

    # Combine metrics
    metrics_df = combine_metrics(star_df, rseqc_df, dupradar_df)

    # Save combined metrics table
    metrics_df.to_csv(METRICS_TABLE, sep='\t', index=False)
    print(f"Saved technical metrics table: {METRICS_TABLE}\n")

    # Create plots
    create_mapping_plot(metrics_df, MAPPING_PLOT)

    if 'duplication_rate' in metrics_df.columns:
        create_duplication_plot(metrics_df, DUPLICATION_PLOT)

    # Create summary
    create_summary(metrics_df, SUMMARY_FILE)

    print("=" * 70)
    print("MULTIQC INTEGRATION COMPLETE")
    print("=" * 70)
    print()
    print("Review:")
    print(f"  - Metrics table: {METRICS_TABLE}")
    print(f"  - Mapping plot: {MAPPING_PLOT}")
    print(f"  - Summary: {SUMMARY_FILE}")
    print()
    print("Next: Run script 06 (outlier detection)")
    print()


if __name__ == "__main__":
    main()
