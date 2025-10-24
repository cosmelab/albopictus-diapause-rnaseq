#!/usr/bin/env python3
"""
Compare Salmon and featureCounts quantification methods

This script validates that Salmon pseudoalignment gives comparable results
to traditional alignment-based counting (featureCounts). This addresses
reviewer concerns about quantification method choice.

Expected correlation: R² > 0.95 for well-expressed genes

Author: L. Cosme
Date: October 22, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import sys

# Set style for publication-quality figures
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

def gather_featurecounts(featurecounts_dir):
    """
    Gather featureCounts output from individual sample files

    featureCounts outputs one file per sample with format:
    Geneid  Chr  Start  End  Strand  Length  Count
    """
    print("Gathering featureCounts data...")

    featurecounts_path = Path(featurecounts_dir)
    fc_files = sorted(featurecounts_path.glob("*.featureCounts.txt"))

    if not fc_files:
        raise FileNotFoundError(f"No featureCounts files found in {featurecounts_dir}")

    print(f"Found {len(fc_files)} featureCounts files")

    # Read first file to get gene IDs
    first_file = pd.read_csv(fc_files[0], sep='\t', comment='#')
    genes = first_file['Geneid'].values

    # Initialize matrix
    count_matrix = pd.DataFrame(index=genes)

    # Read all samples
    for fc_file in fc_files:
        sample_id = fc_file.stem.replace('.featureCounts', '')
        df = pd.read_csv(fc_file, sep='\t', comment='#')

        # The count column is the last column (6th column, 0-indexed position 6)
        count_matrix[sample_id] = df.iloc[:, 6].values

    print(f"featureCounts matrix: {count_matrix.shape[0]} genes × {count_matrix.shape[1]} samples")

    return count_matrix

def load_salmon_counts(salmon_file):
    """
    Load Salmon gene-level counts
    """
    print(f"\nLoading Salmon counts from {salmon_file}...")

    salmon_df = pd.read_csv(salmon_file, sep='\t', index_col=0)

    print(f"Salmon matrix: {salmon_df.shape[0]} genes × {salmon_df.shape[1]} samples")

    return salmon_df

def normalize_gene_ids(gene_ids):
    """
    Normalize gene IDs by removing prefixes like 'gene-', 'rna-'
    """
    return gene_ids.str.replace('gene-', '').str.replace('rna-', '').str.replace('transcript:', '')

def compare_methods(salmon_df, featurecounts_df, output_dir):
    """
    Compare Salmon and featureCounts quantifications
    """
    print("\nComparing methods...")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Normalize gene IDs
    salmon_df.index = normalize_gene_ids(salmon_df.index)
    featurecounts_df.index = normalize_gene_ids(featurecounts_df.index)

    # Find common genes and samples
    common_genes = salmon_df.index.intersection(featurecounts_df.index)
    common_samples = salmon_df.columns.intersection(featurecounts_df.columns)

    print(f"\nCommon genes: {len(common_genes)}")
    print(f"Common samples: {len(common_samples)}")

    if len(common_genes) == 0:
        print("\nERROR: No common genes found!")
        print("Salmon gene IDs (first 5):", salmon_df.index[:5].tolist())
        print("featureCounts gene IDs (first 5):", featurecounts_df.index[:5].tolist())
        sys.exit(1)

    if len(common_samples) == 0:
        print("\nERROR: No common samples found!")
        print("Salmon samples (first 5):", salmon_df.columns[:5].tolist())
        print("featureCounts samples (first 5):", featurecounts_df.columns[:5].tolist())
        sys.exit(1)

    # Subset to common genes and samples
    salmon_subset = salmon_df.loc[common_genes, common_samples]
    fc_subset = featurecounts_df.loc[common_genes, common_samples]

    # Calculate per-sample correlations
    correlations = []
    for sample in common_samples:
        r, p = stats.pearsonr(salmon_subset[sample], fc_subset[sample])
        correlations.append({
            'sample': sample,
            'pearson_r': r,
            'pearson_r2': r**2,
            'p_value': p
        })

    corr_df = pd.DataFrame(correlations)

    print(f"\nPearson correlation statistics:")
    print(f"  Mean R²: {corr_df['pearson_r2'].mean():.4f}")
    print(f"  Median R²: {corr_df['pearson_r2'].median():.4f}")
    print(f"  Min R²: {corr_df['pearson_r2'].min():.4f}")
    print(f"  Max R²: {corr_df['pearson_r2'].max():.4f}")

    # Save correlation statistics
    corr_df.to_csv(output_path / 'correlation_per_sample.tsv', sep='\t', index=False)

    # Create scatter plot for one representative sample
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    # Select 6 samples to plot
    samples_to_plot = common_samples[:6]

    for i, sample in enumerate(samples_to_plot):
        ax = axes[i]

        # Get data for this sample
        x = salmon_subset[sample].values
        y = fc_subset[sample].values

        # Filter for expressed genes (> 10 counts in both methods)
        mask = (x > 10) & (y > 10)
        x_filt = x[mask]
        y_filt = y[mask]

        # Log transform
        x_log = np.log10(x_filt + 1)
        y_log = np.log10(y_filt + 1)

        # Calculate correlation
        r, p = stats.pearsonr(x_log, y_log)

        # Scatter plot
        ax.scatter(x_log, y_log, alpha=0.5, s=10)
        ax.plot([0, x_log.max()], [0, x_log.max()], 'r--', lw=2, label='y=x')

        # Labels
        ax.set_xlabel('Salmon log10(counts + 1)', fontsize=10)
        ax.set_ylabel('featureCounts log10(counts + 1)', fontsize=10)
        ax.set_title(f'{sample}\nR² = {r**2:.4f}', fontsize=11)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path / 'salmon_vs_featurecounts_scatter.pdf', bbox_inches='tight')
    plt.savefig(output_path / 'salmon_vs_featurecounts_scatter.png', bbox_inches='tight', dpi=300)
    plt.close()

    print(f"\nScatter plots saved to {output_path}")

    # Create overall correlation plot (all genes, all samples)
    fig, ax = plt.subplots(figsize=(8, 8))

    # Flatten all values
    all_salmon = salmon_subset.values.flatten()
    all_fc = fc_subset.values.flatten()

    # Filter expressed genes
    mask = (all_salmon > 10) & (all_fc > 10)
    x_filt = all_salmon[mask]
    y_filt = all_fc[mask]

    # Log transform
    x_log = np.log10(x_filt + 1)
    y_log = np.log10(y_filt + 1)

    # Calculate correlation
    r, p = stats.pearsonr(x_log, y_log)

    # Hexbin plot for better visualization with many points
    hb = ax.hexbin(x_log, y_log, gridsize=100, cmap='Blues', mincnt=1)
    ax.plot([0, x_log.max()], [0, x_log.max()], 'r--', lw=2, label='y=x')

    # Labels
    ax.set_xlabel('Salmon log10(counts + 1)', fontsize=12)
    ax.set_ylabel('featureCounts log10(counts + 1)', fontsize=12)
    ax.set_title(f'Salmon vs featureCounts - All Samples\nR² = {r**2:.4f}, n={len(x_log):,} genes×samples', fontsize=14)

    # Colorbar
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Count', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path / 'salmon_vs_featurecounts_overall.pdf', bbox_inches='tight')
    plt.savefig(output_path / 'salmon_vs_featurecounts_overall.png', bbox_inches='tight', dpi=300)
    plt.close()

    print(f"Overall correlation plot saved to {output_path}")

    # Save summary statistics
    summary = {
        'common_genes': len(common_genes),
        'common_samples': len(common_samples),
        'mean_r2': corr_df['pearson_r2'].mean(),
        'median_r2': corr_df['pearson_r2'].median(),
        'min_r2': corr_df['pearson_r2'].min(),
        'max_r2': corr_df['pearson_r2'].max(),
        'overall_r2': r**2,
        'overall_p_value': p,
        'genes_expressed_both': mask.sum()
    }

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(output_path / 'comparison_summary.tsv', sep='\t', index=False)

    print(f"\nSummary statistics saved to {output_path / 'comparison_summary.tsv'}")

    return summary

def main():
    """Main analysis workflow"""

    # Define paths
    project_base = Path("/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq")

    featurecounts_dir = project_base / "results/featurecounts_genelevel"
    salmon_file = project_base / "results/02_count_matrices/salmon_gene_counts.tsv"
    output_dir = project_base / "results/07_method_comparison"

    print("="*70)
    print("Salmon vs featureCounts Comparison")
    print("="*70)
    print(f"\nProject base: {project_base}")
    print(f"featureCounts dir: {featurecounts_dir}")
    print(f"Salmon file: {salmon_file}")
    print(f"Output dir: {output_dir}")

    # Load data
    fc_df = gather_featurecounts(featurecounts_dir)
    salmon_df = load_salmon_counts(salmon_file)

    # Compare methods
    summary = compare_methods(salmon_df, fc_df, output_dir)

    print("\n" + "="*70)
    print("VALIDATION COMPLETE!")
    print("="*70)
    print(f"\nOverall R² = {summary['overall_r2']:.4f}")
    print(f"Mean per-sample R² = {summary['mean_r2']:.4f}")

    if summary['overall_r2'] > 0.95:
        print("\n✅ VALIDATION PASSED: R² > 0.95")
        print("Salmon and featureCounts show excellent agreement.")
        print("Either method is suitable for downstream analysis.")
    elif summary['overall_r2'] > 0.90:
        print("\n⚠️  VALIDATION ACCEPTABLE: 0.90 < R² < 0.95")
        print("Methods show good agreement but some discrepancies exist.")
    else:
        print("\n❌ VALIDATION FAILED: R² < 0.90")
        print("Significant discrepancies between methods - investigate!")

    print(f"\nResults saved to: {output_dir}")
    print("="*70)

if __name__ == "__main__":
    main()
