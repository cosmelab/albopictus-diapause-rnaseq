#!/usr/bin/env python3
"""
Generate publication-ready QC figures from extracted metrics
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

# Set up plotting style
plt.style.use('seaborn-v0_8-white')
sns.set_palette("husl")

# Define project colors (ordered: Adults → Embryos → Larvae)
PROJECT_COLORS = {
    'PRJNA268379': '#1f77b4',  # Blue - Adult females
    'PRJNA158021': '#ff7f0e',  # Orange - Embryos
    'PRJNA187045': '#2ca02c'   # Green - Pharate larvae
}

# Define project labels (ordered: Adults → Embryos → Larvae)
PROJECT_LABELS = {
    'PRJNA268379': 'Adult Females',
    'PRJNA158021': 'Embryos',
    'PRJNA187045': 'Pharate Larvae'
}

# Define consistent stage order
STAGE_ORDER = ['PRJNA268379', 'PRJNA158021', 'PRJNA187045']

def thousands_formatter(x, pos):
    """Format large numbers with thousand separators"""
    return f'{int(x):,}'

def abbreviated_formatter(value):
    """Format numbers with k abbreviation (e.g., 19.8k)"""
    if value >= 1000:
        return f'{value/1000:.1f}k'
    return f'{int(value)}'

def create_figure_1_overview(df, output_dir):
    """Create Figure 1: QC Overview (4 panels)"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Figure 1: RNA-seq Quality Control Overview', fontsize=16, y=0.98)

    # Sort by stage order: Adults → Embryos → Larvae
    df['stage_order'] = df['project_id'].map({p: i for i, p in enumerate(STAGE_ORDER)})
    df_sorted = df.sort_values(['stage_order', 'sample_id'])

    # Panel A: Sequencing depth (grouped by stage)
    ax = axes[0, 0]
    x = range(len(df_sorted))
    bars = ax.bar(x, df_sorted['total_reads_millions'],
                   color=[PROJECT_COLORS[p] for p in df_sorted['project_id']])
    ax.axhline(y=20, color='red', linestyle='--', alpha=0.5, label='20M threshold')

    # Add vertical lines to separate stages
    stage_boundaries = []
    for i, stage in enumerate(STAGE_ORDER[:-1]):
        n_samples = (df_sorted['project_id'] == stage).sum()
        if i == 0:
            boundary = n_samples - 0.5
        else:
            boundary = stage_boundaries[-1] + n_samples
        stage_boundaries.append(boundary)
        ax.axvline(x=boundary, color='black', linestyle='-', linewidth=1.5, alpha=0.3)

    # Add stage labels at bottom
    stage_positions = []
    cumsum = 0
    for stage in STAGE_ORDER:
        n = (df_sorted['project_id'] == stage).sum()
        stage_positions.append(cumsum + n/2 - 0.5)
        cumsum += n
    ax.set_xticks(stage_positions)
    ax.set_xticklabels([PROJECT_LABELS[s] for s in STAGE_ORDER], fontsize=10)

    # Format y-axis with thousand separators (for readability in millions)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))

    ax.set_ylabel('Total Reads (millions)')
    ax.set_title('A. Sequencing Depth by Life Stage')
    ax.legend()

    # Panel B: Mapping rates (grouped stacked bars)
    ax = axes[0, 1]
    width = 0.8

    # Calculate unmapped percentage
    df_sorted['unmapped_percent'] = 100 - df_sorted['reads_mapped_percent']

    # Create stacked bars
    p1 = ax.bar(x, df_sorted['uniquely_mapped_percent'], width,
                label='Uniquely Mapped', color='#2ecc71')
    p2 = ax.bar(x, df_sorted['multi_mapped_percent'], width,
                bottom=df_sorted['uniquely_mapped_percent'],
                label='Multi-mapped', color='#f39c12')
    p3 = ax.bar(x, df_sorted['unmapped_percent'], width,
                bottom=df_sorted['uniquely_mapped_percent'] + df_sorted['multi_mapped_percent'],
                label='Unmapped', color='#e74c3c')

    # Add stage separators
    for boundary in stage_boundaries:
        ax.axvline(x=boundary, color='black', linestyle='-', linewidth=1.5, alpha=0.3)

    ax.set_xticks(stage_positions)
    ax.set_xticklabels([PROJECT_LABELS[s] for s in STAGE_ORDER], fontsize=10)
    ax.set_ylabel('Percentage (%)')
    ax.set_title('B. Read Mapping Distribution by Life Stage')
    ax.legend()
    ax.set_ylim(0, 100)

    # Panel C: Read distribution (genomic features, grouped)
    ax = axes[1, 0]

    # Create stacked bars for genomic features
    p1 = ax.bar(x, df_sorted['exonic_percent'], width,
                label='Exonic', color='#3498db')
    p2 = ax.bar(x, df_sorted['intronic_percent'], width,
                bottom=df_sorted['exonic_percent'],
                label='Intronic', color='#9b59b6')
    p3 = ax.bar(x, df_sorted['intergenic_percent'], width,
                bottom=df_sorted['exonic_percent'] + df_sorted['intronic_percent'],
                label='Intergenic', color='#95a5a6')

    # Add stage separators
    for boundary in stage_boundaries:
        ax.axvline(x=boundary, color='black', linestyle='-', linewidth=1.5, alpha=0.3)

    ax.set_xticks(stage_positions)
    ax.set_xticklabels([PROJECT_LABELS[s] for s in STAGE_ORDER], fontsize=10)
    ax.set_ylabel('Percentage (%)')
    ax.set_title('C. Genomic Feature Distribution by Life Stage')
    ax.legend()

    # Panel D: Gene detection (boxplot by stage)
    ax = axes[1, 1]

    # Box plot by project (in correct order)
    gene_det_data = []
    gene_means = []
    for project in STAGE_ORDER:
        project_data = df_sorted[df_sorted['project_id'] == project]['genes_detected_gt10'].values
        gene_det_data.append(project_data)
        gene_means.append(project_data.mean())

    bp = ax.boxplot(gene_det_data, labels=[PROJECT_LABELS[p] for p in STAGE_ORDER],
                    patch_artist=True)

    # Color boxes by project
    for patch, project in zip(bp['boxes'], STAGE_ORDER):
        patch.set_facecolor(PROJECT_COLORS[project])

    # Add mean annotations inside boxes
    for i, (mean_val, project) in enumerate(zip(gene_means, STAGE_ORDER)):
        # Abbreviated format
        mean_text = abbreviated_formatter(mean_val)
        # Position at median height, center of box
        ax.text(i + 1, mean_val, mean_text,
                ha='center', va='center', fontsize=10, fontweight='bold',
                color='white',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=PROJECT_COLORS[project],
                         edgecolor='none', alpha=0.9))

    # Format y-axis with thousand separators
    ax.yaxis.set_major_formatter(FuncFormatter(thousands_formatter))

    ax.set_ylabel('Number of Genes (counts > 10)')
    ax.set_title('D. Gene Detection by Life Stage')
    ax.tick_params(axis='x', rotation=45)

    # Add project legend
    legend_elements = [Patch(facecolor=PROJECT_COLORS[p], label=PROJECT_LABELS[p])
                      for p in STAGE_ORDER]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, -0.02), ncol=3)

    plt.tight_layout()
    plt.savefig(output_dir / 'figure_1_qc_overview.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(output_dir / 'figure_1_qc_overview.png', bbox_inches='tight', dpi=300)
    plt.close()

def create_figure_2_detailed_qc(df, output_dir):
    """Create Figure 2: Detailed QC Analysis"""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 2: Detailed Quality Control Analysis', fontsize=16, y=0.98)

    # Panel A: 5'/3' bias by project (in stage order)
    ax = axes[0, 0]
    for project in STAGE_ORDER:
        project_data = df[df['project_id'] == project]
        ax.scatter(range(len(project_data)), project_data['five_three_bias'],
                  label=PROJECT_LABELS[project], color=PROJECT_COLORS[project],
                  s=50, alpha=0.7)
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='No bias')
    ax.set_xlabel('Samples')
    ax.set_ylabel("5'/3' Bias Ratio")
    ax.set_title("A. Gene Coverage Bias")
    ax.legend()

    # Panel B: rRNA contamination (in stage order)
    ax = axes[0, 1]
    rrna_data = []
    for project in STAGE_ORDER:
        project_data = df[df['project_id'] == project]['rrna_percent'].values
        rrna_data.append(project_data)

    bp = ax.boxplot(rrna_data, labels=[PROJECT_LABELS[p] for p in STAGE_ORDER],
                    patch_artist=True)

    for patch, project in zip(bp['boxes'], STAGE_ORDER):
        patch.set_facecolor(PROJECT_COLORS[project])

    ax.set_ylabel('rRNA Contamination (%)')
    ax.set_title('B. rRNA Content by Life Stage\n(STAR default: 0% all samples)')
    ax.tick_params(axis='x', rotation=45)

    # Panel C: Insert size distribution (in stage order)
    ax = axes[1, 0]
    for i, project in enumerate(STAGE_ORDER):
        project_data = df[df['project_id'] == project]
        ax.violinplot([project_data['insert_size_mean'].values],
                     positions=[i],
                     widths=0.7, showmeans=True, showextrema=True)

    ax.set_xticks(range(len(STAGE_ORDER)))
    ax.set_xticklabels([PROJECT_LABELS[p] for p in STAGE_ORDER], rotation=45)

    # Format y-axis with thousand separators
    ax.yaxis.set_major_formatter(FuncFormatter(thousands_formatter))

    ax.set_ylabel('Insert Size (bp)')
    ax.set_title('C. Library Insert Size Distribution')

    # Panel D: Duplication rate vs library size (in stage order)
    ax = axes[1, 1]
    for project in STAGE_ORDER:
        project_data = df[df['project_id'] == project]
        ax.scatter(project_data['total_reads_millions'], project_data['duplication_percent'],
                  label=PROJECT_LABELS[project], color=PROJECT_COLORS[project],
                  s=50, alpha=0.7)

    ax.set_xlabel('Library Size (million reads)')
    ax.set_ylabel('Duplication Rate (%)')
    ax.set_title('D. Duplication vs Library Size')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'figure_2_detailed_qc.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(output_dir / 'figure_2_detailed_qc.png', bbox_inches='tight', dpi=300)
    plt.close()

def create_correlation_heatmap(df, output_dir):
    """Create correlation heatmap of QC metrics"""

    # Select numeric columns for correlation
    qc_metrics = [
        'total_reads_millions', 'uniquely_mapped_percent', 'multi_mapped_percent',
        'exonic_percent', 'intronic_percent', 'intergenic_percent',
        'rrna_percent', 'five_three_bias', 'duplication_percent',
        'gc_percent', 'error_rate_percent', 'genes_detected_gt10'
    ]

    # Filter to available columns
    available_metrics = [m for m in qc_metrics if m in df.columns]

    # Remove columns with zero variance (can't compute correlation)
    metrics_with_variance = []
    for metric in available_metrics:
        if df[metric].std() > 0:
            metrics_with_variance.append(metric)
        else:
            print(f"  Note: Excluding {metric} from correlation (zero variance: all values = {df[metric].iloc[0]:.1f})")

    # Calculate correlation matrix
    corr_matrix = df[metrics_with_variance].corr()

    # Create nicer labels
    label_map = {
        'total_reads_millions': 'Total Reads (M)',
        'uniquely_mapped_percent': 'Uniquely Mapped (%)',
        'multi_mapped_percent': 'Multi-mapped (%)',
        'exonic_percent': 'Exonic (%)',
        'intronic_percent': 'Intronic (%)',
        'intergenic_percent': 'Intergenic (%)',
        'rrna_percent': 'rRNA (%)',
        'five_three_bias': "5'/3' Bias",
        'duplication_percent': 'Duplication (%)',
        'gc_percent': 'GC Content (%)',
        'error_rate_percent': 'Error Rate (%)',
        'genes_detected_gt10': 'Genes Detected (>10)'
    }

    # Rename for display
    corr_matrix_labeled = corr_matrix.copy()
    corr_matrix_labeled.columns = [label_map.get(c, c) for c in corr_matrix.columns]
    corr_matrix_labeled.index = [label_map.get(i, i) for i in corr_matrix.index]

    # Create figure
    plt.figure(figsize=(11, 9))

    # Create heatmap
    sns.heatmap(corr_matrix_labeled, annot=True, fmt='.2f', cmap='coolwarm',
                center=0, vmin=-1, vmax=1, square=True,
                cbar_kws={"shrink": 0.8},
                annot_kws={"size": 8})

    plt.title('QC Metrics Correlation Matrix',
              fontsize=12, pad=15)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()

    plt.savefig(output_dir / 'figure_s1_correlation_heatmap.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(output_dir / 'figure_s1_correlation_heatmap.png', bbox_inches='tight', dpi=300)
    plt.close()

def create_summary_table(df, output_dir):
    """Create publication-ready summary table"""
    
    # Define metrics to include in the table
    metrics = {
        'total_reads_millions': 'Total Reads (M)',
        'uniquely_mapped_percent': 'Uniquely Mapped (%)',
        'multi_mapped_percent': 'Multi-mapped (%)',
        'exonic_percent': 'Exonic Reads (%)',
        'intronic_percent': 'Intronic Reads (%)',
        'intergenic_percent': 'Intergenic Reads (%)',
        'rrna_percent': 'rRNA (%)',
        'genes_detected_gt10': 'Genes Detected (>10 counts)',
        'five_three_bias': "5'/3' Bias"
    }
    
    # Create summary by project (in stage order: Adults → Embryos → Larvae)
    summary_data = []

    for project in STAGE_ORDER:
        project_data = df[df['project_id'] == project]

        row = {
            'Life Stage': PROJECT_LABELS[project],
            'n': len(project_data)
        }

        for metric, label in metrics.items():
            if metric in project_data.columns:
                values = project_data[metric]
                row[label] = f"{values.mean():.1f} ± {values.std():.1f}"

        summary_data.append(row)
    
    # Create DataFrame
    summary_df = pd.DataFrame(summary_data)
    
    # Save as CSV
    summary_df.to_csv(output_dir / 'table_1_qc_summary.csv', index=False)

    # Also save as LaTeX (if jinja2 available)
    try:
        latex_table = summary_df.to_latex(index=False, escape=False,
                                         caption="Summary of RNA-seq quality control metrics by life stage",
                                         label="tab:qc_summary")

        with open(output_dir / 'table_1_qc_summary.tex', 'w') as f:
            f.write(latex_table)
    except ImportError:
        print("Note: LaTeX table generation skipped (jinja2 not available)")

    return summary_df

def main():
    """Main function"""
    
    # Check if QC metrics file exists
    qc_file = Path("results/01_quality_control/metrics/all_samples_qc_metrics_fixed.csv")
    if not qc_file.exists():
        print(f"Error: QC metrics file not found at {qc_file}")
        print("Please run 01_extract_qc_metrics.py first")
        return
    
    # Load QC data
    df = pd.read_csv(qc_file)
    print(f"Loaded {len(df)} samples")
    
    # Create output directory for figures
    figures_dir = Path("results/01_quality_control/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    # Create tables directory
    tables_dir = Path("results/01_quality_control/tables")
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate figures
    print("Creating Figure 1: QC Overview...")
    create_figure_1_overview(df, figures_dir)
    
    print("Creating Figure 2: Detailed QC...")
    create_figure_2_detailed_qc(df, figures_dir)
    
    print("Creating Correlation Heatmap...")
    create_correlation_heatmap(df, figures_dir)
    
    print("Creating Summary Table...")
    summary_df = create_summary_table(df, tables_dir)
    print("\nSummary Table:")
    print(summary_df)
    
    print(f"\nAll figures saved to {figures_dir}")
    print(f"Tables saved to {tables_dir}")

if __name__ == "__main__":
    main()