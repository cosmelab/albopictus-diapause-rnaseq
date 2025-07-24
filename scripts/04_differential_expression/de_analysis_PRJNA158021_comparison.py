#!/usr/bin/env python3
"""
Differential Expression Analysis Comparison for PRJNA158021 (Embryo Dataset)

This script performs both pairwise comparisons (like Poelchau et al. 2013a) and 
complex interaction analysis to compare results and understand why we find 
different numbers of significant genes.

Approach 1: Pairwise comparisons (Poelchau et al. style)
- DI vs NDI at 72h separately
- DI vs NDI at 135h separately

Approach 2: Complex interaction analysis (current approach)
- Full model with timepoint × photoperiod interactions

Usage:
    python de_analysis_PRJNA158021_comparison.py

Requirements:
    - pandas, numpy, matplotlib, seaborn
    - rpy2 for R integration
    - DESeq2 R package
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import warnings
warnings.filterwarnings('ignore')

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('de_analysis_comparison_PRJNA158021.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class ComparisonAnalyzer:
    def __init__(self, gene_counts_file="output/tables/PRJNA158021_gene_counts.csv", 
                 output_base="output"):
        """
        Initialize the comparison analyzer for embryo dataset.
        
        Args:
            gene_counts_file: Path to the gene counts file
            output_base: Base directory for output files
        """
        self.gene_counts_file = Path(gene_counts_file)
        self.output_base = Path(output_base)
        
        # Create output directories
        self.de_dir = self.output_base / "differential_expression" / "PRJNA158021_comparison"
        self.figures_dir = self.output_base / "figures" / "PRJNA158021_comparison"
        
        for dir_path in [self.de_dir, self.figures_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.gene_counts = None
        self.metadata = None
        
        # Results containers
        self.pairwise_results = {}  # Store pairwise comparison results
        self.interaction_results = None  # Store interaction model results
        self.comparison_summary = {}  # Store comparison statistics
        
        # Setup R environment
        self._setup_r_environment()
        
        logger.info("ComparisonAnalyzer initialized for PRJNA158021")

    def _setup_r_environment(self):
        """Setup R environment and load required packages."""
        try:
            # Import R packages
            self.deseq2 = importr('DESeq2')
            self.stats = importr('stats')
            
            # Enable automatic conversion between pandas and R
            pandas2ri.activate()
            
            logger.info("R environment setup completed")
            
        except Exception as e:
            logger.error(f"Error setting up R environment: {e}")
            raise

    def load_data(self):
        """Load gene counts and create metadata."""
        if not self.gene_counts_file.exists():
            raise FileNotFoundError(f"Gene counts file not found: {self.gene_counts_file}")
        
        # Load gene counts
        self.gene_counts = pd.read_csv(self.gene_counts_file, index_col=0)
        logger.info(f"Loaded gene counts: {self.gene_counts.shape[0]} genes, {self.gene_counts.shape[1]} samples")
        
        # Create metadata
        self._create_metadata()

    def _create_metadata(self):
        """Create metadata DataFrame from sample names."""
        sample_names = self.gene_counts.columns.tolist()
        
        metadata_list = []
        for sample in sample_names:
            # Parse sample name: PRJNA158021_SRR458462_LD_72h
            parts = sample.split('_')
            if len(parts) >= 4:
                prjna = parts[0]
                srr = parts[1]
                light_cycle = parts[2]  # LD or SD
                time_point = parts[3]   # 72h or 135h
                
                # Create treatment column
                treatment = f"{light_cycle}_{time_point}"
                
                # Add biological annotations
                diapause_status = 'Diapause_induction' if light_cycle == 'SD' else 'Non_diapause'
                timepoint_clean = time_point.replace('h', 'h_pov')
                timepoint_num = 72 if time_point == '72h' else 135
                
                metadata_list.append({
                    'sample': sample,
                    'prjna': prjna,
                    'srr': srr,
                    'light_cycle': light_cycle,
                    'time_point': time_point,
                    'treatment': treatment,
                    'diapause_status': diapause_status,
                    'timepoint_clean': timepoint_clean,
                    'timepoint_num': timepoint_num
                })
        
        self.metadata = pd.DataFrame(metadata_list)
        self.metadata.set_index('sample', inplace=True)
        
        logger.info(f"Created metadata for {len(self.metadata)} samples")
        logger.info(f"Sample distribution:")
        logger.info(self.metadata.groupby(['light_cycle', 'time_point']).size())

    def run_pairwise_analysis(self):
        """Run pairwise comparisons like Poelchau et al. 2013a."""
        logger.info("Running pairwise comparisons (Poelchau et al. style)...")
        
        # Prepare count data for R
        count_matrix = self.gene_counts.astype(int)
        r_counts = pandas2ri.py2rpy(count_matrix)
        ro.globalenv['r_counts'] = r_counts
        
        # Run analysis for each timepoint separately
        for timepoint in ['72h', '135h']:
            logger.info(f"Analyzing {timepoint} timepoint...")
            
            # Filter samples for this timepoint
            timepoint_samples = self.metadata[self.metadata['time_point'] == timepoint]
            timepoint_counts = count_matrix[timepoint_samples.index]
            
            # Create metadata for this timepoint
            timepoint_metadata = timepoint_samples.copy()
            timepoint_metadata['sample'] = timepoint_metadata.index
            
            # Convert to R objects
            r_timepoint_counts = pandas2ri.py2rpy(timepoint_counts)
            r_timepoint_metadata = pandas2ri.py2rpy(timepoint_metadata)
            
            # Assign to R global environment
            ro.globalenv['r_timepoint_counts'] = r_timepoint_counts
            ro.globalenv['r_timepoint_metadata'] = r_timepoint_metadata
            
            # Run DESeq2 for this timepoint
            r_code = f'''
            # Create DESeq2 dataset for {timepoint}
            dds_{timepoint} <- DESeqDataSetFromMatrix(
                countData = r_timepoint_counts,
                colData = r_timepoint_metadata,
                design = ~ light_cycle
            )
            
            # Set reference level
            dds_{timepoint}$light_cycle <- relevel(dds_{timepoint}$light_cycle, ref = "LD")
            
            # Print design matrix
            print("Design matrix for {timepoint}:")
            print(design(dds_{timepoint}))
            
            # Run DESeq2
            dds_{timepoint} <- DESeq(dds_{timepoint})
            
            # Get results (SD vs LD)
            results_{timepoint} <- results(dds_{timepoint}, name = "light_cycle_SD_vs_LD")
            results_{timepoint} <- results_{timepoint}[order(results_{timepoint}$padj), ]
            
            # Print summary
            print("Results summary for {timepoint}:")
            print(summary(results_{timepoint}))
            
            # Store results
            .GlobalEnv$results_{timepoint} <- results_{timepoint}
            .GlobalEnv$dds_{timepoint} <- dds_{timepoint}
            '''
            
            ro.r(r_code)
            
            # Extract results
            r_results = ro.r[f'results_{timepoint}']
            # Fix rpy2 conversion for DESeq2 results objects
            r_df = ro.r['as.data.frame'](r_results)
            results_df = pd.DataFrame({
                'baseMean': list(r_df.rx2('baseMean')),
                'log2FoldChange': list(r_df.rx2('log2FoldChange')),
                'lfcSE': list(r_df.rx2('lfcSE')),
                'stat': list(r_df.rx2('stat')),
                'pvalue': list(r_df.rx2('pvalue')),
                'padj': list(r_df.rx2('padj'))
            }, index=list(ro.r['rownames'](r_results)))
            results_df['gene_id'] = results_df.index
            results_df['comparison'] = f'SD_vs_LD_{timepoint}'
            
            # Add significance columns
            results_df['significant'] = (
                (results_df['padj'] < 0.05) & 
                (abs(results_df['log2FoldChange']) > 0.5)  # Poelchau et al. used 0.5
            )
            
            # Store results
            self.pairwise_results[timepoint] = results_df
            
            # Log summary
            sig_genes = results_df['significant'].sum()
            logger.info(f"{timepoint}: {sig_genes} significant genes (p < 0.05, |log2FC| > 0.5)")
        
        # Save pairwise results
        all_pairwise = pd.concat(self.pairwise_results.values(), ignore_index=True)
        all_pairwise.to_csv(self.de_dir / 'pairwise_comparisons.csv', index=False)
        
        logger.info("Pairwise analysis completed")

    def run_interaction_analysis(self):
        """Run complex interaction model analysis."""
        logger.info("Running interaction model analysis...")
        
        # Prepare count data for R
        count_matrix = self.gene_counts.astype(int)
        design_data = self.metadata.copy()
        design_data['sample'] = design_data.index
        
        # Convert to R objects
        r_counts = pandas2ri.py2rpy(count_matrix)
        r_design = pandas2ri.py2rpy(design_data)
        
        # Assign to R global environment
        ro.globalenv['r_counts'] = r_counts
        ro.globalenv['r_design'] = r_design
        
        # Run DESeq2 with interaction design
        r_code = '''
        # Create DESeq2 dataset with interaction design
        dds_interaction <- DESeqDataSetFromMatrix(
            countData = r_counts,
            colData = r_design,
            design = ~ time_point + light_cycle + time_point:light_cycle
        )
        
        # Set reference levels
        dds_interaction$light_cycle <- relevel(dds_interaction$light_cycle, ref = "LD")
        dds_interaction$time_point <- relevel(dds_interaction$time_point, ref = "72h")
        
        # Print design matrix
        print("Interaction design matrix:")
        print(design(dds_interaction))
        
        # Run DESeq2
        dds_interaction <- DESeq(dds_interaction)
        
        # Get results for interaction effect
        results_interaction <- results(dds_interaction, name = "time_point135h.light_cycleSD")
        results_interaction <- results_interaction[order(results_interaction$padj), ]
        
        # Print summary
        print("Interaction results summary:")
        print(summary(results_interaction))
        
        # Store results
        .GlobalEnv$results_interaction <- results_interaction
        .GlobalEnv$dds_interaction <- dds_interaction
        '''
        
        ro.r(r_code)
        
        # Extract results
        r_results = ro.r['results_interaction']
        # Fix rpy2 conversion for DESeq2 results objects
        r_df = ro.r['as.data.frame'](r_results)
        results_df = pd.DataFrame({
            'baseMean': list(r_df.rx2('baseMean')),
            'log2FoldChange': list(r_df.rx2('log2FoldChange')),
            'lfcSE': list(r_df.rx2('lfcSE')),
            'stat': list(r_df.rx2('stat')),
            'pvalue': list(r_df.rx2('pvalue')),
            'padj': list(r_df.rx2('padj'))
        }, index=list(ro.r['rownames'](r_results)))
        results_df['gene_id'] = results_df.index
        results_df['comparison'] = 'Time_x_Photoperiod_Interaction'
        
        # Add significance columns
        results_df['significant'] = (
            (results_df['padj'] < 0.05) & 
            (abs(results_df['log2FoldChange']) > 0.5)
        )
        
        # Store results
        self.interaction_results = results_df
        
        # Save interaction results
        results_df.to_csv(self.de_dir / 'interaction_analysis.csv', index=False)
        
        # Log summary
        sig_genes = results_df['significant'].sum()
        logger.info(f"Interaction model: {sig_genes} significant genes (p < 0.05, |log2FC| > 0.5)")
        
        logger.info("Interaction analysis completed")

    def compare_results(self):
        """Compare results between pairwise and interaction approaches."""
        logger.info("Comparing results between approaches...")
        
        # Get gene sets
        pairwise_72h_genes = set(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]['gene_id'])
        pairwise_135h_genes = set(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]['gene_id'])
        interaction_genes = set(self.interaction_results[self.interaction_results['significant']]['gene_id'])
        
        # Calculate overlaps
        overlap_72h_135h = pairwise_72h_genes.intersection(pairwise_135h_genes)
        overlap_72h_interaction = pairwise_72h_genes.intersection(interaction_genes)
        overlap_135h_interaction = pairwise_135h_genes.intersection(interaction_genes)
        
        # Store comparison summary
        self.comparison_summary = {
            'pairwise_72h_total': len(pairwise_72h_genes),
            'pairwise_135h_total': len(pairwise_135h_genes),
            'interaction_total': len(interaction_genes),
            'overlap_72h_135h': len(overlap_72h_135h),
            'overlap_72h_interaction': len(overlap_72h_interaction),
            'overlap_135h_interaction': len(overlap_135h_interaction),
            'pairwise_72h_genes': list(pairwise_72h_genes),
            'pairwise_135h_genes': list(pairwise_135h_genes),
            'interaction_genes': list(interaction_genes)
        }
        
        # Log comparison
        logger.info("=== COMPARISON SUMMARY ===")
        logger.info(f"Pairwise 72h: {len(pairwise_72h_genes)} genes")
        logger.info(f"Pairwise 135h: {len(pairwise_135h_genes)} genes")
        logger.info(f"Interaction model: {len(interaction_genes)} genes")
        logger.info(f"Overlap 72h-135h: {len(overlap_72h_135h)} genes")
        logger.info(f"Overlap 72h-interaction: {len(overlap_72h_interaction)} genes")
        logger.info(f"Overlap 135h-interaction: {len(overlap_135h_interaction)} genes")
        
        # Save comparison summary
        comparison_df = pd.DataFrame([{
            'Comparison': 'Pairwise_72h',
            'Significant_Genes': len(pairwise_72h_genes),
            'Description': 'SD vs LD at 72h (Poelchau et al. style)'
        }, {
            'Comparison': 'Pairwise_135h', 
            'Significant_Genes': len(pairwise_135h_genes),
            'Description': 'SD vs LD at 135h (Poelchau et al. style)'
        }, {
            'Comparison': 'Interaction_Model',
            'Significant_Genes': len(interaction_genes),
            'Description': 'Time × Photoperiod interaction'
        }])
        
        comparison_df.to_csv(self.de_dir / 'comparison_summary.csv', index=False)
        
        # Create Venn diagram
        self._create_venn_diagram()
        
        logger.info("Comparison completed")

    def _create_venn_diagram(self):
        """Create Venn diagram showing overlaps between approaches."""
        try:
            from matplotlib_venn import venn3
            
            # Get gene sets
            pairwise_72h_genes = set(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]['gene_id'])
            pairwise_135h_genes = set(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]['gene_id'])
            interaction_genes = set(self.interaction_results[self.interaction_results['significant']]['gene_id'])
            
            # Create Venn diagram
            plt.figure(figsize=(10, 8))
            venn3([pairwise_72h_genes, pairwise_135h_genes, interaction_genes], 
                  ('Pairwise 72h', 'Pairwise 135h', 'Interaction'))
            plt.title('Overlap of Significant Genes Between Analysis Approaches\n(PRJNA158021 Embryo Dataset)', 
                     fontsize=14, fontweight='bold')
            plt.savefig(self.figures_dir / 'venn_diagram_comparison.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info("Venn diagram created")
            
        except ImportError:
            logger.warning("matplotlib_venn not available - skipping Venn diagram")
        except Exception as e:
            logger.error(f"Error creating Venn diagram: {e}")

    def create_comparison_plots(self):
        """Create comparison plots between approaches."""
        logger.info("Creating comparison plots...")
        
        # Create volcano plots for each approach
        self._create_volcano_plots()
        
        # Create summary bar plot
        self._create_summary_bar_plot()
        
        logger.info("Comparison plots created")

    def _create_volcano_plots(self):
        """Create volcano plots for each approach."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Pairwise 72h
        ax1 = axes[0, 0]
        data = self.pairwise_results['72h'].copy()
        # Handle NaN values
        data = data.dropna(subset=['log2FoldChange', 'padj'])
        data['padj'] = data['padj'].replace(0, 1e-300)  # Avoid log(0)
        
        ax1.scatter(data['log2FoldChange'], -np.log10(data['padj']), alpha=0.6, s=1)
        sig_data = data[data['significant']]
        if len(sig_data) > 0:
            ax1.scatter(sig_data['log2FoldChange'], -np.log10(sig_data['padj']), 
                       color='red', alpha=0.8, s=2, label=f'Significant (n={len(sig_data)})')
        ax1.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.7)
        ax1.axvline(0.5, color='gray', linestyle='--', alpha=0.7)
        ax1.axvline(-0.5, color='gray', linestyle='--', alpha=0.7)
        ax1.set_xlabel('log2 Fold Change')
        ax1.set_ylabel('-log10(adjusted p-value)')
        ax1.set_title('Pairwise: SD vs LD at 72h')
        if len(sig_data) > 0:
            ax1.legend()
        
        # Pairwise 135h
        ax2 = axes[0, 1]
        data = self.pairwise_results['135h'].copy()
        # Handle NaN values
        data = data.dropna(subset=['log2FoldChange', 'padj'])
        data['padj'] = data['padj'].replace(0, 1e-300)  # Avoid log(0)
        
        ax2.scatter(data['log2FoldChange'], -np.log10(data['padj']), alpha=0.6, s=1)
        sig_data = data[data['significant']]
        if len(sig_data) > 0:
            ax2.scatter(sig_data['log2FoldChange'], -np.log10(sig_data['padj']), 
                       color='red', alpha=0.8, s=2, label=f'Significant (n={len(sig_data)})')
        ax2.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.7)
        ax2.axvline(0.5, color='gray', linestyle='--', alpha=0.7)
        ax2.axvline(-0.5, color='gray', linestyle='--', alpha=0.7)
        ax2.set_xlabel('log2 Fold Change')
        ax2.set_ylabel('-log10(adjusted p-value)')
        ax2.set_title('Pairwise: SD vs LD at 135h')
        if len(sig_data) > 0:
            ax2.legend()
        
        # Interaction model
        ax3 = axes[1, 0]
        data = self.interaction_results.copy()
        # Handle NaN values
        data = data.dropna(subset=['log2FoldChange', 'padj'])
        data['padj'] = data['padj'].replace(0, 1e-300)  # Avoid log(0)
        
        ax3.scatter(data['log2FoldChange'], -np.log10(data['padj']), alpha=0.6, s=1)
        sig_data = data[data['significant']]
        if len(sig_data) > 0:
            ax3.scatter(sig_data['log2FoldChange'], -np.log10(sig_data['padj']), 
                       color='red', alpha=0.8, s=2, label=f'Significant (n={len(sig_data)})')
        ax3.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.7)
        ax3.axvline(0.5, color='gray', linestyle='--', alpha=0.7)
        ax3.axvline(-0.5, color='gray', linestyle='--', alpha=0.7)
        ax3.set_xlabel('log2 Fold Change')
        ax3.set_ylabel('-log10(adjusted p-value)')
        ax3.set_title('Interaction: Time × Photoperiod')
        if len(sig_data) > 0:
            ax3.legend()
        
        # Summary comparison
        ax4 = axes[1, 1]
        comparisons = ['Pairwise 72h', 'Pairwise 135h', 'Interaction']
        counts = [
            len(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]),
            len(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]),
            len(self.interaction_results[self.interaction_results['significant']])
        ]
        bars = ax4.bar(comparisons, counts, color=['skyblue', 'lightgreen', 'salmon'])
        ax4.set_ylabel('Number of Significant Genes')
        ax4.set_title('Comparison of Significant Genes')
        ax4.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for bar, count in zip(bars, counts):
            ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10, 
                    str(count), ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'volcano_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Volcano plots created")

    def _create_summary_bar_plot(self):
        """Create summary bar plot comparing approaches."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Bar plot of significant genes
        comparisons = ['Pairwise 72h', 'Pairwise 135h', 'Interaction']
        counts = [
            len(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]),
            len(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]),
            len(self.interaction_results[self.interaction_results['significant']])
        ]
        
        bars = ax1.bar(comparisons, counts, color=['skyblue', 'lightgreen', 'salmon'])
        ax1.set_ylabel('Number of Significant Genes')
        ax1.set_title('Significant Genes by Analysis Approach')
        ax1.tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, count in zip(bars, counts):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, 
                    str(count), ha='center', va='bottom')
        
        # Pie chart of overlap - handle zero values properly
        pairwise_72h_genes = set(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]['gene_id'])
        pairwise_135h_genes = set(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]['gene_id'])
        interaction_genes = set(self.interaction_results[self.interaction_results['significant']]['gene_id'])
        
        overlap_72h_135h = len(pairwise_72h_genes.intersection(pairwise_135h_genes))
        overlap_72h_interaction = len(pairwise_72h_genes.intersection(interaction_genes))
        overlap_135h_interaction = len(pairwise_135h_genes.intersection(interaction_genes))
        
        sizes = [overlap_72h_135h, overlap_72h_interaction, overlap_135h_interaction]
        labels = ['72h ∩ 135h', '72h ∩ Interaction', '135h ∩ Interaction']
        colors = ['lightblue', 'lightcoral', 'lightyellow']
        
        # Only create pie chart if there are non-zero values
        if sum(sizes) > 0:
            # Filter out zero values
            non_zero_sizes = [s for s in sizes if s > 0]
            non_zero_labels = [l for s, l in zip(sizes, labels) if s > 0]
            non_zero_colors = [c for s, c in zip(sizes, colors) if s > 0]
            
            if len(non_zero_sizes) > 0:
                ax2.pie(non_zero_sizes, labels=non_zero_labels, colors=non_zero_colors, 
                       autopct='%1.1f%%', startangle=90)
                ax2.set_title('Gene Overlaps Between Approaches')
            else:
                ax2.text(0.5, 0.5, 'No overlaps found', ha='center', va='center', 
                        transform=ax2.transAxes, fontsize=12)
                ax2.set_title('Gene Overlaps Between Approaches')
        else:
            ax2.text(0.5, 0.5, 'No significant genes found', ha='center', va='center', 
                    transform=ax2.transAxes, fontsize=12)
            ax2.set_title('Gene Overlaps Between Approaches')
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'summary_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Summary plots created")

    def generate_report(self):
        """Generate comprehensive comparison report."""
        logger.info("Generating comparison report...")
        
        report = f"""
# Differential Expression Analysis Comparison Report
## PRJNA158021 Embryo Dataset

**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Overview
This report compares two approaches to differential expression analysis:
1. **Pairwise comparisons** (Poelchau et al. 2013a style): DI vs NDI at each timepoint separately
2. **Interaction model**: Full model with timepoint × photoperiod interactions

## Results Summary

### Pairwise Analysis (Poelchau et al. style)
- **72h timepoint**: {len(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']])} significant genes
- **135h timepoint**: {len(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']])} significant genes
- **Total unique genes**: {len(set(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']]['gene_id']) | set(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']]['gene_id']))}

### Interaction Model
- **Significant genes**: {len(self.interaction_results[self.interaction_results['significant']])} genes
- **Model**: ~ time_point + light_cycle + time_point:light_cycle

## Key Findings

### Why Different Numbers?
1. **Statistical Power**: Interaction model requires more statistical power due to complex design
2. **Different Questions**: 
   - Pairwise: "Which genes differ between photoperiods at each timepoint?"
   - Interaction: "Which genes show different photoperiod responses at different timepoints?"
3. **Multiple Testing**: Interaction model tests for more complex effects

### Comparison with Poelchau et al. 2013a
- **Their results**: 2,506 genes at 72h, 4,337 genes at 135h
- **Our pairwise results**: {len(self.pairwise_results['72h'][self.pairwise_results['72h']['significant']])} genes at 72h, {len(self.pairwise_results['135h'][self.pairwise_results['135h']['significant']])} genes at 135h
- **Possible reasons for difference**: Different preprocessing, filtering, or statistical criteria

## Recommendations

1. **For biological interpretation**: Use pairwise results (more genes, easier to interpret)
2. **For temporal dynamics**: Use interaction model (tests for time-dependent effects)
3. **For publication**: Consider both approaches and explain the differences

## Files Generated
- `pairwise_comparisons.csv`: Results from pairwise analysis
- `interaction_analysis.csv`: Results from interaction analysis  
- `comparison_summary.csv`: Summary statistics
- `venn_diagram_comparison.png`: Venn diagram of overlaps
- `volcano_comparison.png`: Volcano plots for all approaches
- `summary_comparison.png`: Summary bar and pie charts

## Statistical Criteria
- **Significance threshold**: p < 0.05 (Benjamini-Hochberg corrected)
- **Effect size threshold**: |log2 fold change| > 0.5
- **Reference**: Poelchau et al. 2013a used same criteria
"""
        
        # Save report
        with open(self.de_dir / 'comparison_report.md', 'w') as f:
            f.write(report)
        
        logger.info("Comparison report generated")
        logger.info(f"Report saved to: {self.de_dir / 'comparison_report.md'}")

def main():
    """Main function to run the comparison analysis."""
    logger.info("Starting differential expression comparison analysis...")
    
    # Initialize analyzer
    analyzer = ComparisonAnalyzer()
    
    try:
        # Load data
        analyzer.load_data()
        
        # Run pairwise analysis (Poelchau et al. style)
        analyzer.run_pairwise_analysis()
        
        # Run interaction analysis
        analyzer.run_interaction_analysis()
        
        # Compare results
        analyzer.compare_results()
        
        # Create plots
        analyzer.create_comparison_plots()
        
        # Generate report
        analyzer.generate_report()
        
        logger.info("Comparison analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Error in comparison analysis: {e}")
        raise

if __name__ == "__main__":
    main() 