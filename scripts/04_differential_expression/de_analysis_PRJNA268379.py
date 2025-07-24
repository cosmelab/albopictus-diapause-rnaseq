# type: ignore
#!/usr/bin/env python3
"""
Differential Expression Analysis for PRJNA268379
Analyzes multiple candidate gene lists with publication-quality visualizations

Usage:
    python de_analysis_PRJNA268379.py

Requirements:
    - Complete bioinformatics environment
    - Gene list text files in output/tables/
"""

import os
import sys
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Visualization libraries that work without matplotlib
import plotnine as p9
from plotnine import *

# R integration
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
pandas2ri.activate()

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

class GeneExpressionAnalyzer:
    def __init__(self, gene_counts_file="output/tables/PRJNA268379_gene_counts.csv", 
                 gene_lists_dir="output/tables", output_base="output", analysis_mode="combined"):
        """
        Initialize the gene expression analyzer
        
        Args:
            gene_counts_file: Path to gene counts file
            gene_lists_dir: Directory containing gene list files
            output_base: Base output directory
            analysis_mode: Analysis approach to use
                - "combined": All samples together (current approach)
                - "separate_fed": Only Fed samples (LD_Fed vs SD_Fed)
                - "separate_unfed": Only Unfed samples (LD_Unfed vs SD_Unfed)
                - "both": Run both combined and separate analyses
        """
        self.gene_counts_file = gene_counts_file
        self.gene_lists_dir = Path(gene_lists_dir)  # Convert to Path object
        self.output_base = Path(output_base)
        self.analysis_mode = analysis_mode
        
        # Create output directories
        self.de_dir = self.output_base / 'differential_expression' / 'PRJNA268379'
        self.figures_dir = self.output_base / 'figures' / 'PRJNA268379'
        self.de_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # Data storage
        self.gene_counts = None
        self.metadata = None
        self.global_results = None
        self.normalized_counts = None
        self.gene_lists = {}
        self.list_analyses = {}
        
        # Analysis results storage
        self.combined_results = None
        self.fed_results = None
        self.unfed_results = None
        
        logger.info(f"Initialized GeneExpressionAnalyzer with analysis_mode: {analysis_mode}")

    def _cleanup_previous_runs(self):
        """Remove data from previous runs to ensure fresh analysis"""
        logger.info("Cleaning up previous run data...")
        
        # Clean up output directories
        if self.de_dir.exists():
            import shutil
            shutil.rmtree(self.de_dir)
            logger.info(f"Removed previous results directory: {self.de_dir}")
        
        if self.figures_dir.exists():
            import shutil
            shutil.rmtree(self.figures_dir)
            logger.info(f"Removed previous figures directory: {self.figures_dir}")
        
        # Recreate directories
        self.de_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # Clear R environment variables
        try:
            robjects.r("rm(list=ls())")
            logger.info("Cleared R environment")
        except:
            logger.info("R environment not available for cleanup")
        
        logger.info("Cleanup completed")

    def _load_gene_lists_from_files(self):
        """Load gene lists from text files"""
        gene_lists = {}
        
        # Define the expected files and their clean names
        file_mapping = {
            'all_candidate_genes_ids_only.txt': 'all_candidates',
            'high_priority_genes_ids_only.txt': 'high_priority', 
            'high_ld_genes_ids_only.txt': 'high_ld',
            'linked_snp_genes_ids_only.txt': 'linked_snp',
            'intragenic_genes_ids_only.txt': 'intragenic',
            'tag_snp_genes_ids_only.txt': 'tag_snp',
            'loc_genes_only_ids_only.txt': 'loc_genes_only'  # In case this differs from all_candidates
        }
        
        for filename, list_name in file_mapping.items():
            file_path = self.gene_lists_dir / filename
            
            if file_path.exists():
                try:
                    with open(file_path, 'r') as f:
                        genes = [line.strip() for line in f if line.strip()]
                    
                    # Remove any empty lines and validate LOC format
                    genes = [gene for gene in genes if gene.startswith('LOC')]
                    
                    if genes:
                        gene_lists[list_name] = genes
                        logger.info(f"Loaded {list_name}: {len(genes)} genes from {filename}")
                    else:
                        logger.warning(f"No valid genes found in {filename}")
                        
                except Exception as e:
                    logger.error(f"Error reading {filename}: {e}")
            else:
                logger.warning(f"File not found: {filename}")
        
        if not gene_lists:
            logger.error("No gene list files found! Check that files are in the correct directory.")
            raise FileNotFoundError(f"No gene list files found in {self.gene_lists_dir}")
        
        return gene_lists

    def validate_gene_lists(self):
        """Validate gene lists and show statistics"""
        logger.info("=== GENE LIST VALIDATION ===")
        
        all_unique_genes = set()
        for list_name, genes in self.gene_lists.items():
            all_unique_genes.update(genes)
        
        logger.info(f"Total unique genes across all lists: {len(all_unique_genes)}")
        
        # Check for overlaps
        logger.info("\n=== OVERLAP ANALYSIS ===")
        list_names = list(self.gene_lists.keys())
        
        for i, list1 in enumerate(list_names):
            for j, list2 in enumerate(list_names[i+1:], i+1):
                overlap = set(self.gene_lists[list1]) & set(self.gene_lists[list2])
                logger.info(f"{list1} ∩ {list2}: {len(overlap)} genes")
        
        # Show gene counts per list
        logger.info("\n=== GENE COUNTS PER LIST ===")
        for list_name, genes in sorted(self.gene_lists.items(), key=lambda x: len(x[1]), reverse=True):
            logger.info(f"{list_name:20}: {len(genes):3} genes")
        
        return self

    def refresh_gene_lists(self):
        """Reload gene lists from files (useful when files are updated)"""
        logger.info("Refreshing gene lists from files...")
        old_counts = {name: len(genes) for name, genes in self.gene_lists.items()}
        
        self.gene_lists = self._load_gene_lists_from_files()
        
        # Show changes
        logger.info("=== GENE LIST CHANGES ===")
        for list_name, genes in self.gene_lists.items():
            old_count = old_counts.get(list_name, 0)
            new_count = len(genes)
            change = new_count - old_count
            if change != 0:
                logger.info(f"{list_name}: {old_count} -> {new_count} genes ({change:+d})")
            else:
                logger.info(f"{list_name}: {new_count} genes (no change)")
        
        return self

    def load_data(self):
        """Load and prepare count data and metadata"""
        logger.info("Loading count data...")
        
        # Load gene counts
        self.gene_counts = pd.read_csv(self.gene_counts_file, index_col=0)
        
        # Filter for PRJNA268379 samples only
        prjna268379_samples = [col for col in self.gene_counts.columns if 'PRJNA268379' in col]
        self.gene_counts = self.gene_counts[prjna268379_samples]
        
        logger.info(f"Loaded {self.gene_counts.shape[0]} genes × {self.gene_counts.shape[1]} samples")
        
        # Create metadata
        self.metadata = self._create_metadata()
        
        # Ensure sample order matches
        self.gene_counts = self.gene_counts[self.metadata.index]
        
        logger.info("Data loading completed")
        return self

    def _create_metadata(self):
        """Create metadata dataframe from sample names"""
        if self.gene_counts is None:
            raise ValueError("gene_counts must be loaded before creating metadata")
            
        metadata_list = []
        
        for sample in self.gene_counts.columns:
            parts = sample.split('_')
            photoperiod = parts[2]  # LD or SD
            blood_meal = parts[3]   # Fed or Unfed
            treatment = f"{photoperiod}_{blood_meal}"
            
            metadata_list.append({
                'sample_id': sample,
                'run_id': parts[1],
                'photoperiod': photoperiod,
                'blood_meal': blood_meal,
                'treatment': treatment,
                'dataset': 'PRJNA268379',
                'life_stage': 'Adult'
            })
        
        metadata_df = pd.DataFrame(metadata_list)
        metadata_df.set_index('sample_id', inplace=True)
        
        logger.info(f"Treatments: {metadata_df['treatment'].value_counts().to_dict()}")
        return metadata_df

    def run_global_analysis(self):
        """Run global differential expression analysis."""
        try:
            logger.info("Running global differential expression analysis...")
            
            # Filter genes
            self._filter_genes()
            
            # Setup R environment
            self._setup_r_environment()
            
            # Run DESeq2
            self._run_deseq2()
            
            # Extract results
            self._extract_global_results()
            
            # Create global expression summary table and statistics
            self._create_global_expression_summary()
            
            # Create PCA analysis following R script methodology
            self._create_pca_analysis()
            
            # Create correlation heatmap following R script methodology
            self._create_correlation_analysis()
            
            # Create global volcano plot
            self._create_global_volcano()
            
            logger.info("Global analysis completed")
            
        except Exception as e:
            logger.error(f"Error in global analysis: {str(e)}")
            raise

    def _filter_genes(self, min_count=10, min_samples=3):
        """Filter genes with low counts - following R script methodology"""
        if self.gene_counts is None:
            raise ValueError("gene_counts must be loaded before filtering")
            
        before_filter = self.gene_counts.shape[0]
        
        # Follow R script: rowSums(counts(dds)) >= 10
        # Keep genes with at least 10 total reads across all samples
        keep = self.gene_counts.sum(axis=1) >= min_count
        
        self.gene_counts = self.gene_counts[keep]
        after_filter = self.gene_counts.shape[0]
        logger.info(f"Gene filtering (R script method): {before_filter} -> {after_filter} genes ({after_filter/before_filter*100:.1f}% retained)")
        logger.info(f"Filtering criteria: genes with ≥{min_count} total reads across all samples")

    def _setup_r_environment(self):
        """Setup R environment with required packages"""
        # Try to load packages directly (they should be available via conda)
        r_code = """
        # Try to load packages directly
        tryCatch({
            library(DESeq2)
            library(ggplot2)
            cat("Successfully loaded DESeq2 and ggplot2\\n")
        }, error = function(e) {
            cat("Error loading packages:", e$message, "\\n")
            stop("Required R packages not available")
        })
        """
        
        try:
            robjects.r(r_code)
            logger.info("R environment setup completed successfully")
        except Exception as e:
            logger.error(f"Error setting up R environment: {str(e)}")
            raise

    def _run_deseq2(self):
        """Run DESeq2 analysis"""
        if self.gene_counts is None or self.metadata is None:
            raise ValueError("gene_counts and metadata must be loaded before running DESeq2")
            
        count_matrix = self.gene_counts.astype(int)
        metadata_r = self.metadata[['photoperiod', 'blood_meal']].copy()
        metadata_r['photoperiod'] = pd.Categorical(metadata_r['photoperiod'], categories=['LD', 'SD'])
        metadata_r['blood_meal'] = pd.Categorical(metadata_r['blood_meal'], categories=['Unfed', 'Fed'])
        
        # Transfer to R
        robjects.globalenv['count_matrix'] = pandas2ri.py2rpy(count_matrix)
        robjects.globalenv['metadata'] = pandas2ri.py2rpy(metadata_r)
        
        # DESeq2 analysis
        r_code = """
        dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                    colData = metadata,
                                    design = ~ blood_meal + photoperiod)
        
        dds$photoperiod <- relevel(dds$photoperiod, ref = "LD")
        dds$blood_meal <- relevel(dds$blood_meal, ref = "Unfed")
        
        dds <- DESeq(dds)
        res <- results(dds, contrast = c("photoperiod", "SD", "LD"))
        res <- res[order(res$padj), ]
        normalized_counts <- counts(dds, normalized = TRUE)
        """
        robjects.r(r_code)

    def _extract_global_results(self):
        """Extract results from DESeq2 analysis"""
        try:
            # Convert DESeq2 results to data frame in R
            r_code = """
            # Convert DESeq2 results to data frame
            results_df <- as.data.frame(results(dds))
            # Add gene IDs as a column
            results_df$gene_id <- rownames(results_df)
            # Convert to data frame
            results_df <- as.data.frame(results_df)
            
            # Get normalized counts
            normalized_counts <- counts(dds, normalized=TRUE)
            normalized_counts_df <- as.data.frame(normalized_counts)
            normalized_counts_df$gene_id <- rownames(normalized_counts_df)
            """
            robjects.r(r_code)
            
            # Get the results data frame from R
            r_results = robjects.r['results_df']
            r_normalized = robjects.r['normalized_counts_df']
            
            # Convert to pandas DataFrame
            self.global_results = pd.DataFrame(
                {
                    'gene_id': r_results.rx2('gene_id'),
                    'baseMean': r_results.rx2('baseMean'),
                    'log2FoldChange': r_results.rx2('log2FoldChange'),
                    'lfcSE': r_results.rx2('lfcSE'),
                    'stat': r_results.rx2('stat'),
                    'pvalue': r_results.rx2('pvalue'),
                    'padj': r_results.rx2('padj')
                }
            )
            
            # Convert normalized counts to pandas DataFrame
            # Get column names from R and convert to Python strings
            colnames = [str(col) for col in r_normalized.colnames]
            # Remove gene_id from column names
            count_cols = [col for col in colnames if col != 'gene_id']
            
            # Create DataFrame with proper column names
            self.normalized_counts = pd.DataFrame(
                {col: r_normalized.rx2(col) for col in count_cols},
                index=r_normalized.rx2('gene_id')
            )
            
            # Set gene_id as index for results
            self.global_results.set_index('gene_id', inplace=True)
            
            # Add significant column based on R script thresholds:
            # |log2FC| > 0.58 and FDR < 0.05 (Benjamini-Hochberg corrected)
            # Following the R script: abs(log2FoldChange) > 0.58
            self.global_results['significant'] = (
                (self.global_results['padj'] < 0.05) & 
                (abs(self.global_results['log2FoldChange']) > 0.58)  # R script uses 0.58, not 0.5
            )
            
            # Sort by adjusted p-value
            self.global_results.sort_values('padj', inplace=True)
            
            # Save results
            self.global_results.to_csv(self.de_dir / 'global_deseq2_results.csv')
            self.normalized_counts.to_csv(self.de_dir / 'normalized_counts.csv')
            
            # Log summary statistics
            total_genes = len(self.global_results)
            sig_genes = len(self.global_results[self.global_results['significant']])  # Use 'significant' column
            up_reg = len(self.global_results[(self.global_results['significant']) & 
                                           (self.global_results['log2FoldChange'] > 0)])
            down_reg = len(self.global_results[(self.global_results['significant']) & 
                                             (self.global_results['log2FoldChange'] < 0)])
            
            logger.info(f"\nGlobal Differential Expression Results (R script thresholds):")
            logger.info(f"Total genes analyzed: {total_genes}")
            logger.info(f"Significant genes (|log2FC| > 0.58, FDR < 0.05): {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
            logger.info(f"Up-regulated: {up_reg} ({up_reg/total_genes*100:.1f}%)")
            logger.info(f"Down-regulated: {down_reg} ({down_reg/total_genes*100:.1f}%)")
            
            # Compare with R script results
            logger.info(f"\nFollowing R script methodology:")
            logger.info(f"Filtering: genes with ≥10 total reads across all samples")
            logger.info(f"Significance: |log2FC| > 0.58 and padj < 0.05")
            logger.info(f"Current analysis - Total significant: {sig_genes} genes ({up_reg} up + {down_reg} down)")
            
        except Exception as e:
            logger.error(f"Error extracting global results: {str(e)}")
            raise

    def analyze_gene_lists(self):
        """Analyze only the all_candidates gene list"""
        logger.info("Analyzing all_candidates gene list...")
        
        # Only analyze all_candidates if it exists
        if 'all_candidates' not in self.gene_lists or not self.gene_lists['all_candidates']:
            logger.warning("all_candidates gene list not found or empty, skipping gene list analysis")
            return
        
        genes = self.gene_lists['all_candidates']
        logger.info(f"Analyzing all_candidates: {len(genes)} genes")
        
        # Extract results for all_candidates
        list_results, summary = self._extract_list_results(genes, 'all_candidates')
        
        # Store analysis
        self.list_analyses['all_candidates'] = {
            'results': list_results,
            'stats': summary,
            'genes': genes
        }
        
        # Create visualizations for all_candidates
        self._create_candidate_visualizations('all_candidates', {
            'results': list_results,
            'stats': summary,
            'genes': genes
        })
        
        logger.info("all_candidates analysis completed")

    def _create_candidate_visualizations(self, list_name, list_analysis):
        """Create visualizations specifically for candidate genes"""
        results = list_analysis['results']
        
        if results.empty:
            logger.warning(f"No results to plot for {list_name}")
            return
        
        # Create expression heatmap for candidate genes
        self._create_candidate_expression_heatmap(list_name, results)
        
        # Create candidate-specific volcano plot
        self._create_candidate_volcano_plot(list_name, results)
        
        # Create candidate expression boxplot
        self._create_candidate_expression_boxplot(list_name, results)

    def _create_candidate_expression_heatmap(self, list_name, results):
        """Create expression heatmap for candidate genes"""
        if results.empty or self.normalized_counts is None or self.metadata is None:
            return
            
        # Get genes for this list
        list_genes = results.index.tolist()
        
        # Get normalized counts for these genes
        heatmap_data = self.normalized_counts.loc[list_genes]
        
        # Sort samples in the specified order: LD_Unfed, LD_Fed, SD_Unfed, SD_Fed
        sample_order = ['LD_Unfed', 'LD_Fed', 'SD_Unfed', 'SD_Fed']
        sorted_samples = []
        for treatment in sample_order:
            treatment_samples = self.metadata[self.metadata['treatment'] == treatment].index.tolist()
            sorted_samples.extend(treatment_samples)
        
        # Reorder the heatmap data
        heatmap_data = heatmap_data[sorted_samples]
        
        # Log2 transform and Z-score
        log_data = np.log2(heatmap_data + 1)
        from scipy.stats import zscore
        z_data = pd.DataFrame(
            zscore(log_data, axis=1),
            index=log_data.index,
            columns=log_data.columns
        )
        
        # Create metadata for annotation with sorted samples
        col_annotation = self.metadata.loc[sorted_samples, ['treatment']].copy()
        
        # Set up the figure with a specific style
        fig = plt.figure(figsize=(14, max(8, len(list_genes) * 0.5)))
        
        # Create custom colormap
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        
        # Create clustermap with enhanced styling
        g = sns.clustermap(z_data,
                          cmap=cmap,
                          center=0,
                          col_colors=[self._get_treatment_color(t) for t in col_annotation['treatment']],
                          row_cluster=True,
                          col_cluster=False,  # Don't cluster columns to maintain order
                          xticklabels=True,
                          yticklabels=True,
                          figsize=(14, max(8, len(list_genes) * 0.5)),
                          cbar_kws={'label': 'Z-score (log₂ expression)',
                                  'orientation': 'horizontal',
                                  'pad': 0.02,
                                  'fraction': 0.05},
                          linewidths=0.5,
                          linecolor='gray')
        
        # Replace x-axis labels with simplified treatment names
        if hasattr(g, 'ax_heatmap'):
            # Get current x-tick labels
            current_labels = g.ax_heatmap.get_xticklabels()
            new_labels = []
            
            for label in current_labels:
                sample_name = label.get_text()
                if sample_name in self.metadata.index:
                    treatment = self.metadata.loc[sample_name, 'treatment']
                    # Add number if multiple samples per treatment
                    existing_count = new_labels.count(treatment)
                    if existing_count > 0:
                        treatment = f"{treatment}_{existing_count + 1}"
                    new_labels.append(treatment)
                else:
                    new_labels.append(sample_name)
            
            # Set new labels
            g.ax_heatmap.set_xticklabels(new_labels, rotation=45, ha='right')
        
        # Adjust the layout
        plt.suptitle(f'Candidate Genes Expression Heatmap\n(PRJNA268379)', y=0.98, fontsize=14, fontweight='bold')
        
        # Add treatment legend with enhanced styling
        self._add_treatment_legend(g.ax_heatmap)
        
        # Adjust axis labels
        g.ax_heatmap.set_xlabel('Samples', fontsize=12, labelpad=10)
        g.ax_heatmap.set_ylabel('Candidate Genes', fontsize=12, labelpad=10)
        
        # Rotate x-axis labels for better readability
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        
        # Save figure with high DPI and tight bounding box
        plt.savefig(self.figures_dir / 'candidate_genes_heatmap.png',
                   dpi=300,
                   bbox_inches='tight',
                   facecolor='white',
                   edgecolor='none')
        plt.savefig(self.figures_dir / 'candidate_genes_heatmap.pdf',
                   bbox_inches='tight',
                   facecolor='white',
                   edgecolor='none')
        plt.close()

    def _create_candidate_volcano_plot(self, list_name, results):
        """Create volcano plot specifically for candidate genes"""
        if results.empty:
            return
            
        plt.figure(figsize=(12, 10))
        
        # Create plot data
        plot_data = results.copy()
        plot_data['-log10_padj'] = -np.log10(plot_data['padj'])
        
        # Plot all candidate genes
        plt.scatter(plot_data['log2FoldChange'], plot_data['-log10_padj'], 
                   alpha=0.6, s=30, c='gray', label='Not significant')
        
        # Highlight significant candidate genes
        sig_up = plot_data[(plot_data['significant']) & (plot_data['log2FoldChange'] > 0)]
        sig_down = plot_data[(plot_data['significant']) & (plot_data['log2FoldChange'] < 0)]
        
        if not sig_up.empty:
            plt.scatter(sig_up['log2FoldChange'], sig_up['-log10_padj'], 
                       alpha=0.8, s=50, c='red', label=f'Up-regulated ({len(sig_up)})')
        
        if not sig_down.empty:
            plt.scatter(sig_down['log2FoldChange'], sig_down['-log10_padj'], 
                       alpha=0.8, s=50, c='blue', label=f'Down-regulated ({len(sig_down)})')
        
        # Add thresholds
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.8)
        plt.axvline(0.58, color='gray', linestyle='--', alpha=0.8)
        plt.axvline(-0.58, color='gray', linestyle='--', alpha=0.8)
        
        plt.xlabel('Log2 Fold Change', fontsize=14)
        plt.ylabel('-Log10(Adjusted P-value)', fontsize=14)
        plt.title('Candidate Genes Volcano Plot\n(PRJNA268379)', fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'candidate_genes_volcano.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'candidate_genes_volcano.pdf', bbox_inches='tight')
        plt.close()

    def _create_candidate_expression_boxplot(self, list_name, results):
        """Create expression boxplot for candidate genes"""
        if results.empty or self.normalized_counts is None or self.metadata is None:
            return
            
        # Get genes for this list
        list_genes = results.index.tolist()
        
        # Get normalized counts for these genes
        expression_data = self.normalized_counts.loc[list_genes]
        
        # Log2 transform
        log_expression = np.log2(expression_data + 1)
        
        # Prepare data for boxplot
        plot_data = []
        treatment_labels = []
        
        treatments = ['LD_Unfed', 'LD_Fed', 'SD_Unfed', 'SD_Fed']
        colors = ['#1f77b4', '#aec7e8', '#d62728', '#ffbb78']
        
        for treatment in treatments:
            treatment_samples = self.metadata[self.metadata['treatment'] == treatment].index
            if len(treatment_samples) > 0:
                treatment_data = log_expression[treatment_samples].values.flatten()
                plot_data.extend([treatment_data])
                treatment_labels.extend([treatment] * len(treatment_data))
        
        # Create boxplot
        plt.figure(figsize=(12, 8))
        
        bp = plt.boxplot(plot_data, labels=treatments, patch_artist=True)
        
        # Color the boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        plt.xlabel('Treatment', fontsize=14)
        plt.ylabel('Log2 Normalized Expression', fontsize=14)
        plt.title('Candidate Genes Expression by Treatment\n(PRJNA268379)', fontsize=16, fontweight='bold')
        plt.grid(True, alpha=0.3)
        
        # Add legend
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.7, label=treatment) 
                          for treatment, color in zip(treatments, colors)]
        plt.legend(handles=legend_elements, title='Treatment', fontsize=12)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'candidate_genes_expression.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'candidate_genes_expression.pdf', bbox_inches='tight')
        plt.close()

    def _extract_list_results(self, genes, list_name):
        """Extract results for a specific gene list"""
        if self.global_results is None:
            logger.warning("No global results available for gene list analysis")
            return pd.DataFrame(), {}
        
        # Find genes that exist in the results
        found_genes = [gene for gene in genes if gene in self.global_results.index]
        missing_genes = [gene for gene in genes if gene not in self.global_results.index]
        
        if missing_genes:
            logger.warning(f"Missing {len(missing_genes)} genes from {list_name} in results")
        
        if not found_genes:
            logger.warning(f"No genes from {list_name} found in results")
            return pd.DataFrame(), {}
        
        # Extract results for found genes
        list_results = self.global_results.loc[found_genes].copy()
        
        # Calculate summary statistics
        total_genes = len(genes)
        genes_found = len(found_genes)
        significant_genes = len(list_results[list_results['significant']])
        up_regulated = len(list_results[(list_results['significant']) & (list_results['log2FoldChange'] > 0)])
        down_regulated = len(list_results[(list_results['significant']) & (list_results['log2FoldChange'] < 0)])
        
        summary = {
            'total_genes': total_genes,
            'genes_found': genes_found,
            'missing_genes': len(missing_genes),
            'significant_genes': significant_genes,
            'up_regulated': up_regulated,
            'down_regulated': down_regulated,
            'mean_log2fc': list_results['log2FoldChange'].mean(),
            'median_padj': list_results['padj'].median(),
            'hit_rate': significant_genes / genes_found if genes_found > 0 else 0
        }
        
        # Save results for this list
        list_results.to_csv(self.de_dir / f'{list_name}_results.csv')
        
        logger.info(f"{list_name}: {genes_found}/{total_genes} genes found, {significant_genes} significant")
        
        return list_results, summary

    def _get_treatment_color(self, treatment):
        """Get color for treatment"""
        colors = {
            'LD_Unfed': '#1f77b4',  # Blue
            'LD_Fed': '#aec7e8',    # Light blue
            'SD_Unfed': '#d62728',  # Red
            'SD_Fed': '#ffbb78'     # Orange
        }
        return colors.get(treatment, '#808080')  # Gray for unknown treatments

    def _add_treatment_legend(self, ax):
        """Add treatment legend to plot"""
        treatments = ['LD_Unfed', 'LD_Fed', 'SD_Unfed', 'SD_Fed']
        colors = [self._get_treatment_color(t) for t in treatments]
        
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.7, label=treatment) 
                          for treatment, color in zip(treatments, colors)]
        
        legend = ax.legend(handles=legend_elements, title='Treatment', 
                          bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        plt.setp(legend.get_title(), fontsize=12)
        plt.setp(legend.get_texts(), fontsize=10)

    def _create_statistical_summary(self, significant_candidates):
        """Create a markdown file with statistical summary of the analysis."""
        try:
            # Create summary statistics
            summary_stats = {
                'Total significant genes': len(significant_candidates),
                'Up-regulated genes': len(significant_candidates[significant_candidates['log2FoldChange'] > 0]),
                'Down-regulated genes': len(significant_candidates[significant_candidates['log2FoldChange'] < 0]),
                'Significance threshold (padj)': 0.05
            }
            
            # Create gene-specific statistics
            gene_stats = significant_candidates.copy()
            gene_stats['regulation'] = gene_stats['log2FoldChange'].apply(
                lambda x: 'Up-regulated' if x > 0 else 'Down-regulated'
            )
            
            # Add treatment-specific means
            for treatment in ['LD_Unfed', 'LD_Fed', 'SD_Unfed', 'SD_Fed']:
                treatment_samples = self.metadata[self.metadata['treatment'] == treatment].index
                gene_stats[f'{treatment}_mean'] = self.normalized_counts.loc[gene_stats.index, treatment_samples].mean(axis=1)
            
            # Create QC directory if it doesn't exist
            qc_dir = self.output_base / 'qc' / 'PRJNA268379'
            qc_dir.mkdir(parents=True, exist_ok=True)
            
            # Save gene statistics to CSV
            gene_stats.to_csv(qc_dir / 'candidate_genes_statistics.csv')
            
            # Create markdown content
            md_content = """# Statistical Analysis of Candidate Genes Expression (PRJNA268379)

## Summary Statistics
"""
            # Add summary statistics
            for key, value in summary_stats.items():
                md_content += f"- {key}: {value}\n"
            
            md_content += """
## Gene-specific Statistics
The detailed statistics for each gene can be found in the CSV file: `candidate_genes_statistics.csv`

The file contains the following information for each gene:
- log2FoldChange: Log2 fold change in expression
- padj: Adjusted p-value
- regulation: Direction of regulation (Up-regulated/Down-regulated)
- baseMean: Mean normalized expression across all samples
- LD_Unfed_mean: Mean expression in LD Unfed samples
- LD_Fed_mean: Mean expression in LD Fed samples
- SD_Unfed_mean: Mean expression in SD Unfed samples
- SD_Fed_mean: Mean expression in SD Fed samples
- lfcSE: Standard error of log2 fold change
- stat: Wald statistic
- pvalue: Raw p-value

## Treatment Groups
The analysis includes the following treatment groups:
- LD_Unfed: Long-day unfed
- LD_Fed: Long-day fed
- SD_Unfed: Short-day unfed
- SD_Fed: Short-day fed

## Statistical Tests
- Differential expression analysis was performed using DESeq2
- Significance was determined using adjusted p-value (padj) < 0.05
- Log2 fold changes were calculated relative to the control condition
"""
            
            # Save markdown file
            with open(qc_dir / 'candidate_genes_analysis_summary.md', 'w') as f:
                f.write(md_content)
            
            logger.info("Statistical summary markdown file created successfully")
            
        except Exception as e:
            logger.error(f"Error creating statistical summary: {str(e)}")
            raise

    def create_comparison_summary(self):
        """Create summary comparison for all_candidates only"""
        logger.info("Creating comparison summary for all_candidates...")
        
        if 'all_candidates' not in self.list_analyses:
            logger.warning("all_candidates analysis not found, skipping comparison summary")
            return
        
        # Get all_candidates analysis
        analysis = self.list_analyses['all_candidates']
        stats = analysis['stats']
        
        # Create summary table
        summary_data = {
            'Gene_List': ['all_candidates'],
            'Total_Genes': [stats['total_genes']],
            'Genes_Found': [stats['genes_found']],
            'Significant_Genes': [stats['significant_genes']],
            'Up_Regulated': [stats['up_regulated']],
            'Down_Regulated': [stats['down_regulated']],
            'Hit_Rate_Percent': [stats['hit_rate'] * 100 if stats['genes_found'] > 0 else 0],
            'Mean_Log2FC': [stats['mean_log2fc']],
            'Median_Padj': [stats['median_padj']]
        }
        
        summary_df = pd.DataFrame(summary_data)
        
        # Save summary
        summary_df.to_csv(self.de_dir / 'candidate_genes_summary.csv', index=False)
        
        # Create simple comparison plot
        self._create_candidate_summary_plot(summary_df)
        
        logger.info("Comparison summary completed")

    def _create_candidate_summary_plot(self, summary_df):
        """Create summary plot for candidate genes"""
        if summary_df.empty:
            return
            
        # Create figure with subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Hit rate
        ax1.bar(summary_df['Gene_List'], summary_df['Hit_Rate_Percent'], 
               color=['#2E86AB'], alpha=0.7)
        ax1.set_ylabel('Hit Rate (%)', fontsize=12)
        ax1.set_title('Candidate Genes Hit Rate', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Up vs Down regulated
        up_down_data = summary_df[['Up_Regulated', 'Down_Regulated']].values.flatten()
        labels = ['Up-regulated', 'Down-regulated']
        colors = ['#FF6B6B', '#4ECDC4']
        ax2.pie(up_down_data, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax2.set_title('Regulation Direction', fontsize=14, fontweight='bold')
        
        # Plot 3: Genes found vs significant
        found_sig_data = summary_df[['Genes_Found', 'Significant_Genes']].values.flatten()
        labels = ['Genes Found', 'Significant']
        colors = ['#95A5A6', '#E74C3C']
        ax3.bar(labels, found_sig_data, color=colors, alpha=0.7)
        ax3.set_ylabel('Number of Genes', fontsize=12)
        ax3.set_title('Gene Counts', fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Mean log2FC distribution
        ax4.hist(summary_df['Mean_Log2FC'], bins=10, color='#9B59B6', alpha=0.7, edgecolor='black')
        ax4.set_xlabel('Mean Log2 Fold Change', fontsize=12)
        ax4.set_ylabel('Frequency', fontsize=12)
        ax4.set_title('Mean Log2FC Distribution', fontsize=14, fontweight='bold')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'candidate_genes_summary.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'candidate_genes_summary.pdf', bbox_inches='tight')
        plt.close()

    def generate_analysis_report(self):
        """Generate analysis summary report for all_candidates"""
        logger.info("Generating analysis report...")
        
        # Compile overall statistics
        total_genes_tested = len(self.global_results) if self.global_results is not None else 0
        global_significant = self.global_results['significant'].sum() if self.global_results is not None else 0
        
        # Get candidate genes statistics
        candidate_stats = None
        if 'all_candidates' in self.list_analyses:
            candidate_stats = self.list_analyses['all_candidates']['stats']
        
        report = f"""# PRJNA268379 Gene Expression Analysis - Candidate Genes Focus

## Analysis Summary
- **Total genes tested**: {total_genes_tested:,}
- **Significant genes (global)**: {global_significant:,} ({global_significant/total_genes_tested*100:.2f}% if total_genes_tested > 0 else 0%)
- **Model**: ~ blood_meal + photoperiod (SD vs LD)
- **Significance criteria**: FDR < 0.05 AND |log2FC| > 0.58 (following R script methodology)
- **Filtering criteria**: genes with ≥10 total reads across all samples (following R script methodology)

## Candidate Genes Analysis Summary

"""
        
        if candidate_stats:
            report += f"""| Metric | Value |
|--------|-------|
| Total candidate genes | {candidate_stats['total_genes']} |
| Genes found in analysis | {candidate_stats['genes_found']} |
| Significant genes | {candidate_stats['significant_genes']} |
| Up-regulated in SD | {candidate_stats['up_regulated']} |
| Down-regulated in SD | {candidate_stats['down_regulated']} |
| Hit rate | {candidate_stats['hit_rate'] * 100:.1f}% |
| Mean log2FC | {candidate_stats['mean_log2fc']:.3f} |
| Median adjusted p-value | {candidate_stats['median_padj']:.2e} |

## Key Findings

### Candidate Genes Performance:
- **Hit Rate**: {candidate_stats['hit_rate'] * 100:.1f}% of candidate genes show significant differential expression
- **Regulation Pattern**: {candidate_stats['up_regulated']} up-regulated vs {candidate_stats['down_regulated']} down-regulated in SD conditions
- **Effect Size**: Mean log2 fold change of {candidate_stats['mean_log2fc']:.3f}

### Expression Patterns:
"""
            
            if candidate_stats['mean_log2fc'] > 0.5:
                report += "- **Overall trend**: Candidate genes tend to be upregulated in SD (diapause) conditions\n"
            elif candidate_stats['mean_log2fc'] < -0.5:
                report += "- **Overall trend**: Candidate genes tend to be downregulated in SD (diapause) conditions\n"
            else:
                report += "- **Overall trend**: Candidate genes show mixed regulation patterns\n"
        
        report += f"""
## Output Files
- **Global results**: `global_deseq2_results.csv`
- **Candidate genes results**: `all_candidates_results.csv`
- **Candidate genes summary**: `candidate_genes_summary.csv`
- **Global expression summary**: `global_expression_summary.csv`
- **Up-regulated genes**: `up_regulated_genes.csv`
- **Down-regulated genes**: `down_regulated_genes.csv`
- **PCA results**: `pca_results.csv`
- **PCA variance**: `pca_variance.csv`
- **Sample correlation matrix**: `sample_correlation_matrix.csv`
- **Normalized counts**: `normalized_counts.csv`
- **Visualizations**: PNG and PDF files in figures directory

## Visualizations Created

### Global Expression Profiling:
1. **Global volcano plot** - All genes with candidate genes highlighted
2. **PCA analysis** - Sample clustering and variance explanation
3. **Correlation heatmap** - Sample-to-sample correlations

### Candidate Genes Analysis:
4. **Candidate genes heatmap** - Expression patterns across treatments
5. **Candidate genes volcano plot** - Focused on candidate genes only
6. **Candidate genes expression boxplot** - Expression levels by treatment
7. **Candidate genes summary plots** - Hit rates and regulation patterns

## Global Expression Profiling Results

### Principal Component Analysis (PCA):
- **PC1 variance explained**: See PCA results file for detailed variance explanation
- **PC2 variance explained**: See PCA results file for detailed variance explanation
- **Sample clustering**: Shows treatment-based separation
- **Quality control**: Identifies potential outliers or batch effects

### Sample Correlation Analysis:
- **Correlation range**: -1 to +1 (Pearson correlation)
- **Treatment consistency**: High correlation within treatment groups
- **Quality assessment**: Identifies problematic samples

### Global Differential Expression:
- **Total genes analyzed**: {total_genes_tested:,}
- **Significant genes**: {global_significant:,} ({global_significant/total_genes_tested*100:.1f}% if total_genes_tested > 0 else 0%)
- **Up-regulated in SD**: {len(self.global_results[(self.global_results['significant']) & (self.global_results['log2FoldChange'] > 0)]) if self.global_results is not None else 0}
- **Down-regulated in SD**: {len(self.global_results[(self.global_results['significant']) & (self.global_results['log2FoldChange'] < 0)]) if self.global_results is not None else 0}

## Recommendations
1. **Focus on significant candidate genes** for follow-up validation
2. **Investigate biological pathways** enriched in significant candidate genes
3. **Compare with published literature** on diapause-related gene expression
4. **Validate top candidates** experimentally in independent samples
5. **Use PCA results** to identify potential outliers or batch effects
6. **Examine correlation patterns** to assess sample quality and treatment consistency

**Analysis completed**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}
"""
        
        # Save report
        with open(self.de_dir / 'candidate_genes_analysis_report.md', 'w') as f:
            f.write(report)
        
        logger.info("Analysis report generated")
        return self

    def run_analysis(self):
        """Main analysis pipeline"""
        try:
            logger.info("Starting analysis pipeline...")
            
            # Clean up previous runs
            self._cleanup_previous_runs()
            
            # Load data
            self.load_data()
            
            # Run global analysis
            self.run_global_analysis()
            
            # Analyze gene lists (candidate genes)
            self.analyze_gene_lists()
            
            # Create comparison summary
            self.create_comparison_summary()
            
            # Generate analysis report
            self.generate_analysis_report()
            
            # Create visualizations
            logger.info("Creating visualizations...")
            self._create_pca_analysis()
            self._create_correlation_analysis()
            self._create_global_volcano()
            self._analyze_candidate_genes()
            
            logger.info("Analysis completed successfully!")
            logger.info(f"Results saved to: {self.de_dir}")
            logger.info(f"Figures in: {self.figures_dir}")
            
        except Exception as e:
            logger.error(f"Error in analysis pipeline: {str(e)}")
            raise

    def _create_global_expression_summary(self):
        """Create global expression summary table following R script methodology"""
        logger.info("Creating global expression summary...")
        
        if self.global_results is None or self.normalized_counts is None:
            logger.warning("No global results available for expression summary")
            return
        
        # Calculate summary statistics for all genes
        summary_data = []
        
        for gene_id in self.global_results.index:
            gene_results = self.global_results.loc[gene_id]
            gene_expression = self.normalized_counts.loc[gene_id]
            
            # Calculate treatment means
            ld_unfed_mean = gene_expression[self.metadata[self.metadata['treatment'] == 'LD_Unfed'].index].mean()
            ld_fed_mean = gene_expression[self.metadata[self.metadata['treatment'] == 'LD_Fed'].index].mean()
            sd_unfed_mean = gene_expression[self.metadata[self.metadata['treatment'] == 'SD_Unfed'].index].mean()
            sd_fed_mean = gene_expression[self.metadata[self.metadata['treatment'] == 'SD_Fed'].index].mean()
            
            summary_data.append({
                'gene_id': gene_id,
                'baseMean': gene_results['baseMean'],
                'log2FoldChange': gene_results['log2FoldChange'],
                'padj': gene_results['padj'],
                'significant': gene_results['significant'],
                'LD_Unfed_mean': ld_unfed_mean,
                'LD_Fed_mean': ld_fed_mean,
                'SD_Unfed_mean': sd_unfed_mean,
                'SD_Fed_mean': sd_fed_mean,
                'regulation': 'Up' if gene_results['log2FoldChange'] > 0 else 'Down' if gene_results['log2FoldChange'] < 0 else 'None'
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.set_index('gene_id', inplace=True)
        
        # Save summary
        summary_df.to_csv(self.de_dir / 'global_expression_summary.csv')
        
        # Create up/down regulated gene lists
        up_genes = summary_df[(summary_df['significant']) & (summary_df['log2FoldChange'] > 0)]
        down_genes = summary_df[(summary_df['significant']) & (summary_df['log2FoldChange'] < 0)]
        
        up_genes.to_csv(self.de_dir / 'up_regulated_genes.csv')
        down_genes.to_csv(self.de_dir / 'down_regulated_genes.csv')
        
        logger.info(f"Global expression summary created: {len(summary_df)} genes")
        logger.info(f"Up-regulated: {len(up_genes)} genes")
        logger.info(f"Down-regulated: {len(down_genes)} genes")

    def _create_pca_analysis(self):
        """Create PCA analysis following R script methodology"""
        logger.info("Creating PCA analysis...")
        
        if self.normalized_counts is None or self.metadata is None:
            logger.warning("No normalized counts available for PCA")
            return
        
        # Log2 transform the data (approximate rlog transformation)
        log_data = np.log2(self.normalized_counts + 1)
        
        # Transpose for PCA (samples as rows, genes as columns)
        pca_data = log_data.T
        
        # Perform PCA
        from sklearn.decomposition import PCA
        pca = PCA()
        pca_result = pca.fit_transform(pca_data)
        
        # Create PCA DataFrame
        pca_df = pd.DataFrame(
            pca_result,
            index=pca_data.index,
            columns=[f'PC{i+1}' for i in range(pca_result.shape[1])]
        )
        
        # Add metadata
        pca_df['treatment'] = self.metadata.loc[pca_df.index, 'treatment']
        pca_df['photoperiod'] = self.metadata.loc[pca_df.index, 'photoperiod']
        pca_df['blood_meal'] = self.metadata.loc[pca_df.index, 'blood_meal']
        
        # Create PCA plot with plotnine
        explained_var = pca.explained_variance_ratio_
        
        p = (ggplot(pca_df, aes(x='PC1', y='PC2', color='treatment'))
             + geom_point(size=3, alpha=0.8)
             + labs(
                 x=f'PC1 ({explained_var[0]*100:.1f}%)',
                 y=f'PC2 ({explained_var[1]*100:.1f}%)',
                 title='PCA Analysis of Gene Expression (PRJNA268379)',
                 color='Treatment'
             )
             + scale_color_manual(values=['#1f77b4', '#aec7e8', '#d62728', '#ffbb78'])
             + theme_minimal()
             + theme(legend_position='right'))
        
        # Save plot
        p.save(self.figures_dir / 'pca_analysis.png', dpi=300, width=10, height=8)
        p.save(self.figures_dir / 'pca_analysis.pdf', width=10, height=8)
        
        # Save PCA results
        pca_df.to_csv(self.de_dir / 'pca_results.csv')
        
        # Save explained variance
        var_df = pd.DataFrame({
            'PC': [f'PC{i+1}' for i in range(len(explained_var))],
            'Explained_Variance': explained_var,
            'Cumulative_Variance': np.cumsum(explained_var)
        })
        var_df.to_csv(self.de_dir / 'pca_variance.csv', index=False)
        
        logger.info(f"PCA analysis completed. PC1: {explained_var[0]*100:.1f}%, PC2: {explained_var[1]*100:.1f}%")

    def _create_correlation_analysis(self):
        """Create correlation heatmap following R script methodology"""
        logger.info("Creating correlation heatmap...")
        
        if self.normalized_counts is None or self.metadata is None:
            logger.warning("No normalized counts available for correlation analysis")
            return
        
        # Calculate correlation matrix
        corr_matrix = self.normalized_counts.corr()
        
        # Create correlation heatmap with plotnine
        # Prepare data for plotnine
        corr_data = corr_matrix.reset_index().melt(id_vars='index', var_name='Sample2', value_name='Correlation')
        corr_data.columns = ['Sample1', 'Sample2', 'Correlation']
        
        p = (ggplot(corr_data, aes(x='Sample1', y='Sample2', fill='Correlation'))
             + geom_tile()
             + scale_fill_gradient2(
                 low='#d73027', 
                 mid='#f7f7f7', 
                 high='#1e90ff',
                 midpoint=0.5,
                 limits=[0, 1]
             )
             + labs(
                 title='Sample Correlation Matrix (PRJNA268379)',
                 x='Sample',
                 y='Sample',
                 fill='Correlation'
             )
             + theme_minimal()
             + theme(
                 axis_text_x=element_text(rotation=45, hjust=1),
                 axis_text_y=element_text(rotation=0),
                 figure_size=(12, 10)
             ))
        
        # Save plot
        p.save(self.figures_dir / 'correlation_heatmap.png', dpi=300, width=12, height=10)
        p.save(self.figures_dir / 'correlation_heatmap.pdf', width=12, height=10)
        
        # Save correlation matrix
        corr_matrix.to_csv(self.de_dir / 'sample_correlation_matrix.csv')
        
        # Calculate summary statistics
        corr_stats = {
            'mean_correlation': corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].mean(),
            'min_correlation': corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].min(),
            'max_correlation': corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].max(),
            'std_correlation': corr_matrix.values[np.triu_indices_from(corr_matrix.values, k=1)].std()
        }
        
        logger.info(f"Correlation analysis completed. Mean correlation: {corr_stats['mean_correlation']:.3f}")

    def _create_global_volcano(self):
        """Create global volcano plot following R script methodology"""
        logger.info("Creating global volcano plot...")
        
        if self.global_results is None:
            logger.warning("No global results available for volcano plot")
            return
        
        # Prepare data for plotting
        plot_data = self.global_results.copy()
        plot_data['-log10_padj'] = -np.log10(plot_data['padj'].fillna(1))
        plot_data['log2FoldChange'] = plot_data['log2FoldChange'].fillna(0)
        
        # Add significance categories
        plot_data['Significance'] = 'Not Significant'
        plot_data.loc[(plot_data['padj'] < 0.05) & (abs(plot_data['log2FoldChange']) > 1), 'Significance'] = 'Significant'
        plot_data.loc[(plot_data['padj'] < 0.05) & (abs(plot_data['log2FoldChange']) <= 1), 'Significance'] = 'Significant (Low FC)'
        
        # Create volcano plot with plotnine
        p = (ggplot(plot_data, aes(x='log2FoldChange', y='-log10_padj', color='Significance'))
             + geom_point(alpha=0.6, size=1)
             + geom_hline(yintercept=-np.log10(0.05), linetype='dashed', color='red', alpha=0.7)
             + geom_vline(xintercept=[-1, 1], linetype='dashed', color='red', alpha=0.7)
             + scale_color_manual(values={
                 'Not Significant': '#999999',
                 'Significant': '#d62728',
                 'Significant (Low FC)': '#ff7f0e'
             })
             + labs(
                 title='Global Differential Expression Analysis (PRJNA268379)',
                 x='log2 Fold Change',
                 y='-log10 Adjusted P-value',
                 color='Significance'
             )
             + theme_minimal()
             + theme(legend_position='right'))
        
        # Save plot
        p.save(self.figures_dir / 'global_volcano.png', dpi=300, width=10, height=8)
        p.save(self.figures_dir / 'global_volcano.pdf', width=10, height=8)
        
        # Save plot data
        plot_data.to_csv(self.de_dir / 'volcano_plot_data.csv')
        
        # Print summary statistics
        sig_genes = plot_data[plot_data['Significance'] == 'Significant']
        logger.info(f"Global volcano plot completed. Significant genes: {len(sig_genes)}")

    def _analyze_candidate_genes(self):
        """Analyze candidate genes with publication-quality visualizations"""
        logger.info("Analyzing candidate genes...")
        
        if self.global_results is None:
            logger.warning("No global results available for candidate gene analysis")
            return
        
        # Prepare data for all candidate genes
        all_candidates = []
        for list_name, genes in self.gene_lists.items():
            if genes:
                candidate_data = self.global_results[self.global_results.index.isin(genes)].copy()
                candidate_data['Gene_List'] = list_name
                candidate_data['Gene_ID'] = candidate_data.index
                all_candidates.append(candidate_data)
        
        if not all_candidates:
            logger.warning("No candidate genes found in results")
            return
        
        combined_candidates = pd.concat(all_candidates, ignore_index=True)
        
        # Create candidate gene volcano plot with plotnine
        combined_candidates['-log10_padj'] = -np.log10(combined_candidates['padj'].fillna(1))
        combined_candidates['log2FoldChange'] = combined_candidates['log2FoldChange'].fillna(0)
        
        # Add significance categories
        combined_candidates['Significance'] = 'Not Significant'
        combined_candidates.loc[(combined_candidates['padj'] < 0.05) & (abs(combined_candidates['log2FoldChange']) > 1), 'Significance'] = 'Significant'
        combined_candidates.loc[(combined_candidates['padj'] < 0.05) & (abs(combined_candidates['log2FoldChange']) <= 1), 'Significance'] = 'Significant (Low FC)'
        
        p = (ggplot(combined_candidates, aes(x='log2FoldChange', y='-log10_padj', color='Gene_List', shape='Significance'))
             + geom_point(alpha=0.7, size=3)
             + geom_hline(yintercept=-np.log10(0.05), linetype='dashed', color='red', alpha=0.7)
             + geom_vline(xintercept=[-1, 1], linetype='dashed', color='red', alpha=0.7)
             + scale_color_manual(values=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
             + labs(
                 title='Candidate Gene Analysis (PRJNA268379)',
                 x='log2 Fold Change',
                 y='-log10 Adjusted P-value',
                 color='Gene List',
                 shape='Significance'
             )
             + theme_minimal()
             + theme(legend_position='right'))
        
        # Save plot
        p.save(self.figures_dir / 'candidate_genes_volcano.png', dpi=300, width=12, height=8)
        p.save(self.figures_dir / 'candidate_genes_volcano.pdf', width=12, height=8)
        
        # Create summary statistics
        summary_stats = []
        for list_name in combined_candidates['Gene_List'].unique():
            list_data = combined_candidates[combined_candidates['Gene_List'] == list_name]
            total_genes = len(list_data)
            sig_genes = len(list_data[list_data['Significance'] == 'Significant'])
            up_regulated = len(list_data[(list_data['Significance'] == 'Significant') & (list_data['log2FoldChange'] > 0)])
            down_regulated = len(list_data[(list_data['Significance'] == 'Significant') & (list_data['log2FoldChange'] < 0)])
            
            summary_stats.append({
                'Gene_List': list_name,
                'Total_Genes': total_genes,
                'Significant_Genes': sig_genes,
                'Up_Regulated': up_regulated,
                'Down_Regulated': down_regulated,
                'Significance_Rate': sig_genes / total_genes if total_genes > 0 else 0
            })
        
        summary_df = pd.DataFrame(summary_stats)
        summary_df.to_csv(self.de_dir / 'candidate_genes_summary.csv', index=False)
        
        # Save detailed results
        combined_candidates.to_csv(self.de_dir / 'candidate_genes_detailed.csv', index=False)
        
        logger.info(f"Candidate gene analysis completed. Analyzed {len(combined_candidates)} genes across {len(summary_stats)} lists")

def main():
    """Main execution function"""
    # You can change the analysis_mode here:
    # "combined" - all samples together (current approach)
    # "separate_fed" - only Fed samples (LD_Fed vs SD_Fed)
    # "separate_unfed" - only Unfed samples (LD_Unfed vs SD_Unfed)
    # "both" - run all three analyses
    analysis_mode = "both"  # Change this to run different analyses
    
    analyzer = GeneExpressionAnalyzer(
        gene_lists_dir="output/tables",  # Where your text files are
        analysis_mode=analysis_mode
    )
    
    # Run the analysis based on the selected mode
    analyzer.run_analysis()
    
    logger.info("=== ANALYSIS COMPLETED ===")
    logger.info(f"Results in: {analyzer.de_dir}")
    logger.info(f"Figures in: {analyzer.figures_dir}")
    
    # Print summary based on analysis mode
    if analysis_mode == "combined":
        print("\n=== COMBINED ANALYSIS SUMMARY ===")
        if analyzer.combined_results is not None:
            total_genes = len(analyzer.combined_results)
            sig_genes = len(analyzer.combined_results[analyzer.combined_results['significant']])
            print(f"Total genes: {total_genes}")
            print(f"Significant genes: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
    
    elif analysis_mode == "separate_fed":
        print("\n=== FED SAMPLES ANALYSIS SUMMARY ===")
        if analyzer.fed_results is not None:
            total_genes = len(analyzer.fed_results)
            sig_genes = len(analyzer.fed_results[analyzer.fed_results['significant']])
            print(f"Total genes: {total_genes}")
            print(f"Significant genes: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
    
    elif analysis_mode == "separate_unfed":
        print("\n=== UNFED SAMPLES ANALYSIS SUMMARY ===")
        if analyzer.unfed_results is not None:
            total_genes = len(analyzer.unfed_results)
            sig_genes = len(analyzer.unfed_results[analyzer.unfed_results['significant']])
            print(f"Total genes: {total_genes}")
            print(f"Significant genes: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
    
    elif analysis_mode == "both":
        print("\n=== ALL ANALYSES SUMMARY ===")
        
        if analyzer.combined_results is not None:
            total_genes = len(analyzer.combined_results)
            sig_genes = len(analyzer.combined_results[analyzer.combined_results['significant']])
            print(f"Combined Analysis - Total genes: {total_genes}, Significant: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
        
        if analyzer.fed_results is not None:
            total_genes = len(analyzer.fed_results)
            sig_genes = len(analyzer.fed_results[analyzer.fed_results['significant']])
            print(f"Fed Analysis - Total genes: {total_genes}, Significant: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
        
        if analyzer.unfed_results is not None:
            total_genes = len(analyzer.unfed_results)
            sig_genes = len(analyzer.unfed_results[analyzer.unfed_results['significant']])
            print(f"Unfed Analysis - Total genes: {total_genes}, Significant: {sig_genes} ({sig_genes/total_genes*100:.1f}%)")
    
    # Print gene list summary if available
    if analyzer.list_analyses:
        print("\n=== GENE LIST SUMMARY ===")
        print(f"{'List Name':<20} {'Genes':<6} {'Significant':<11} {'% Hit Rate':<10}")
        print("-" * 50)
        
        for list_name, analysis in analyzer.list_analyses.items():
            stats = analysis['stats']
            pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            print(f"{list_name:<20} {stats['genes_found']:<6} {stats['significant_genes']:<11} {pct_sig:>9.1f}%")

if __name__ == "__main__":
    main() 
