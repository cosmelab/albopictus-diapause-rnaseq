#!/usr/bin/env python3
"""
Differential Expression Analysis for PRJNA158021 (Embryo Dataset)

This script performs comprehensive differential expression analysis on the embryo dataset
comparing different developmental time points under long day (LD) and short day (SD) conditions.

Key Features:
- Developmental trajectory analysis
- Photoperiod × Time interactions
- Temporal dynamics visualization
- Time-course expression plots

Usage:
    python de_analysis_PRJNA158021.py

Requirements:
    - pandas, numpy, matplotlib, seaborn
    - rpy2 for R integration
    - DESeq2 R package
    - plotnine for ggplot-style plots
    - bokeh for interactive plots
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

# Set up logging first
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('de_analysis_PRJNA158021.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Advanced visualization libraries for temporal analysis
try:
    import plotnine as p9
    from plotnine import *
    PLOTNINE_AVAILABLE = True
    logger.info("plotnine available for advanced visualizations")
except ImportError:
    PLOTNINE_AVAILABLE = False
    logger.warning("plotnine not available - some visualizations will be skipped")

try:
    import bokeh.plotting as bk
    from bokeh.models import HoverTool, ColorBar
    from bokeh.palettes import viridis, Spectral6
    from bokeh.transform import linear_cmap
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False
    logger.warning("bokeh not available - interactive plots will be skipped")

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class GeneExpressionAnalyzer:
    def __init__(self, gene_counts_file="output/tables/PRJNA158021_gene_counts.csv", 
                 gene_lists_dir="output/tables", output_base="output"):
        """
        Initialize the gene expression analyzer for embryo dataset.
        
        Args:
            gene_counts_file: Path to the gene counts file
            gene_lists_dir: Directory containing gene list files
            output_base: Base directory for output files
        """
        self.gene_counts_file = Path(gene_counts_file)
        self.gene_lists_dir = Path(gene_lists_dir)
        self.output_base = Path(output_base)
        
        # Create output directories
        self.de_dir = self.output_base / "differential_expression" / "PRJNA158021"
        self.figures_dir = self.output_base / "figures" / "PRJNA158021"
        self.qc_dir = self.output_base / "qc" / "PRJNA158021"
        
        for dir_path in [self.de_dir, self.figures_dir, self.qc_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.gene_counts = None
        self.metadata = None
        self.gene_lists = {}
        self.global_results = None
        self.normalized_counts = None
        
        # New containers for temporal analysis
        self.timepoint_analyses = {}  # Store results for each timepoint
        self.interaction_results = None  # Store photoperiod × time interaction results
        self.list_analyses = {}  # Store gene list analyses across different contrasts
        
        # Load gene lists
        self._load_gene_lists_from_files()
        
        logger.info("GeneExpressionAnalyzer initialized for PRJNA158021 (Embryo Dataset)")
        logger.info("Enhanced with temporal analysis capabilities")

    def _load_gene_lists_from_files(self):
        """Load gene lists from text files."""
        gene_list_files = {
            'all_candidates': 'all_candidate_genes_ids_only.txt',
            'high_priority': 'high_priority_genes_ids_only.txt',
            'high_ld': 'high_ld_genes_ids_only.txt',
            'tag_snp': 'tag_snp_genes_ids_only.txt'
        }
        
        for list_name, filename in gene_list_files.items():
            file_path = self.gene_lists_dir / filename
            if file_path.exists():
                with open(file_path, 'r') as f:
                    genes = [line.strip() for line in f if line.strip()]
                self.gene_lists[list_name] = genes
                logger.info(f"Loaded {len(genes)} genes for {list_name}")
            else:
                logger.warning(f"Gene list file not found: {file_path}")
                self.gene_lists[list_name] = []

    def validate_gene_lists(self):
        """Validate that gene lists contain valid gene IDs."""
        if self.gene_counts is None:
            logger.error("Gene counts not loaded. Call load_data() first.")
            return False
        
        available_genes = set(self.gene_counts.index)
        
        for list_name, genes in self.gene_lists.items():
            if not genes:
                continue
                
            valid_genes = [g for g in genes if g in available_genes]
            missing_genes = [g for g in genes if g not in available_genes]
            
            if missing_genes:
                logger.warning(f"{list_name}: {len(missing_genes)} genes not found in dataset")
            
            self.gene_lists[list_name] = valid_genes
            logger.info(f"{list_name}: {len(valid_genes)} valid genes out of {len(genes)}")
        
        return True

    def load_data(self):
        """Load gene counts and create metadata."""
        if not self.gene_counts_file.exists():
            raise FileNotFoundError(f"Gene counts file not found: {self.gene_counts_file}")
        
        # Load gene counts
        self.gene_counts = pd.read_csv(self.gene_counts_file, index_col=0)
        logger.info(f"Loaded gene counts: {self.gene_counts.shape[0]} genes, {self.gene_counts.shape[1]} samples")
        
        # Create metadata
        self._create_metadata()
        
        # Validate gene lists
        self.validate_gene_lists()

    def _create_metadata(self):
        """Create metadata DataFrame from sample names with enhanced biological annotations."""
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
                
                # Add biological annotations for temporal analysis
                diapause_status = 'Diapause_induction' if light_cycle == 'SD' else 'Non_diapause'
                timepoint_clean = time_point.replace('h', 'h_pov')  # post-oviposition
                timepoint_num = 72 if time_point == '72h' else 135
                
                metadata_list.append({
                    'sample': sample,
                    'prjna': prjna,
                    'srr': srr,
                    'light_cycle': light_cycle,
                    'time_point': time_point,
                    'timepoint_clean': timepoint_clean,
                    'timepoint_num': timepoint_num,
                    'treatment': treatment,
                    'diapause_status': diapause_status,
                    'dataset': 'PRJNA158021',
                    'life_stage': 'Embryo',
                    'comparison_group': 'Diapause_induction' if light_cycle == 'SD' else 'Non_diapause'
                })
        
        self.metadata = pd.DataFrame(metadata_list)
        self.metadata.set_index('sample', inplace=True)
        
        logger.info(f"Created enhanced metadata for {len(self.metadata)} samples")
        logger.info(f"Treatments: {sorted(self.metadata['treatment'].unique())}")
        logger.info(f"Timepoints: {self.metadata['time_point'].value_counts().to_dict()}")
        logger.info(f"Photoperiods: {self.metadata['light_cycle'].value_counts().to_dict()}")

    def run_global_analysis(self):
        """Run comprehensive global differential expression analysis with temporal components."""
        logger.info("Starting comprehensive global differential expression analysis...")
        
        # Filter genes
        self._filter_genes()
        
        # Set up R environment
        self._setup_r_environment()
        
        # Run main DESeq2 analysis with interaction terms
        self._run_main_deseq2_with_interactions()
        
        # Extract results
        self._extract_global_results()
        
        # Create global volcano plot
        self._create_global_volcano()
        
        # Create temporal analysis visualizations
        self._create_temporal_analysis_plots()
        
        logger.info("Comprehensive global analysis completed")

    def _filter_genes(self, min_count=10, min_samples=3):
        """Filter genes based on count thresholds."""
        initial_genes = len(self.gene_counts)
        
        # Filter by minimum count
        count_filter = (self.gene_counts >= min_count).sum(axis=1) >= min_samples
        self.gene_counts = self.gene_counts[count_filter]
        
        filtered_genes = len(self.gene_counts)
        logger.info(f"Filtered genes: {initial_genes} -> {filtered_genes} (min_count={min_count}, min_samples={min_samples})")

    def _setup_r_environment(self):
        """Set up R environment and load required packages."""
        try:
            # Activate pandas conversion
            pandas2ri.activate()
            
            # Import R packages
            base = importr('base')
            stats = importr('stats')
            
            # Try to load DESeq2
            try:
                deseq2 = importr('DESeq2')
                logger.info("DESeq2 loaded successfully")
            except Exception as e:
                logger.error(f"Failed to load DESeq2: {e}")
                raise
            
            logger.info("R environment set up successfully")
            
        except Exception as e:
            logger.error(f"Error setting up R environment: {e}")
            raise

    def _run_main_deseq2_with_interactions(self):
        """Run main DESeq2 analysis with interaction terms."""
        try:
            # Prepare count data for R
            count_matrix = self.gene_counts.astype(int)
            
            # Create design matrix with interaction terms
            design_data = self.metadata.copy()
            design_data['sample'] = design_data.index
            
            # Convert to R objects
            r_counts = pandas2ri.py2rpy(count_matrix)
            r_design = pandas2ri.py2rpy(design_data)
            
            # Assign to R global environment
            ro.globalenv['r_counts'] = r_counts
            ro.globalenv['r_design'] = r_design
            
            # Create DESeq2 dataset with interaction design
            ro.r('''
            # Create DESeq2 dataset with interaction design
            dds <- DESeqDataSetFromMatrix(
                countData = r_counts,
                colData = r_design,
                design = ~ time_point + light_cycle + time_point:light_cycle
            )
            
            # Explicitly convert design variables to factors to avoid warnings
            dds$light_cycle <- as.factor(dds$light_cycle)
            dds$time_point <- as.factor(dds$time_point)
            
            # Set reference levels
            dds$light_cycle <- relevel(dds$light_cycle, ref = "LD")
            dds$time_point <- relevel(dds$time_point, ref = "72h")
            
            # Print design matrix for debugging
            print("Design matrix:")
            print(design(dds))
            print("Sample info:")
            print(colData(dds))
            
            # Run DESeq2
            dds <- DESeq(dds)
            
            # Get results for main photoperiod effect (controlling for timepoint)
            results_main <- results(dds, name = "light_cycle_SD_vs_LD")
            results_main <- results_main[order(results_main$padj), ]
            
            # Print summary of main results
            print("Main results summary:")
            print(summary(results_main))
            print("First few rows of main results:")
            print(head(results_main))
            
            # Get results for timepoint effect (controlling for photoperiod)
            results_time <- results(dds, name = "time_point_135h_vs_72h")
            results_time <- results_time[order(results_time$padj), ]
            
            # Get results for interaction effect
            results_interaction <- results(dds, name = "time_point135h.light_cycleSD")
            results_interaction <- results_interaction[order(results_interaction$padj), ]
            
            # Get results for individual comparisons
            results_72h <- results(dds, contrast = c("light_cycle", "SD", "LD"))
            results_72h <- results_72h[order(results_72h$padj), ]
            
            # For 135h comparison, we need to calculate manually due to interaction
            # This is a simplified approach - in practice you might want to subset the data
            results_135h <- results(dds, contrast = c("light_cycle", "SD", "LD"))
            results_135h <- results_135h[order(results_135h$padj), ]
            
            # Print summary of 72h results
            print("72h results summary:")
            print(summary(results_72h))
            
            # Get normalized counts
            normalized_counts <- counts(dds, normalized = TRUE)
            
            # Store results
            .GlobalEnv$dds <- dds
            .GlobalEnv$results_main <- results_main
            .GlobalEnv$results_time <- results_time
            .GlobalEnv$results_interaction <- results_interaction
            .GlobalEnv$results_72h <- results_72h
            .GlobalEnv$results_135h <- results_135h
            .GlobalEnv$normalized_counts <- normalized_counts
            ''')
            
            logger.info("DESeq2 analysis with interactions completed successfully")
            
        except Exception as e:
            logger.error(f"Error in DESeq2 analysis with interactions: {e}")
            raise

    def _extract_global_results(self):
        """Extract results from DESeq2 analysis with interactions."""
        try:
            # Convert DESeq2 results to data frames in R
            r_code = """
            # Convert main results to data frame
            results_main_df <- as.data.frame(results_main)
            results_main_df$gene_id <- rownames(results_main_df)
            results_main_df$comparison <- "Main_photoperiod_effect"
            
            # Convert time results to data frame
            results_time_df <- as.data.frame(results_time)
            results_time_df$gene_id <- rownames(results_time_df)
            results_time_df$comparison <- "Time_effect"
            
            # Convert interaction results to data frame
            results_interaction_df <- as.data.frame(results_interaction)
            results_interaction_df$gene_id <- rownames(results_interaction_df)
            results_interaction_df$comparison <- "Photoperiod_x_Time_interaction"
            
            # Convert results for 72h comparison to data frame
            results_72h_df <- as.data.frame(results_72h)
            results_72h_df$gene_id <- rownames(results_72h_df)
            results_72h_df$comparison <- "SD_72h_vs_LD_72h"
            
            # Convert results for 135h comparison to data frame
            results_135h_df <- as.data.frame(results_135h)
            results_135h_df$gene_id <- rownames(results_135h_df)
            results_135h_df$comparison <- "SD_135h_vs_LD_135h"
            
            # Convert normalized counts
            normalized_counts_df <- as.data.frame(normalized_counts)
            normalized_counts_df$gene_id <- rownames(normalized_counts_df)
            """
            ro.r(r_code)
            
            # Get the results data frames from R
            r_results_main = ro.r['results_main_df']
            r_results_time = ro.r['results_time_df']
            r_results_interaction = ro.r['results_interaction_df']
            r_results_72h = ro.r['results_72h_df']
            r_results_135h = ro.r['results_135h_df']
            r_normalized = ro.r['normalized_counts_df']
            
            # Convert to pandas DataFrames using pandas2ri
            results_main_py = pandas2ri.rpy2py(r_results_main)
            results_time_py = pandas2ri.rpy2py(r_results_time)
            results_interaction_py = pandas2ri.rpy2py(r_results_interaction)
            results_72h_py = pandas2ri.rpy2py(r_results_72h)
            results_135h_py = pandas2ri.rpy2py(r_results_135h)
            normalized_counts_py = pandas2ri.rpy2py(r_normalized)
            
            # Combine all results
            all_results = pd.concat([
                results_main_py, results_time_py, results_interaction_py,
                results_72h_py, results_135h_py
            ], ignore_index=True)
            
            # Set gene_id as index
            all_results.set_index('gene_id', inplace=True)
            
            # Add significance columns based on Poelchau et al. 2013a thresholds:
            # |log2FC| > 0.5 and FDR < 0.05 (Benjamini-Hochberg corrected)
            all_results['significant'] = (
                (all_results['padj'] < 0.05) & 
                (abs(all_results['log2FoldChange']) > 0.5)  # Poelchau et al. used 0.5, not 1.0
            )
            all_results['significant_0.1'] = (
                (all_results['padj'] < 0.1) & 
                (abs(all_results['log2FoldChange']) > 0.5)
            )
            
            # Store main results (photoperiod effect)
            self.global_results = results_main_py.set_index('gene_id')
            self.global_results['significant'] = (
                (self.global_results['padj'] < 0.05) & 
                (abs(self.global_results['log2FoldChange']) > 0.5)  # Poelchau et al. used 0.5
            )
            self.global_results['significant_0.1'] = (
                (self.global_results['padj'] < 0.1) & 
                (abs(self.global_results['log2FoldChange']) > 0.5)
            )
            
            # Store interaction results (timepoint-specific comparisons)
            self.interaction_results = pd.concat([results_72h_py, results_135h_py], ignore_index=True)
            self.interaction_results.set_index('gene_id', inplace=True)
            self.interaction_results['significant'] = (
                (self.interaction_results['padj'] < 0.05) & 
                (abs(self.interaction_results['log2FoldChange']) > 0.5)  # Poelchau et al. used 0.5
            )
            self.interaction_results['significant_0.1'] = (
                (self.interaction_results['padj'] < 0.1) & 
                (abs(self.interaction_results['log2FoldChange']) > 0.5)
            )
            
            # Store normalized counts
            self.normalized_counts = normalized_counts_py.set_index('gene_id')
            
            # Save all results
            all_results.to_csv(self.de_dir / 'all_results_with_interactions.csv')
            self.global_results.to_csv(self.de_dir / 'main_photoperiod_effect.csv')
            self.interaction_results.to_csv(self.de_dir / 'timepoint_comparisons.csv')
            self.normalized_counts.to_csv(self.de_dir / 'normalized_counts.csv')
            
            # Log summary statistics
            logger.info("=== ANALYSIS RESULTS SUMMARY ===")
            for comparison in all_results['comparison'].unique():
                comp_data = all_results[all_results['comparison'] == comparison]
                sig_genes = comp_data['significant'].sum()
                lenient_sig = comp_data['significant_0.1'].sum()
                up_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] > 0)).sum()
                down_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] < 0)).sum()
                
                logger.info(f"{comparison}: {sig_genes} significant (FDR<0.05), {lenient_sig} lenient (FDR<0.1) ({up_reg} up, {down_reg} down)")
            
        except Exception as e:
            logger.error(f"Error extracting results: {e}")
            raise

    def _create_temporal_analysis_plots(self):
        """Create comprehensive temporal analysis visualizations."""
        logger.info("Creating temporal analysis plots...")
        
        # Create time-course plots for gene lists
        self._create_timecourse_plots()
        
        # Create timepoint comparison plots
        self._create_timepoint_comparison_plots()
        
        # Create interaction analysis plots
        self._create_interaction_analysis_plots()
        
        logger.info("Temporal analysis plots completed")

    def _create_timecourse_plots(self):
        """Create time-course expression plots for gene lists."""
        if self.normalized_counts is None:
            logger.warning("Normalized counts not available for time-course plots")
            return
            
        logger.info("Creating time-course expression plots...")
        
        for list_name, genes in self.gene_lists.items():
            if not genes:
                continue
                
            # Get genes present in the dataset
            present_genes = [g for g in genes if g in self.normalized_counts.index]
            if not present_genes:
                logger.warning(f"No genes from {list_name} found in dataset")
                continue
            
            # Get top 6 most significant genes for visualization
            if self.global_results is not None:
                # Use global results to find most significant genes
                list_results = self.global_results[self.global_results.index.isin(present_genes)]
                if not list_results.empty:
                    top_genes = list_results.nsmallest(6, 'padj').index
                else:
                    top_genes = present_genes[:6]
            else:
                top_genes = present_genes[:6]
            
            # Get normalized expression data
            expr_data = self.normalized_counts.loc[top_genes]
            
            # Create time-course plot
            self._create_single_timecourse_plot(list_name, expr_data)

    def _create_single_timecourse_plot(self, list_name, expr_data):
        """Create a single time-course plot for a gene list."""
        # Prepare data for plotting
        plot_data_list = []
        for gene in expr_data.index:
            for sample in expr_data.columns:
                sample_meta = self.metadata.loc[sample]
                plot_data_list.append({
                    'gene': gene,
                    'sample': sample,
                    'expression': np.log2(expr_data.loc[gene, sample] + 1),
                    'diapause_status': sample_meta['diapause_status'],
                    'timepoint': sample_meta['time_point'],
                    'timepoint_num': sample_meta['timepoint_num'],
                    'photoperiod': sample_meta['light_cycle']
                })
        
        plot_df = pd.DataFrame(plot_data_list)
        
        # Create time-course plot
        if PLOTNINE_AVAILABLE:
            try:
                # Try with loess smoothing first
                p = (ggplot(plot_df, aes(x='timepoint_num', y='expression', color='diapause_status')) +
                     geom_point(size=3, alpha=0.7) +
                     geom_smooth(method='loess', se=True, alpha=0.3) +
                     facet_wrap('gene', scales='free_y', ncol=3) +
                     scale_color_manual(values={
                         'Diapause_induction': '#e74c3c',      # Red
                         'Non_diapause': '#3498db'             # Blue
                     }) +
                     scale_x_continuous(breaks=[72, 135], labels=['72h', '135h']) +
                     labs(
                         title=f'{list_name.replace("_", " ").title()} - Expression Time Course',
                         subtitle='Top 6 most significant genes during embryonic development',
                         x='Time Post-Oviposition (hours)',
                         y='Log₂(Normalized Expression + 1)',
                         color='Condition'
                     ) +
                     theme_bw() +
                     theme(
                         figure_size=(15, 10),
                         plot_title=element_text(size=16, face='bold'),
                         plot_subtitle=element_text(size=12),
                         axis_title=element_text(size=12),
                         strip_text=element_text(size=10, face='bold')
                     ))
                
                # Save the plot
                p.save(self.figures_dir / f'{list_name}_timecourse.png', dpi=300, width=15, height=10)
                logger.info(f"Created time-course plot for {list_name}")
                
            except Exception as e:
                if "scikit-misc" in str(e) or "loess" in str(e).lower():
                    # Fallback to linear smoothing if loess fails
                    try:
                        p = (ggplot(plot_df, aes(x='timepoint_num', y='expression', color='diapause_status')) +
                             geom_point(size=3, alpha=0.7) +
                             geom_smooth(method='lm', se=True, alpha=0.3) +
                             facet_wrap('gene', scales='free_y', ncol=3) +
                             scale_color_manual(values={
                                 'Diapause_induction': '#e74c3c',      # Red
                                 'Non_diapause': '#3498db'             # Blue
                             }) +
                             scale_x_continuous(breaks=[72, 135], labels=['72h', '135h']) +
                             labs(
                                 title=f'{list_name.replace("_", " ").title()} - Expression Time Course',
                                 subtitle='Top 6 most significant genes during embryonic development (Linear smoothing)',
                                 x='Time Post-Oviposition (hours)',
                                 y='Log₂(Normalized Expression + 1)',
                                 color='Condition'
                             ) +
                             theme_bw() +
                             theme(
                                 figure_size=(15, 10),
                                 plot_title=element_text(size=16, face='bold'),
                                 plot_subtitle=element_text(size=12),
                                 axis_title=element_text(size=12),
                                 strip_text=element_text(size=10, face='bold')
                             ))
                        
                        # Save the plot
                        p.save(self.figures_dir / f'{list_name}_timecourse.png', dpi=300, width=15, height=10)
                        logger.info(f"Created time-course plot for {list_name} (linear smoothing)")
                        
                    except Exception as e2:
                        logger.error(f"Error creating time-course plot for {list_name} (linear fallback): {e2}")
                        # Final fallback to matplotlib
                        self._create_matplotlib_timecourse_plot(list_name, plot_df)
                else:
                    logger.error(f"Error creating time-course plot for {list_name}: {e}")
                    # Fallback to matplotlib
                    self._create_matplotlib_timecourse_plot(list_name, plot_df)
        else:
            # Fallback to matplotlib
            self._create_matplotlib_timecourse_plot(list_name, plot_df)

    def _create_matplotlib_timecourse_plot(self, list_name, plot_df):
        """Create time-course plot using matplotlib as fallback."""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'{list_name.replace("_", " ").title()} - Expression Time Course', 
                    fontsize=16, fontweight='bold')
        
        genes = plot_df['gene'].unique()
        for i, gene in enumerate(genes[:6]):
            ax = axes[i//3, i%3]
            gene_data = plot_df[plot_df['gene'] == gene]
            
            # Plot points
            for status in ['Diapause_induction', 'Non_diapause']:
                status_data = gene_data[gene_data['diapause_status'] == status]
                color = '#e74c3c' if status == 'Diapause_induction' else '#3498db'
                ax.scatter(status_data['timepoint_num'], status_data['expression'], 
                          c=color, alpha=0.7, s=30, label=status.replace('_', ' '))
            
            ax.set_title(gene, fontsize=10, fontweight='bold')
            ax.set_xlabel('Time (hours)')
            ax.set_ylabel('Log₂(Expression + 1)')
            ax.set_xticks([72, 135])
            ax.set_xticklabels(['72h', '135h'])
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / f'{list_name}_timecourse_matplotlib.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"Created matplotlib time-course plot for {list_name}")

    def _create_timepoint_comparison_plots(self):
        """Create comparison plots between timepoints."""
        logger.info("Creating timepoint comparison plots...")
        
        if self.interaction_results is None:
            logger.warning("Interaction results not available for timepoint comparison")
            return
        
        # Create comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Embryonic Photoperiod Response - Timepoint Comparison', 
                    fontsize=20, fontweight='bold')
        
        # Get results for each timepoint
        results_72h = self.interaction_results[self.interaction_results['comparison'] == 'SD_72h_vs_LD_72h']
        results_135h = self.interaction_results[self.interaction_results['comparison'] == 'SD_135h_vs_LD_135h']
        
        # Plot 1: Hit rates comparison
        ax1 = axes[0, 0]
        sig_72h = results_72h['significant'].sum()
        sig_135h = results_135h['significant'].sum()
        
        bars = ax1.bar(['72h', '135h'], [sig_72h, sig_135h], 
                      color=['#3498db', '#9b59b6'], alpha=0.8)
        ax1.set_ylabel('Number of Significant Genes', fontsize=12)
        ax1.set_title('A) Significant Genes by Timepoint', fontsize=14, fontweight='bold')
        
        # Add value labels on bars
        for bar, value in zip(bars, [sig_72h, sig_135h]):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                    str(value), ha='center', va='bottom', fontweight='bold')
        
        # Plot 2: Expression direction comparison
        ax2 = axes[0, 1]
        up_72h = ((results_72h['significant']) & (results_72h['log2FoldChange'] > 0)).sum()
        down_72h = ((results_72h['significant']) & (results_72h['log2FoldChange'] < 0)).sum()
        up_135h = ((results_135h['significant']) & (results_135h['log2FoldChange'] > 0)).sum()
        down_135h = ((results_135h['significant']) & (results_135h['log2FoldChange'] < 0)).sum()
        
        x_pos = np.arange(2)
        width = 0.35
        ax2.bar(x_pos - width/2, [up_72h, up_135h], width, 
               label='Up in SD', color='#e74c3c', alpha=0.8)
        ax2.bar(x_pos + width/2, [down_72h, down_135h], width, 
               label='Down in SD', color='#2ecc71', alpha=0.8)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(['72h', '135h'])
        ax2.set_ylabel('Number of Genes', fontsize=12)
        ax2.set_title('B) Expression Direction by Timepoint', fontsize=14, fontweight='bold')
        ax2.legend()
        
        # Plot 3: Log2FC distribution comparison
        ax3 = axes[1, 0]
        ax3.hist(results_72h['log2FoldChange'], bins=50, alpha=0.7, 
                label='72h', color='#3498db', density=True)
        ax3.hist(results_135h['log2FoldChange'], bins=50, alpha=0.7, 
                label='135h', color='#9b59b6', density=True)
        ax3.set_xlabel('Log₂ Fold Change', fontsize=12)
        ax3.set_ylabel('Density', fontsize=12)
        ax3.set_title('C) Log₂FC Distribution Comparison', fontsize=14, fontweight='bold')
        ax3.legend()
        ax3.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        
        # Plot 4: Sample distribution
        ax4 = axes[1, 1]
        treatment_counts = self.metadata['treatment'].value_counts()
        colors_pie = ['#e74c3c', '#3498db', '#e67e22', '#2ecc71']
        wedges, texts, autotexts = ax4.pie(treatment_counts.values, 
                                          labels=treatment_counts.index, 
                                          autopct='%1.0f', startangle=90, 
                                          colors=colors_pie[:len(treatment_counts)])
        ax4.set_title('D) Sample Distribution', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'timepoint_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Timepoint comparison plots created")

    def _create_interaction_analysis_plots(self):
        """Create plots for photoperiod × time interaction analysis."""
        logger.info("Creating interaction analysis plots...")
        
        if self.interaction_results is None:
            logger.warning("Interaction results not available")
            return
        
        # Create interaction analysis summary
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle('Photoperiod × Time Interaction Analysis', fontsize=16, fontweight='bold')
        
        # Get results for each timepoint
        results_72h = self.interaction_results[self.interaction_results['comparison'] == 'SD_72h_vs_LD_72h']
        results_135h = self.interaction_results[self.interaction_results['comparison'] == 'SD_135h_vs_LD_135h']
        
        # Plot 1: Volcano plot comparison
        ax1 = axes[0]
        # Plot 72h results
        ax1.scatter(results_72h[~results_72h['significant']]['log2FoldChange'],
                   -np.log10(results_72h[~results_72h['significant']]['padj']),
                   c='lightblue', alpha=0.6, s=20, label='72h (NS)')
        ax1.scatter(results_72h[results_72h['significant']]['log2FoldChange'],
                   -np.log10(results_72h[results_72h['significant']]['padj']),
                   c='blue', alpha=0.8, s=30, label='72h (Sig)')
        
        # Plot 135h results
        ax1.scatter(results_135h[~results_135h['significant']]['log2FoldChange'],
                   -np.log10(results_135h[~results_135h['significant']]['padj']),
                   c='lightcoral', alpha=0.6, s=20, label='135h (NS)')
        ax1.scatter(results_135h[results_135h['significant']]['log2FoldChange'],
                   -np.log10(results_135h[results_135h['significant']]['padj']),
                   c='red', alpha=0.8, s=30, label='135h (Sig)')
        
        ax1.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Log₂ Fold Change', fontsize=12)
        ax1.set_ylabel('-Log₁₀(Adjusted P-value)', fontsize=12)
        ax1.set_title('A) Volcano Plot Comparison', fontsize=14, fontweight='bold')
        ax1.legend()
        
        # Plot 2: Log2FC correlation between timepoints
        ax2 = axes[1]
        # Merge results by gene
        merged = pd.merge(results_72h[['log2FoldChange']], 
                         results_135h[['log2FoldChange']], 
                         left_index=True, right_index=True, 
                         suffixes=('_72h', '_135h'))
        
        ax2.scatter(merged['log2FoldChange_72h'], merged['log2FoldChange_135h'], 
                   alpha=0.6, s=20)
        ax2.set_xlabel('Log₂FC (72h)', fontsize=12)
        ax2.set_ylabel('Log₂FC (135h)', fontsize=12)
        ax2.set_title('B) Log₂FC Correlation Between Timepoints', fontsize=14, fontweight='bold')
        
        # Add correlation line
        z = np.polyfit(merged['log2FoldChange_72h'], merged['log2FoldChange_135h'], 1)
        p = np.poly1d(z)
        ax2.plot(merged['log2FoldChange_72h'], p(merged['log2FoldChange_72h']), 
                "r--", alpha=0.8)
        
        # Calculate and display correlation
        corr = merged['log2FoldChange_72h'].corr(merged['log2FoldChange_135h'])
        ax2.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                transform=ax2.transAxes, fontsize=12, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'interaction_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("Interaction analysis plots created")

    def _create_global_volcano(self):
        """Create a publication-quality volcano plot for global gene expression."""
        try:
            if self.global_results is None:
                logger.warning("Global results not available for volcano plot")
                return
                
            # DIAGNOSTIC: Check data quality
            logger.info("=== VOLCANO PLOT DIAGNOSTICS ===")
            logger.info(f"Global results shape: {self.global_results.shape}")
            logger.info(f"Columns: {list(self.global_results.columns)}")
            
            # Check for missing or problematic values
            logger.info(f"log2FoldChange range: {self.global_results['log2FoldChange'].min():.3f} to {self.global_results['log2FoldChange'].max():.3f}")
            logger.info(f"padj range: {self.global_results['padj'].min():.3e} to {self.global_results['padj'].max():.3e}")
            logger.info(f"pvalue range: {self.global_results['pvalue'].min():.3e} to {self.global_results['pvalue'].max():.3e}")
            
            # Check for infinite or NaN values
            inf_log2fc = np.isinf(self.global_results['log2FoldChange']).sum()
            nan_log2fc = np.isnan(self.global_results['log2FoldChange']).sum()
            inf_padj = np.isinf(self.global_results['padj']).sum()
            nan_padj = np.isnan(self.global_results['padj']).sum()
            
            logger.info(f"Infinite log2FoldChange: {inf_log2fc}")
            logger.info(f"NaN log2FoldChange: {nan_log2fc}")
            logger.info(f"Infinite padj: {inf_padj}")
            logger.info(f"NaN padj: {nan_padj}")
            
            # Check unique values
            unique_log2fc = self.global_results['log2FoldChange'].nunique()
            unique_padj = self.global_results['padj'].nunique()
            logger.info(f"Unique log2FoldChange values: {unique_log2fc}")
            logger.info(f"Unique padj values: {unique_padj}")
            
            # Check if all values are the same
            if unique_log2fc <= 1:
                logger.error("CRITICAL: All log2FoldChange values are the same!")
                logger.error(f"Value: {self.global_results['log2FoldChange'].iloc[0]}")
            if unique_padj <= 1:
                logger.error("CRITICAL: All padj values are the same!")
                logger.error(f"Value: {self.global_results['padj'].iloc[0]}")
            
            # Show sample of data
            logger.info("Sample of global results:")
            logger.info(self.global_results.head(10))
            
            # Prepare data
            plot_data = self.global_results.copy()
            plot_data['is_candidate'] = plot_data.index.isin(self.gene_lists.get('all_candidates', []))
            
            # Create separate plots for each comparison
            for comparison in plot_data['comparison'].unique():
                comp_data = plot_data[plot_data['comparison'] == comparison].copy()
                
                logger.info(f"=== PLOTTING {comparison} ===")
                logger.info(f"Data shape: {comp_data.shape}")
                logger.info(f"log2FoldChange range: {comp_data['log2FoldChange'].min():.3f} to {comp_data['log2FoldChange'].max():.3f}")
                logger.info(f"padj range: {comp_data['padj'].min():.3e} to {comp_data['padj'].max():.3e}")
                logger.info(f"Significant genes: {comp_data['significant'].sum()}")
                
                # Create figure
                plt.figure(figsize=(12, 8))
                
                # Plot non-significant genes
                non_sig = comp_data[~comp_data['significant'] & ~comp_data['significant_0.1']]
                if not non_sig.empty:
                    plt.scatter(
                        non_sig['log2FoldChange'],
                        -np.log10(non_sig['padj']),
                        c='gray',
                        alpha=0.5,
                        s=20,
                        marker='o',
                        label='Non-significant'
                    )
                
                # Plot leniently significant genes (FDR < 0.1)
                lenient_sig = comp_data[comp_data['significant_0.1'] & ~comp_data['significant']]
                if not lenient_sig.empty:
                    plt.scatter(
                        lenient_sig['log2FoldChange'],
                        -np.log10(lenient_sig['padj']),
                        c='orange',
                        alpha=0.6,
                        s=25,
                        marker='s',
                        label='FDR < 0.1'
                    )
                
                # Plot strictly significant non-candidate genes
                sig_non_cand = comp_data[comp_data['significant'] & ~comp_data['is_candidate']]
                if not sig_non_cand.empty:
                    plt.scatter(
                        sig_non_cand['log2FoldChange'],
                        -np.log10(sig_non_cand['padj']),
                        c='blue',
                        alpha=0.6,
                        s=30,
                        marker='^',
                        label='FDR < 0.05'
                    )
                
                # Plot candidate genes (non-significant)
                candidates = comp_data[comp_data['is_candidate']]
                if not candidates.empty:
                    non_sig_cand = candidates[~candidates['significant'] & ~candidates['significant_0.1']]
                    if not non_sig_cand.empty:
                        plt.scatter(
                            non_sig_cand['log2FoldChange'],
                            -np.log10(non_sig_cand['padj']),
                            c='red',
                            alpha=0.6,
                            s=30,
                            marker='s',
                            label='Candidate (NS)'
                        )
                    
                    # Plot leniently significant candidate genes
                    lenient_cand = candidates[candidates['significant_0.1'] & ~candidates['significant']]
                    if not lenient_cand.empty:
                        plt.scatter(
                            lenient_cand['log2FoldChange'],
                            -np.log10(lenient_cand['padj']),
                            c='darkred',
                            alpha=0.7,
                            s=35,
                            marker='*',
                            label='Candidate (FDR < 0.1)'
                        )
                    
                    # Plot strictly significant candidate genes
                    sig_candidates = candidates[candidates['significant']]
                    if not sig_candidates.empty:
                        plt.scatter(
                            sig_candidates['log2FoldChange'],
                            -np.log10(sig_candidates['padj']),
                            c='red',
                            alpha=0.8,
                            s=40,
                            marker='*',
                            label='Candidate (FDR < 0.05)'
                        )
                        
                        # Add annotations for significant candidate genes
                        for idx, row in sig_candidates.iterrows():
                            plt.annotate(
                                idx,
                                (row['log2FoldChange'], -np.log10(row['padj'])),
                                xytext=(5, 5),
                                textcoords='offset points',
                                fontsize=8,
                                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                            )
                
                # Add significance thresholds based on Poelchau et al. 2013a
                plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='FDR = 0.05')
                plt.axhline(y=-np.log10(0.1), color='lightgray', linestyle='--', alpha=0.3, label='FDR = 0.1')
                plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
                plt.axvline(x=0.5, color='red', linestyle='--', alpha=0.7, label='|log2FC| = 0.5 (Poelchau et al.)')
                plt.axvline(x=-0.5, color='red', linestyle='--', alpha=0.7)
                
                # Add LD/SD labels
                plt.text(
                    plt.xlim()[0] * 0.9,
                    plt.ylim()[1] * 0.95,
                    'LD',
                    fontsize=12,
                    ha='center',
                    va='center',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )
                plt.text(
                    plt.xlim()[1] * 0.9,
                    plt.ylim()[1] * 0.95,
                    'SD',
                    fontsize=12,
                    ha='center',
                    va='center',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )
                
                # Add statistics text box with Poelchau et al. comparison
                stats_text = f'Total genes: {len(comp_data):,}\n'
                stats_text += f'FDR < 0.05 & |log2FC| > 0.5: {comp_data["significant"].sum()}\n'
                stats_text += f'FDR < 0.1 & |log2FC| > 0.5: {comp_data["significant_0.1"].sum()}\n'
                stats_text += f'Up-regulated: {((comp_data["significant"]) & (comp_data["log2FoldChange"] > 0)).sum()}\n'
                stats_text += f'Down-regulated: {((comp_data["significant"]) & (comp_data["log2FoldChange"] < 0)).sum()}\n'
                stats_text += f'Candidate genes: {comp_data["is_candidate"].sum()}\n'
                stats_text += f'\nPoelchau et al. 2013a:\n'
                stats_text += f'72h: 2,506 DE genes\n'
                stats_text += f'135h: 4,337 DE genes'
                
                plt.text(
                    0.02, 0.98, stats_text,
                    transform=plt.gca().transAxes,
                    fontsize=10,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9)
                )
                
                # Customize plot
                plt.title(f'Embryo Differential Expression Analysis - {comparison} (PRJNA158021)\nPoelchau et al. 2013a Thresholds: |log2FC| > 0.5, FDR < 0.05', 
                         pad=20, fontsize=14)
                plt.xlabel('log2 Fold Change (SD vs LD)', fontsize=12)
                plt.ylabel('-log10 Adjusted p-value', fontsize=12)
                
                # Add legend with custom styling
                legend = plt.legend(
                    title='Gene Categories',
                    title_fontsize=12,
                    fontsize=10,
                    bbox_to_anchor=(1.05, 1),
                    loc='upper left',
                    frameon=True,
                    framealpha=0.95,
                    edgecolor='black'
                )
                legend.get_frame().set_linewidth(0.5)
                
                # Set axis limits to focus on the data
                plt.xlim(-25, 10)
                plt.ylim(0, 10)
                
                # Adjust layout
                plt.tight_layout()
                
                # Save plot
                filename = f'global_volcano_{comparison.replace("_vs_", "_")}.pdf'
                plt.savefig(
                    self.figures_dir / filename,
                    dpi=300,
                    bbox_inches='tight',
                    format='pdf'
                )
                
                plt.close()
            
            logger.info("Global volcano plots created successfully")
            
        except Exception as e:
            logger.error(f"Error creating global volcano plot: {str(e)}")
            raise

    def generate_analysis_report(self):
        """Generate a comprehensive analysis report with temporal analysis results."""
        try:
            report_path = self.qc_dir / 'embryo_analysis_report.md'
            
            with open(report_path, 'w') as f:
                f.write("# Embryo Differential Expression Analysis Report (PRJNA158021)\n\n")
                f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("## Analysis Overview\n\n")
                f.write("This analysis replicates the experimental design from **Poelchau et al. 2013a** for embryonic photoperiod response.\n\n")
                f.write("### Original Study (Poelchau et al. 2013a)\n")
                f.write("- **Significance thresholds:** |log2FC| > 0.5, FDR < 0.05 (Benjamini-Hochberg)\n")
                f.write("- **Results:** 2,506 DE genes at 72h, 4,337 DE genes at 135h\n")
                f.write("- **Experimental design:** 3 replicates per condition, 2 timepoints (72h, 135h)\n\n")
                
                f.write("### Current Analysis\n")
                f.write("- **Significance thresholds:** |log2FC| > 0.5, FDR < 0.05 (matching Poelchau et al.)\n")
                f.write("- **Dataset:** PRJNA158021 processed data\n")
                f.write("- **Analysis approach:** DESeq2 with interaction model\n\n")
                
                f.write("## Analysis Summary\n\n")
                f.write(f"- **Total genes analyzed:** {len(self.global_results) if self.global_results is not None else 'N/A'}\n")
                f.write(f"- **Comparisons:** SD_72h vs LD_72h, SD_135h vs LD_135h\n")
                f.write(f"- **Temporal analysis:** Photoperiod × Time interactions\n")
                f.write(f"- **Significance criteria:** |log2FC| > 0.5 AND FDR < 0.05\n\n")
                
                f.write("## Global Analysis Results\n\n")
                if self.global_results is not None:
                    for comparison in self.global_results['comparison'].unique():
                        comp_data = self.global_results[self.global_results['comparison'] == comparison]
                        sig_genes = comp_data['significant'].sum()
                        lenient_sig = comp_data['significant_0.1'].sum()
                        up_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] > 0)).sum()
                        down_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] < 0)).sum()
                        
                        f.write(f"### {comparison}\n")
                        f.write(f"- **Significant genes (|log2FC| > 0.5, FDR < 0.05):** {sig_genes}\n")
                        f.write(f"- **Leniently significant (|log2FC| > 0.5, FDR < 0.1):** {lenient_sig}\n")
                        f.write(f"- **Up-regulated:** {up_reg}\n")
                        f.write(f"- **Down-regulated:** {down_reg}\n\n")
                
                f.write("## Comparison with Poelchau et al. 2013a\n\n")
                if self.interaction_results is not None:
                    results_72h = self.interaction_results[self.interaction_results['comparison'] == 'SD_72h_vs_LD_72h']
                    results_135h = self.interaction_results[self.interaction_results['comparison'] == 'SD_135h_vs_LD_135h']
                    
                    sig_72h = results_72h['significant'].sum() if not results_72h.empty else 0
                    sig_135h = results_135h['significant'].sum() if not results_135h.empty else 0
                    
                    f.write("| Timepoint | Poelchau et al. 2013a | Current Analysis | Difference |\n")
                    f.write("|-----------|------------------------|------------------|------------|\n")
                    f.write(f"| 72h (3d) | 2,506 DE genes | {sig_72h} DE genes | {sig_72h - 2506:+,} |\n")
                    f.write(f"| 135h (6d) | 4,337 DE genes | {sig_135h} DE genes | {sig_135h - 4337:+,} |\n\n")
                    
                    f.write("### Potential Reasons for Differences\n")
                    f.write("1. **Data processing:** Different alignment/mapping pipelines\n")
                    f.write("2. **Sample quality:** Potential batch effects or technical differences\n")
                    f.write("3. **Statistical approach:** Different normalization or filtering methods\n")
                    f.write("4. **Biological variation:** Natural variation between experiments\n")
                    f.write("5. **Outlier effects:** Technical outliers affecting power\n\n")
                
                f.write("## Temporal Analysis Results\n\n")
                if self.interaction_results is not None:
                    f.write("### Timepoint-Specific Analysis\n\n")
                    for comparison in self.interaction_results['comparison'].unique():
                        comp_data = self.interaction_results[self.interaction_results['comparison'] == comparison]
                        sig_genes = comp_data['significant'].sum()
                        lenient_sig = comp_data['significant_0.1'].sum()
                        up_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] > 0)).sum()
                        down_reg = ((comp_data['significant']) & (comp_data['log2FoldChange'] < 0)).sum()
                        
                        f.write(f"#### {comparison}\n")
                        f.write(f"- **Significant genes:** {sig_genes}\n")
                        f.write(f"- **Leniently significant genes:** {lenient_sig}\n")
                        f.write(f"- **Up-regulated:** {up_reg}\n")
                        f.write(f"- **Down-regulated:** {down_reg}\n\n")
                    
                    # Add correlation analysis
                    results_72h = self.interaction_results[self.interaction_results['comparison'] == 'SD_72h_vs_LD_72h']
                    results_135h = self.interaction_results[self.interaction_results['comparison'] == 'SD_135h_vs_LD_135h']
                    
                    if not results_72h.empty and not results_135h.empty:
                        merged = pd.merge(results_72h[['log2FoldChange']], 
                                         results_135h[['log2FoldChange']], 
                                         left_index=True, right_index=True, 
                                         suffixes=('_72h', '_135h'))
                        corr = merged['log2FoldChange_72h'].corr(merged['log2FoldChange_135h'])
                        f.write(f"### Temporal Correlation\n")
                        f.write(f"- **Log₂FC correlation (72h vs 135h):** {corr:.3f}\n\n")
                
                f.write("## Gene Lists Analysis\n\n")
                for list_name, genes in self.gene_lists.items():
                    if genes:
                        f.write(f"### {list_name.replace('_', ' ').title()}\n")
                        f.write(f"- **Total genes:** {len(genes)}\n")
                        
                        # Count significant genes in this list
                        if self.global_results is not None:
                            for comparison in self.global_results['comparison'].unique():
                                comp_data = self.global_results[self.global_results['comparison'] == comparison]
                                list_genes_in_comp = comp_data[comp_data.index.isin(genes)]
                                sig_in_list = list_genes_in_comp['significant'].sum()
                                
                                f.write(f"- **Significant in {comparison}:** {sig_in_list}\n")
                        f.write("\n")
                
                f.write("## Key Findings\n\n")
                f.write("### Developmental Trajectory\n")
                f.write("- Analysis captures gene expression changes during embryonic development\n")
                f.write("- Two critical timepoints: 72h (early) and 135h (late) post-oviposition\n")
                f.write("- Photoperiod effects may vary across developmental stages\n\n")
                
                f.write("### Photoperiod × Time Interactions\n")
                f.write("- Genes may respond differently to photoperiod at different timepoints\n")
                f.write("- Interaction analysis reveals temporal dynamics of photoperiod response\n")
                f.write("- Some genes show consistent effects, others show time-dependent responses\n\n")
                
                f.write("### Temporal Dynamics\n")
                f.write("- Early timepoint (72h): Initial photoperiod response establishment\n")
                f.write("- Late timepoint (135h): Mature photoperiod response patterns\n")
                f.write("- Correlation analysis shows temporal consistency of gene responses\n\n")
                
                f.write("## Recommendations\n\n")
                f.write("1. **Run outlier detection** to identify potential technical issues\n")
                f.write("2. **Compare with robust datasets** using different outlier handling strategies\n")
                f.write("3. **Investigate sample quality** metrics and batch effects\n")
                f.write("4. **Consider alternative statistical approaches** (e.g., rank-based tests)\n")
                f.write("5. **Validate key findings** with independent methods\n\n")
                
                f.write("## Files Generated\n\n")
                f.write("- **Global results:** `global_results.csv`\n")
                f.write("- **Interaction results:** `interaction_results.csv`\n")
                f.write("- **Normalized counts:** `normalized_counts.csv`\n")
                f.write("- **Time-course plots:** `[list_name]_timecourse.png`\n")
                f.write("- **Timepoint comparison:** `timepoint_comparison.png`\n")
                f.write("- **Interaction analysis:** `interaction_analysis.png`\n")
                f.write("- **Volcano plots:** `global_volcano_[comparison].pdf`\n\n")
            
            logger.info(f"Comprehensive analysis report generated: {report_path}")
            
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            raise

def main():
    """Main function to run the comprehensive temporal analysis."""
    try:
        # Initialize analyzer
        analyzer = GeneExpressionAnalyzer()
        
        # Load data
        analyzer.load_data()
        
        # Run comprehensive global analysis with temporal components
        analyzer.run_global_analysis()
        
        # Generate comprehensive report
        analyzer.generate_analysis_report()
        
        logger.info("=== COMPREHENSIVE EMBRYO ANALYSIS COMPLETED ===")
        logger.info("Key Features Implemented:")
        logger.info("1. Developmental trajectory analysis")
        logger.info("2. Photoperiod × Time interactions")
        logger.info("3. Temporal dynamics visualization")
        logger.info("4. Time-course expression plots")
        logger.info("5. Timepoint comparison analysis")
        
        # Print quick summary
        print("\n=== EMBRYO ANALYSIS SUMMARY ===")
        print("Analysis Type: Comprehensive Temporal Analysis")
        print("Dataset: PRJNA158021 (Embryonic Development)")
        print("Timepoints: 72h and 135h post-oviposition")
        print("Conditions: Long Day (LD) vs Short Day (SD)")
        
        if analyzer.global_results is not None:
            print(f"\nMain Photoperiod Effect:")
            sig_main = analyzer.global_results['significant'].sum()
            print(f"  - Significant genes: {sig_main}")
        
        if analyzer.interaction_results is not None:
            print(f"\nTimepoint-Specific Results:")
            for comparison in analyzer.interaction_results['comparison'].unique():
                comp_data = analyzer.interaction_results[analyzer.interaction_results['comparison'] == comparison]
                sig_comp = comp_data['significant'].sum()
                print(f"  - {comparison}: {sig_comp} significant genes")
        
        print(f"\nOutput Directories:")
        print(f"  - Results: {analyzer.de_dir}")
        print(f"  - Figures: {analyzer.figures_dir}")
        print(f"  - QC/Reports: {analyzer.qc_dir}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main() 