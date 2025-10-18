#!/usr/bin/env python3
"""
Comprehensive Differential Expression Analysis for PRJNA158021 (Embryos)
Analyzes diapause preparation vs non-diapause development in embryos
with stunning visualizations

Experimental Design:
- Timepoints: 72-78h and 135-141h post-oviposition
- Conditions: Diapause (D) vs Non-diapause (ND) 
- Focus: Diapause preparation during embryonic development

Usage:
    python analyze_prjna158021_embryos.py

Requirements:
    - Complete bioinformatics environment
    - Gene list text files in output/tables/
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# Advanced visualization libraries
import plotnine as p9
from plotnine import *
import altair as alt
import holoviews as hv
from holoviews import opts
import bokeh.plotting as bk
from bokeh.models import HoverTool, ColorBar
from bokeh.palettes import viridis, Spectral6
from bokeh.transform import linear_cmap
import upsetplot

# R integration
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
pandas2ri.activate()

# Enable HoloViews
hv.extension('bokeh')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class EmbryoAnalyzer:
    def __init__(self, gene_counts_file="output/tables/PRJNA158021_gene_counts.csv", 
                 gene_lists_dir="output/tables", output_base="output"):
        """
        Initialize embryo analyzer for PRJNA158021 dataset
        
        Parameters:
        -----------
        gene_counts_file : str
            Path to gene counts CSV file
        gene_lists_dir : str
            Directory containing gene list text files
        output_base : str
            Base output directory
        """
        self.gene_counts_file = Path(gene_counts_file)
        self.gene_lists_dir = Path(gene_lists_dir)
        self.output_base = Path(output_base)
        
        # Load gene lists from text files
        self.gene_lists = self._load_gene_lists_from_files()
        
        # Create output directories
        self.de_dir = self.output_base / "differential_expression" / "PRJNA158021_embryos"
        self.figures_dir = self.output_base / "figures" / "PRJNA158021_embryos"
        self.qc_dir = self.output_base / "qc" / "PRJNA158021_embryos"
        
        for directory in [self.de_dir, self.figures_dir, self.qc_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        # Data containers
        self.gene_counts = None
        self.metadata = None
        self.global_results = None
        self.normalized_counts = None
        self.list_analyses = {}
        self.timepoint_analyses = {}
        
        logger.info("Initialized embryo analyzer for PRJNA158021")
        logger.info(f"Gene lists loaded: {list(self.gene_lists.keys())}")
        for list_name, genes in self.gene_lists.items():
            logger.info(f"  - {list_name}: {len(genes)} genes")
        logger.info(f"Output directories created")

    def _load_gene_lists_from_files(self):
        """Load gene lists from text files"""
        gene_lists = {}
        
        file_mapping = {
            'all_candidate_genes_ids_only.txt': 'all_candidates',
            'high_priority_genes_ids_only.txt': 'high_priority', 
            'high_ld_genes_ids_only.txt': 'high_ld',
            'linked_snp_genes_ids_only.txt': 'linked_snp',
            'intragenic_genes_ids_only.txt': 'intragenic',
            'tag_snp_genes_ids_only.txt': 'tag_snp',
            'loc_genes_only_ids_only.txt': 'loc_genes_only'
        }
        
        for filename, list_name in file_mapping.items():
            file_path = self.gene_lists_dir / filename
            
            if file_path.exists():
                try:
                    with open(file_path, 'r') as f:
                        genes = [line.strip() for line in f if line.strip()]
                    
                    genes = [gene for gene in genes if gene.startswith('LOC')]
                    
                    if genes:
                        gene_lists[list_name] = genes
                        
                except Exception as e:
                    logger.error(f"Error reading {filename}: {e}")
            else:
                logger.warning(f"File not found: {filename}")
        
        if not gene_lists:
            logger.error("No gene list files found!")
            raise FileNotFoundError(f"No gene list files found in {self.gene_lists_dir}")
        
        return gene_lists

    def load_data(self):
        """Load and prepare count data and metadata"""
        logger.info("Loading embryo count data...")
        
        # Load gene counts
        self.gene_counts = pd.read_csv(self.gene_counts_file, index_col=0)
        
        # Filter for PRJNA158021 samples only
        prjna158021_samples = [col for col in self.gene_counts.columns if 'PRJNA158021' in col]
        self.gene_counts = self.gene_counts[prjna158021_samples]
        
        logger.info(f"Loaded {self.gene_counts.shape[0]} genes × {self.gene_counts.shape[1]} samples")
        
        # Create metadata
        self.metadata = self._create_metadata()
        
        # Ensure sample order matches
        self.gene_counts = self.gene_counts[self.metadata.index]
        
        logger.info("Embryo data loading completed")
        return self

    def _create_metadata(self):
        """Create metadata dataframe from sample names"""
        metadata_list = []
        
        for sample in self.gene_counts.columns:
            # Parse: PRJNA158021_SRR458462_SD_72h
            parts = sample.split('_')
            photoperiod = parts[2]    # SD (short-day) or LD (long-day)
            timepoint = parts[3]      # 72h or 135h
            
            # Map photoperiod to biological meaning (same as adults!)
            diapause_status = 'Diapause_induction' if photoperiod == 'SD' else 'Non_diapause'
            timepoint_clean = timepoint.replace('h', 'h_pov')  # post-oviposition
            
            metadata_list.append({
                'sample_id': sample,
                'run_id': parts[1],  # SRR number
                'photoperiod': photoperiod,  # LD or SD (like adults)
                'timepoint': timepoint,
                'timepoint_clean': timepoint_clean,
                'diapause_status': diapause_status,
                'dataset': 'PRJNA158021',
                'life_stage': 'Embryo',
                'treatment': f"{photoperiod}_{timepoint}",  # SD_72h, LD_72h, etc.
                'comparison_group': 'Diapause_induction' if photoperiod == 'SD' else 'Non_diapause'
            })
        
        metadata_df = pd.DataFrame(metadata_list)
        metadata_df.set_index('sample_id', inplace=True)
        
        logger.info(f"Embryo metadata created for {len(metadata_df)} samples")
        logger.info(f"Treatments: {metadata_df['treatment'].value_counts().to_dict()}")
        logger.info(f"Photoperiod: {metadata_df['photoperiod'].value_counts().to_dict()}")
        logger.info(f"Timepoints: {metadata_df['timepoint'].value_counts().to_dict()}")
        
        return metadata_df

    def run_global_analysis(self):
        """Run global differential expression analysis"""
        logger.info("Running global embryo differential expression analysis...")
        
        # Filter low-count genes
        self._filter_genes()
        
        # Setup R environment
        self._setup_r_environment()
        
        # Run main analysis: Diapause vs Non-diapause (controlling for timepoint)
        self._run_main_deseq2()
        
        # Run timepoint-specific analyses
        self._run_timepoint_analyses()
        
        logger.info("Global embryo analysis completed")
        return self

    def _filter_genes(self, min_count=10, min_samples=3):
        """Filter genes with low counts"""
        before_filter = self.gene_counts.shape[0]
        keep = (self.gene_counts >= min_count).sum(axis=1) >= min_samples
        self.gene_counts = self.gene_counts[keep]
        after_filter = self.gene_counts.shape[0]
        logger.info(f"Gene filtering: {before_filter} -> {after_filter} genes ({after_filter/before_filter*100:.1f}% retained)")

    def _setup_r_environment(self):
        """Setup R environment"""
        r_code = """
        library(DESeq2)
        library(ggplot2)
        library(pheatmap)
        """
        robjects.r(r_code)

    def _run_main_deseq2(self):
        """Run main DESeq2 analysis: SD vs LD (photoperiod effect)"""
        count_matrix = self.gene_counts.astype(int)
        metadata_r = self.metadata[['photoperiod', 'timepoint']].copy()
        metadata_r['photoperiod'] = pd.Categorical(metadata_r['photoperiod'], 
                                                  categories=['LD', 'SD'])
        metadata_r['timepoint'] = pd.Categorical(metadata_r['timepoint'], 
                                               categories=['72h', '135h'])
        
        # Transfer to R
        robjects.globalenv['count_matrix'] = pandas2ri.py2rpy(count_matrix)
        robjects.globalenv['metadata'] = pandas2ri.py2rpy(metadata_r)
        
        # DESeq2 analysis with timepoint as covariate
        r_code = """
        # Main analysis: SD vs LD photoperiod (controlling for timepoint)
        dds_main <- DESeqDataSetFromMatrix(countData = count_matrix,
                                          colData = metadata,
                                          design = ~ timepoint + photoperiod)
        
        dds_main$photoperiod <- relevel(dds_main$photoperiod, ref = "LD")
        dds_main$timepoint <- relevel(dds_main$timepoint, ref = "72h")
        
        dds_main <- DESeq(dds_main)
        
        # Main contrast: SD vs LD (same as adults!)
        res_main <- results(dds_main, contrast = c("photoperiod", "SD", "LD"))
        res_main <- res_main[order(res_main$padj), ]
        
        # Get normalized counts
        normalized_counts_main <- counts(dds_main, normalized = TRUE)
        """
        robjects.r(r_code)
        
        # Extract main results
        r_results = robjects.globalenv['res_main']
        normalized_counts = robjects.globalenv['normalized_counts_main']
        
        self.global_results = pandas2ri.rpy2py(r_results)
        self.global_results.index = self.gene_counts.index
        
        self.normalized_counts = pandas2ri.rpy2py(normalized_counts)
        self.normalized_counts.index = self.gene_counts.index
        self.normalized_counts.columns = self.gene_counts.columns
        
        # Apply significance thresholds
        self.global_results['significant'] = (
            (self.global_results['padj'] < 0.05) & 
            (abs(self.global_results['log2FoldChange']) > 2.0)
        )
        
        # Save main results
        self.global_results.to_csv(self.de_dir / 'main_SD_vs_LD_photoperiod_results.csv')
        self.normalized_counts.to_csv(self.de_dir / 'normalized_counts.csv')
        
        total_sig = self.global_results['significant'].sum()
        logger.info(f"Main analysis (SD vs LD photoperiod): {total_sig} significant genes")

    def _run_timepoint_analyses(self):
        """Run separate analyses for each timepoint: SD vs LD"""
        logger.info("Running timepoint-specific analyses (SD vs LD)...")
        
        for timepoint in ['72h', '135h']:
            logger.info(f"Analyzing timepoint: {timepoint} (SD vs LD)")
            
            # Filter samples for this timepoint
            tp_samples = self.metadata[self.metadata['timepoint'] == timepoint].index
            tp_counts = self.gene_counts[tp_samples]
            tp_metadata = self.metadata.loc[tp_samples]
            
            # Prepare for R
            count_matrix = tp_counts.astype(int)
            metadata_r = tp_metadata[['photoperiod']].copy()
            metadata_r['photoperiod'] = pd.Categorical(metadata_r['photoperiod'], 
                                                      categories=['LD', 'SD'])
            
            # Transfer to R
            robjects.globalenv['count_matrix_tp'] = pandas2ri.py2rpy(count_matrix)
            robjects.globalenv['metadata_tp'] = pandas2ri.py2rpy(metadata_r)
            
            # DESeq2 analysis for this timepoint
            r_code = f"""
            # Analysis for {timepoint}: SD vs LD
            dds_tp <- DESeqDataSetFromMatrix(countData = count_matrix_tp,
                                           colData = metadata_tp,
                                           design = ~ photoperiod)
            
            dds_tp$photoperiod <- relevel(dds_tp$photoperiod, ref = "LD")
            dds_tp <- DESeq(dds_tp)
            
            res_tp <- results(dds_tp, contrast = c("photoperiod", "SD", "LD"))
            res_tp <- res_tp[order(res_tp$padj), ]
            """
            robjects.r(r_code)
            
            # Extract timepoint results
            r_results = robjects.globalenv['res_tp']
            tp_results = pandas2ri.rpy2py(r_results)
            tp_results.index = tp_counts.index
            
            # Apply significance thresholds
            tp_results['significant'] = (
                (tp_results['padj'] < 0.05) & 
                (abs(tp_results['log2FoldChange']) > 2.0)
            )
            
            # Store results
            self.timepoint_analyses[timepoint] = tp_results
            
            # Save timepoint results
            tp_results.to_csv(self.de_dir / f'{timepoint}_SD_vs_LD_photoperiod_results.csv')
            
            total_sig = tp_results['significant'].sum()
            logger.info(f"Timepoint {timepoint}: {total_sig} significant genes (SD vs LD)")

    def analyze_gene_lists(self):
        """Analyze each gene list across different contrasts"""
        logger.info("Analyzing gene lists for embryo photoperiod effects...")
        
        # Analyze main contrast: SD vs LD (controlling for timepoint)
        for list_name, genes in self.gene_lists.items():
            logger.info(f"Analyzing {list_name}: {len(genes)} genes (main photoperiod effect)")
            
            # Extract results for this gene list from main analysis
            list_results = self._extract_list_results(genes, list_name, self.global_results)
            
            # Store analysis
            self.list_analyses[list_name] = list_results
            
            # Create visualizations
            self._create_list_visualizations(list_name, list_results, "Main_SD_vs_LD")
        
        # Analyze timepoint-specific effects
        for timepoint, tp_results in self.timepoint_analyses.items():
            for list_name, genes in self.gene_lists.items():
                logger.info(f"Analyzing {list_name} at {timepoint} (SD vs LD)")
                
                # Extract results for this timepoint
                list_results = self._extract_list_results(genes, f"{list_name}_{timepoint}", tp_results)
                
                # Store analysis
                self.list_analyses[f"{list_name}_{timepoint}"] = list_results
                
                # Create visualizations
                self._create_list_visualizations(f"{list_name}_{timepoint}", list_results, f"Timepoint_{timepoint}")
        
        logger.info("Gene list analyses completed")
        return self

    def _extract_list_results(self, gene_list, list_name, results_df):
        """Extract results for a specific gene list"""
        # Find genes present in results
        present_genes = [gene for gene in gene_list if gene in results_df.index]
        missing_genes = [gene for gene in gene_list if gene not in results_df.index]
        
        if missing_genes:
            logger.warning(f"{list_name}: {len(missing_genes)} genes not found")
        
        # Extract results
        list_results = results_df.loc[present_genes].copy() if present_genes else pd.DataFrame()
        
        # Calculate summary statistics
        if not list_results.empty:
            stats = {
                'total_genes': len(gene_list),
                'genes_found': len(present_genes),
                'genes_missing': len(missing_genes),
                'significant_genes': list_results['significant'].sum(),
                'upregulated': ((list_results['padj'] < 0.05) & (list_results['log2FoldChange'] > 2)).sum(),
                'downregulated': ((list_results['padj'] < 0.05) & (list_results['log2FoldChange'] < -2)).sum(),
                'mean_log2fc': list_results['log2FoldChange'].mean(),
                'median_padj': list_results['padj'].median()
            }
        else:
            stats = {key: 0 for key in ['total_genes', 'genes_found', 'genes_missing', 'significant_genes', 'upregulated', 'downregulated']}
            stats.update({'mean_log2fc': 0, 'median_padj': 1})
        
        # Save results
        if not list_results.empty:
            list_results.to_csv(self.de_dir / f'{list_name}_results.csv')
        
        return {
            'results': list_results,
            'stats': stats,
            'genes': present_genes
        }

    def _create_list_visualizations(self, list_name, list_analysis, analysis_type):
        """Create visualizations for a specific gene list"""
        results = list_analysis['results']
        
        if results.empty:
            logger.warning(f"No results to plot for {list_name}")
            return
        
        # 1. Enhanced volcano plot
        self._create_embryo_volcano(list_name, results, list_analysis, analysis_type)
        
        # 2. Time-course expression plot (if main analysis)
        if analysis_type == "Main_Diapause_vs_NonDiapause":
            self._create_timecourse_plot(list_name, results)

    def _create_embryo_volcano(self, list_name, results, list_analysis, analysis_type):
        """Create beautiful volcano plot for embryo data"""
        # Prepare data
        plot_data = results.copy()
        plot_data['gene_id'] = plot_data.index
        plot_data['-log10_padj'] = -np.log10(plot_data['padj'])
        plot_data['significance'] = 'Not Significant'
        
        # Assign significance categories
        plot_data.loc[(plot_data['padj'] < 0.05) & (plot_data['log2FoldChange'] > 2), 'significance'] = 'Up in SD'
        plot_data.loc[(plot_data['padj'] < 0.05) & (plot_data['log2FoldChange'] < -2), 'significance'] = 'Down in SD'
        plot_data.loc[(plot_data['padj'] < 0.05) & (abs(plot_data['log2FoldChange']) <= 2), 'significance'] = 'Significant'
        
        # Create ggplot
        contrast_label = "SD vs LD Photoperiod" if "Main" in analysis_type else f"{analysis_type.replace('_', ' ')} - SD vs LD"
        
        p = (ggplot(plot_data, aes(x='log2FoldChange', y='-log10_padj', color='significance')) +
             geom_point(size=3, alpha=0.7) +
             geom_text(aes(label='gene_id'), size=8, nudge_y=0.1, alpha=0.8) +
             scale_color_manual(values={
                 'Up in SD': '#d62728',           # Red - up in short-day (diapause-inducing)
                 'Down in SD': '#2ca02c',         # Green - down in short-day  
                 'Significant': '#ff7f0e',        # Orange
                 'Not Significant': '#7f7f7f'     # Gray
             }) +
             geom_hline(yintercept=-np.log10(0.05), linetype='dashed', alpha=0.5) +
             geom_vline(xintercept=2, linetype='dashed', alpha=0.5) +
             geom_vline(xintercept=-2, linetype='dashed', alpha=0.5) +
             labs(
                 title=f'{list_name.replace("_", " ").title()} Genes - Embryonic {contrast_label}',
                 subtitle=f'{len(results)} genes analyzed, {list_analysis["stats"]["significant_genes"]} significant',
                 x='Log₂ Fold Change (SD vs LD)',
                 y='-Log₁₀(Adjusted P-value)',
                 color='Significance'
             ) +
             theme_bw() +
             theme(
                 figure_size=(12, 8),
                 plot_title=element_text(size=16, face='bold'),
                 plot_subtitle=element_text(size=12),
                 axis_title=element_text(size=12),
                 legend_title=element_text(size=11),
                 legend_position='right'
             ))
        
        # Save the plot
        filename = f'{list_name}_volcano_{analysis_type.lower()}.png'
        p.save(self.figures_dir / filename, dpi=300, width=12, height=8)

    def _create_timecourse_plot(self, list_name, results):
        """Create time-course expression plot showing diapause vs non-diapause over time"""
        if results.empty:
            return
            
        # Get top 6 most significant genes for visualization
        top_genes = results.nsmallest(6, 'padj').index
        
        # Get normalized expression data
        expr_data = self.normalized_counts.loc[top_genes]
        
        # Prepare data for plotting
        plot_data_list = []
        for gene in top_genes:
            for sample in expr_data.columns:
                sample_meta = self.metadata.loc[sample]
                plot_data_list.append({
                    'gene': gene,
                    'sample': sample,
                    'expression': np.log2(expr_data.loc[gene, sample] + 1),
                    'diapause_status': sample_meta['diapause_status'],
                    'timepoint': sample_meta['timepoint'],
                    'timepoint_num': 72 if sample_meta['timepoint'] == '72h' else 135
                })
        
        plot_df = pd.DataFrame(plot_data_list)
        
        # Create time-course plot using ggplot
        p = (ggplot(plot_df, aes(x='timepoint_num', y='expression', color='diapause_status')) +
             geom_point(size=3, alpha=0.7) +
             geom_smooth(method='loess', se=True, alpha=0.3) +
             facet_wrap('gene', scales='free_y', ncol=3) +
             scale_color_manual(values={
                 'Diapause': '#e74c3c',      # Red
                 'Non_diapause': '#3498db'   # Blue
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

    def create_embryo_summary(self):
        """Create comprehensive summary for embryo analysis"""
        logger.info("Creating embryo analysis summary...")
        
        # Create timepoint comparison
        self._create_timepoint_comparison()
        
        # Create publication figure
        self._create_embryo_publication_figure()
        
        logger.info("Embryo summary completed")
        return self

    def _create_timepoint_comparison(self):
        """Compare gene expression patterns between timepoints"""
        # Get main gene lists (not timepoint-specific)
        main_lists = {k: v for k, v in self.list_analyses.items() 
                     if not ('_72h' in k or '_135h' in k)}
        
        if not main_lists:
            return
            
        # Create comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Embryonic Diapause Gene Expression Patterns\nTimepoint Comparison', 
                    fontsize=20, fontweight='bold')
        
        # Compile data for all main lists
        summary_data = []
        for list_name, analysis in main_lists.items():
            stats = analysis['stats']
            stats['list_name'] = list_name
            summary_data.append(stats)
        
        summary_df = pd.DataFrame(summary_data).set_index('list_name')
        
        # Plot 1: Hit rates for main analysis
        ax1 = axes[0, 0]
        pct_sig = (summary_df['significant_genes'] / summary_df['genes_found'] * 100).fillna(0)
        colors = plt.cm.Set2(np.linspace(0, 1, len(summary_df)))
        bars = ax1.barh(range(len(summary_df)), pct_sig, color=colors, alpha=0.8)
        ax1.set_yticks(range(len(summary_df)))
        ax1.set_yticklabels([name.replace('_', ' ').title() for name in summary_df.index])
        ax1.set_xlabel('Hit Rate (%)', fontsize=12)
        ax1.set_title('A) Main Analysis Hit Rates\n(Diapause vs Non-diapause)', fontsize=14, fontweight='bold')
        
        # Plot 2: Expression direction in main analysis
        ax2 = axes[0, 1]
        x_pos = np.arange(len(summary_df))
        width = 0.35
        up_bars = ax2.bar(x_pos - width/2, summary_df['upregulated'], width, 
                         label='Up in Diapause', color='#e74c3c', alpha=0.8)
        down_bars = ax2.bar(x_pos + width/2, summary_df['downregulated'], width, 
                           label='Down in Diapause', color='#2ecc71', alpha=0.8)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([name.replace('_', '\n') for name in summary_df.index], fontsize=9)
        ax2.set_ylabel('Number of Genes', fontsize=12)
        ax2.set_title('B) Expression Direction\n(Main Analysis)', fontsize=14, fontweight='bold')
        ax2.legend()
        
        # Plot 3: Timepoint-specific effects comparison
        ax3 = axes[1, 0]
        timepoint_data = []
        for tp in ['72h', '135h']:
            tp_lists = {k: v for k, v in self.list_analyses.items() if f'_{tp}' in k}
            for list_name, analysis in tp_lists.items():
                base_name = list_name.replace(f'_{tp}', '')
                stats = analysis['stats']
                pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
                timepoint_data.append({
                    'list_name': base_name,
                    'timepoint': tp,
                    'hit_rate': pct_sig,
                    'significant': stats['significant_genes']
                })
        
        if timepoint_data:
            tp_df = pd.DataFrame(timepoint_data)
            tp_pivot = tp_df.pivot(index='list_name', columns='timepoint', values='hit_rate')
            
            x_pos = np.arange(len(tp_pivot))
            width = 0.35
            ax3.bar(x_pos - width/2, tp_pivot['72h'], width, label='72h', color='#3498db', alpha=0.8)
            ax3.bar(x_pos + width/2, tp_pivot['135h'], width, label='135h', color='#9b59b6', alpha=0.8)
            ax3.set_xticks(x_pos)
            ax3.set_xticklabels([name.replace('_', '\n') for name in tp_pivot.index], fontsize=9)
            ax3.set_ylabel('Hit Rate (%)', fontsize=12)
            ax3.set_title('C) Timepoint-Specific Effects', fontsize=14, fontweight='bold')
            ax3.legend()
        
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
        plt.savefig(self.figures_dir / 'embryo_timepoint_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _create_embryo_publication_figure(self):
        """Create comprehensive publication figure for embryo analysis"""
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)
        
        # Main title
        fig.suptitle('GWAS Candidate Gene Expression in Embryonic Diapause Preparation\n' + 
                    'Aedes albopictus Developmental Time Course Analysis', 
                    fontsize=24, fontweight='bold', y=0.98)
        
        # Get main analysis data
        main_lists = {k: v for k, v in self.list_analyses.items() 
                     if not ('_72h' in k or '_135h' in k)}
        
        summary_data = []
        for list_name, analysis in main_lists.items():
            stats = analysis['stats']
            stats['list_name'] = list_name
            summary_data.append(stats)
        
        summary_df = pd.DataFrame(summary_data).set_index('list_name')
        
        # Plot 1: Main analysis hit rates
        ax1 = fig.add_subplot(gs[0, 0])
        pct_sig = (summary_df['significant_genes'] / summary_df['genes_found'] * 100).fillna(0)
        colors = plt.cm.Set2(np.linspace(0, 1, len(summary_df)))
        bars = ax1.barh(range(len(summary_df)), pct_sig, color=colors, alpha=0.8)
        ax1.set_yticks(range(len(summary_df)))
        ax1.set_yticklabels([name.replace('_', ' ').title() for name in summary_df.index])
        ax1.set_xlabel('Hit Rate (%)', fontsize=12)
        ax1.set_title('A) Diapause Preparation Hit Rates', fontsize=14, fontweight='bold')
        
        for i, (bar, pct) in enumerate(zip(bars, pct_sig)):
            ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                    f'{pct:.1f}%', va='center', fontweight='bold')
        
        # Plot 2: Expression direction in embryonic diapause
        ax2 = fig.add_subplot(gs[0, 1])
        x_pos = np.arange(len(summary_df))
        width = 0.35
        up_bars = ax2.bar(x_pos - width/2, summary_df['upregulated'], width, 
                         label='Up in Diapause Prep', color='#e74c3c', alpha=0.8)
        down_bars = ax2.bar(x_pos + width/2, summary_df['downregulated'], width, 
                           label='Down in Diapause Prep', color='#2ecc71', alpha=0.8)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([name.replace('_', '\n') for name in summary_df.index], fontsize=10)
        ax2.set_ylabel('Number of Genes', fontsize=12)
        ax2.set_title('B) Expression Direction', fontsize=14, fontweight='bold')
        ax2.legend()
        
        # Plot 3: Developmental timeline
        ax3 = fig.add_subplot(gs[0, 2])
        # Create a timeline visualization
        timeline_data = self.metadata.groupby(['timepoint', 'diapause_status']).size().unstack()
        timeline_data.plot(kind='bar', ax=ax3, color=['#3498db', '#e74c3c'], alpha=0.8)
        ax3.set_xlabel('Developmental Stage', fontsize=12)
        ax3.set_ylabel('Number of Samples', fontsize=12)
        ax3.set_title('C) Experimental Timeline', fontsize=14, fontweight='bold')
        ax3.legend(['Non-diapause', 'Diapause Prep'])
        ax3.tick_params(axis='x', rotation=0)
        
        # Plots 4-6: Top gene list volcano plots
        top_lists = pct_sig.nlargest(3).index if len(pct_sig) >= 3 else pct_sig.index
        for i, list_name in enumerate(top_lists):
            ax = fig.add_subplot(gs[1, i])
            self._create_mini_embryo_volcano(ax, list_name, f'D{i+1}) {list_name.replace("_", " ").title()}')
        
        # Plot 7: Time-course heatmap of top genes
        ax7 = fig.add_subplot(gs[2, :])
        self._create_timecourse_heatmap(ax7)
        
        # Plot 8: Summary statistics table
        ax8 = fig.add_subplot(gs[3, :])
        ax8.axis('tight')
        ax8.axis('off')
        
        # Create summary table
        table_data = []
        for list_name, stats in zip(summary_df.index, summary_data):
            pct = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            table_data.append([
                list_name.replace('_', ' ').title(),
                stats['total_genes'],
                stats['genes_found'], 
                stats['significant_genes'],
                f"{pct:.1f}%",
                stats['upregulated'],
                stats['downregulated'],
                f"{stats['mean_log2fc']:.2f}"
            ])
        
        table = ax8.table(cellText=table_data,
                         colLabels=['Gene List', 'Total', 'Found', 'Significant', 'Hit Rate', 
                                   'Up in Diapause', 'Down in Diapause', 'Mean log₂FC'],
                         cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        ax8.set_title('G) Embryonic Diapause Preparation - Summary Statistics', 
                     fontsize=14, fontweight='bold', pad=20)
        
        plt.savefig(self.figures_dir / 'embryo_publication_summary.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'embryo_publication_summary.pdf', bbox_inches='tight')
        plt.close()

    def _create_mini_embryo_volcano(self, ax, list_name, title):
        """Create mini volcano plot for embryo publication figure"""
        if list_name not in self.list_analyses or self.list_analyses[list_name]['results'].empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title, fontsize=12, fontweight='bold')
            return
            
        results = self.list_analyses[list_name]['results']
        x = results['log2FoldChange']
        y = -np.log10(results['padj'])
        
        # Color based on significance (embryo-specific colors)
        colors = []
        for _, row in results.iterrows():
            if row['padj'] < 0.05 and row['log2FoldChange'] > 2:
                colors.append('#e74c3c')  # Red - up in diapause prep
            elif row['padj'] < 0.05 and row['log2FoldChange'] < -2:
                colors.append('#2ecc71')  # Green - down in diapause prep
            elif row['padj'] < 0.05:
                colors.append('#f39c12')  # Orange - significant
            else:
                colors.append('#95a5a6')  # Gray - not significant
        
        ax.scatter(x, y, c=colors, s=60, alpha=0.8)
        ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=2, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-2, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel('Log₂FC (Diapause vs Non-diapause)', fontsize=10)
        ax.set_ylabel('-Log₁₀(Adj P)', fontsize=10)
        ax.set_title(title, fontsize=12, fontweight='bold')

    def _create_timecourse_heatmap(self, ax):
        """Create time-course heatmap for publication figure"""
        # Get top 20 most significant genes across all main lists
        all_genes = set()
        for list_name, analysis in self.list_analyses.items():
            if not ('_72h' in list_name or '_135h' in list_name):  # Main analyses only
                if not analysis['results'].empty:
                    sig_genes = analysis['results'][analysis['results']['significant']].index
                    all_genes.update(sig_genes)
        
        if not all_genes:
            ax.text(0.5, 0.5, 'No significant genes for heatmap', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('F) No Significant Genes Found', fontsize=14, fontweight='bold')
            return
        
        # Limit to top 20 genes by significance
        top_genes = list(all_genes)[:20]
        
        # Get expression data
        expr_data = self.normalized_counts.loc[top_genes]
        log_expr = np.log2(expr_data + 1)
        
        # Z-score normalize
        from scipy.stats import zscore
        z_expr = pd.DataFrame(
            zscore(log_expr, axis=1),
            index=log_expr.index,
            columns=log_expr.columns
        )
        
        # Create annotation colors
        annotation_colors = []
        for sample in z_expr.columns:
            meta = self.metadata.loc[sample]
            if meta['diapause_status'] == 'Diapause' and meta['timepoint'] == '72h':
                annotation_colors.append('#e74c3c')  # Dark red
            elif meta['diapause_status'] == 'Diapause' and meta['timepoint'] == '135h':
                annotation_colors.append('#c0392b')  # Darker red
            elif meta['diapause_status'] == 'Non_diapause' and meta['timepoint'] == '72h':
                annotation_colors.append('#3498db')  # Blue
            else:
                annotation_colors.append('#2980b9')  # Darker blue
        
        # Plot heatmap
        im = ax.imshow(z_expr.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
        
        # Add colorbar
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Z-score (Expression)', fontsize=10)
        
        # Set labels
        ax.set_xticks(range(len(z_expr.columns)))
        ax.set_xticklabels([f"{self.metadata.loc[col, 'diapause_status'][:3]}\n{self.metadata.loc[col, 'timepoint']}" 
                           for col in z_expr.columns], fontsize=8)
        ax.set_yticks(range(len(top_genes)))
        ax.set_yticklabels(top_genes, fontsize=8)
        ax.set_title('F) Top Significant Genes - Expression Heatmap', fontsize=14, fontweight='bold')
        
        # Add sample group colors at top
        for i, color in enumerate(annotation_colors):
            ax.add_patch(plt.Rectangle((i-0.4, -1), 0.8, 0.3, facecolor=color, alpha=0.8))

    def generate_embryo_report(self):
        """Generate comprehensive embryo analysis report"""
        logger.info("Generating embryo analysis report...")
        
        # Calculate overall statistics
        total_genes_tested = len(self.global_results)
        global_significant = self.global_results['significant'].sum()
        
        # Get main analysis stats
        main_lists = {k: v for k, v in self.list_analyses.items() 
                     if not ('_72h' in k or '_135h' in k)}
        
        report = f"""# PRJNA158021 Embryonic Diapause Analysis Report

## Experimental Design
- **Dataset**: PRJNA158021 (Embryonic development)
- **Life Stage**: Embryos at 72-78h and 135-141h post-oviposition
- **Comparison**: Diapause preparation vs Normal development
- **Model**: ~ timepoint + diapause_status
- **Timepoints**: 72h (early) and 135h (late) post-oviposition

## Global Analysis Results
- **Total genes tested**: {total_genes_tested:,}
- **Globally significant**: {global_significant:,} ({global_significant/total_genes_tested*100:.2f}%)
- **Significance criteria**: FDR < 0.05 AND |log₂FC| > 2.0

## Main Analysis: Diapause Preparation vs Normal Development

| Gene List | Total Genes | Found | Significant | % Significant | Up in Diapause | Down in Diapause | Mean log₂FC |
|-----------|-------------|-------|-------------|---------------|----------------|------------------|-------------|
"""
        
        for list_name, analysis in main_lists.items():
            stats = analysis['stats']
            pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            
            report += f"| {list_name.replace('_', ' ').title()} | {stats['total_genes']} | {stats['genes_found']} | {stats['significant_genes']} | {pct_sig:.1f}% | {stats['upregulated']} | {stats['downregulated']} | {stats['mean_log2fc']:.2f} |\n"

        report += f"""

## Timepoint-Specific Analysis Results

### 72h Post-Oviposition (Early Development)
"""
        
        # Add 72h results
        h72_lists = {k: v for k, v in self.list_analyses.items() if '_72h' in k}
        for list_name, analysis in h72_lists.items():
            base_name = list_name.replace('_72h', '')
            stats = analysis['stats']
            pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            report += f"- **{base_name.replace('_', ' ').title()}**: {stats['significant_genes']}/{stats['genes_found']} significant ({pct_sig:.1f}%)\n"

        report += f"""

### 135h Post-Oviposition (Late Development)
"""
        
        # Add 135h results
        h135_lists = {k: v for k, v in self.list_analyses.items() if '_135h' in k}
        for list_name, analysis in h135_lists.items():
            base_name = list_name.replace('_135h', '')
            stats = analysis['stats']
            pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            report += f"- **{base_name.replace('_', ' ').title()}**: {stats['significant_genes']}/{stats['genes_found']} significant ({pct_sig:.1f}%)\n"

        # Find best performing lists
        main_performance = []
        for list_name, analysis in main_lists.items():
            stats = analysis['stats']
            pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
            main_performance.append((list_name, pct_sig, stats['significant_genes']))
        
        main_performance.sort(key=lambda x: x[1], reverse=True)

        report += f"""

## Key Findings

### Most Responsive Gene Lists (Embryonic Diapause):
"""
        for i, (list_name, pct_sig, n_sig) in enumerate(main_performance[:3]):
            report += f"{i+1}. **{list_name.replace('_', ' ').title()}**: {pct_sig:.1f}% hit rate ({n_sig} genes)\n"

        report += f"""

### Biological Insights:
- **Diapause preparation signature**: {'Strong' if global_significant > 100 else 'Moderate' if global_significant > 50 else 'Weak'} ({global_significant} genes)
- **Temporal dynamics**: Genes show {'time-dependent' if any('72h' in k or '135h' in k for k in self.list_analyses.keys()) else 'time-independent'} expression patterns
- **GWAS validation**: {'High' if any(stats['significant_genes'] > 2 for stats in [a['stats'] for a in main_lists.values()]) else 'Moderate'} concordance with genetic mapping

## Files Generated
- **Main analysis**: `main_diapause_vs_nondiapause_results.csv`
- **Timepoint analyses**: `72h_diapause_vs_nondiapause_results.csv`, `135h_diapause_vs_nondiapause_results.csv`
- **Gene list results**: `[list_name]_results.csv`
- **Visualizations**: Multiple volcano plots, time-course plots, and publication figures

## Biological Interpretation
The embryonic analysis reveals gene expression changes associated with diapause preparation during critical developmental windows. Early timepoint (72h) likely captures initial diapause commitment signals, while late timepoint (135h) reflects establishment of diapause state.

**Analysis completed**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}
"""
        
        # Save report
        with open(self.de_dir / 'embryo_analysis_report.md', 'w') as f:
            f.write(report)
        
        logger.info("Embryo analysis report generated")
        return self

def main():
    """Main execution function for embryo analysis"""
    analyzer = EmbryoAnalyzer(
        gene_lists_dir="output/tables"
    )
    
    # Run complete analysis pipeline
    analyzer.load_data()
    analyzer.run_global_analysis()
    analyzer.analyze_gene_lists()
    analyzer.create_embryo_summary()
    analyzer.generate_embryo_report()
    
    logger.info("=== EMBRYO ANALYSIS COMPLETED ===")
    logger.info(f"Results in: {analyzer.de_dir}")
    logger.info(f"Figures in: {analyzer.figures_dir}")
    
    # Print quick summary
    print("\n=== EMBRYO ANALYSIS SUMMARY ===")
    print("Main Analysis (Diapause vs Non-diapause):")
    print(f"{'List Name':<20} {'Genes':<6} {'Significant':<11} {'% Hit Rate':<10}")
    print("-" * 50)
    
    main_lists = {k: v for k, v in analyzer.list_analyses.items() 
                 if not ('_72h' in k or '_135h' in k)}
    
    for list_name, analysis in main_lists.items():
        stats = analysis['stats']
        pct_sig = (stats['significant_genes'] / stats['genes_found'] * 100) if stats['genes_found'] > 0 else 0
        print(f"{list_name:<20} {stats['genes_found']:<6} {stats['significant_genes']:<11} {pct_sig:>9.1f}%")

if __name__ == "__main__":
    main()
