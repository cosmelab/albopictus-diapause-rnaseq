#!/usr/bin/env python3
"""
Outlier Detection and Handling for RNA-seq Analysis (PRJNA158021)

This script provides comprehensive outlier detection and handling methods for RNA-seq data,
including sample-level and gene-level outlier detection, visualization, and handling strategies.

Key Features:
- Sample-level outlier detection (PCA, hierarchical clustering, sample distances)
- Gene-level outlier detection (expression patterns, Cook's distance)
- Quality control metrics (library size, gene detection, GC content)
- Multiple handling strategies (removal, winsorization, robust methods)
- Interactive visualizations for outlier assessment

Usage:
    python outlier_detection_PRJNA158021.py

Requirements:
    - pandas, numpy, matplotlib, seaborn
    - scikit-learn for clustering and PCA
    - scipy for statistical tests
    - plotly for interactive plots
    - rpy2 for R integration (DESeq2 Cook's distance)
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
import warnings
warnings.filterwarnings('ignore')

# Statistical and ML libraries
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.metrics import pairwise_distances
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, iqr

# Interactive plotting
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("plotly not available - interactive plots will be skipped")

# R integration for DESeq2 Cook's distance
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    RPY2_AVAILABLE = True
except ImportError:
    RPY2_AVAILABLE = False
    print("rpy2 not available - Cook's distance analysis will be skipped")

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('outlier_detection_PRJNA158021.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class RNAseqOutlierDetector:
    def __init__(self, gene_counts_file="output/tables/PRJNA158021_gene_counts.csv", 
                 output_base="output"):
        """
        Initialize the RNA-seq outlier detector.
        
        Args:
            gene_counts_file: Path to the gene counts file
            output_base: Base directory for output files
        """
        self.gene_counts_file = Path(gene_counts_file)
        self.output_base = Path(output_base)
        
        # Create output directories
        self.outlier_dir = self.output_base / "outlier_analysis" / "PRJNA158021"
        self.figures_dir = self.output_base / "figures" / "PRJNA158021"
        self.qc_dir = self.output_base / "qc" / "PRJNA158021"
        
        for dir_path in [self.outlier_dir, self.figures_dir, self.qc_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.gene_counts = None
        self.metadata = None
        self.normalized_counts = None
        
        # Outlier detection results
        self.sample_outliers = {}
        self.gene_outliers = {}
        self.quality_metrics = {}
        
        logger.info("RNAseqOutlierDetector initialized for PRJNA158021")

    def load_data(self):
        """Load gene counts and create metadata."""
        if not self.gene_counts_file.exists():
            raise FileNotFoundError(f"Gene counts file not found: {self.gene_counts_file}")
        
        # Load gene counts
        self.gene_counts = pd.read_csv(self.gene_counts_file, index_col=0)
        logger.info(f"Loaded gene counts: {self.gene_counts.shape[0]} genes, {self.gene_counts.shape[1]} samples")
        
        # Create metadata
        self._create_metadata()
        
        # Calculate basic quality metrics
        self._calculate_quality_metrics()

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
                
                treatment = f"{light_cycle}_{time_point}"
                
                metadata_list.append({
                    'sample': sample,
                    'prjna': prjna,
                    'srr': srr,
                    'light_cycle': light_cycle,
                    'time_point': time_point,
                    'treatment': treatment
                })
        
        self.metadata = pd.DataFrame(metadata_list)
        self.metadata.set_index('sample', inplace=True)
        logger.info(f"Created metadata for {len(self.metadata)} samples")

    def _calculate_quality_metrics(self):
        """Calculate basic quality metrics for each sample."""
        logger.info("Calculating quality metrics...")
        
        # Library size (total counts per sample)
        library_sizes = self.gene_counts.sum()
        
        # Number of detected genes (genes with > 0 counts)
        detected_genes = (self.gene_counts > 0).sum()
        
        # Percentage of detected genes
        gene_detection_rate = (detected_genes / len(self.gene_counts)) * 100
        
        # Create quality metrics DataFrame
        self.quality_metrics = pd.DataFrame({
            'library_size': library_sizes,
            'detected_genes': detected_genes,
            'gene_detection_rate': gene_detection_rate,
            'light_cycle': self.metadata['light_cycle'],
            'time_point': self.metadata['time_point'],
            'treatment': self.metadata['treatment']
        })
        
        logger.info("Quality metrics calculated")

    def detect_sample_outliers(self, methods=['pca', 'clustering', 'distance', 'quality']):
        """
        Detect sample-level outliers using multiple methods.
        
        Args:
            methods: List of methods to use for outlier detection
        """
        logger.info("Starting sample-level outlier detection...")
        
        # Filter genes for analysis (remove low-count genes)
        min_count = 10
        min_samples = 3
        keep_genes = (self.gene_counts >= min_count).sum(axis=1) >= min_samples
        filtered_counts = self.gene_counts[keep_genes]
        
        logger.info(f"Using {filtered_counts.shape[0]} genes for outlier detection")
        
        # 1. PCA-based outlier detection
        if 'pca' in methods:
            self._detect_pca_outliers(filtered_counts)
        
        # 2. Hierarchical clustering-based outlier detection
        if 'clustering' in methods:
            self._detect_clustering_outliers(filtered_counts)
        
        # 3. Sample distance-based outlier detection
        if 'distance' in methods:
            self._detect_distance_outliers(filtered_counts)
        
        # 4. Quality metrics-based outlier detection
        if 'quality' in methods:
            self._detect_quality_outliers()
        
        # Combine results
        self._combine_sample_outlier_results()
        
        logger.info("Sample-level outlier detection completed")

    def _detect_pca_outliers(self, counts, n_components=2, threshold=2.0):
        """Detect outliers using PCA."""
        logger.info("Running PCA-based outlier detection...")
        
        # Log-transform and center the data
        log_counts = np.log2(counts + 1)
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(log_counts.T)
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(scaled_data)
        
        # Calculate Mahalanobis distance for outlier detection
        from scipy.spatial.distance import mahalanobis
        from scipy.linalg import inv
        
        # Calculate covariance matrix
        cov_matrix = np.cov(pca_result.T)
        inv_cov_matrix = inv(cov_matrix)
        
        # Calculate Mahalanobis distance for each sample
        mahal_distances = []
        for i in range(len(pca_result)):
            mahal_dist = mahalanobis(pca_result[i], np.mean(pca_result, axis=0), inv_cov_matrix)
            mahal_distances.append(mahal_dist)
        
        # Identify outliers (samples with high Mahalanobis distance)
        mahal_threshold = np.percentile(mahal_distances, 95)  # 95th percentile
        pca_outliers = [self.gene_counts.columns[i] for i, dist in enumerate(mahal_distances) if dist > mahal_threshold]
        
        self.sample_outliers['pca'] = {
            'outliers': pca_outliers,
            'distances': dict(zip(self.gene_counts.columns, mahal_distances)),
            'threshold': mahal_threshold,
            'pca_result': pca_result,
            'explained_variance': pca.explained_variance_ratio_
        }
        
        logger.info(f"PCA detected {len(pca_outliers)} potential outliers")

    def _detect_clustering_outliers(self, counts, n_clusters=4, distance_threshold=0.8):
        """Detect outliers using hierarchical clustering."""
        logger.info("Running clustering-based outlier detection...")
        
        # Log-transform the data
        log_counts = np.log2(counts + 1)
        
        # Calculate correlation distance matrix
        corr_matrix = log_counts.corr()
        distance_matrix = 1 - np.abs(corr_matrix)
        
        # Perform hierarchical clustering
        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage='ward',
            metric='precomputed'
        )
        
        # Convert distance matrix to condensed form
        condensed_distances = squareform(distance_matrix.values)
        cluster_labels = clustering.fit_predict(distance_matrix)
        
        # Identify samples that are far from their cluster center
        cluster_outliers = []
        for cluster_id in range(n_clusters):
            cluster_samples = [self.gene_counts.columns[i] for i, label in enumerate(cluster_labels) if label == cluster_id]
            if len(cluster_samples) > 1:
                cluster_distances = distance_matrix.loc[cluster_samples, cluster_samples]
                avg_distances = cluster_distances.mean()
                
                # Find samples with high average distance to cluster members
                for sample in cluster_samples:
                    if avg_distances[sample] > distance_threshold:
                        cluster_outliers.append(sample)
        
        self.sample_outliers['clustering'] = {
            'outliers': cluster_outliers,
            'cluster_labels': dict(zip(self.gene_counts.columns, cluster_labels)),
            'distance_matrix': distance_matrix
        }
        
        logger.info(f"Clustering detected {len(cluster_outliers)} potential outliers")

    def _detect_distance_outliers(self, counts, threshold_percentile=95):
        """Detect outliers based on sample-to-sample distances."""
        logger.info("Running distance-based outlier detection...")
        
        # Log-transform the data
        log_counts = np.log2(counts + 1)
        
        # Calculate Euclidean distance matrix
        distance_matrix = pairwise_distances(log_counts.T, metric='euclidean')
        distance_df = pd.DataFrame(distance_matrix, 
                                 index=self.gene_counts.columns, 
                                 columns=self.gene_counts.columns)
        
        # Calculate average distance for each sample
        avg_distances = distance_df.mean()
        
        # Identify outliers (samples with high average distance)
        distance_threshold = np.percentile(avg_distances, threshold_percentile)
        distance_outliers = avg_distances[avg_distances > distance_threshold].index.tolist()
        
        self.sample_outliers['distance'] = {
            'outliers': distance_outliers,
            'avg_distances': avg_distances.to_dict(),
            'threshold': distance_threshold,
            'distance_matrix': distance_df
        }
        
        logger.info(f"Distance-based detection found {len(distance_outliers)} potential outliers")

    def _detect_quality_outliers(self, library_size_threshold=0.1, detection_rate_threshold=0.1):
        """Detect outliers based on quality metrics."""
        logger.info("Running quality-based outlier detection...")
        
        quality_outliers = []
        
        # Library size outliers (samples with very low or very high library size)
        library_sizes = self.quality_metrics['library_size']
        lib_size_lower = np.percentile(library_sizes, library_size_threshold * 100)
        lib_size_upper = np.percentile(library_sizes, (1 - library_size_threshold) * 100)
        
        lib_size_outliers = library_sizes[
            (library_sizes < lib_size_lower) | (library_sizes > lib_size_upper)
        ].index.tolist()
        
        # Gene detection rate outliers
        detection_rates = self.quality_metrics['gene_detection_rate']
        det_rate_lower = np.percentile(detection_rates, detection_rate_threshold * 100)
        det_rate_upper = np.percentile(detection_rates, (1 - detection_rate_threshold) * 100)
        
        det_rate_outliers = detection_rates[
            (detection_rates < det_rate_lower) | (detection_rates > det_rate_upper)
        ].index.tolist()
        
        quality_outliers = list(set(lib_size_outliers + det_rate_outliers))
        
        self.sample_outliers['quality'] = {
            'outliers': quality_outliers,
            'library_size_outliers': lib_size_outliers,
            'detection_rate_outliers': det_rate_outliers,
            'library_size_thresholds': (lib_size_lower, lib_size_upper),
            'detection_rate_thresholds': (det_rate_lower, det_rate_upper)
        }
        
        logger.info(f"Quality-based detection found {len(quality_outliers)} potential outliers")

    def _combine_sample_outlier_results(self):
        """Combine results from different outlier detection methods."""
        logger.info("Combining sample outlier results...")
        
        # Collect all outliers
        all_outliers = set()
        method_counts = {}
        
        for method, results in self.sample_outliers.items():
            outliers = results['outliers']
            all_outliers.update(outliers)
            method_counts[method] = len(outliers)
        
        # Calculate how many methods flagged each sample as outlier
        outlier_scores = {}
        for sample in all_outliers:
            score = sum(1 for method, results in self.sample_outliers.items() 
                       if sample in results['outliers'])
            outlier_scores[sample] = score
        
        # Classify samples based on consensus
        consensus_outliers = [sample for sample, score in outlier_scores.items() if score >= 2]
        potential_outliers = [sample for sample, score in outlier_scores.items() if score == 1]
        
        self.sample_outliers['combined'] = {
            'consensus_outliers': consensus_outliers,
            'potential_outliers': potential_outliers,
            'outlier_scores': outlier_scores,
            'method_counts': method_counts
        }
        
        logger.info(f"Combined results: {len(consensus_outliers)} consensus outliers, "
                   f"{len(potential_outliers)} potential outliers")

    def detect_gene_outliers(self, methods=['expression', 'variance', 'cooks_distance']):
        """
        Detect gene-level outliers using multiple methods.
        
        Args:
            methods: List of methods to use for gene outlier detection
        """
        logger.info("Starting gene-level outlier detection...")
        
        # Filter genes for analysis
        min_count = 10
        min_samples = 3
        keep_genes = (self.gene_counts >= min_count).sum(axis=1) >= min_samples
        filtered_counts = self.gene_counts[keep_genes]
        
        # 1. Expression pattern outliers
        if 'expression' in methods:
            self._detect_expression_outliers(filtered_counts)
        
        # 2. Variance outliers
        if 'variance' in methods:
            self._detect_variance_outliers(filtered_counts)
        
        # 3. Cook's distance outliers (if R is available)
        if 'cooks_distance' in methods and RPY2_AVAILABLE:
            self._detect_cooks_distance_outliers(filtered_counts)
        
        logger.info("Gene-level outlier detection completed")

    def _detect_expression_outliers(self, counts, threshold=3.0):
        """Detect genes with unusual expression patterns."""
        logger.info("Running expression pattern outlier detection...")
        
        # Log-transform the data
        log_counts = np.log2(counts + 1)
        
        # Calculate z-scores for each gene across samples
        z_scores = log_counts.apply(lambda x: zscore(x), axis=1)
        
        # Find genes with extreme z-scores
        extreme_genes = []
        for gene in log_counts.index:
            gene_z_scores = z_scores.loc[gene]
            if any(abs(score) > threshold for score in gene_z_scores):
                extreme_genes.append(gene)
        
        self.gene_outliers['expression'] = {
            'outliers': extreme_genes,
            'z_scores': z_scores,
            'threshold': threshold
        }
        
        logger.info(f"Expression pattern detection found {len(extreme_genes)} potential outliers")

    def _detect_variance_outliers(self, counts, threshold_percentile=95):
        """Detect genes with unusually high or low variance."""
        logger.info("Running variance-based outlier detection...")
        
        # Log-transform the data
        log_counts = np.log2(counts + 1)
        
        # Calculate variance for each gene
        gene_variances = log_counts.var(axis=1)
        
        # Identify high variance outliers
        high_var_threshold = np.percentile(gene_variances, threshold_percentile)
        high_var_outliers = gene_variances[gene_variances > high_var_threshold].index.tolist()
        
        # Identify low variance outliers (genes with very consistent expression)
        low_var_threshold = np.percentile(gene_variances, 100 - threshold_percentile)
        low_var_outliers = gene_variances[gene_variances < low_var_threshold].index.tolist()
        
        variance_outliers = high_var_outliers + low_var_outliers
        
        self.gene_outliers['variance'] = {
            'outliers': variance_outliers,
            'high_variance_outliers': high_var_outliers,
            'low_variance_outliers': low_var_outliers,
            'variances': gene_variances.to_dict(),
            'thresholds': (low_var_threshold, high_var_threshold)
        }
        
        logger.info(f"Variance-based detection found {len(variance_outliers)} potential outliers")

    def _detect_cooks_distance_outliers(self, counts):
        """Detect outliers using Cook's distance from DESeq2."""
        logger.info("Running Cook's distance outlier detection...")
        
        try:
            # Set up R environment
            pandas2ri.activate()
            
            # Import R packages
            base = importr('base')
            deseq2 = importr('DESeq2')
            
            # Prepare data for R
            count_matrix = counts.astype(int)
            metadata_r = self.metadata[['light_cycle', 'time_point']].copy()
            
            # Convert to R objects
            r_counts = pandas2ri.py2rpy(count_matrix)
            r_metadata = pandas2ri.py2rpy(metadata_r)
            
            # Assign to R global environment
            ro.globalenv['r_counts'] = r_counts
            ro.globalenv['r_metadata'] = r_metadata
            
            # Run DESeq2 and calculate Cook's distance
            r_code = """
            # Create DESeq2 dataset
            dds <- DESeqDataSetFromMatrix(
                countData = r_counts,
                colData = r_metadata,
                design = ~ time_point + light_cycle
            )
            
            # Set reference levels
            dds$light_cycle <- relevel(dds$light_cycle, ref = "LD")
            dds$time_point <- relevel(dds$time_point, ref = "72h")
            
            # Run DESeq2
            dds <- DESeq(dds)
            
            # Calculate Cook's distance
            cooks_dist <- assays(dds)[["cooks"]]
            
            # Store results
            .GlobalEnv$cooks_dist <- cooks_dist
            .GlobalEnv$dds <- dds
            """
            
            ro.r(r_code)
            
            # Extract Cook's distance
            cooks_dist_r = ro.r['cooks_dist']
            cooks_dist_py = pandas2ri.rpy2py(cooks_dist_r)
            
            # Identify genes with high Cook's distance
            cooks_threshold = 4.0 / len(counts.columns)  # Standard threshold
            high_cooks_genes = []
            
            for gene in counts.index:
                if gene in cooks_dist_py.index:
                    gene_cooks = cooks_dist_py.loc[gene]
                    if any(cook_dist > cooks_threshold for cook_dist in gene_cooks):
                        high_cooks_genes.append(gene)
            
            self.gene_outliers['cooks_distance'] = {
                'outliers': high_cooks_genes,
                'cooks_distances': cooks_dist_py,
                'threshold': cooks_threshold
            }
            
            logger.info(f"Cook's distance detection found {len(high_cooks_genes)} potential outliers")
            
        except Exception as e:
            logger.error(f"Error in Cook's distance analysis: {e}")
            self.gene_outliers['cooks_distance'] = {
                'outliers': [],
                'error': str(e)
            }

    def create_outlier_visualizations(self):
        """Create comprehensive visualizations for outlier detection results."""
        logger.info("Creating outlier detection visualizations...")
        
        # 1. Quality metrics plots
        self._create_quality_plots()
        
        # 2. Sample outlier plots
        self._create_sample_outlier_plots()
        
        # 3. Gene outlier plots
        self._create_gene_outlier_plots()
        
        # 4. Summary plots
        self._create_summary_plots()
        
        logger.info("Outlier detection visualizations completed")

    def _create_quality_plots(self):
        """Create quality metrics plots."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Sample Quality Metrics', fontsize=16, fontweight='bold')
        
        # Library size distribution
        axes[0, 0].hist(self.quality_metrics['library_size'], bins=20, alpha=0.7, color='skyblue')
        axes[0, 0].set_xlabel('Library Size')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Library Size Distribution')
        
        # Gene detection rate
        axes[0, 1].hist(self.quality_metrics['gene_detection_rate'], bins=20, alpha=0.7, color='lightgreen')
        axes[0, 1].set_xlabel('Gene Detection Rate (%)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Gene Detection Rate Distribution')
        
        # Library size by treatment
        sns.boxplot(data=self.quality_metrics, x='treatment', y='library_size', ax=axes[1, 0])
        axes[1, 0].set_title('Library Size by Treatment')
        axes[1, 0].tick_params(axis='x', rotation=45)
        
        # Gene detection rate by treatment
        sns.boxplot(data=self.quality_metrics, x='treatment', y='gene_detection_rate', ax=axes[1, 1])
        axes[1, 1].set_title('Gene Detection Rate by Treatment')
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'quality_metrics.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _create_sample_outlier_plots(self):
        """Create sample outlier detection plots."""
        if 'pca' in self.sample_outliers:
            # PCA plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            pca_result = self.sample_outliers['pca']['pca_result']
            distances = self.sample_outliers['pca']['distances']
            
            # Color samples by outlier status
            colors = []
            for sample in self.gene_counts.columns:
                if sample in self.sample_outliers['pca']['outliers']:
                    colors.append('red')
                else:
                    colors.append('blue')
            
            scatter = ax.scatter(pca_result[:, 0], pca_result[:, 1], c=colors, alpha=0.7, s=100)
            ax.set_xlabel(f'PC1 ({self.sample_outliers["pca"]["explained_variance"][0]:.1%} variance)')
            ax.set_ylabel(f'PC2 ({self.sample_outliers["pca"]["explained_variance"][1]:.1%} variance)')
            ax.set_title('PCA-based Sample Outlier Detection')
            
            # Add sample labels
            for i, sample in enumerate(self.gene_counts.columns):
                ax.annotate(sample.split('_')[-2:][0], (pca_result[i, 0], pca_result[i, 1]), 
                           xytext=(5, 5), textcoords='offset points', fontsize=8)
            
            # Add legend
            from matplotlib.lines import Line2D
            legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                                    markersize=10, label='Outlier'),
                              Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
                                    markersize=10, label='Normal')]
            ax.legend(handles=legend_elements)
            
            plt.tight_layout()
            plt.savefig(self.figures_dir / 'pca_outliers.png', dpi=300, bbox_inches='tight')
            plt.close()

    def _create_gene_outlier_plots(self):
        """Create gene outlier detection plots."""
        if 'variance' in self.gene_outliers:
            # Variance distribution plot
            fig, ax = plt.subplots(figsize=(10, 6))
            
            variances = list(self.gene_outliers['variance']['variances'].values())
            thresholds = self.gene_outliers['variance']['thresholds']
            
            ax.hist(variances, bins=50, alpha=0.7, color='lightblue')
            ax.axvline(thresholds[0], color='red', linestyle='--', label='Low variance threshold')
            ax.axvline(thresholds[1], color='red', linestyle='--', label='High variance threshold')
            ax.set_xlabel('Gene Variance (log2 counts)')
            ax.set_ylabel('Frequency')
            ax.set_title('Gene Variance Distribution')
            ax.legend()
            
            plt.tight_layout()
            plt.savefig(self.figures_dir / 'gene_variance_outliers.png', dpi=300, bbox_inches='tight')
            plt.close()

    def _create_summary_plots(self):
        """Create summary plots of outlier detection results."""
        # Method comparison plot
        if 'combined' in self.sample_outliers:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            method_counts = self.sample_outliers['combined']['method_counts']
            methods = list(method_counts.keys())
            counts = list(method_counts.values())
            
            bars = ax.bar(methods, counts, color=['skyblue', 'lightgreen', 'lightcoral', 'gold'])
            ax.set_xlabel('Detection Method')
            ax.set_ylabel('Number of Outliers Detected')
            ax.set_title('Sample Outliers Detected by Method')
            
            # Add value labels on bars
            for bar, count in zip(bars, counts):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                       str(count), ha='center', va='bottom')
            
            plt.tight_layout()
            plt.savefig(self.figures_dir / 'outlier_method_comparison.png', dpi=300, bbox_inches='tight')
            plt.close()

    def generate_outlier_report(self):
        """Generate a comprehensive outlier detection report."""
        logger.info("Generating outlier detection report...")
        
        report_path = self.outlier_dir / 'outlier_detection_report.md'
        
        with open(report_path, 'w') as f:
            f.write("# RNA-seq Outlier Detection Report\n\n")
            f.write(f"**Dataset:** PRJNA158021 (Embryo Dataset)\n")
            f.write(f"**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Summary statistics
            f.write("## Summary Statistics\n\n")
            f.write(f"- **Total samples:** {len(self.metadata)}\n")
            f.write(f"- **Total genes:** {len(self.gene_counts)}\n")
            f.write(f"- **Treatment groups:** {', '.join(self.metadata['treatment'].unique())}\n\n")
            
            # Quality metrics summary
            f.write("## Quality Metrics Summary\n\n")
            f.write("| Metric | Mean | Median | Min | Max |\n")
            f.write("|--------|------|--------|-----|-----|\n")
            f.write(f"| Library Size | {self.quality_metrics['library_size'].mean():.0f} | "
                   f"{self.quality_metrics['library_size'].median():.0f} | "
                   f"{self.quality_metrics['library_size'].min():.0f} | "
                   f"{self.quality_metrics['library_size'].max():.0f} |\n")
            f.write(f"| Gene Detection Rate (%) | {self.quality_metrics['gene_detection_rate'].mean():.1f} | "
                   f"{self.quality_metrics['gene_detection_rate'].median():.1f} | "
                   f"{self.quality_metrics['gene_detection_rate'].min():.1f} | "
                   f"{self.quality_metrics['gene_detection_rate'].max():.1f} |\n\n")
            
            # Sample outlier results
            if 'combined' in self.sample_outliers:
                f.write("## Sample Outlier Detection Results\n\n")
                f.write(f"- **Consensus outliers:** {len(self.sample_outliers['combined']['consensus_outliers'])}\n")
                f.write(f"- **Potential outliers:** {len(self.sample_outliers['combined']['potential_outliers'])}\n\n")
                
                if self.sample_outliers['combined']['consensus_outliers']:
                    f.write("### Consensus Outliers\n\n")
                    for sample in self.sample_outliers['combined']['consensus_outliers']:
                        score = self.sample_outliers['combined']['outlier_scores'][sample]
                        f.write(f"- **{sample}** (flagged by {score} methods)\n")
                    f.write("\n")
                
                if self.sample_outliers['combined']['potential_outliers']:
                    f.write("### Potential Outliers\n\n")
                    for sample in self.sample_outliers['combined']['potential_outliers']:
                        score = self.sample_outliers['combined']['outlier_scores'][sample]
                        f.write(f"- **{sample}** (flagged by {score} method)\n")
                    f.write("\n")
            
            # Gene outlier results
            f.write("## Gene Outlier Detection Results\n\n")
            for method, results in self.gene_outliers.items():
                if 'outliers' in results:
                    f.write(f"### {method.replace('_', ' ').title()}\n\n")
                    f.write(f"- **Outliers detected:** {len(results['outliers'])}\n")
                    if results['outliers']:
                        f.write(f"- **Percentage:** {len(results['outliers'])/len(self.gene_counts)*100:.1f}%\n")
                    f.write("\n")
            
            # Recommendations
            f.write("## Recommendations\n\n")
            f.write("### For Sample Outliers:\n\n")
            f.write("1. **Investigate consensus outliers:** Samples flagged by multiple methods should be carefully examined\n")
            f.write("2. **Check experimental conditions:** Verify sample collection, processing, and sequencing quality\n")
            f.write("3. **Consider removal:** If outliers are due to technical issues, consider removing them from analysis\n")
            f.write("4. **Robust analysis:** If outliers are biological, consider using robust statistical methods\n\n")
            
            f.write("### For Gene Outliers:\n\n")
            f.write("1. **High variance genes:** May represent true biological variation or technical noise\n")
            f.write("2. **Expression pattern outliers:** Check for biological relevance and experimental artifacts\n")
            f.write("3. **Cook's distance outliers:** May indicate influential observations in statistical models\n")
            f.write("4. **Filtering strategy:** Consider removing extreme outliers that are likely technical artifacts\n\n")
            
            f.write("### Handling Strategies:\n\n")
            f.write("1. **Conservative approach:** Remove only consensus outliers with clear technical issues\n")
            f.write("2. **Robust methods:** Use methods that are less sensitive to outliers (e.g., rank-based tests)\n")
            f.write("3. **Winsorization:** Replace extreme values with less extreme ones\n")
            f.write("4. **Multiple analyses:** Run analyses with and without outliers to assess impact\n\n")
        
        logger.info(f"Outlier detection report generated: {report_path}")

    def get_outlier_recommendations(self):
        """Get specific recommendations for handling outliers."""
        recommendations = {
            'sample_outliers': [],
            'gene_outliers': [],
            'general_recommendations': []
        }
        
        # Sample outlier recommendations
        if 'combined' in self.sample_outliers:
            consensus_outliers = self.sample_outliers['combined']['consensus_outliers']
            if consensus_outliers:
                recommendations['sample_outliers'].append(
                    f"Remove {len(consensus_outliers)} consensus outliers: {', '.join(consensus_outliers)}"
                )
            else:
                recommendations['sample_outliers'].append("No consensus outliers detected - all samples appear acceptable")
        
        # Gene outlier recommendations
        total_gene_outliers = 0
        for method, results in self.gene_outliers.items():
            if 'outliers' in results:
                total_gene_outliers += len(results['outliers'])
        
        if total_gene_outliers > 0:
            recommendations['gene_outliers'].append(
                f"Consider filtering {total_gene_outliers} gene outliers across all detection methods"
            )
        else:
            recommendations['gene_outliers'].append("No gene outliers detected")
        
        # General recommendations
        recommendations['general_recommendations'].extend([
            "Run differential expression analysis with and without outliers to assess impact",
            "Use robust statistical methods that are less sensitive to outliers",
            "Consider winsorization for extreme values rather than complete removal",
            "Document all outlier handling decisions for reproducibility"
        ])
        
        return recommendations

    def handle_sample_outliers(self, strategy='conservative'):
        """
        Handle sample outliers using different strategies.
        
        Args:
            strategy: 'conservative' (remove consensus outliers), 'aggressive' (remove all flagged), 
                     'robust' (keep all, use robust methods), 'manual' (return recommendations)
        
        Returns:
            dict: Cleaned data and metadata, or recommendations
        """
        logger.info(f"Handling sample outliers using {strategy} strategy...")
        
        if 'combined' not in self.sample_outliers:
            logger.warning("No outlier detection results available. Run detect_sample_outliers() first.")
            return None
        
        consensus_outliers = self.sample_outliers['combined']['consensus_outliers']
        potential_outliers = self.sample_outliers['combined']['potential_outliers']
        
        if strategy == 'conservative':
            # Remove only consensus outliers (flagged by 2+ methods)
            samples_to_remove = consensus_outliers
            logger.info(f"Conservative strategy: removing {len(samples_to_remove)} consensus outliers")
            
        elif strategy == 'aggressive':
            # Remove all flagged outliers
            samples_to_remove = consensus_outliers + potential_outliers
            logger.info(f"Aggressive strategy: removing {len(samples_to_remove)} total outliers")
            
        elif strategy == 'robust':
            # Keep all samples, recommend robust analysis
            samples_to_remove = []
            logger.info("Robust strategy: keeping all samples, recommend robust statistical methods")
            
        elif strategy == 'manual':
            # Return recommendations without removing anything
            return {
                'recommendations': self.get_outlier_recommendations(),
                'consensus_outliers': consensus_outliers,
                'potential_outliers': potential_outliers,
                'outlier_scores': self.sample_outliers['combined']['outlier_scores']
            }
        
        # Create cleaned dataset
        if samples_to_remove:
            cleaned_counts = self.gene_counts.drop(columns=samples_to_remove)
            cleaned_metadata = self.metadata.drop(index=samples_to_remove)
            
            logger.info(f"Created cleaned dataset: {cleaned_counts.shape[1]} samples remaining")
            
            return {
                'cleaned_counts': cleaned_counts,
                'cleaned_metadata': cleaned_metadata,
                'removed_samples': samples_to_remove,
                'strategy': strategy
            }
        else:
            return {
                'cleaned_counts': self.gene_counts,
                'cleaned_metadata': self.metadata,
                'removed_samples': [],
                'strategy': strategy
            }

    def handle_gene_outliers(self, strategy='conservative', threshold_percentile=95):
        """
        Handle gene outliers using different strategies.
        
        Args:
            strategy: 'conservative' (remove extreme outliers), 'aggressive' (remove all flagged),
                     'winsorize' (cap extreme values), 'robust' (keep all, use robust methods)
            threshold_percentile: Percentile threshold for outlier detection
        
        Returns:
            dict: Cleaned gene counts or recommendations
        """
        logger.info(f"Handling gene outliers using {strategy} strategy...")
        
        if not self.gene_outliers:
            logger.warning("No gene outlier detection results available. Run detect_gene_outliers() first.")
            return None
        
        # Collect all gene outliers
        all_gene_outliers = set()
        for method, results in self.gene_outliers.items():
            if 'outliers' in results:
                all_gene_outliers.update(results['outliers'])
        
        if strategy == 'conservative':
            # Remove only extreme outliers (high variance genes)
            if 'variance' in self.gene_outliers:
                high_var_outliers = self.gene_outliers['variance']['high_variance_outliers']
                genes_to_remove = high_var_outliers
                logger.info(f"Conservative strategy: removing {len(genes_to_remove)} high variance outliers")
            else:
                genes_to_remove = []
                
        elif strategy == 'aggressive':
            # Remove all flagged gene outliers
            genes_to_remove = list(all_gene_outliers)
            logger.info(f"Aggressive strategy: removing {len(genes_to_remove)} total gene outliers")
            
        elif strategy == 'winsorize':
            # Winsorize extreme values instead of removing genes
            return self._winsorize_gene_counts(threshold_percentile)
            
        elif strategy == 'robust':
            # Keep all genes, recommend robust analysis
            genes_to_remove = []
            logger.info("Robust strategy: keeping all genes, recommend robust statistical methods")
        
        # Create cleaned dataset
        if genes_to_remove:
            cleaned_counts = self.gene_counts.drop(index=genes_to_remove)
            logger.info(f"Created cleaned dataset: {cleaned_counts.shape[0]} genes remaining")
            
            return {
                'cleaned_counts': cleaned_counts,
                'removed_genes': genes_to_remove,
                'strategy': strategy
            }
        else:
            return {
                'cleaned_counts': self.gene_counts,
                'removed_genes': [],
                'strategy': strategy
            }

    def _winsorize_gene_counts(self, percentile=95):
        """
        Winsorize gene counts to handle outliers by capping extreme values.
        
        Args:
            percentile: Percentile to use for winsorization (e.g., 95 = cap at 95th percentile)
        
        Returns:
            dict: Winsorized counts and information
        """
        logger.info(f"Winsorizing gene counts at {percentile}th percentile...")
        
        # Log-transform the data
        log_counts = np.log2(self.gene_counts + 1)
        
        # Calculate percentiles for each gene
        lower_percentile = 100 - percentile
        upper_percentile = percentile
        
        winsorized_counts = log_counts.copy()
        
        # Track which genes were winsorized
        winsorized_genes = []
        
        for gene in log_counts.index:
            gene_values = log_counts.loc[gene]
            lower_bound = np.percentile(gene_values, lower_percentile)
            upper_bound = np.percentile(gene_values, upper_percentile)
            
            # Check if any values were winsorized
            original_min = gene_values.min()
            original_max = gene_values.max()
            
            # Winsorize
            winsorized_values = np.clip(gene_values, lower_bound, upper_bound)
            winsorized_counts.loc[gene] = winsorized_values
            
            # Check if winsorization occurred
            if (winsorized_values != gene_values).any():
                winsorized_genes.append(gene)
        
        # Convert back to counts (approximate)
        winsorized_counts_rounded = np.round(2**winsorized_counts - 1)
        winsorized_counts_rounded = np.maximum(winsorized_counts_rounded, 0).astype(int)
        
        logger.info(f"Winsorized {len(winsorized_genes)} genes")
        
        return {
            'winsorized_counts': winsorized_counts_rounded,
            'winsorized_genes': winsorized_genes,
            'percentile': percentile,
            'strategy': 'winsorize'
        }

    def create_robust_analysis_datasets(self):
        """
        Create multiple datasets for robust analysis with different outlier handling strategies.
        
        Returns:
            dict: Multiple datasets with different outlier handling approaches
        """
        logger.info("Creating robust analysis datasets...")
        
        datasets = {}
        
        # 1. Original dataset (no outlier handling)
        datasets['original'] = {
            'counts': self.gene_counts,
            'metadata': self.metadata,
            'description': 'Original dataset with no outlier handling'
        }
        
        # 2. Conservative sample outlier removal
        sample_handling = self.handle_sample_outliers(strategy='conservative')
        if sample_handling and 'cleaned_counts' in sample_handling:
            datasets['conservative_samples'] = {
                'counts': sample_handling['cleaned_counts'],
                'metadata': sample_handling['cleaned_metadata'],
                'description': f"Conservative sample outlier removal ({len(sample_handling['removed_samples'])} samples removed)"
            }
        
        # 3. Aggressive sample outlier removal
        sample_handling_agg = self.handle_sample_outliers(strategy='aggressive')
        if sample_handling_agg and 'cleaned_counts' in sample_handling_agg:
            datasets['aggressive_samples'] = {
                'counts': sample_handling_agg['cleaned_counts'],
                'metadata': sample_handling_agg['cleaned_metadata'],
                'description': f"Aggressive sample outlier removal ({len(sample_handling_agg['removed_samples'])} samples removed)"
            }
        
        # 4. Winsorized gene counts
        gene_winsorized = self.handle_gene_outliers(strategy='winsorize')
        if gene_winsorized and 'winsorized_counts' in gene_winsorized:
            datasets['winsorized_genes'] = {
                'counts': gene_winsorized['winsorized_counts'],
                'metadata': self.metadata,
                'description': f"Winsorized gene counts ({len(gene_winsorized['winsorized_genes'])} genes affected)"
            }
        
        # 5. Conservative gene outlier removal
        gene_handling = self.handle_gene_outliers(strategy='conservative')
        if gene_handling and 'cleaned_counts' in gene_handling:
            datasets['conservative_genes'] = {
                'counts': gene_handling['cleaned_counts'],
                'metadata': self.metadata,
                'description': f"Conservative gene outlier removal ({len(gene_handling['removed_genes'])} genes removed)"
            }
        
        # 6. Combined conservative approach
        if 'conservative_samples' in datasets and 'conservative_genes' in datasets:
            # Apply gene filtering to sample-cleaned data
            sample_cleaned_counts = datasets['conservative_samples']['counts']
            gene_cleaned_counts = datasets['conservative_genes']['counts']
            
            # Find intersection of genes
            common_genes = sample_cleaned_counts.index.intersection(gene_cleaned_counts.index)
            combined_counts = sample_cleaned_counts.loc[common_genes]
            combined_metadata = datasets['conservative_samples']['metadata']
            
            datasets['combined_conservative'] = {
                'counts': combined_counts,
                'metadata': combined_metadata,
                'description': f"Combined conservative approach: {combined_counts.shape[1]} samples, {combined_counts.shape[0]} genes"
            }
        
        logger.info(f"Created {len(datasets)} analysis datasets")
        
        return datasets

    def save_robust_datasets(self, datasets, output_dir=None):
        """
        Save robust analysis datasets to files.
        
        Args:
            datasets: Dictionary of datasets from create_robust_analysis_datasets()
            output_dir: Directory to save datasets (default: outlier_dir)
        """
        if output_dir is None:
            output_dir = self.outlier_dir / "robust_datasets"
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Saving robust datasets to {output_dir}")
        
        # Save each dataset
        for name, dataset in datasets.items():
            counts_file = output_dir / f"{name}_counts.csv"
            metadata_file = output_dir / f"{name}_metadata.csv"
            info_file = output_dir / f"{name}_info.txt"
            
            # Save counts
            dataset['counts'].to_csv(counts_file)
            
            # Save metadata
            dataset['metadata'].to_csv(metadata_file)
            
            # Save info
            with open(info_file, 'w') as f:
                f.write(f"Dataset: {name}\n")
                f.write(f"Description: {dataset['description']}\n")
                f.write(f"Shape: {dataset['counts'].shape[0]} genes x {dataset['counts'].shape[1]} samples\n")
                f.write(f"Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # Create summary file
        summary_file = output_dir / "datasets_summary.md"
        with open(summary_file, 'w') as f:
            f.write("# Robust Analysis Datasets Summary\n\n")
            f.write(f"**Created:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Available Datasets\n\n")
            for name, dataset in datasets.items():
                f.write(f"### {name}\n")
                f.write(f"- **Description:** {dataset['description']}\n")
                f.write(f"- **Shape:** {dataset['counts'].shape[0]} genes x {dataset['counts'].shape[1]} samples\n")
                f.write(f"- **Files:** {name}_counts.csv, {name}_metadata.csv, {name}_info.txt\n\n")
        
        logger.info(f"Saved {len(datasets)} datasets to {output_dir}")

    def generate_outlier_handling_report(self):
        """Generate a comprehensive report on outlier handling strategies."""
        logger.info("Generating outlier handling report...")
        
        report_path = self.outlier_dir / 'outlier_handling_strategies.md'
        
        with open(report_path, 'w') as f:
            f.write("# RNA-seq Outlier Handling Strategies\n\n")
            f.write(f"**Dataset:** PRJNA158021 (Embryo Dataset)\n")
            f.write(f"**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Overview\n\n")
            f.write("This report provides practical strategies for handling outliers in RNA-seq data.\n\n")
            
            f.write("## Sample-Level Outlier Handling\n\n")
            f.write("### Conservative Approach\n")
            f.write("- Remove only samples flagged as outliers by multiple detection methods\n")
            f.write("- Safest approach, minimizes false positives\n")
            f.write("- Recommended for most analyses\n\n")
            
            f.write("### Aggressive Approach\n")
            f.write("- Remove all samples flagged by any detection method\n")
            f.write("- May remove biologically relevant samples\n")
            f.write("- Use only when technical issues are suspected\n\n")
            
            f.write("### Robust Approach\n")
            f.write("- Keep all samples, use robust statistical methods\n")
            f.write("- Methods: rank-based tests, median-based normalization\n")
            f.write("- Preserves biological variation\n\n")
            
            f.write("## Gene-Level Outlier Handling\n\n")
            f.write("### Conservative Approach\n")
            f.write("- Remove only genes with extremely high variance\n")
            f.write("- Focus on likely technical artifacts\n")
            f.write("- Preserves most biological signal\n\n")
            
            f.write("### Winsorization\n")
            f.write("- Cap extreme values rather than removing genes\n")
            f.write("- Preserves gene presence in analysis\n")
            f.write("- Reduces impact of outliers on statistics\n\n")
            
            f.write("### Aggressive Approach\n")
            f.write("- Remove all genes flagged as outliers\n")
            f.write("- May remove biologically interesting genes\n")
            f.write("- Use with caution\n\n")
            
            f.write("## Recommended Workflow\n\n")
            f.write("1. **Run outlier detection** to identify potential issues\n")
            f.write("2. **Investigate consensus outliers** for technical problems\n")
            f.write("3. **Use conservative approach** for sample outliers\n")
            f.write("4. **Consider winsorization** for gene outliers\n")
            f.write("5. **Run multiple analyses** with different strategies\n")
            f.write("6. **Compare results** to assess robustness\n")
            f.write("7. **Document decisions** for reproducibility\n\n")
            
            f.write("## Implementation\n\n")
            f.write("Use the `create_robust_analysis_datasets()` method to generate multiple datasets:\n\n")
            f.write("```python\n")
            f.write("detector = RNAseqOutlierDetector()\n")
            f.write("detector.load_data()\n")
            f.write("detector.detect_sample_outliers()\n")
            f.write("detector.detect_gene_outliers()\n")
            f.write("datasets = detector.create_robust_analysis_datasets()\n")
            f.write("detector.save_robust_datasets(datasets)\n")
            f.write("```\n\n")
            
            f.write("## Quality Control\n\n")
            f.write("- Always compare results across different outlier handling strategies\n")
            f.write("- Check that key biological findings are robust\n")
            f.write("- Document any samples or genes removed and the rationale\n")
            f.write("- Consider the impact on statistical power\n\n")
        
        logger.info(f"Outlier handling report generated: {report_path}")

def main():
    """Main function to run outlier detection analysis."""
    logger.info("=== RNA-seq Outlier Detection Analysis ===")
    
    # Initialize detector
    detector = RNAseqOutlierDetector()
    
    try:
        # Load data
        detector.load_data()
        
        # Run sample-level outlier detection
        detector.detect_sample_outliers()
        
        # Run gene-level outlier detection
        detector.detect_gene_outliers()
        
        # Create visualizations
        detector.create_outlier_visualizations()
        
        # Generate reports
        detector.generate_outlier_report()
        detector.generate_outlier_handling_report()
        
        # Create robust analysis datasets
        logger.info("Creating robust analysis datasets...")
        datasets = detector.create_robust_analysis_datasets()
        detector.save_robust_datasets(datasets)
        
        # Get recommendations
        recommendations = detector.get_outlier_recommendations()
        
        # Print summary
        logger.info("=== OUTLIER DETECTION COMPLETED ===")
        logger.info("Key Results:")
        
        if 'combined' in detector.sample_outliers:
            consensus_count = len(detector.sample_outliers['combined']['consensus_outliers'])
            potential_count = len(detector.sample_outliers['combined']['potential_outliers'])
            logger.info(f"- Sample outliers: {consensus_count} consensus, {potential_count} potential")
        
        total_gene_outliers = sum(len(results.get('outliers', [])) for results in detector.gene_outliers.values())
        logger.info(f"- Gene outliers: {total_gene_outliers} total across all methods")
        
        logger.info(f"\nCreated {len(datasets)} robust analysis datasets:")
        for name, dataset in datasets.items():
            logger.info(f"  - {name}: {dataset['description']}")
        
        logger.info("\nRecommendations:")
        for category, recs in recommendations.items():
            logger.info(f"\n{category.replace('_', ' ').title()}:")
            for rec in recs:
                logger.info(f"  - {rec}")
        
        logger.info(f"\nReports generated:")
        logger.info(f"  - Outlier detection: {detector.outlier_dir / 'outlier_detection_report.md'}")
        logger.info(f"  - Handling strategies: {detector.outlier_dir / 'outlier_handling_strategies.md'}")
        logger.info(f"  - Robust datasets: {detector.outlier_dir / 'robust_datasets/datasets_summary.md'}")
        
        logger.info("\n=== OUTLIER HANDLING STRATEGIES ===")
        logger.info("1. Conservative: Remove only consensus outliers (recommended)")
        logger.info("2. Aggressive: Remove all flagged outliers")
        logger.info("3. Winsorization: Cap extreme values instead of removal")
        logger.info("4. Robust: Keep all data, use robust statistical methods")
        logger.info("\nUse the created datasets for differential expression analysis with different outlier handling approaches.")
        
    except Exception as e:
        logger.error(f"Error in outlier detection analysis: {e}")
        raise

if __name__ == "__main__":
    main() 