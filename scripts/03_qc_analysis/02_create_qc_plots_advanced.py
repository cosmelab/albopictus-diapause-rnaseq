import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy import stats

# Set style and create output directory
plt.style.use('default')
sns.set_palette("husl")
Path("output/figures").mkdir(parents=True, exist_ok=True)

# Read the QC metrics
df = pd.read_csv('output/tables/Table_S2_QC_Metrics.csv')

# 1. Sample Quality Overview Heatmap
def create_quality_heatmap():
    # Select numeric columns for heatmap
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    heatmap_data = df[numeric_cols].copy()
    
    # Z-score normalization
    scaler = StandardScaler()
    heatmap_normalized = pd.DataFrame(
        scaler.fit_transform(heatmap_data),
        columns=heatmap_data.columns,
        index=df['Sample ID']
    )
    
    # Create figure
    plt.figure(figsize=(20, 12))
    
    # Create heatmap
    sns.heatmap(heatmap_normalized, 
                cmap='RdYlGn_r', 
                center=0,
                cbar_kws={'label': 'Z-score'},
                xticklabels=True,
                yticklabels=True)
    
    plt.title('Sample Quality Overview Heatmap', pad=20, fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save in multiple formats
    plt.savefig('output/figures/QC_Heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/QC_Heatmap.eps', format='eps', bbox_inches='tight')
    plt.savefig('output/figures/QC_Heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

# 2. Technical Replicates Correlation
def create_replicate_correlations():
    # Define replicate pairs
    replicate_pairs = [
        ('SD_72h', 'LD_72h'),
        ('SD_135h', 'LD_135h'),
        ('SD_Unfed', 'LD_Unfed'),
        ('SD_Fed', 'LD_Fed'),
        ('SD_11d', 'LD_11d'),
        ('SD_21d', 'LD_21d'),
        ('SD_40d', 'LD_40d')
    ]
    
    # Create subplots
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    # Metrics to plot
    metrics = ['Total_Reads (M)', 'Mapping_Rate (%)', 'Uniquely_Mapped (%)', 
              'Multi_Mapped (%)', 'Duplication_Rate (%)', 'Adapter_Contamination (%)',
              'Error_Rate (%)', 'Insert_Size (bp)']
    
    for idx, (sd, ld) in enumerate(replicate_pairs):
        if idx < len(axes):
            # Get data for the pair
            sd_data = df[df['Sample ID'].str.contains(sd)][metrics].mean()
            ld_data = df[df['Sample ID'].str.contains(ld)][metrics].mean()
            
            # Calculate correlation
            corr = stats.pearsonr(sd_data, ld_data)[0]
            
            # Create scatter plot
            axes[idx].scatter(sd_data, ld_data, alpha=0.6)
            
            # Add diagonal line
            min_val = min(min(sd_data), min(ld_data))
            max_val = max(max(sd_data), max(ld_data))
            axes[idx].plot([min_val, max_val], [min_val, max_val], 
                         'r--', alpha=0.3)
            
            # Add labels and title
            axes[idx].set_xlabel(f'{sd} Values')
            axes[idx].set_ylabel(f'{ld} Values')
            axes[idx].set_title(f'{sd} vs {ld}\nCorrelation: {corr:.2f}')
            
            # Add grid
            axes[idx].grid(True, alpha=0.3)
            
            # Add metric labels
            for i, metric in enumerate(metrics):
                axes[idx].annotate(metric.replace('_', ' ').replace(' (%)', ''),
                                 (sd_data[i], ld_data[i]),
                                 fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('output/figures/QC_Replicate_Correlations.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/QC_Replicate_Correlations.eps', format='eps', bbox_inches='tight')
    plt.savefig('output/figures/QC_Replicate_Correlations.png', dpi=300, bbox_inches='tight')
    plt.close()

# 3. Batch Effect Assessment (PCA)
def create_pca_plot():
    # Select key metrics for PCA
    metrics = [
        'Total_Reads (M)', 'Read_Length (bp)', 'Duplication_Rate (%)',
        'Adapter_Contamination (%)', 'Mapping_Rate (%)', 'Uniquely_Mapped (%)',
        'Multi_Mapped (%)', 'Exonic_Reads (%)', 'Intronic_Reads (%)',
        'Intergenic_Reads (%)', 'Insert_Size (bp)', 'Properly_Paired (%)',
        'Error_Rate (%)', 'Five_Three_Bias (ratio)'
    ]
    
    # Prepare data for PCA
    pca_data = df[metrics].copy()
    
    # Standardize the data
    scaler = StandardScaler()
    pca_scaled = scaler.fit_transform(pca_data)
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(pca_scaled)
    
    # Create DataFrame with PCA results
    pca_df = pd.DataFrame(
        data=pca_result,
        columns=['PC1', 'PC2']
    )
    pca_df['Project'] = df['Project']
    pca_df['Sample ID'] = df['Sample ID']
    
    # Calculate variance explained
    var_explained = pca.explained_variance_ratio_ * 100
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Plot points
    for project in pca_df['Project'].unique():
        project_data = pca_df[pca_df['Project'] == project]
        plt.scatter(project_data['PC1'], project_data['PC2'], 
                   label=project, alpha=0.7)
    
    # Add labels and title
    plt.xlabel(f'PC1 ({var_explained[0]:.1f}% variance explained)')
    plt.ylabel(f'PC2 ({var_explained[1]:.1f}% variance explained)')
    plt.title('PCA of QC Metrics by Project')
    
    # Add legend inside plot at top
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3)
    plt.grid(True, alpha=0.3)
    
    # Add sample labels
    for i, txt in enumerate(pca_df['Sample ID']):
        plt.annotate(txt, (pca_df['PC1'].iloc[i], pca_df['PC2'].iloc[i]),
                    fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('output/figures/QC_PCA.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/QC_PCA.eps', format='eps', bbox_inches='tight')
    plt.savefig('output/figures/QC_PCA.png', dpi=300, bbox_inches='tight')
    plt.close()

# 4. Project-wise Boxplot Comparisons
def create_project_boxplots():
    # Select key metrics for comparison
    metrics = ['Mapping_Rate (%)', 'Uniquely_Mapped (%)', 'Multi_Mapped (%)',
              'Duplication_Rate (%)', 'Adapter_Contamination (%)', 'Error_Rate (%)']
    
    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        # Create boxplot
        sns.boxplot(data=df, x='Project', y=metric, ax=axes[idx])
        
        # Add individual points
        sns.stripplot(data=df, x='Project', y=metric, ax=axes[idx],
                     color='black', alpha=0.3)
        
        # Customize plot
        axes[idx].set_title(metric.replace('_', ' '))
        axes[idx].set_xticklabels(axes[idx].get_xticklabels(), rotation=45)
        axes[idx].grid(True, axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('output/figures/QC_Project_Boxplots.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/QC_Project_Boxplots.eps', format='eps', bbox_inches='tight')
    plt.savefig('output/figures/QC_Project_Boxplots.png', dpi=300, bbox_inches='tight')
    plt.close()

# Create all plots
def create_all_plots():
    create_quality_heatmap()
    create_replicate_correlations()
    create_pca_plot()
    create_project_boxplots()
    print("All plots have been generated in output/figures/")

if __name__ == "__main__":
    create_all_plots() 