import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import re
import ast
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.decomposition import PCA
from matplotlib import gridspec
from sklearn.preprocessing import StandardScaler

# Set style and create output directory
plt.style.use('default')
sns.set_palette("colorblind")  # Colorblind-friendly palette
Path("output/figures/publication").mkdir(parents=True, exist_ok=True)

# Set font sizes
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10

# Read the data
df = pd.read_csv("output/tables/Table_S2_QC_Metrics.csv")

# Load sample metadata for matching
meta = pd.read_csv('data/metadata/samples.csv')

# Merge metadata with QC table for matching
qc_annot = df.merge(meta, left_on=['Project', 'SRR'], right_on=['BioProject', 'Run'], how='left')

# Define sample order and grouping
sample_order = [
    # PRJNA268379 (Fed/Unfed)
    'SD_Unfed', 'LD_Unfed', 'SD_Fed', 'LD_Fed',
    # PRJNA158021 (72h/135h)
    'SD_72h', 'LD_72h', 'SD_135h', 'LD_135h',
    # PRJNA187045 (11d/21d/40d)
    'SD_11d', 'LD_11d', 'SD_21d', 'LD_21d', 'SD_40d', 'LD_40d'
]

# Create a mapping for sorting
def get_sample_order(sample_id):
    for i, pattern in enumerate(sample_order):
        if pattern in sample_id:
            return i
    return len(sample_order)

# Sort the dataframe
df['sort_order'] = df['Sample ID'].apply(get_sample_order)
df = df.sort_values('sort_order')

# Panel A: Boxplots of key metrics
plt.figure(figsize=(10, 6))
metrics = ['Mapping_Rate (%)', 'Uniquely_Mapped (%)', 'Duplication_Rate (%)']
data = pd.melt(df, id_vars=['Project'], value_vars=metrics)
sns.boxplot(data=data, x='Project', y='value', hue='variable')
plt.title('A', fontsize=14, fontweight='bold', pad=20)
plt.xlabel('Project')
plt.ylabel('Percentage (%)')
plt.legend(title='Metric', bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=3)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('output/figures/publication/Figure1A_QC_Metrics.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.savefig('output/figures/publication/Figure1A_QC_Metrics.png', format='png', bbox_inches='tight', dpi=300)
plt.close()

# Panel B: PCA plot
plt.figure(figsize=(10, 6))
# Select all QC metrics for PCA
qc_metrics = [col for col in df.columns if col not in ['Project', 'SRR', 'Sample ID', 'sort_order']]
X = df[qc_metrics].values
X_scaled = StandardScaler().fit_transform(X)
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Create color map for projects
project_colors = sns.color_palette("Set2", n_colors=len(df['Project'].unique()))
project_color_map = dict(zip(df['Project'].unique(), project_colors))

# Plot PCA with filled circles and colored borders
for project in df['Project'].unique():
    mask = df['Project'] == project
    plt.scatter(X_pca[mask, 0], X_pca[mask, 1], 
               c=[project_color_map[project]], 
               edgecolor='black',
               s=100,
               label=project)

plt.title('B', fontsize=14, fontweight='bold', pad=20)
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=3)

# Add sample labels
for i, txt in enumerate(df['Sample ID']):
    plt.annotate(txt, (X_pca[i, 0], X_pca[i, 1]), fontsize=8)

plt.tight_layout()
plt.savefig('output/figures/publication/Figure1B_PCA.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.savefig('output/figures/publication/Figure1B_PCA.png', format='png', bbox_inches='tight', dpi=300)
plt.close()

# Panel C: Stacked bar chart of RNA-seq feature distribution
plt.figure(figsize=(10, 6))
metrics = ['Exonic_Reads (%)', 'Intronic_Reads (%)', 'Intergenic_Reads (%)']
bottom = np.zeros(len(df))
for metric in metrics:
    plt.bar(range(len(df)), df[metric], bottom=bottom, label=metric.replace('_', ' '))
    bottom += df[metric]

plt.title('C', fontsize=14, fontweight='bold', pad=20)
plt.xlabel('Samples')
plt.ylabel('Percentage (%)')
plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=3)
plt.tight_layout()
plt.savefig('output/figures/publication/Figure1C_Feature_Distribution.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.savefig('output/figures/publication/Figure1C_Feature_Distribution.png', format='png', bbox_inches='tight', dpi=300)
plt.close()

# Supplementary Figure S2: Gene coverage profiles
plt.figure(figsize=(18, 7))  # Wider figure
for project in df['Project'].unique():
    mask = df['Project'] == project
    plt.scatter(
        df[mask].index, df[mask]['Five_Three_Bias (ratio)'],
        label=project,
        s=80,
        c=[project_color_map[project]],
        edgecolor='black',
        zorder=3
    )

plt.axhline(y=1.0, color='r', linestyle='--', label='No bias')
plt.title('5\' to 3\' Coverage Bias', pad=20, fontsize=14)
plt.xlabel('Samples', labelpad=10)
plt.ylabel('5\'/3\' Ratio', labelpad=10)
plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=len(df['Project'].unique())+1)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('output/figures/publication/FigureS2_Coverage_Bias.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.savefig('output/figures/publication/FigureS2_Coverage_Bias.png', format='png', bbox_inches='tight', dpi=300)
plt.close()

# Supplementary Figure S3: Technical replicates correlation matrix
# Now, match SD and LD samples by biological condition (using 'New abreviation' and 'Replicate' columns)
def create_correlation_matrix_by_condition(project):
    project_qc = qc_annot[qc_annot['Project'] == project]
    # Use only rows with both SD and LD for the same biological condition (New abreviation minus SD/LD prefix)
    # Extract base condition (e.g., 72h, 135h, Fed, Unfed, 11d, 21d, 40d) and replicate
    def base_cond(row):
        abbr = row['New abreviation']
        if abbr.startswith('SD_'):
            return abbr[3:]
        elif abbr.startswith('LD_'):
            return abbr[3:]
        return abbr
    project_qc['base_cond'] = project_qc.apply(base_cond, axis=1)
    # Pair by base_cond and Replicate
    pairs = []
    for cond in project_qc['base_cond'].unique():
        for rep in project_qc['Replicate'].unique():
            sd = project_qc[(project_qc['New abreviation'].str.startswith('SD_')) & (project_qc['base_cond'] == cond) & (project_qc['Replicate'] == rep)]
            ld = project_qc[(project_qc['New abreviation'].str.startswith('LD_')) & (project_qc['base_cond'] == cond) & (project_qc['Replicate'] == rep)]
            if len(sd) == 1 and len(ld) == 1:
                pairs.append((sd.iloc[0], ld.iloc[0]))
    if not pairs:
        return None, None
    # Build matrices for correlation
    metrics = ['Mapping_Rate (%)', 'Uniquely_Mapped (%)', 'Duplication_Rate (%)',
               'Exonic_Reads (%)', 'Intronic_Reads (%)', 'Intergenic_Reads (%)']
    sd_mat = np.array([[pair[0][m] for m in metrics] for pair in pairs])
    ld_mat = np.array([[pair[1][m] for m in metrics] for pair in pairs])
    corr_matrix = np.corrcoef(sd_mat.T, ld_mat.T)[:len(metrics), len(metrics):]
    return corr_matrix, metrics

# Create correlation matrices for each project using improved matching
projects = df['Project'].unique()
fig, axes = plt.subplots(1, len(projects), figsize=(8 * len(projects), 7))  # Wider and taller

for i, project in enumerate(projects):
    corr_matrix, metrics = create_correlation_matrix_by_condition(project)
    if corr_matrix is not None and metrics is not None:
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        sns.heatmap(
            corr_matrix,
            annot=True,
            cmap=cmap,
            center=0,
            xticklabels=metrics,
            yticklabels=metrics,
            ax=axes[i],
            square=True,
            fmt='.2f',
            annot_kws={'size': 16},  # Larger annotation font
            cbar_kws={'label': 'Pearson correlation'}
        )
        axes[i].set_title(f'{project}\nSD vs LD Correlation', fontsize=18)
        axes[i].set_xlabel('LD Metrics', fontsize=16)
        axes[i].set_ylabel('SD Metrics', fontsize=16)
        axes[i].tick_params(axis='x', labelsize=14, rotation=45)
        axes[i].tick_params(axis='y', labelsize=14)
        # Add gridlines for clarity
        for _, spine in axes[i].spines.items():
            spine.set_visible(True)
    else:
        axes[i].text(0.5, 0.5, f'No matching SD/LD pairs\nfor {project}',
                     horizontalalignment='center', verticalalignment='center',
                     transform=axes[i].transAxes, fontsize=16)
        axes[i].set_title(f'{project}\nSD vs LD Correlation', fontsize=18)
plt.tight_layout()
plt.savefig('output/figures/publication/FigureS3_Technical_Replicates.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.savefig('output/figures/publication/FigureS3_Technical_Replicates.png', format='png', bbox_inches='tight', dpi=300)
plt.close()

def create_supplementary_figures():
    """Create supplementary figures"""
    # S1: Junction saturation curves
    plt.figure(figsize=(15, 10))
    junction_files = glob.glob('results/zenodo/*/multiqc/*_rseqc_junction_saturation_plot_All_Junctions.txt')
    
    for file in junction_files:
        project = re.search(r'PRJNA\d+', file).group()
        srr_id = re.search(r'SRR\d+', file).group()
        
        # Get the corresponding Sample ID from the dataframe
        matching_samples = df[df['SRR'] == srr_id]
        if len(matching_samples) > 0:
            sample_id = matching_samples['Sample ID'].iloc[0]
        else:
            print(f"Warning: Could not find matching Sample ID for {srr_id}")
            continue
        
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    data_line = lines[1].strip().split('\t')
                    reads = []
                    junctions = []
                    for value in data_line[1:]:
                        if value.startswith('('):
                            read, junction = ast.literal_eval(value)
                            reads.append(read)
                            junctions.append(junction)
                    
                    plt.plot(reads, junctions, label=sample_id, alpha=0.5)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    plt.title('Supplementary Figure S1: Junction Saturation Analysis', pad=20, fontsize=20)
    plt.xlabel('Number of Reads', fontsize=16)
    plt.ylabel('Number of Junctions', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('output/figures/publication/FigureS1_Junction_Saturation.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/publication/FigureS1_Junction_Saturation.eps', format='eps', bbox_inches='tight')
    plt.close()
    
    # S2: Gene coverage profiles
    plt.figure(figsize=(20, 10))
    coverage_files = glob.glob('results/zenodo/*/multiqc/*_qualimap_gene_coverage_profile_Normalised.txt')
    
    for file in coverage_files:
        project = re.search(r'PRJNA\d+', file).group()
        srr_id = re.search(r'SRR\d+', file).group()
        
        # Get the corresponding Sample ID from the dataframe
        matching_samples = df[df['SRR'] == srr_id]
        if len(matching_samples) > 0:
            sample_id = matching_samples['Sample ID'].iloc[0]
        else:
            print(f"Warning: Could not find matching Sample ID for {srr_id}")
            continue
        
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    data_line = lines[1].strip().split('\t')
                    positions = []
                    coverages = []
                    for value in data_line[1:]:
                        if value.startswith('('):
                            pos, cov = ast.literal_eval(value)
                            positions.append(pos)
                            coverages.append(cov)
                    
                    plt.plot(positions, coverages, alpha=0.5)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    plt.title('Supplementary Figure S2: Gene Coverage Profiles', pad=20, fontsize=20)
    plt.xlabel('Gene Position (5\' to 3\')', fontsize=16)
    plt.ylabel('Normalized Coverage', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('output/figures/publication/FigureS2_Gene_Coverage.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/publication/FigureS2_Gene_Coverage.eps', format='eps', bbox_inches='tight')
    plt.close()
    
    # S3: Technical replicates correlation
    plt.figure(figsize=(15, 10))
    
    # Group by project
    for project in df['Project'].unique():
        project_data = df[df['Project'] == project]
        
        # Get SD and LD samples for this project
        sd_samples = project_data[project_data['Sample ID'].str.startswith('SD')]
        ld_samples = project_data[project_data['Sample ID'].str.startswith('LD')]
        
        if len(sd_samples) > 0 and len(ld_samples) > 0:
            # Get common metrics
            metrics = ['Total_Reads (M)', 'Duplication_Rate (%)', 'Insert_Size (bp)']
            
            for metric in metrics:
                sd_data = sd_samples[metric].values
                ld_data = ld_samples[metric].values
                
                # Ensure arrays have same length by using the minimum length
                min_len = min(len(sd_data), len(ld_data))
                if min_len > 0:
                    sd_data = sd_data[:min_len]
                    ld_data = ld_data[:min_len]
                    
                    corr = stats.pearsonr(sd_data, ld_data)[0]
                    
                    plt.scatter(sd_data, ld_data, 
                              label=f'{project} - {metric} (r={corr:.2f})',
                              alpha=0.6)
    
    plt.title('Supplementary Figure S3: Technical Replicates Correlation', pad=20, fontsize=20)
    plt.xlabel('SD Values', fontsize=16)
    plt.ylabel('LD Values', fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('output/figures/publication/FigureS3_Technical_Replicates.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('output/figures/publication/FigureS3_Technical_Replicates.eps', format='eps', bbox_inches='tight')
    plt.close()

def create_all_figures():
    create_supplementary_figures()
    print("All publication-ready figures have been generated in output/figures/publication/")

if __name__ == "__main__":
    create_all_figures() 