import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import re
import ast

# Set style and create output directory
plt.style.use('default')
sns.set_palette("husl")
Path("output/figures").mkdir(parents=True, exist_ok=True)

# Read the QC metrics
df = pd.read_csv('output/tables/Table_S2_QC_Metrics.csv')

def extract_condition(sample_id):
    """Extract biological condition from sample ID"""
    if 'Unfed' in sample_id:
        return 'Unfed'
    elif 'Fed' in sample_id:
        return 'Fed'
    elif any(t in sample_id for t in ['72h', '135h', '11d', '21d', '40d']):
        return 'Timepoint'
    return 'Other'

# Add condition column
df['Condition'] = df['Sample ID'].apply(extract_condition)

# 1. Gene Coverage Profile
def create_gene_coverage_plot():
    # Read gene coverage files
    coverage_files = glob.glob('results/zenodo/*/multiqc/*_qualimap_gene_coverage_profile_Normalised.txt')
    
    coverage_data = []
    for file in coverage_files:
        project = re.search(r'PRJNA\d+', file).group()
        sample = re.search(r'SRR\d+', file).group()
        
        try:
            # Read the file
            with open(file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Check if file has content
                    # Skip header and get data line
                    data_line = lines[1].strip().split('\t')
                    sample_name = data_line[0]
                    
                    # Extract position and coverage values from the tuple strings
                    positions = []
                    coverages = []
                    for value in data_line[1:]:  # Skip the Sample column
                        if value.startswith('('):
                            pos, cov = ast.literal_eval(value)
                            positions.append(pos)
                            coverages.append(cov)
                    
                    coverage_data.append({
                        'Project': project,
                        'Sample': sample,
                        'Position': positions,
                        'Coverage': coverages
                    })
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if coverage_data:
        # Create plot
        plt.figure(figsize=(15, 8))
        
        for data in coverage_data:
            plt.plot(data['Position'], data['Coverage'], 
                    label=f"{data['Project']} - {data['Sample']}", alpha=0.5)
        
        plt.title('Gene Coverage Profile Across Samples')
        plt.xlabel('Gene Position (5\' to 3\')')
        plt.ylabel('Normalized Coverage')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plots
        plt.savefig('output/figures/QC_Gene_Coverage.pdf', format='pdf', bbox_inches='tight')
        plt.savefig('output/figures/QC_Gene_Coverage.eps', format='eps', bbox_inches='tight')
        plt.savefig('output/figures/QC_Gene_Coverage.png', dpi=300, bbox_inches='tight')
        plt.close()

# 2. Junction Analysis
def create_junction_saturation_plot():
    # Read junction saturation files
    junction_files = glob.glob('results/zenodo/*/multiqc/*_rseqc_junction_saturation_plot_All_Junctions.txt')
    
    junction_data = []
    for file in junction_files:
        project = re.search(r'PRJNA\d+', file).group()
        sample = re.search(r'SRR\d+', file).group()
        
        try:
            # Read the file
            with open(file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Check if file has content
                    # Skip header and get data line
                    data_line = lines[1].strip().split('\t')
                    sample_name = data_line[0]
                    
                    # Extract reads and junctions from the tuple strings
                    reads = []
                    junctions = []
                    for value in data_line[1:]:  # Skip the Sample column
                        if value.startswith('('):
                            read, junction = ast.literal_eval(value)
                            reads.append(read)
                            junctions.append(junction)
                    
                    junction_data.append({
                        'Project': project,
                        'Sample': sample,
                        'Reads': reads,
                        'Junctions': junctions
                    })
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if junction_data:
        # Create plot
        plt.figure(figsize=(12, 8))
        
        for data in junction_data:
            plt.plot(data['Reads'], data['Junctions'], 
                    label=f"{data['Project']} - {data['Sample']}", alpha=0.5)
        
        plt.title('Junction Saturation Analysis')
        plt.xlabel('Number of Reads')
        plt.ylabel('Number of Junctions')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plots
        plt.savefig('output/figures/QC_Junction_Saturation.pdf', format='pdf', bbox_inches='tight')
        plt.savefig('output/figures/QC_Junction_Saturation.eps', format='eps', bbox_inches='tight')
        plt.savefig('output/figures/QC_Junction_Saturation.png', dpi=300, bbox_inches='tight')
        plt.close()

# 3. Generate Summary Statistics Table
def create_summary_statistics():
    # Select numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    
    # Calculate statistics
    stats = pd.DataFrame({
        'Mean': df[numeric_cols].mean(),
        'Std': df[numeric_cols].std(),
        'Min': df[numeric_cols].min(),
        'Max': df[numeric_cols].max()
    })
    
    # Format statistics
    stats['Mean ± Std'] = stats.apply(lambda x: f"{x['Mean']:.2f} ± {x['Std']:.2f}", axis=1)
    stats['Range'] = stats.apply(lambda x: f"{x['Min']:.2f} - {x['Max']:.2f}", axis=1)
    
    # Save to CSV
    stats[['Mean ± Std', 'Range']].to_csv('output/tables/QC_Summary_Statistics.csv')
    
    # Create a LaTeX table
    latex_table = stats[['Mean ± Std', 'Range']].to_latex(
        float_format=lambda x: f"{x:.2f}",
        escape=False,
        index=True,
        caption="Summary Statistics of QC Metrics",
        label="tab:qc_summary"
    )
    
    # Save LaTeX table
    with open('output/tables/QC_Summary_Statistics.tex', 'w') as f:
        f.write(latex_table)

def create_all_plots():
    create_gene_coverage_plot()
    create_junction_saturation_plot()
    create_summary_statistics()
    print("All supplementary plots and tables have been generated in output/figures/ and output/tables/")

if __name__ == "__main__":
    create_all_plots() 