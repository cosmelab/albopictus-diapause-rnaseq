#!/usr/bin/env python3

import pandas as pd
import glob
import os

def read_general_stats(file_path):
    """Read and process general stats file."""
    df = pd.read_csv(file_path, sep='\t')
    return df

def read_read_distribution(file_path):
    """Read and process read distribution file."""
    df = pd.read_csv(file_path, sep='\t')
    return df

def create_qc_table():
    # Read sample information table
    sample_info = pd.read_csv('output/tables/Table_S1_Detailed_Sample_Information.csv')
    # Create a dictionary mapping SRR to Sample ID
    srr_to_sample = dict(zip(sample_info['Run'], sample_info['Sample ID']))
    
    # Initialize empty list to store all data
    all_data = []
    
    # Process each project directory
    for project_dir in glob.glob('results/zenodo/PRJNA*/multiqc'):
        project_id = os.path.basename(os.path.dirname(project_dir))
        
        # Process each sample's general stats
        for stats_file in glob.glob(f'{project_dir}/*_multiqc_general_stats.txt'):
            full_sample = os.path.basename(stats_file).split('_multiqc_')[0]
            srr = full_sample.split('_')[1]  # Extract SRR number
            
            # Read general stats
            stats_df = read_general_stats(stats_file)
            
            # Get main sample row (not _1 or _2)
            main_row = stats_df[~stats_df['Sample'].str.contains('_1|_2')].iloc[0]
            
            # Get metrics from _1 row (first read)
            read1_row = stats_df[stats_df['Sample'].str.endswith('_1')].iloc[0]
            
            # Get metrics from _2 row (second read)
            read2_row = stats_df[stats_df['Sample'].str.endswith('_2')].iloc[0]
            
            # Calculate average metrics from read1 and read2
            avg_duplication = (read1_row['fastqc_raw-percent_duplicates'] + read2_row['fastqc_raw-percent_duplicates']) / 2
            avg_adapter = (read1_row['cutadapt-percent_trimmed'] + read2_row['cutadapt-percent_trimmed']) / 2
            avg_read_length = (read1_row['fastqc_raw-avg_sequence_length'] + read2_row['fastqc_raw-avg_sequence_length']) / 2
            
            # Calculate multi-mapped percentage
            multi_mapped_percent = main_row['star-mapped_percent'] - main_row['star-uniquely_mapped_percent']
            
            # Extract metrics
            sample_data = {
                'Project': project_id,
                'SRR': srr,
                'Sample ID': srr_to_sample.get(srr, 'NA'),  # Get Sample ID from mapping, use 'NA' if not found
                # Sequencing Quality
                'Total_Reads (M)': main_row['star-total_reads'],
                'Read_Length (bp)': avg_read_length,
                'Duplication_Rate (%)': avg_duplication,
                'Adapter_Contamination (%)': avg_adapter,
                
                # Alignment Quality
                'Mapping_Rate (%)': main_row['star-mapped_percent'],
                'Uniquely_Mapped (%)': main_row['star-uniquely_mapped_percent'],
                'Multi_Mapped (%)': multi_mapped_percent,
                
                # RNA-seq Specific Metrics
                'Exonic_Reads (%)': main_row['qualimap_rnaseq-reads_aligned_exonic'],
                'Intronic_Reads (%)': main_row['qualimap_rnaseq-reads_aligned_intronic'],
                'Intergenic_Reads (%)': main_row['qualimap_rnaseq-reads_aligned_intergenic'],
                'Insert_Size (bp)': main_row['samtools_stats-insert_size_average'],
                
                # Technical Assessment
                'Properly_Paired (%)': main_row['samtools_stats-reads_properly_paired_percent'],
                'Error_Rate (%)': main_row['samtools_stats-error_rate'],
                'Five_Three_Bias (ratio)': main_row['qualimap_rnaseq-5_3_bias']
            }
            
            all_data.append(sample_data)
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(all_data)
    
    # Round numeric columns to 2 decimal places
    numeric_cols = df.select_dtypes(include=['float64']).columns
    df[numeric_cols] = df[numeric_cols].round(2)
    
    # Save to CSV in output/tables directory
    output_file = 'output/tables/Table_S2_QC_Metrics.csv'
    os.makedirs('output/tables', exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"QC summary table saved to {output_file}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(df.groupby('Project')[['Total_Reads (M)', 'Mapping_Rate (%)', 'Uniquely_Mapped (%)', 'Multi_Mapped (%)']].mean().round(2))

if __name__ == "__main__":
    create_qc_table() 