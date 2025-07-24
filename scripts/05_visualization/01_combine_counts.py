#!/usr/bin/env python3

import pandas as pd
import glob
import os
from pathlib import Path

def read_metadata():
    """Read and process metadata files to create sample mapping."""
    # Read the metadata file
    metadata = pd.read_csv('data/metadata/samples.csv')
    
    # Create a mapping dictionary for sample IDs
    sample_mapping = {}
    for _, row in metadata.iterrows():
        key = f"{row['BioProject']}_{row['Run']}"
        sample_mapping[key] = f"{row['BioProject']}_{row['Run']}_{row['New abreviation']}"
    
    return sample_mapping

def process_counts(count_type='gene'):
    """
    Process count files for either genes or transcripts.
    
    Args:
        count_type (str): Either 'gene' or 'transcript'
    """
    # Create output directory if it doesn't exist
    output_dir = Path('output/tables')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get sample mapping
    sample_mapping = read_metadata()
    
    # Dictionary to store dataframes for each project
    project_dfs = {}
    all_counts = {}
    
    # Process each project directory
    for project_dir in glob.glob(f'results/zenodo/PRJNA*/{count_type}_counts'):
        project_id = os.path.basename(os.path.dirname(project_dir))
        project_dfs[project_id] = {}
        
        # Process each count file
        for count_file in glob.glob(f'{project_dir}/*_{count_type}_counts.tsv'):
            # Extract sample identifier
            sample_id = os.path.basename(count_file).replace(f'_{count_type}_counts.tsv', '')
            
            # Read the count file
            df = pd.read_csv(count_file, sep='\t')
            
            # Get the mapping key
            mapping_key = '_'.join(sample_id.split('_')[:2])
            
            # Get the column name for this sample
            column_name = sample_mapping.get(mapping_key, 'NA')
            
            # Store the counts (using the last column which contains the counts)
            project_dfs[project_id][column_name] = df.iloc[:, -1]
            all_counts[column_name] = df.iloc[:, -1]
            
            # Store the IDs from the first file
            if 'gene_id' not in project_dfs[project_id]:
                project_dfs[project_id]['gene_id'] = df['gene_id']
                # For transcripts, also store tx
                if count_type == 'transcript':
                    project_dfs[project_id]['tx'] = df['tx']
            if 'gene_id' not in all_counts:
                all_counts['gene_id'] = df['gene_id']
                # For transcripts, also store tx
                if count_type == 'transcript':
                    all_counts['tx'] = df['tx']
    
    # Create and save project-specific tables
    for project_id, counts in project_dfs.items():
        # Create DataFrame
        project_df = pd.DataFrame(counts)
        
        # Set index based on count type
        if count_type == 'transcript':
            project_df.set_index(['tx', 'gene_id'], inplace=True)
        else:
            project_df.set_index('gene_id', inplace=True)
        
        # Save to CSV
        output_file = output_dir / f'{project_id}_{count_type}_counts.csv'
        project_df.to_csv(output_file)
        print(f"Saved {project_id} {count_type} counts to {output_file}")
    
    # Create and save combined table
    combined_df = pd.DataFrame(all_counts)
    
    # Set index based on count type
    if count_type == 'transcript':
        combined_df.set_index(['tx', 'gene_id'], inplace=True)
    else:
        combined_df.set_index('gene_id', inplace=True)
    
    output_file = output_dir / f'combined_{count_type}_counts.csv'
    combined_df.to_csv(output_file)
    print(f"Saved combined {count_type} counts to {output_file}")
    
    # Print summary
    print(f"\nSummary of {count_type} counts:")
    print(f"Total number of samples: {len(all_counts) - (2 if count_type == 'transcript' else 1)}")
    print(f"Total number of {count_type}s: {len(combined_df)}")
    print("\nSamples included:")
    for sample in sorted(all_counts.keys()):
        if sample not in ['gene_id', 'tx']:
            print(f"- {sample}")

def main():
    # Process gene counts
    print("Processing gene counts...")
    process_counts(count_type='gene')
    
    # Process transcript counts
    print("\nProcessing transcript counts...")
    process_counts(count_type='transcript')

if __name__ == "__main__":
    main() 