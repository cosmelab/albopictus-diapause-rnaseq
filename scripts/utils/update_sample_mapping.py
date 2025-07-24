#!/usr/bin/env python3
"""
Update sample mapping files to include phenotype information and ensure consistent columns
"""

import pandas as pd
import os

def add_phenotype_column(file_path):
    """Add phenotype column to sample mapping file"""
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Add phenotype column based on condition
    if 'PRJNA158021' in file_path:  # Embryos
        df['phenotype'] = df['condition'].apply(lambda x: 'diapause-inducing' if 'Diapause inducing' in x else 'non-diapause-inducing')
        # Add experimental details
        df['photoperiod'] = '16L:8D'  # From PMC3797908
        df['temperature'] = '21°C'
        df['humidity'] = '80%'
        df['time_point'] = df['life_stage'].apply(lambda x: x.split(',')[0].strip())
    elif 'PRJNA187045' in file_path:  # Pharate Larvae
        df['phenotype'] = df['condition'].apply(lambda x: 'diapause-inducing' if x.startswith('D_') else 'non-diapause-inducing')
        # Already has experimental details
    elif 'PRJNA268379' in file_path:  # Adult Females
        df['phenotype'] = df['condition'].apply(lambda x: 'diapause-inducing' if 'BM' in x else 'non-diapause-inducing')
        # Add experimental details
        df['photoperiod'] = df['condition'].apply(lambda x: '16L:8D' if 'LD' in x else '8L:16D')
        df['temperature'] = '25°C'  # Standard lab conditions
        df['humidity'] = '80%'
        df['time_point'] = 'Adult'
    
    # Ensure consistent column order
    columns = ['accession', 'dataset', 'condition', 'life_stage', 'description', 
              'fastq_files', 'sequencing_type', 'phenotype', 'photoperiod', 
              'temperature', 'humidity', 'time_point']
    
    # Reorder columns and fill missing values with NA
    df = df.reindex(columns=columns, fill_value='NA')
    
    # Save back to the same file
    df.to_csv(file_path, index=False)
    print(f"Updated file: {file_path}")

def main():
    # Update all sample mapping files
    mapping_files = [
        'data/sample_mapping_PRJNA158021.csv',
        'data/sample_mapping_PRJNA187045.csv',
        'data/sample_mapping_PRJNA268379.csv'
    ]
    
    for file_path in mapping_files:
        if os.path.exists(file_path):
            print(f"\nProcessing {file_path}...")
            add_phenotype_column(file_path)
        else:
            print(f"File not found: {file_path}")

if __name__ == "__main__":
    main() 