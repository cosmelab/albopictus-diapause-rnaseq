#!/usr/bin/env python3
"""
Create nf-core/rnaseq samplesheet from our existing table data
"""
import pandas as pd
from pathlib import Path

def create_samplesheet(detailed_info_file, output_file, data_dir):
    """Create nf-core/rnaseq compatible samplesheet from our detailed sample information"""
    
    # Read detailed sample information
    df = pd.read_csv(detailed_info_file)
    
    # Create rows for samplesheet
    rows = []
    for _, row in df.iterrows():
        sample = row['Run']  # Use SRA Run ID as sample name
        
        # Determine strandedness (assume unstranded for now)
        strandedness = "unstranded"
        
        # Check for fastq files
        fastq_dir = Path(data_dir) / "raw" / row['BioProject']
        
        # Look for paired-end files first
        fastq_1 = fastq_dir / f"{sample}_1.fastq.gz"
        fastq_2 = fastq_dir / f"{sample}_2.fastq.gz"
        
        if fastq_2.exists():
            # Paired-end
            rows.append([sample, fastq_1, fastq_2, strandedness])
        else:
            # Single-end
            fastq_1 = fastq_dir / f"{sample}.fastq.gz"
            if fastq_1.exists():
                rows.append([sample, fastq_1, "", strandedness])
            else:
                print(f"Warning: No fastq files found for sample {sample}")
    
    # Create DataFrame
    samplesheet_df = pd.DataFrame(rows, columns=['sample', 'fastq_1', 'fastq_2', 'strandedness'])
    
    # Save samplesheet
    samplesheet_df.to_csv(output_file, index=False)
    print(f"Samplesheet saved to {output_file}")
    
    # Print summary
    print("\nSamplesheet Summary:")
    print(f"Total samples: {len(samplesheet_df)}")
    print(f"Paired-end samples: {len(samplesheet_df[samplesheet_df['fastq_2'] != ''])}")
    print(f"Single-end samples: {len(samplesheet_df[samplesheet_df['fastq_2'] == ''])}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python create_samplesheet.py <detailed_sample_info.csv> <output.csv> <data_dir>")
        sys.exit(1)
    
    create_samplesheet(sys.argv[1], sys.argv[2], sys.argv[3])
