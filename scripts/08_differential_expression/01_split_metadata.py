#!/usr/bin/env python3
"""
Split metadata by developmental stage for stage-specific analyses

This script splits the combined metadata file into three stage-specific files
to enable independent differential expression analysis for each developmental stage.
This is necessary due to complete confounding between platform and stage.

Author: L. Cosme
Date: October 2025
"""

import pandas as pd
import os
from pathlib import Path

# Set paths
BASE_DIR = Path("/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq")
METADATA_DIR = BASE_DIR / "data" / "metadata"
COUNT_DIR = BASE_DIR / "results" / "count_matrices"

def load_and_check_metadata():
    """Load the main metadata file and perform quality checks"""

    metadata_file = METADATA_DIR / "samples.csv"

    if not metadata_file.exists():
        print(f"Error: Metadata file not found at {metadata_file}")
        print("Creating from available information...")
        create_metadata_from_counts()

    # Load metadata
    metadata = pd.read_csv(metadata_file)

    print(f"Loaded metadata for {len(metadata)} samples")
    print(f"Columns: {list(metadata.columns)}")

    return metadata

def create_metadata_from_counts():
    """Create metadata file from count matrix column names if it doesn't exist"""

    # Try to get sample names from count matrix
    count_file = COUNT_DIR / "salmon_gene_counts.tsv"

    if not count_file.exists():
        print("Error: No count matrix found to infer metadata from")
        print("Please create data/metadata/samples.csv manually with columns:")
        print("Sample, BioProject, Stage, Platform, Condition, Run")
        return

    # Read count matrix header
    counts = pd.read_csv(count_file, sep='\t', nrows=1)
    samples = [col for col in counts.columns if col not in ['Name', 'Gene_ID', 'gene_id']]

    # Parse sample names to extract metadata
    metadata_list = []

    for sample in samples:
        # Expected format: PRJNAXXXXXX_SRRXXXXXXX or similar
        parts = sample.split('_')

        if len(parts) >= 2:
            bioproject = parts[0]
            run = parts[1] if parts[1].startswith('SRR') else parts[-1]
        else:
            bioproject = "Unknown"
            run = sample

        # Infer stage from BioProject
        stage_map = {
            'PRJNA268379': 'Adult',
            'PRJNA158021': 'Embryo',
            'PRJNA187045': 'Pharate_Larvae'
        }

        # Infer platform from BioProject
        platform_map = {
            'PRJNA268379': 'Illumina HiSeq 2000',
            'PRJNA158021': 'Illumina HiSeq',
            'PRJNA187045': 'Illumina Genome Analyzer IIx'
        }

        stage = stage_map.get(bioproject, 'Unknown')
        platform = platform_map.get(bioproject, 'Unknown')

        # Infer condition (this is approximate - should be verified)
        # You'll need to check the original study designs
        condition = 'Unknown'  # This needs to be filled in based on study design

        metadata_list.append({
            'Sample': sample,
            'BioProject': bioproject,
            'Stage': stage,
            'Platform': platform,
            'Condition': condition,
            'Run': run
        })

    # Create DataFrame
    metadata = pd.DataFrame(metadata_list)

    # Save
    output_file = METADATA_DIR / "samples.csv"
    metadata.to_csv(output_file, index=False)
    print(f"Created metadata file: {output_file}")
    print("WARNING: Condition column needs to be filled in manually!")
    print("Check original publications for diapause vs non-diapause assignments")

    return metadata

def split_metadata_by_stage(metadata):
    """Split metadata into stage-specific files"""

    # Check if Stage column exists
    if 'Stage' not in metadata.columns:
        print("Error: 'Stage' column not found in metadata")
        print("Available columns:", list(metadata.columns))
        return

    # Get unique stages
    stages = metadata['Stage'].unique()
    print(f"\nFound {len(stages)} developmental stages: {stages}")

    # Split by stage
    for stage in stages:
        # Filter for this stage
        stage_data = metadata[metadata['Stage'] == stage].copy()

        # Clean stage name for filename
        stage_clean = stage.lower().replace(' ', '_').replace('pharate_larvae', 'larvae')

        # Output filename
        output_file = METADATA_DIR / f"samples_{stage_clean}.csv"

        # Save
        stage_data.to_csv(output_file, index=False)

        print(f"\n{stage}:")
        print(f"  Samples: {len(stage_data)}")
        print(f"  BioProject: {stage_data['BioProject'].iloc[0] if len(stage_data) > 0 else 'N/A'}")
        print(f"  Platform: {stage_data['Platform'].iloc[0] if len(stage_data) > 0 else 'N/A'}")

        # Check condition distribution
        if 'Condition' in stage_data.columns:
            condition_counts = stage_data['Condition'].value_counts()
            print(f"  Conditions:")
            for condition, count in condition_counts.items():
                print(f"    - {condition}: {count} samples")

        print(f"  Saved to: {output_file}")

    return stages

def split_count_matrices(metadata):
    """Split count matrices by developmental stage"""

    print("\n" + "="*50)
    print("Splitting count matrices by stage...")
    print("="*50)

    # Load count matrix
    count_files = [
        COUNT_DIR / "salmon_gene_counts.tsv",
        COUNT_DIR.parent / "02_count_matrices" / "raw" / "salmon_counts_all.tsv"
    ]

    count_file = None
    for cf in count_files:
        if cf.exists():
            count_file = cf
            break

    if not count_file:
        print("Warning: No count matrix found to split")
        print("Looked for:", count_files)
        return

    print(f"Loading counts from: {count_file}")
    counts = pd.read_csv(count_file, sep='\t', index_col=0)
    print(f"Loaded {counts.shape[0]} genes × {counts.shape[1]} samples")

    # Get stages
    stages = metadata['Stage'].unique()

    for stage in stages:
        # Get samples for this stage
        stage_samples = metadata[metadata['Stage'] == stage]['Sample'].tolist()

        # Find matching columns in count matrix
        matching_cols = []
        for col in counts.columns:
            # Check if any stage sample matches this column
            for sample in stage_samples:
                if sample in col or col in sample:
                    matching_cols.append(col)
                    break

        if not matching_cols:
            print(f"Warning: No matching samples found for {stage}")
            continue

        # Filter count matrix
        stage_counts = counts[matching_cols]

        # Clean stage name for filename
        stage_clean = stage.lower().replace(' ', '_').replace('pharate_larvae', 'larvae')

        # Save filtered counts
        output_file = COUNT_DIR.parent / "02_count_matrices" / "raw" / f"{stage_clean}_counts.tsv"
        output_file.parent.mkdir(parents=True, exist_ok=True)

        stage_counts.to_csv(output_file, sep='\t')

        print(f"\n{stage} counts:")
        print(f"  Samples: {len(matching_cols)}")
        print(f"  Genes: {stage_counts.shape[0]}")
        print(f"  Saved to: {output_file}")

def create_sample_condition_mapping():
    """Create a reference file mapping samples to conditions based on original studies"""

    print("\n" + "="*50)
    print("Creating sample-condition mapping reference...")
    print("="*50)

    # This information comes from the original publications
    # You'll need to verify these assignments

    condition_map = {
        # Adults (PRJNA268379) - Huang et al. 2015
        'PRJNA268379': {
            'diapause': ['SRR1663527', 'SRR1663528', 'SRR1663529',  # Add actual SRR numbers
                        'SRR1663530', 'SRR1663531', 'SRR1663532',
                        'SRR1663533', 'SRR1663534', 'SRR1663535'],
            'non_diapause': ['SRR1663536', 'SRR1663537', 'SRR1663538',
                            'SRR1663539', 'SRR1663540', 'SRR1663541',
                            'SRR1663542']
        },
        # Embryos (PRJNA158021) - Poelchau et al. 2013a
        'PRJNA158021': {
            'diapause': [],  # Add actual SRR numbers
            'non_diapause': []
        },
        # Larvae (PRJNA187045) - Poelchau et al. 2013b
        'PRJNA187045': {
            'diapause': [],  # Add actual SRR numbers
            'non_diapause': []
        }
    }

    # Save as reference
    reference_file = METADATA_DIR / "sample_condition_reference.txt"

    with open(reference_file, 'w') as f:
        f.write("# Sample-Condition Mapping Reference\n")
        f.write("# This file needs to be completed with actual SRR accessions\n")
        f.write("# Check original publications for correct assignments\n\n")

        for bioproject, conditions in condition_map.items():
            f.write(f"\n## {bioproject}\n")
            for condition, samples in conditions.items():
                f.write(f"\n{condition}:\n")
                for sample in samples:
                    f.write(f"  - {sample}\n")

    print(f"Created reference file: {reference_file}")
    print("NOTE: This file needs to be completed with actual sample-condition mappings")

def main():
    """Main execution function"""

    print("="*50)
    print("Metadata Splitting for Stage-Specific Analysis")
    print("="*50)

    # Create metadata directory if it doesn't exist
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    # Load metadata
    metadata = load_and_check_metadata()

    if metadata is None:
        print("\nError: Could not load or create metadata")
        print("Please create data/metadata/samples.csv manually")
        return

    # Split metadata by stage
    stages = split_metadata_by_stage(metadata)

    # Split count matrices
    split_count_matrices(metadata)

    # Create reference file
    create_sample_condition_mapping()

    print("\n" + "="*50)
    print("Summary")
    print("="*50)

    print("\nFiles created:")
    for stage in ['adults', 'embryos', 'larvae']:
        meta_file = METADATA_DIR / f"samples_{stage}.csv"
        count_file = COUNT_DIR.parent / "02_count_matrices" / "raw" / f"{stage}_counts.tsv"

        print(f"\n{stage.title()}:")
        print(f"  Metadata: {meta_file.exists()} - {meta_file if meta_file.exists() else 'NOT FOUND'}")
        print(f"  Counts: {count_file.exists()} - {count_file if count_file.exists() else 'NOT FOUND'}")

    print("\n✓ Metadata splitting complete")
    print("\nNext steps:")
    print("1. Verify sample-condition assignments in metadata files")
    print("2. Run DESeq2 on adults: Rscript scripts/06_differential_expression/03_adults_deseq2.R")
    print("3. Check results in results/03_differential_expression/adults/")

if __name__ == "__main__":
    main()