#!/usr/bin/env python3
"""
Combine and Organize Count Matrices by Stage and Comparison

This script reads featureCounts and Salmon outputs from the nf-core/rnaseq pipeline
and creates organized count matrices split by:
  - Life stage (adults, embryos, larvae)
  - Biological comparison (unfed/fed for adults, timepoints for embryos/larvae)

Input:
  - featureCounts: output/star_salmon/featurecounts/*.featureCounts.tsv (45 files)
  - Salmon: output/star_salmon/*/quant.sf (45 directories)
  - Metadata: data/metadata/samples.csv

Output:
  - Organized matrices in output/count_matrices/featurecounts/ and output/count_matrices/salmon/
  - Metadata files for each comparison group
  - Summary log file

Created: October 27, 2025
Author: Automated analysis pipeline
"""

import pandas as pd
import glob
import os
from pathlib import Path
import sys

# =============================================================================
# Configuration
# =============================================================================

METADATA_FILE = "data/metadata/samples.csv"
FEATURECOUNTS_DIR = "output/star_salmon/featurecounts"
SALMON_DIR = "output/star_salmon"
OUTPUT_BASE = "output/count_matrices"

# Sample group definitions based on experimental design
SAMPLE_GROUPS = {
    'adults': {
        'bioproject': 'PRJNA268379',
        'unfed': ['SD_Unfed', 'LD_Unfed'],
        'fed': ['SD_Fed', 'LD_Fed'],
        'all': ['SD_Unfed', 'LD_Unfed', 'SD_Fed', 'LD_Fed']
    },
    'embryos': {
        'bioproject': 'PRJNA158021',
        '72h': ['SD_72h', 'LD_72h'],
        '135h': ['SD_135h', 'LD_135h'],
        'SD_timeseries': ['SD_72h', 'SD_135h'],
        'LD_timeseries': ['LD_72h', 'LD_135h'],
        'all': ['SD_72h', 'LD_72h', 'SD_135h', 'LD_135h']
    },
    'larvae': {
        'bioproject': 'PRJNA187045',
        '11d': ['SD_11d', 'LD_11d'],
        'SD_timeseries': ['SD_11d', 'SD_21d', 'SD_40d'],
        'all': ['SD_11d', 'LD_11d', 'SD_21d', 'LD_21d', 'SD_40d', 'LD_40d']
    }
}

# =============================================================================
# Helper Functions
# =============================================================================

def load_metadata():
    """Load and validate metadata file."""
    print(f"\n{'='*70}")
    print("Loading metadata...")
    print(f"{'='*70}")

    if not os.path.exists(METADATA_FILE):
        print(f"ERROR: Metadata file not found: {METADATA_FILE}")
        sys.exit(1)

    metadata = pd.read_csv(METADATA_FILE)

    # Validate required columns
    required_cols = ['BioProject', 'Run', 'New abreviation', 'Sample collection', 'Photoperiod']
    missing = [col for col in required_cols if col not in metadata.columns]
    if missing:
        print(f"ERROR: Missing required columns in metadata: {missing}")
        sys.exit(1)

    print(f"Loaded metadata: {len(metadata)} samples")
    print(f"BioProjects: {', '.join(metadata['BioProject'].unique())}")
    print(f"Sample IDs: {', '.join(sorted(metadata['New abreviation'].unique()))}")

    return metadata


def read_featurecounts_file(filepath):
    """Read a featureCounts output file and extract gene counts."""
    # featureCounts format: Geneid, Chr, Start, End, Strand, Length, counts
    df = pd.read_csv(filepath, sep='\t', comment='#')

    # Extract just Geneid and the count column (last column)
    gene_id = df.iloc[:, 0]  # First column is Geneid
    counts = df.iloc[:, -1]  # Last column is counts

    return pd.DataFrame({'gene_id': gene_id, 'counts': counts})


def read_salmon_file(quant_dir):
    """Read a Salmon quant.sf file and aggregate to gene level."""
    quant_file = os.path.join(quant_dir, 'quant.sf')

    if not os.path.exists(quant_file):
        print(f"WARNING: Salmon file not found: {quant_file}")
        return None

    # Salmon format: Name, Length, EffectiveLength, TPM, NumReads
    df = pd.read_csv(quant_file, sep='\t')

    # For now, return transcript-level counts (NumReads)
    # We'll aggregate to gene level using tximport in R later
    # For this script, we'll sum transcripts to genes using gene_id from transcript name

    # Extract gene_id from transcript name (format: transcript_id|gene_id)
    # If no pipe, use the whole name as gene_id
    if '|' in df['Name'].iloc[0]:
        df['gene_id'] = df['Name'].str.split('|').str[1]
    else:
        # Try to match with featureCounts gene IDs
        # For AalbF3, transcripts are named like gene-LOC109400916-RA
        # Gene IDs are like LOC109400916
        df['gene_id'] = df['Name'].str.extract(r'(LOC\d+)')[0]

        # If no match, just use the Name
        df['gene_id'] = df['gene_id'].fillna(df['Name'])

    # Sum counts per gene
    gene_counts = df.groupby('gene_id')['NumReads'].sum().reset_index()
    gene_counts.columns = ['gene_id', 'counts']

    return gene_counts


def create_count_matrix(metadata, sample_ids, count_files, count_type='featurecounts'):
    """
    Create a count matrix for specified samples.

    Args:
        metadata: Full metadata DataFrame
        sample_ids: List of sample IDs (e.g., ['SD_Unfed', 'LD_Unfed'])
        count_files: Dict mapping BioProject_Run to file path
        count_type: 'featurecounts' or 'salmon'

    Returns:
        count_matrix: DataFrame with genes as rows, samples as columns
        sample_metadata: Metadata for included samples
    """
    # Filter metadata to samples of interest
    sample_meta = metadata[metadata['New abreviation'].isin(sample_ids)].copy()

    if len(sample_meta) == 0:
        print(f"WARNING: No samples found for {sample_ids}")
        return None, None

    print(f"  Creating matrix for {len(sample_meta)} samples: {', '.join(sample_meta['New abreviation'].tolist())}")

    # Initialize count matrix
    all_counts = {}
    gene_ids = None

    for _, row in sample_meta.iterrows():
        sample_key = f"{row['BioProject']}_{row['Run']}"

        if sample_key not in count_files:
            print(f"  WARNING: Count file not found for {sample_key}")
            continue

        # Read count file
        if count_type == 'featurecounts':
            count_data = read_featurecounts_file(count_files[sample_key])
        else:  # salmon
            count_data = read_salmon_file(count_files[sample_key])

        if count_data is None:
            continue

        # Store counts with unique sample ID (BioProject_Run)
        sample_name = f"{row['BioProject']}_{row['Run']}"
        all_counts[sample_name] = count_data.set_index('gene_id')['counts']

        # Store gene IDs from first sample
        if gene_ids is None:
            gene_ids = count_data['gene_id'].tolist()

    if len(all_counts) == 0:
        print(f"  ERROR: No count data loaded")
        return None, None

    # Create DataFrame
    count_matrix = pd.DataFrame(all_counts)

    # Sort columns by sample name
    count_matrix = count_matrix[sorted(count_matrix.columns)]

    print(f"  Matrix shape: {count_matrix.shape} (genes x samples)")

    return count_matrix, sample_meta


def save_outputs(count_matrix, sample_meta, output_dir, comparison_name, count_type):
    """Save count matrix and metadata to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save count matrix
    matrix_file = output_dir / f"{comparison_name}_counts.tsv"
    count_matrix.to_csv(matrix_file, sep='\t')
    print(f"  Saved: {matrix_file}")

    # Save metadata
    meta_file = output_dir / f"{comparison_name}_metadata.tsv"
    sample_meta.to_csv(meta_file, sep='\t', index=False)
    print(f"  Saved: {meta_file}")

    return matrix_file, meta_file


# =============================================================================
# Main Processing Functions
# =============================================================================

def process_featurecounts(metadata):
    """Process all featureCounts files and create organized matrices."""
    print(f"\n{'='*70}")
    print("Processing featureCounts")
    print(f"{'='*70}")

    output_dir = Path(OUTPUT_BASE) / 'featurecounts'

    # Find all featureCounts files
    count_files = {}
    for filepath in glob.glob(f"{FEATURECOUNTS_DIR}/*.featureCounts.tsv"):
        filename = os.path.basename(filepath)
        # Extract BioProject_Run (e.g., PRJNA268379_SRR1663685)
        sample_key = filename.replace('.featureCounts.tsv', '')
        count_files[sample_key] = filepath

    print(f"Found {len(count_files)} featureCounts files")

    if len(count_files) == 0:
        print("ERROR: No featureCounts files found")
        return

    # Process each life stage
    for stage, groups in SAMPLE_GROUPS.items():
        print(f"\n{'-'*70}")
        print(f"Processing {stage.upper()}")
        print(f"{'-'*70}")

        stage_dir = output_dir / stage

        # Create matrices for each comparison group
        for group_name, sample_ids in groups.items():
            if group_name == 'bioproject':
                continue

            print(f"\nComparison: {group_name}")

            # Create count matrix
            count_matrix, sample_meta = create_count_matrix(
                metadata, sample_ids, count_files, count_type='featurecounts'
            )

            if count_matrix is not None:
                save_outputs(count_matrix, sample_meta, stage_dir, group_name, 'featurecounts')


def process_salmon(metadata):
    """Process all Salmon quantification files and create organized matrices."""
    print(f"\n{'='*70}")
    print("Processing Salmon")
    print(f"{'='*70}")

    output_dir = Path(OUTPUT_BASE) / 'salmon'

    # Find all Salmon directories
    count_files = {}
    for quant_dir in glob.glob(f"{SALMON_DIR}/PRJNA*_SRR*/"):
        dirname = os.path.basename(quant_dir.rstrip('/'))
        # dirname is like PRJNA268379_SRR1663685
        count_files[dirname] = quant_dir

    print(f"Found {len(count_files)} Salmon quantification directories")

    if len(count_files) == 0:
        print("ERROR: No Salmon directories found")
        return

    # Process each life stage
    for stage, groups in SAMPLE_GROUPS.items():
        print(f"\n{'-'*70}")
        print(f"Processing {stage.upper()}")
        print(f"{'-'*70}")

        stage_dir = output_dir / stage

        # Create matrices for each comparison group
        for group_name, sample_ids in groups.items():
            if group_name == 'bioproject':
                continue

            print(f"\nComparison: {group_name}")

            # Create count matrix
            count_matrix, sample_meta = create_count_matrix(
                metadata, sample_ids, count_files, count_type='salmon'
            )

            if count_matrix is not None:
                save_outputs(count_matrix, sample_meta, stage_dir, group_name, 'salmon')


def print_summary(metadata):
    """Print final summary of created files."""
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    # Count created files
    fc_files = list(Path(OUTPUT_BASE).glob('featurecounts/**/*.tsv'))
    salmon_files = list(Path(OUTPUT_BASE).glob('salmon/**/*.tsv'))

    print(f"\nfeatureCounts matrices created: {len([f for f in fc_files if 'counts.tsv' in f.name])}")
    print(f"featureCounts metadata files: {len([f for f in fc_files if 'metadata.tsv' in f.name])}")

    print(f"\nSalmon matrices created: {len([f for f in salmon_files if 'counts.tsv' in f.name])}")
    print(f"Salmon metadata files: {len([f for f in salmon_files if 'metadata.tsv' in f.name])}")

    print(f"\nOutput directory: {OUTPUT_BASE}")
    print("\nNext step: Run 03_sanity_check_counts.py to validate outputs")


# =============================================================================
# Main
# =============================================================================

def main():
    """Main execution function."""
    print("\n" + "="*70)
    print("COMBINE AND ORGANIZE COUNT MATRICES")
    print("="*70)
    print("\nThis script creates organized count matrices split by:")
    print("  - Life stage (adults, embryos, larvae)")
    print("  - Biological comparison (unfed/fed, timepoints)")
    print("\nBoth featureCounts and Salmon counts will be processed.")

    # Load metadata
    metadata = load_metadata()

    # Process featureCounts
    process_featurecounts(metadata)

    # Process Salmon
    process_salmon(metadata)

    # Print summary
    print_summary(metadata)

    print("\nDone!")


if __name__ == "__main__":
    main()
