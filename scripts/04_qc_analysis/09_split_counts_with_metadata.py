#!/usr/bin/env python3
"""
Split Count Matrices by Developmental Stage (Using Proper Metadata)

Purpose: Split Salmon count matrix into stage-specific matrices with correct condition info

Input:
    - results/count_matrices/salmon_gene_counts.tsv
    - data/metadata/samples_adults.csv
    - data/metadata/samples_embryos.csv
    - data/metadata/samples_larvae.csv

Output:
    - results/count_matrices/adults_counts.tsv
    - results/count_matrices/embryos_counts.tsv
    - results/count_matrices/larvae_counts.tsv
    - results/count_matrices/adults_metadata.csv
    - results/count_matrices/embryos_metadata.csv
    - results/count_matrices/larvae_metadata.csv

Author: LCosme
Date: October 21, 2025
"""

import pandas as pd
from pathlib import Path

# Project paths
PROJECT_BASE = Path("/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq")
COUNT_DIR = PROJECT_BASE / "results" / "count_matrices"
METADATA_DIR = PROJECT_BASE / "data" / "metadata"

print("=" * 70)
print("Splitting Count Matrices by Developmental Stage")
print("(Using Proper Metadata Files)")
print("=" * 70)
print()

# Load count matrix
print("Loading Salmon gene count matrix...")
counts = pd.read_csv(COUNT_DIR / "salmon_gene_counts.tsv", sep="\t", index_col=0)
print(f"  Shape: {counts.shape}")
print()

# Process each stage
stages = {
    'adults': 'samples_adults.csv',
    'embryos': 'samples_embryos.csv',
    'larvae': 'samples_larvae.csv'
}

for stage_name, metadata_file in stages.items():
    print(f"Processing {stage_name.capitalize()}...")

    # Load stage-specific metadata
    metadata_path = METADATA_DIR / metadata_file
    if not metadata_path.exists():
        print(f"  WARNING: Metadata file not found: {metadata_file}")
        continue

    metadata = pd.read_csv(metadata_path)
    print(f"  Loaded metadata: {len(metadata)} samples")

    # Check required columns
    if 'Sample' not in metadata.columns:
        print(f"  ERROR: 'Sample' column not found in {metadata_file}")
        continue

    if 'Condition' not in metadata.columns:
        print(f"  ERROR: 'Condition' column not found in {metadata_file}")
        continue

    # Get sample IDs
    sample_ids = metadata['Sample'].tolist()

    # Verify samples exist in count matrix
    missing_samples = [s for s in sample_ids if s not in counts.columns]
    if missing_samples:
        print(f"  WARNING: {len(missing_samples)} samples not found in count matrix:")
        for s in missing_samples[:5]:
            print(f"    - {s}")

    available_samples = [s for s in sample_ids if s in counts.columns]
    print(f"  Available samples: {len(available_samples)}")

    if len(available_samples) == 0:
        print(f"  ERROR: No samples found for {stage_name}")
        continue

    # Extract counts for these samples
    stage_counts = counts[available_samples]

    # Save stage-specific count matrix
    out_file = COUNT_DIR / f"{stage_name}_counts.tsv"
    stage_counts.to_csv(out_file, sep="\t")
    print(f"  Saved counts: {out_file.name}")
    print(f"    Shape: {stage_counts.shape}")

    # Prepare metadata for DESeq2
    # Keep relevant columns and rename for consistency
    deseq_meta = metadata[metadata['Sample'].isin(available_samples)].copy()

    # Simplify condition labels for DESeq2
    # Map from long descriptive names to simple labels
    def simplify_condition(cond):
        if pd.isna(cond):
            return 'unknown'
        cond_lower = str(cond).lower()
        # Check for "non-diapause" FIRST (before checking for "diapause")
        if 'non-diapause' in cond_lower or 'non_diapause' in cond_lower:
            return 'non_diapause'
        elif 'diapause' in cond_lower:
            return 'diapause'
        else:
            return cond

    deseq_meta['condition'] = deseq_meta['Condition'].apply(simplify_condition)

    # Keep useful columns
    useful_cols = ['Sample', 'condition']

    # Add other metadata if available
    optional_cols = ['Run', 'Photoperiod', 'Temperature (oC)', 'Replicate',
                     'BioProject', 'Plattaform', 'Read length']

    for col in optional_cols:
        if col in deseq_meta.columns:
            useful_cols.append(col)

    # Create final metadata
    final_meta = deseq_meta[useful_cols].copy()
    final_meta.rename(columns={'Sample': 'sample'}, inplace=True)

    # Save metadata
    out_file = COUNT_DIR / f"{stage_name}_metadata.csv"
    final_meta.to_csv(out_file, index=False)
    print(f"  Saved metadata: {out_file.name}")

    # Show condition summary
    print(f"  Condition summary:")
    for cond, count in final_meta['condition'].value_counts().items():
        print(f"    {cond}: {count}")
    print()

# Summary
print("=" * 70)
print("Summary")
print("=" * 70)
print(f"Count matrices split into {len(stages)} developmental stages")
print(f"Output directory: {COUNT_DIR}")
print()
print("Stage-specific files created:")
for stage_name in stages.keys():
    count_file = COUNT_DIR / f"{stage_name}_counts.tsv"
    meta_file = COUNT_DIR / f"{stage_name}_metadata.csv"
    if count_file.exists():
        print(f"  {stage_name.capitalize()}:")
        print(f"    Counts: {count_file.name}")
        print(f"    Metadata: {meta_file.name}")
print("=" * 70)