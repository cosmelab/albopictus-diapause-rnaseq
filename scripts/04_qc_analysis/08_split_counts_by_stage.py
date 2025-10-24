#!/usr/bin/env python3
"""
Split Count Matrices by Developmental Stage

Purpose: Split the combined Salmon count matrix into stage-specific matrices
         (adults, embryos, larvae) for stage-specific differential expression

Input:
    - results/count_matrices/salmon_gene_counts.tsv
    - data/metadata/samplesheet.csv

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
print("=" * 70)
print()

# Load sample metadata
print("Loading sample metadata...")
samplesheet = pd.read_csv(METADATA_DIR / "samplesheet.csv")
print(f"  Total samples: {len(samplesheet)}")

# Map BioProject to stage
# Based on the project documentation:
# PRJNA158021 = Embryos (12 samples)
# PRJNA187045 = Larvae (17 samples)
# PRJNA268379 = Adults (16 samples)

def get_stage(bioproject):
    """Map BioProject to developmental stage"""
    stage_map = {
        'PRJNA158021': 'embryo',
        'PRJNA187045': 'larvae',
        'PRJNA268379': 'adult'
    }
    return stage_map.get(bioproject, 'unknown')

# Add stage column
if 'stage' not in samplesheet.columns:
    samplesheet['bioproject'] = samplesheet['sample'].str.split('_').str[0]
    samplesheet['stage'] = samplesheet['bioproject'].apply(get_stage)

# Verify stages
print("\nSamples per stage:")
for stage, count in samplesheet['stage'].value_counts().items():
    print(f"  {stage.capitalize()}: {count}")
print()

# Load count matrix
print("Loading Salmon gene count matrix...")
counts = pd.read_csv(COUNT_DIR / "salmon_gene_counts.tsv", sep="\t", index_col=0)
print(f"  Shape: {counts.shape}")
print()

# Split by stage
stages = ['adult', 'embryo', 'larvae']

for stage in stages:
    print(f"Processing {stage.capitalize()} samples...")

    # Get samples for this stage
    stage_samples = samplesheet[samplesheet['stage'] == stage]
    sample_ids = stage_samples['sample'].tolist()

    print(f"  Samples: {len(sample_ids)}")

    if len(sample_ids) == 0:
        print(f"  WARNING: No samples found for stage: {stage}")
        continue

    # Extract counts for these samples
    stage_counts = counts[sample_ids]

    # Save stage-specific count matrix
    out_file = COUNT_DIR / f"{stage}s_counts.tsv"
    stage_counts.to_csv(out_file, sep="\t")
    print(f"  Saved counts: {out_file}")
    print(f"    Shape: {stage_counts.shape}")

    # Prepare metadata for DESeq2
    # Need columns: sample, condition (diapause/non-diapause)
    stage_meta = stage_samples[['sample']].copy()

    # Parse condition from sample name or use existing column
    if 'condition' in stage_samples.columns:
        stage_meta['condition'] = stage_samples['condition'].values
    else:
        # Infer from sample name if not present
        # This needs to be verified with actual data
        stage_meta['condition'] = 'diapause'  # Placeholder

    # Add other useful columns
    if 'strandedness' in stage_samples.columns:
        stage_meta['strandedness'] = stage_samples['strandedness'].values

    stage_meta['stage'] = stage
    stage_meta['bioproject'] = stage_samples['bioproject'].values

    # Save metadata
    out_file = COUNT_DIR / f"{stage}s_metadata.csv"
    stage_meta.to_csv(out_file, index=False)
    print(f"  Saved metadata: {out_file}")
    print(f"    Columns: {', '.join(stage_meta.columns)}")
    print()

# Summary
print("=" * 70)
print("Summary")
print("=" * 70)
print(f"Count matrices split into {len(stages)} developmental stages")
print(f"Output directory: {COUNT_DIR}")
print()
print("Stage-specific files created:")
for stage in stages:
    count_file = COUNT_DIR / f"{stage}s_counts.tsv"
    meta_file = COUNT_DIR / f"{stage}s_metadata.csv"
    if count_file.exists():
        print(f"  {stage.capitalize()}:")
        print(f"    Counts: {count_file}")
        print(f"    Metadata: {meta_file}")
print("=" * 70)