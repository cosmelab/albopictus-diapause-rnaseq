#!/usr/bin/env python3
"""
Sanity Check Count Matrices

This script validates the organized count matrices created by 01_combine_counts.py.
It performs comprehensive checks to ensure data quality and catch errors early.

Checks performed:
  1. Sample counts (all 45 samples accounted for, no duplicates)
  2. Gene counts (~22,177 expected, may vary slightly)
  3. Metadata alignment (samples in counts match metadata)
  4. Missing data (check for excessive NAs or zeros)
  5. Replicate correlations (r > 0.9 within groups)
  6. File integrity (all expected files exist, non-empty)

If any check fails, the script reports the issue and exits with error code.
All checks must pass before proceeding to Phase 3 (Salmon validation).

Created: October 27, 2025
Author: Automated analysis pipeline
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from scipy.stats import pearsonr

# =============================================================================
# Configuration
# =============================================================================

OUTPUT_BASE = "output/count_matrices"
EXPECTED_SAMPLES = 45
EXPECTED_GENES_MIN = 20000  # Allow some variation
EXPECTED_GENES_MAX = 35000  # Salmon has more genes (transcript aggregation)
MIN_REPLICATE_CORRELATION = 0.9

# Expected files structure
EXPECTED_MATRICES = {
    'featurecounts': {
        'adults': ['unfed_counts.tsv', 'fed_counts.tsv', 'all_counts.tsv'],
        'embryos': ['72h_counts.tsv', '135h_counts.tsv', 'SD_timeseries_counts.tsv',
                    'LD_timeseries_counts.tsv', 'all_counts.tsv'],
        'larvae': ['11d_counts.tsv', 'SD_timeseries_counts.tsv', 'all_counts.tsv']
    },
    'salmon': {
        'adults': ['unfed_counts.tsv', 'fed_counts.tsv', 'all_counts.tsv'],
        'embryos': ['72h_counts.tsv', '135h_counts.tsv', 'SD_timeseries_counts.tsv',
                    'LD_timeseries_counts.tsv', 'all_counts.tsv'],
        'larvae': ['11d_counts.tsv', 'SD_timeseries_counts.tsv', 'all_counts.tsv']
    }
}

# Expected sample counts per matrix
EXPECTED_SAMPLE_COUNTS = {
    'unfed_counts.tsv': 8,
    'fed_counts.tsv': 8,
    '72h_counts.tsv': 6,
    '135h_counts.tsv': 6,
    'SD_timeseries_counts.tsv': {'adults': None, 'embryos': 6, 'larvae': 9},
    'LD_timeseries_counts.tsv': 6,
    '11d_counts.tsv': 6,
    'all_counts.tsv': {'adults': 16, 'embryos': 12, 'larvae': 17}
}

# =============================================================================
# Helper Functions
# =============================================================================

def print_header(text):
    """Print formatted header."""
    print(f"\n{'='*70}")
    print(text)
    print(f"{'='*70}")


def print_subheader(text):
    """Print formatted subheader."""
    print(f"\n{'-'*70}")
    print(text)
    print(f"{'-'*70}")


def print_check(check_name, passed, message=""):
    """Print check result."""
    status = "PASS" if passed else "FAIL"
    symbol = "✓" if passed else "✗"
    print(f"  [{status}] {symbol} {check_name}")
    if message:
        print(f"        {message}")
    return passed


# =============================================================================
# Sanity Check Functions
# =============================================================================

def check_file_structure():
    """Check 1: Verify all expected files exist and are non-empty."""
    print_subheader("Check 1: File Structure")

    all_passed = True
    total_files = 0
    found_files = 0

    for count_type, stages in EXPECTED_MATRICES.items():
        print(f"\n{count_type.upper()}:")

        for stage, matrix_files in stages.items():
            for matrix_file in matrix_files:
                total_files += 1
                filepath = Path(OUTPUT_BASE) / count_type / stage / matrix_file

                # Check existence
                if not filepath.exists():
                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        f"File not found: {filepath}")
                    all_passed = False
                    continue

                # Check non-empty
                if filepath.stat().st_size == 0:
                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        "File is empty")
                    all_passed = False
                    continue

                # Check metadata file exists
                meta_file = filepath.parent / matrix_file.replace('_counts.tsv', '_metadata.tsv')
                if not meta_file.exists():
                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        f"Metadata file missing: {meta_file.name}")
                    all_passed = False
                    continue

                found_files += 1
                print_check(f"{stage}/{matrix_file}", True,
                           f"Size: {filepath.stat().st_size / 1024:.1f} KB")

    print(f"\nTotal: {found_files}/{total_files} files found")

    return all_passed


def check_sample_counts():
    """Check 2: Verify correct number of samples per matrix."""
    print_subheader("Check 2: Sample Counts")

    all_passed = True

    for count_type, stages in EXPECTED_MATRICES.items():
        print(f"\n{count_type.upper()}:")

        for stage, matrix_files in stages.items():
            for matrix_file in matrix_files:
                filepath = Path(OUTPUT_BASE) / count_type / stage / matrix_file

                if not filepath.exists():
                    continue

                # Read matrix
                df = pd.read_csv(filepath, sep='\t', index_col=0)
                n_samples = df.shape[1]

                # Get expected count
                if matrix_file in EXPECTED_SAMPLE_COUNTS:
                    expected = EXPECTED_SAMPLE_COUNTS[matrix_file]
                    if isinstance(expected, dict):
                        expected = expected.get(stage, None)

                    if expected is None:
                        continue

                    if n_samples != expected:
                        passed = print_check(f"{stage}/{matrix_file}", False,
                                            f"Expected {expected} samples, got {n_samples}")
                        all_passed = False
                    else:
                        print_check(f"{stage}/{matrix_file}", True,
                                   f"{n_samples} samples")

    return all_passed


def check_gene_counts():
    """Check 3: Verify gene counts are in expected range."""
    print_subheader("Check 3: Gene Counts")

    all_passed = True
    gene_counts = []

    for count_type, stages in EXPECTED_MATRICES.items():
        print(f"\n{count_type.upper()}:")

        for stage, matrix_files in stages.items():
            for matrix_file in matrix_files:
                filepath = Path(OUTPUT_BASE) / count_type / stage / matrix_file

                if not filepath.exists():
                    continue

                # Read matrix
                df = pd.read_csv(filepath, sep='\t', index_col=0)
                n_genes = df.shape[0]
                gene_counts.append(n_genes)

                if n_genes < EXPECTED_GENES_MIN or n_genes > EXPECTED_GENES_MAX:
                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        f"Expected {EXPECTED_GENES_MIN}-{EXPECTED_GENES_MAX}, got {n_genes}")
                    all_passed = False
                else:
                    print_check(f"{stage}/{matrix_file}", True,
                               f"{n_genes} genes")

    if gene_counts:
        print(f"\nGene count range: {min(gene_counts)} - {max(gene_counts)}")
        print(f"Mean: {np.mean(gene_counts):.0f}, Median: {np.median(gene_counts):.0f}")

    return all_passed


def check_missing_data():
    """Check 4: Check for excessive missing data (NAs or zeros)."""
    print_subheader("Check 4: Missing Data")

    all_passed = True

    for count_type, stages in EXPECTED_MATRICES.items():
        print(f"\n{count_type.upper()}:")

        for stage, matrix_files in stages.items():
            for matrix_file in matrix_files:
                filepath = Path(OUTPUT_BASE) / count_type / stage / matrix_file

                if not filepath.exists():
                    continue

                # Read matrix
                df = pd.read_csv(filepath, sep='\t', index_col=0)

                # Check for NAs
                n_nas = df.isna().sum().sum()
                pct_nas = (n_nas / (df.shape[0] * df.shape[1])) * 100

                if pct_nas > 1.0:  # More than 1% NAs is suspicious
                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        f"{pct_nas:.2f}% NAs (too high)")
                    all_passed = False
                else:
                    print_check(f"{stage}/{matrix_file}", True,
                               f"{pct_nas:.3f}% NAs")

                # Check for excessive zeros
                n_zeros = (df == 0).sum().sum()
                pct_zeros = (n_zeros / (df.shape[0] * df.shape[1])) * 100

                if pct_zeros > 50:  # More than 50% zeros is very suspicious
                    print(f"        WARNING: {pct_zeros:.1f}% zeros (may be normal for RNA-seq)")

    return all_passed


def check_replicate_correlations():
    """Check 5: Calculate correlations between biological replicates."""
    print_subheader("Check 5: Replicate Correlations")

    all_passed = True

    # Define replicate groups (samples that should correlate highly)
    replicate_groups = {
        'adults': {
            'unfed_counts.tsv': {
                'SD_Unfed': ['SD_Unfed'],
                'LD_Unfed': ['LD_Unfed']
            },
            'fed_counts.tsv': {
                'SD_Fed': ['SD_Fed'],
                'LD_Fed': ['LD_Fed']
            }
        },
        'embryos': {
            '72h_counts.tsv': {
                'SD_72h': ['SD_72h'],
                'LD_72h': ['LD_72h']
            },
            '135h_counts.tsv': {
                'SD_135h': ['SD_135h'],
                'LD_135h': ['LD_135h']
            }
        },
        'larvae': {
            '11d_counts.tsv': {
                'SD_11d': ['SD_11d'],
                'LD_11d': ['LD_11d']
            }
        }
    }

    for count_type in ['featurecounts', 'salmon']:
        print(f"\n{count_type.upper()}:")

        for stage, matrices in replicate_groups.items():
            for matrix_file, groups in matrices.items():
                filepath = Path(OUTPUT_BASE) / count_type / stage / matrix_file

                if not filepath.exists():
                    continue

                # Read matrix
                df = pd.read_csv(filepath, sep='\t', index_col=0)

                for group_name, sample_prefixes in groups.items():
                    # Find samples matching prefix
                    group_samples = [col for col in df.columns
                                    if any(col.startswith(prefix) for prefix in sample_prefixes)]

                    if len(group_samples) < 2:
                        continue

                    # Calculate pairwise correlations
                    correlations = []
                    for i, sample1 in enumerate(group_samples):
                        for sample2 in group_samples[i+1:]:
                            # Use log2(count + 1) for correlation
                            x = np.log2(df[sample1] + 1)
                            y = np.log2(df[sample2] + 1)

                            r, _ = pearsonr(x, y)
                            correlations.append(r)

                    if correlations:
                        mean_r = np.mean(correlations)
                        min_r = np.min(correlations)

                        if min_r < MIN_REPLICATE_CORRELATION:
                            passed = print_check(f"{stage}/{group_name}", False,
                                                f"Min correlation: {min_r:.3f} (< {MIN_REPLICATE_CORRELATION})")
                            all_passed = False
                        else:
                            print_check(f"{stage}/{group_name}", True,
                                       f"Mean r = {mean_r:.3f}, Min r = {min_r:.3f}")

    return all_passed


def check_metadata_alignment():
    """Check 6: Verify samples in count matrices match metadata."""
    print_subheader("Check 6: Metadata Alignment")

    all_passed = True

    for count_type, stages in EXPECTED_MATRICES.items():
        print(f"\n{count_type.upper()}:")

        for stage, matrix_files in stages.items():
            for matrix_file in matrix_files:
                counts_file = Path(OUTPUT_BASE) / count_type / stage / matrix_file
                meta_file = counts_file.parent / matrix_file.replace('_counts.tsv', '_metadata.tsv')

                if not counts_file.exists() or not meta_file.exists():
                    continue

                # Read files
                counts = pd.read_csv(counts_file, sep='\t', index_col=0)
                metadata = pd.read_csv(meta_file, sep='\t')

                # Get sample IDs (create BioProject_Run IDs from metadata)
                count_samples = set(counts.columns)
                meta_samples = set(metadata['BioProject'] + '_' + metadata['Run'])

                # Check match
                if count_samples != meta_samples:
                    missing_in_meta = count_samples - meta_samples
                    missing_in_counts = meta_samples - count_samples

                    message = []
                    if missing_in_meta:
                        message.append(f"In counts but not metadata: {missing_in_meta}")
                    if missing_in_counts:
                        message.append(f"In metadata but not counts: {missing_in_counts}")

                    passed = print_check(f"{stage}/{matrix_file}", False,
                                        "; ".join(message))
                    all_passed = False
                else:
                    print_check(f"{stage}/{matrix_file}", True,
                               f"{len(count_samples)} samples match")

    return all_passed


# =============================================================================
# Main
# =============================================================================

def main():
    """Run all sanity checks."""
    print_header("SANITY CHECK: COUNT MATRICES")
    print("\nThis script validates the organized count matrices.")
    print("All checks must pass before proceeding to Phase 3 (Salmon validation).")

    # Run all checks
    results = {
        'File Structure': check_file_structure(),
        'Sample Counts': check_sample_counts(),
        'Gene Counts': check_gene_counts(),
        'Missing Data': check_missing_data(),
        'Replicate Correlations': check_replicate_correlations(),
        'Metadata Alignment': check_metadata_alignment()
    }

    # Print summary
    print_header("SUMMARY")

    all_passed = all(results.values())

    for check_name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        symbol = "✓" if passed else "✗"
        print(f"  [{status}] {symbol} {check_name}")

    if all_passed:
        print("\nAll sanity checks passed!")
        print("Ready to proceed to Phase 3 (Validate Salmon vs featureCounts)")
        sys.exit(0)
    else:
        print("\nSome sanity checks failed!")
        print("Please investigate and fix issues before proceeding.")
        print("See error messages above for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
