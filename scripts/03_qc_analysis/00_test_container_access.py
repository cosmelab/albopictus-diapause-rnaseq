#!/usr/bin/env python3
"""
Test script to verify container can access files correctly
"""
import os
from pathlib import Path

print("Testing file access in container...")
print(f"Current working directory: {os.getcwd()}")
print()

# Test access to output directory
output_dir = Path("output")
if output_dir.exists():
    print("✓ Found 'output' directory")
    projects = list(output_dir.glob("PRJNA*"))
    print(f"  Found {len(projects)} projects: {[p.name for p in projects]}")
    
    # Check first sample in first project
    if projects:
        samples = list(projects[0].glob("SRR*"))
        print(f"  Found {len(samples)} samples in {projects[0].name}")
        
        if samples:
            # Check for key files
            test_files = [
                "multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt",
                "star_salmon/salmon.merged.gene_counts.tsv",
                "star_salmon/qualimap/*/rnaseq_qc_results.txt"
            ]
            
            print(f"\n  Checking files in {samples[0].name}:")
            for test_file in test_files:
                if list(samples[0].glob(test_file)):
                    print(f"    ✓ Found: {test_file}")
                else:
                    print(f"    ✗ Missing: {test_file}")
else:
    print("✗ Cannot find 'output' directory")

print("\nFile system check complete.")