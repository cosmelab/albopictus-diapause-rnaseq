#!/usr/bin/env python3
"""
Filter problematic transcripts from GTF file
Author: L. Cosme
Date: October 2025

This script identifies and removes transcripts that span multiple chromosomes,
which can cause issues with RNA-seq quantification tools.
"""

import sys
import os
from collections import defaultdict

def check_and_filter_transcripts(input_gtf, output_gtf):
    """
    Check for transcripts spanning multiple chromosomes and filter them out
    """

    print(f"Checking GTF file: {input_gtf}")

    # First pass: identify problematic transcripts
    transcript_chromosomes = defaultdict(set)
    gene_chromosomes = defaultdict(set)

    with open(input_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature = parts[2]
            attributes = parts[8]

            # Extract transcript_id and gene_id
            transcript_id = None
            gene_id = None

            for attr in attributes.split(';'):
                attr = attr.strip()
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]
                if 'gene_id' in attr:
                    gene_id = attr.split('"')[1]

            if transcript_id:
                transcript_chromosomes[transcript_id].add(chrom)
            if gene_id and feature == 'gene':
                gene_chromosomes[gene_id].add(chrom)

    # Identify problematic transcripts
    problematic_transcripts = set()
    problematic_genes = set()

    for transcript_id, chroms in transcript_chromosomes.items():
        if len(chroms) > 1:
            problematic_transcripts.add(transcript_id)
            print(f"  Problematic transcript: {transcript_id} on chromosomes: {', '.join(sorted(chroms))}")

    for gene_id, chroms in gene_chromosomes.items():
        if len(chroms) > 1:
            problematic_genes.add(gene_id)
            print(f"  Problematic gene: {gene_id} on chromosomes: {', '.join(sorted(chroms))}")

    print(f"\nFound {len(problematic_transcripts)} problematic transcripts")
    print(f"Found {len(problematic_genes)} problematic genes")

    if len(problematic_transcripts) == 0 and len(problematic_genes) == 0:
        print("No problematic features found - GTF is clean!")
        # Just copy the file
        with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
            outfile.write(infile.read())
        return output_gtf

    # Second pass: filter out problematic transcripts
    print(f"\nFiltering problematic features to: {output_gtf}")

    lines_removed = 0
    lines_kept = 0

    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                outfile.write(line)
                continue

            attributes = parts[8]

            # Extract transcript_id and gene_id
            transcript_id = None
            gene_id = None

            for attr in attributes.split(';'):
                attr = attr.strip()
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]
                if 'gene_id' in attr:
                    gene_id = attr.split('"')[1]

            # Skip if transcript or gene is problematic
            if transcript_id and transcript_id in problematic_transcripts:
                lines_removed += 1
                continue
            if gene_id and gene_id in problematic_genes:
                lines_removed += 1
                continue

            outfile.write(line)
            lines_kept += 1

    print(f"Lines removed: {lines_removed}")
    print(f"Lines kept: {lines_kept}")

    return output_gtf

def main():
    """Main function"""

    # Set paths
    base_dir = "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
    input_gtf = os.path.join(base_dir, "data/references/AalbF3_processed/AalbF3_raw.gtf")
    output_gtf = os.path.join(base_dir, "data/references/AalbF3_processed/AalbF3_filtered.gtf")

    if not os.path.exists(input_gtf):
        print(f"Error: Input GTF file not found: {input_gtf}")
        sys.exit(1)

    # Check and filter transcripts
    check_and_filter_transcripts(input_gtf, output_gtf)

    # Check the output
    if os.path.exists(output_gtf):
        file_size = os.path.getsize(output_gtf) / (1024 * 1024)  # Size in MB
        with open(output_gtf, 'r') as f:
            line_count = sum(1 for _ in f)
        print(f"\nOutput GTF file:")
        print(f"  Path: {output_gtf}")
        print(f"  Size: {file_size:.1f} MB")
        print(f"  Lines: {line_count:,}")

if __name__ == "__main__":
    main()