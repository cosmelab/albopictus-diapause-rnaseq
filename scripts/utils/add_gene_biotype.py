#!/usr/bin/env python3
"""
Add gene_biotype attribute to all GTF lines based on gbkey.

This fixes the issue where featureCounts expects gene_biotype but the GTF
only has it on some lines. We map gbkey values to standard gene_biotype terms.
"""

import sys
import gzip
import re

def map_gbkey_to_biotype(gbkey):
    """Map NCBI gbkey to standard gene_biotype."""
    mapping = {
        'mRNA': 'protein_coding',
        'CDS': 'protein_coding',
        'ncRNA': 'ncRNA',
        'misc_RNA': 'misc_RNA',
        'tRNA': 'tRNA',
        'rRNA': 'rRNA',
        'exon': 'protein_coding',  # exons are typically from protein-coding genes
        'Gene': 'protein_coding',  # generic gene features
    }
    return mapping.get(gbkey, 'unknown')

def add_gene_biotype(input_file, output_file):
    """Add gene_biotype attribute to all GTF lines."""

    print(f"Reading input file: {input_file}")
    lines_processed = 0
    lines_with_biotype = 0
    lines_added_biotype = 0

    # Open input and output files (handle gzipped)
    if input_file.endswith('.gz'):
        infile = gzip.open(input_file, 'rt')
    else:
        infile = open(input_file, 'r')

    if output_file.endswith('.gz'):
        outfile = gzip.open(output_file, 'wt')
    else:
        outfile = open(output_file, 'w')

    with infile, outfile:
        for line in infile:
            lines_processed += 1

            # Pass through comment lines unchanged
            if line.startswith('#'):
                outfile.write(line)
                continue

            # Check if line already has gene_biotype
            if 'gene_biotype' in line:
                lines_with_biotype += 1
                outfile.write(line)
                continue

            # Extract gbkey value
            gbkey_match = re.search(r'gbkey "([^"]+)"', line)
            if not gbkey_match:
                # No gbkey, write line unchanged
                outfile.write(line)
                continue

            gbkey = gbkey_match.group(1)
            biotype = map_gbkey_to_biotype(gbkey)

            # Add gene_biotype before the trailing newline
            # Strip trailing whitespace/newline, remove trailing semicolon if present
            line = line.rstrip().rstrip(';')
            line = line + f'; gene_biotype "{biotype}";\n'

            outfile.write(line)
            lines_added_biotype += 1

    print(f"\nProcessing Summary:")
    print(f"Total lines processed: {lines_processed}")
    print(f"Lines already with gene_biotype: {lines_with_biotype}")
    print(f"Lines with gene_biotype added: {lines_added_biotype}")
    print(f"\nOutput written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_gene_biotype.py input.gtf.gz output.gtf.gz")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    add_gene_biotype(input_file, output_file)
