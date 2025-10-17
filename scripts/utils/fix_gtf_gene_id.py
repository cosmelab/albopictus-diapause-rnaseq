#!/usr/bin/env python3
import gzip
import re
import sys

def fix_gtf_gene_id(input_file, output_file):
    """Add missing gene_id attributes to GTF file."""

    lines_processed = 0
    lines_fixed = 0

    with open(input_file, 'r') as infile, \
         gzip.open(output_file, 'wt') as outfile:

        for line in infile:
            lines_processed += 1

            if line.startswith('#'):
                outfile.write(line)
                continue

            # Check if gene_id is missing
            if 'gene_id' not in line:
                # Extract gene name from gene_name or gene attribute
                gene_match = re.search(r'gene[_ ]"([^"]+)"', line)
                if gene_match:
                    gene_id = gene_match.group(1)
                    # Remove "gene-" prefix if present
                    gene_id = gene_id.replace('gene-', '')
                    # Strip trailing whitespace and semicolons
                    line = line.rstrip().rstrip(';')
                    # Add gene_id
                    line = line + f'; gene_id "{gene_id}";\n'
                    lines_fixed += 1

            outfile.write(line)

    print(f"Processed {lines_processed} lines")
    print(f"Fixed {lines_fixed} lines by adding gene_id")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_gtf_gene_id.py input.gtf output.gtf.gz")
        sys.exit(1)

    fix_gtf_gene_id(sys.argv[1], sys.argv[2])
