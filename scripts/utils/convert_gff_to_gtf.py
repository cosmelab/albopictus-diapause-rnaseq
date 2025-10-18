#!/usr/bin/env python3
"""
Convert GFF3 to GTF format.

This script converts NCBI GFF3 format to GTF format, properly handling:
- Attribute format conversion (key=value -> key "value")
- Adding transcript_id and gene_id to all features
- Proper GTF formatting with semicolon-space separation
"""

import sys
import gzip
import re

def parse_gff_attributes(attr_string):
    """Parse GFF3 attributes (key=value;key=value)."""
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs

def format_gtf_attributes(attrs_dict, feature_type):
    """Format attributes for GTF (key "value"; key "value")."""
    # Always include transcript_id and gene_id first if available
    gtf_attrs = []

    # Extract transcript_id (from Parent or transcript_id)
    transcript_id = None
    if 'Parent' in attrs_dict:
        transcript_id = attrs_dict['Parent']
    elif 'transcript_id' in attrs_dict:
        transcript_id = attrs_dict['transcript_id']
    elif 'ID' in attrs_dict and feature_type in ['mRNA', 'tRNA', 'rRNA', 'ncRNA', 'transcript']:
        transcript_id = attrs_dict['ID']

    # Extract gene_id (from gene or ID for gene features)
    gene_id = None
    if 'gene' in attrs_dict:
        gene_id = attrs_dict['gene']
    elif feature_type == 'gene' and 'ID' in attrs_dict:
        gene_id = attrs_dict['ID']
    elif 'Parent' in attrs_dict:
        # For exons/CDS, extract gene from parent
        parent = attrs_dict['Parent']
        if 'gene' in attrs_dict:
            gene_id = attrs_dict['gene']

    # Add transcript_id and gene_id
    if transcript_id:
        gtf_attrs.append(f'transcript_id "{transcript_id}"')
    if gene_id:
        gtf_attrs.append(f'gene_id "{gene_id}"')

    # Add other relevant attributes
    for key in ['gene_name', 'gene', 'Dbxref', 'Name', 'gbkey', 'product', 'Note',
                'partial', 'start_range', 'end_range', 'protein_id', 'inference',
                'anticodon', 'model_evidence', 'gene_biotype']:
        if key in attrs_dict:
            value = attrs_dict[key]
            # Don't duplicate transcript_id or gene in attributes
            if key == 'gene' and gene_id:
                continue
            gtf_attrs.append(f'{key} "{value}"')

    return '; '.join(gtf_attrs) + ';'

def convert_gff_to_gtf(input_file, output_file):
    """Convert GFF3 to GTF format."""
    print(f"Reading GFF3 file: {input_file}")

    lines_processed = 0
    lines_written = 0

    # Open files
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

            # Skip comment lines
            if line.startswith('#'):
                continue

            # Parse GFF3 line
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            chrom, source, feature_type, start, end, score, strand, frame, attributes = fields

            # Skip gene features (we only need transcript-level and below)
            if feature_type == 'gene':
                continue

            # Convert feature types
            # mRNA/tRNA/rRNA/ncRNA -> transcript in GTF (but keep original for some tools)
            if feature_type in ['mRNA', 'tRNA', 'rRNA', 'ncRNA']:
                gtf_feature = 'transcript'
            else:
                gtf_feature = feature_type

            # Parse and convert attributes
            gff_attrs = parse_gff_attributes(attributes)
            gtf_attrs = format_gtf_attributes(gff_attrs, feature_type)

            # Write GTF line
            gtf_line = f"{chrom}\t{source}\t{gtf_feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attrs}\n"
            outfile.write(gtf_line)
            lines_written += 1

    print(f"\nConversion Summary:")
    print(f"Lines processed: {lines_processed}")
    print(f"Features written: {lines_written}")
    print(f"\nGTF file created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_gff_to_gtf.py input.gff.gz output.gtf.gz")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    convert_gff_to_gtf(input_file, output_file)
