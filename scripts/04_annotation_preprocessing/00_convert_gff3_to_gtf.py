#!/usr/bin/env python3
"""
Convert GFF3 to GTF format for RNA-seq analysis
Author: L. Cosme
Date: October 2025

This script converts the Aedes albopictus GFF3 annotation to GTF format.
"""

import sys
import os
import re
from collections import defaultdict

def parse_attributes(attr_string, gff_format=True):
    """Parse GFF3 or GTF attributes"""
    attrs = {}
    if gff_format:
        # GFF3 format: key=value pairs separated by semicolons
        for item in attr_string.strip().split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                attrs[key] = value
    else:
        # GTF format: key "value" pairs
        pattern = r'(\w+)\s+"([^"]+)"'
        for match in re.finditer(pattern, attr_string):
            attrs[match.group(1)] = match.group(2)
    return attrs

def format_gtf_attributes(gene_id, transcript_id=None, **extra_attrs):
    """Format attributes in GTF style"""
    parts = []

    # Required attributes first
    if gene_id:
        parts.append(f'gene_id "{gene_id}"')
    if transcript_id:
        parts.append(f'transcript_id "{transcript_id}"')

    # Add any extra attributes
    for key, value in extra_attrs.items():
        if value and key not in ['gene_id', 'transcript_id']:
            parts.append(f'{key} "{value}"')

    return '; '.join(parts) + ';'

def gff3_to_gtf(gff3_file, gtf_file):
    """Convert GFF3 to GTF format"""

    print(f"Converting {gff3_file} to {gtf_file}")

    # Track gene and transcript relationships
    gene_to_transcripts = defaultdict(set)
    transcript_to_gene = {}
    feature_counts = defaultdict(int)

    with open(gff3_file, 'r') as infile, open(gtf_file, 'w') as outfile:
        # Write header
        outfile.write('##gtf-version 2.2\n')

        for line_num, line in enumerate(infile, 1):
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = parts
            attrs = parse_attributes(attributes, gff_format=True)

            feature_counts[feature] += 1

            # Skip features we don't need in GTF
            if feature in ['region', 'sequence_feature', 'match', 'cDNA_match']:
                continue

            # Handle different feature types
            if feature == 'gene':
                gene_id = attrs.get('ID', attrs.get('Name', f'gene_{line_num}'))
                gene_name = attrs.get('Name', gene_id)

                # Write gene feature
                gtf_attrs = format_gtf_attributes(
                    gene_id=gene_id,
                    gene_name=gene_name,
                    gene_biotype=attrs.get('biotype', 'protein_coding')
                )
                outfile.write('\t'.join([chrom, source, 'gene', start, end, score, strand, frame, gtf_attrs]) + '\n')

            elif feature in ['mRNA', 'transcript', 'rRNA', 'tRNA', 'ncRNA', 'snoRNA', 'snRNA']:
                transcript_id = attrs.get('ID', f'transcript_{line_num}')
                parent = attrs.get('Parent', '').replace('gene:', '')

                if parent:
                    gene_id = parent
                    transcript_to_gene[transcript_id] = gene_id
                    gene_to_transcripts[gene_id].add(transcript_id)
                else:
                    # If no parent, use transcript ID as gene ID
                    gene_id = transcript_id

                # Write transcript feature
                gtf_attrs = format_gtf_attributes(
                    gene_id=gene_id,
                    transcript_id=transcript_id,
                    transcript_biotype=attrs.get('biotype', 'protein_coding')
                )
                outfile.write('\t'.join([chrom, source, 'transcript', start, end, score, strand, frame, gtf_attrs]) + '\n')

            elif feature == 'exon':
                parent = attrs.get('Parent', '')

                # Handle multiple parents (exons can belong to multiple transcripts)
                parents = parent.split(',') if parent else []

                for parent_id in parents:
                    parent_id = parent_id.replace('transcript:', '').replace('mRNA:', '')

                    if parent_id in transcript_to_gene:
                        gene_id = transcript_to_gene[parent_id]
                        transcript_id = parent_id
                    else:
                        # If parent not found, try to extract gene from parent ID
                        gene_id = parent_id.split('.')[0] if '.' in parent_id else parent_id
                        transcript_id = parent_id

                    # Calculate exon number (simplified - would need more logic for accuracy)
                    exon_number = attrs.get('exon_number', attrs.get('rank', '1'))

                    # Write exon feature
                    gtf_attrs = format_gtf_attributes(
                        gene_id=gene_id,
                        transcript_id=transcript_id,
                        exon_number=exon_number,
                        exon_id=attrs.get('ID', f'{transcript_id}_exon_{exon_number}')
                    )
                    outfile.write('\t'.join([chrom, source, 'exon', start, end, score, strand, frame, gtf_attrs]) + '\n')

            elif feature == 'CDS':
                parent = attrs.get('Parent', '')
                parents = parent.split(',') if parent else []

                for parent_id in parents:
                    parent_id = parent_id.replace('transcript:', '').replace('mRNA:', '')

                    if parent_id in transcript_to_gene:
                        gene_id = transcript_to_gene[parent_id]
                        transcript_id = parent_id
                    else:
                        gene_id = parent_id.split('.')[0] if '.' in parent_id else parent_id
                        transcript_id = parent_id

                    # Write CDS feature
                    gtf_attrs = format_gtf_attributes(
                        gene_id=gene_id,
                        transcript_id=transcript_id,
                        protein_id=attrs.get('protein_id', transcript_id)
                    )
                    outfile.write('\t'.join([chrom, source, 'CDS', start, end, score, strand, frame, gtf_attrs]) + '\n')

            elif feature in ['five_prime_UTR', 'three_prime_UTR']:
                parent = attrs.get('Parent', '')
                parents = parent.split(',') if parent else []

                for parent_id in parents:
                    parent_id = parent_id.replace('transcript:', '').replace('mRNA:', '')

                    if parent_id in transcript_to_gene:
                        gene_id = transcript_to_gene[parent_id]
                        transcript_id = parent_id
                    else:
                        gene_id = parent_id.split('.')[0] if '.' in parent_id else parent_id
                        transcript_id = parent_id

                    # Convert feature name for GTF
                    gtf_feature = feature.replace('_prime_', '_prime_')

                    # Write UTR feature
                    gtf_attrs = format_gtf_attributes(
                        gene_id=gene_id,
                        transcript_id=transcript_id
                    )
                    outfile.write('\t'.join([chrom, source, gtf_feature, start, end, score, strand, frame, gtf_attrs]) + '\n')

    # Print summary statistics
    print("\nConversion complete!")
    print("Feature counts in GFF3:")
    for feature, count in sorted(feature_counts.items()):
        print(f"  {feature}: {count}")
    print(f"\nTotal genes with transcripts: {len(gene_to_transcripts)}")
    print(f"Total transcripts: {len(transcript_to_gene)}")

    return gtf_file

def main():
    """Main function"""

    # Set paths
    base_dir = "/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
    gff3_file = os.path.join(base_dir, "data/references/AalbF3_processed/AalbF3.gff3")
    gtf_file = os.path.join(base_dir, "data/references/AalbF3_processed/AalbF3_raw.gtf")

    if not os.path.exists(gff3_file):
        print(f"Error: GFF3 file not found: {gff3_file}")
        sys.exit(1)

    # Convert GFF3 to GTF
    gff3_to_gtf(gff3_file, gtf_file)

    # Check the output
    if os.path.exists(gtf_file):
        file_size = os.path.getsize(gtf_file) / (1024 * 1024)  # Size in MB
        with open(gtf_file, 'r') as f:
            line_count = sum(1 for _ in f)
        print(f"\nOutput GTF file:")
        print(f"  Path: {gtf_file}")
        print(f"  Size: {file_size:.1f} MB")
        print(f"  Lines: {line_count:,}")

        # Check for gene_id presence
        with open(gtf_file, 'r') as f:
            sample_lines = [next(f) for _ in range(min(100, line_count))]

        lines_with_gene_id = sum(1 for line in sample_lines if 'gene_id' in line and not line.startswith('#'))
        print(f"  Sample check: {lines_with_gene_id}/100 lines have gene_id")

if __name__ == "__main__":
    main()