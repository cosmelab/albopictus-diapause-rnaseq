#!/bin/bash
#
# Strip 'chr' prefix from Dryad reference files
#
# Problem: Dryad files use "chr" prefix (chr1.1, chr2.206)
#          SNP chip uses no prefix (1.1, 2.206)
# Solution: Strip "chr" prefix from FASTA headers and GFF3 chromosome column
#
# Run from project root:
# bash scripts/01_download_data/strip_chr_prefix.sh

set -euo pipefail

echo "========================================"
echo "Stripping chr prefix from reference files"
echo "========================================"
echo ""

# Define paths
DRYAD_DIR="data/references/dryad_download/doi_10_5061_dryad_mgqnk98z4__v20210127"
OUTPUT_DIR="data/references/AalbF3_genome"

# Check input files exist
if [ ! -f "${DRYAD_DIR}/AALBF3_assembly.fa.gz" ]; then
    echo "ERROR: ${DRYAD_DIR}/AALBF3_assembly.fa.gz not found"
    exit 1
fi

if [ ! -f "${DRYAD_DIR}/AALBF3_annotation.gff3.gz" ]; then
    echo "ERROR: ${DRYAD_DIR}/AALBF3_annotation.gff3.gz not found"
    exit 1
fi

# Create output directory if needed
mkdir -p "${OUTPUT_DIR}"

# Process genome FASTA
echo "Processing genome FASTA (this takes 2-3 minutes for 423 MB file)..."
zcat "${DRYAD_DIR}/AALBF3_assembly.fa.gz" | \
  sed 's/^>chr/>/' | \
  gzip > "${OUTPUT_DIR}/AalbF3_genome.fa.gz"

echo "  Created: ${OUTPUT_DIR}/AalbF3_genome.fa.gz"
ls -lh "${OUTPUT_DIR}/AalbF3_genome.fa.gz"
echo ""

# Process GFF3 annotation
echo "Processing GFF3 annotation (this takes ~30 seconds for 6.3 MB file)..."
zcat "${DRYAD_DIR}/AALBF3_annotation.gff3.gz" | \
  sed 's/^chr//' | \
  gzip > "${OUTPUT_DIR}/AalbF3_annotation.gff3.gz"

echo "  Created: ${OUTPUT_DIR}/AalbF3_annotation.gff3.gz"
ls -lh "${OUTPUT_DIR}/AalbF3_annotation.gff3.gz"
echo ""

# Validate results
echo "========================================"
echo "Validating chromosome names"
echo "========================================"
echo ""

echo "Genome FASTA headers (first 5):"
zcat "${OUTPUT_DIR}/AalbF3_genome.fa.gz" | grep "^>" | head -5
echo ""

echo "GFF3 chromosome names (unique, first 10):"
zcat "${OUTPUT_DIR}/AalbF3_annotation.gff3.gz" | grep -v "^#" | cut -f1 | sort -u | head -10
echo ""

echo "========================================"
echo "SUCCESS: chr prefix removed"
echo "========================================"
echo ""
echo "Output files:"
echo "  ${OUTPUT_DIR}/AalbF3_genome.fa.gz"
echo "  ${OUTPUT_DIR}/AalbF3_annotation.gff3.gz"
echo ""
echo "Next step: Convert GFF3 to GTF for nf-core pipeline"
