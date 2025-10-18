#!/bin/bash
#
# Clean rRNA FASTA headers for nf-core compatibility
#
# Problem: Original GenBank headers contain spaces and descriptions:
#   >OM350214.1 Aedes albopictus isolate KH_S1 small subunit ribosomal RNA, partial sequence
#
# Solution: Simplify to just accession + rRNA type:
#   >OM350214.1_18S
#
# Why: nf-core/rnaseq uses --additional_fasta to auto-generate GTF entries
#      Simple headers avoid parsing issues with Salmon quantification
#
# Run from project root:
# bash scripts/01_download_data/clean_rrna_headers.sh

set -euo pipefail

echo "========================================"
echo "Cleaning rRNA FASTA headers"
echo "========================================"
echo ""

# Define paths
INPUT_DIR="data/references/AalbF3_genome"

# Check input files exist in raw_downloads
if [ ! -f "${INPUT_DIR}/raw_downloads/18S_rRNA_sequences.fa" ]; then
    echo "ERROR: ${INPUT_DIR}/raw_downloads/18S_rRNA_sequences.fa not found"
    echo "Run download_rrna_sequences.sh first"
    exit 1
fi

if [ ! -f "${INPUT_DIR}/raw_downloads/28S_rRNA_sequences.fa" ]; then
    echo "ERROR: ${INPUT_DIR}/raw_downloads/28S_rRNA_sequences.fa not found"
    echo "Run download_rrna_sequences.sh first"
    exit 1
fi

# Clean 18S headers: keep only accession, add _18S suffix
echo "Cleaning 18S rRNA headers..."
sed 's/^\(>OM[0-9]*\.[0-9]*\) .*/\1_18S/' \
  "${INPUT_DIR}/raw_downloads/18S_rRNA_sequences.fa" > \
  "${INPUT_DIR}/18S_rRNA_sequences.fa"

echo "  Created: ${INPUT_DIR}/18S_rRNA_sequences.fa"
echo "  Sequences: $(grep -c "^>" ${INPUT_DIR}/18S_rRNA_sequences.fa)"
echo ""

# Clean 28S headers: keep only accession, add _28S suffix
echo "Cleaning 28S rRNA headers..."
sed 's/^\(>OM[0-9]*\.[0-9]*\) .*/\1_28S/' \
  "${INPUT_DIR}/raw_downloads/28S_rRNA_sequences.fa" > \
  "${INPUT_DIR}/28S_rRNA_sequences.fa"

echo "  Created: ${INPUT_DIR}/28S_rRNA_sequences.fa"
echo "  Sequences: $(grep -c "^>" ${INPUT_DIR}/28S_rRNA_sequences.fa)"
echo ""

# Combine cleaned files
echo "Creating combined rRNA file..."
cat "${INPUT_DIR}/18S_rRNA_sequences.fa" \
    "${INPUT_DIR}/28S_rRNA_sequences.fa" > \
    "${INPUT_DIR}/combined_rRNA_sequences.fa"

echo "  Created: ${INPUT_DIR}/combined_rRNA_sequences.fa"
echo "  Total sequences: $(grep -c "^>" ${INPUT_DIR}/combined_rRNA_sequences.fa)"
echo ""

# Show example of cleaned headers
echo "========================================"
echo "Example cleaned headers:"
echo "========================================"
grep "^>" "${INPUT_DIR}/combined_rRNA_sequences.fa" | head -3
echo "..."
grep "^>" "${INPUT_DIR}/combined_rRNA_sequences.fa" | tail -3
echo ""

echo "========================================"
echo "Header cleaning complete!"
echo "========================================"
echo ""
echo "Raw downloads preserved in: ${INPUT_DIR}/raw_downloads/"
echo ""
echo "For nf-core pipeline, use:"
echo "  --additional_fasta ${INPUT_DIR}/combined_rRNA_sequences.fa"
