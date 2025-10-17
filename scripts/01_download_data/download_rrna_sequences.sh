#!/bin/bash
#
# Download rRNA sequences for Aedes albopictus ONLY
#
# Source: Koh et al. 2023, eLife 82762 (Appendix 1)
# Data: https://elifesciences.org/articles/82762
#
# Aedes albopictus accessions (12 sequences total):
#   18S rRNA: OM350214, OM350220, OM350316, OM350217, OM350218, OM350219 (6 sequences)
#   28S rRNA: OM542460, OM542373, OM542374, OM542342, OM542343, OM542344 (6 sequences)
#
# Run from project root:
# bash scripts/01_download_data/download_rrna_sequences.sh

set -euo pipefail

echo "========================================"
echo "Downloading Aedes albopictus rRNA sequences"
echo "========================================"
echo ""

# Define output directory
OUTPUT_DIR="data/references/AalbF3_genome/raw_downloads"
mkdir -p "${OUTPUT_DIR}"

# NCBI E-utilities base URL
EUTILS_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Download 18S rRNA sequences (6 Ae. albopictus sequences)
echo "Downloading 18S rRNA sequences (6 Ae. albopictus)..."
curl -s "${EUTILS_BASE}/efetch.fcgi?db=nucleotide&id=OM350214,OM350220,OM350316,OM350217,OM350218,OM350219&rettype=fasta&retmode=text" > "${OUTPUT_DIR}/18S_rRNA_sequences.fa"
echo "  Saved: ${OUTPUT_DIR}/18S_rRNA_sequences.fa"
echo "  Sequences: $(grep -c "^>" ${OUTPUT_DIR}/18S_rRNA_sequences.fa)"
echo ""

# Download 28S rRNA sequences (6 Ae. albopictus sequences)
echo "Downloading 28S rRNA sequences (6 Ae. albopictus)..."
curl -s "${EUTILS_BASE}/efetch.fcgi?db=nucleotide&id=OM542460,OM542373,OM542374,OM542342,OM542343,OM542344&rettype=fasta&retmode=text" > "${OUTPUT_DIR}/28S_rRNA_sequences.fa"
echo "  Saved: ${OUTPUT_DIR}/28S_rRNA_sequences.fa"
echo "  Sequences: $(grep -c "^>" ${OUTPUT_DIR}/28S_rRNA_sequences.fa)"
echo ""

# Create combined rRNA file
echo "Creating combined rRNA file..."
cat "${OUTPUT_DIR}/18S_rRNA_sequences.fa" "${OUTPUT_DIR}/28S_rRNA_sequences.fa" > "${OUTPUT_DIR}/combined_rRNA_sequences.fa"
echo "  Saved: ${OUTPUT_DIR}/combined_rRNA_sequences.fa"
echo "  Total rRNA sequences: $(grep -c "^>" ${OUTPUT_DIR}/combined_rRNA_sequences.fa)"
echo ""

echo "========================================"
echo "Download complete!"
echo "========================================"
echo ""
echo "Files created in raw_downloads/ (have GenBank descriptions):"
echo "  - ${OUTPUT_DIR}/18S_rRNA_sequences.fa"
echo "  - ${OUTPUT_DIR}/28S_rRNA_sequences.fa"
echo "  - ${OUTPUT_DIR}/combined_rRNA_sequences.fa"
echo ""
echo "Next step: Clean headers for nf-core compatibility"
echo "  bash scripts/01_download_data/clean_rrna_headers.sh"
