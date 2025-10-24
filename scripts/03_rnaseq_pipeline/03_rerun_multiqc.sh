#!/bin/bash
#SBATCH --job-name=multiqc_rerun
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=epyc
#SBATCH -o /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/multiqc_rerun_%j.o.txt
#SBATCH -e /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/multiqc_rerun_%j.e.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lcosme@ucr.edu

#############################################################################
# MultiQC Re-run Script
#
# Purpose: Re-generate MultiQC summary report after pipeline completion
#          Original run failed due to /tmp space exhaustion
#
# Solution: Use project directory for temp files instead of /tmp
#
# Expected runtime: 10-30 minutes
#
# Date: October 19, 2025
#############################################################################

set -euo pipefail

#############################################################################
# Configuration
#############################################################################

PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
OUTPUT_DIR="${PROJECT_BASE}/output"
MULTIQC_DIR="${OUTPUT_DIR}/multiqc"
TEMP_DIR="${OUTPUT_DIR}/multiqc_temp"

# Create temp directory with plenty of space
mkdir -p "${TEMP_DIR}"
mkdir -p "${MULTIQC_DIR}"

# Load required modules FIRST (module purge clears env vars)
module purge
module load singularity/3.9.3

# THEN set environment variables to use our temp directory
export TMPDIR="${TEMP_DIR}"
export TEMP="${TEMP_DIR}"
export TMP="${TEMP_DIR}"

echo "============================================="
echo "MultiQC Re-run"
echo "============================================="
echo "Start time: $(date)"
echo "Temp directory: ${TMPDIR}"
echo "Output directory: ${MULTIQC_DIR}"
echo "---------------------------------------------"

#############################################################################
# Run MultiQC
#############################################################################

# Use the same singularity container that nf-core used
MULTIQC_CONTAINER="${OUTPUT_DIR}/singularity/nf-core-multiqc-1.28.0-4e3d224d7f16bb26.img"

if [[ ! -f "${MULTIQC_CONTAINER}" ]]; then
    echo "ERROR: MultiQC container not found at ${MULTIQC_CONTAINER}"
    echo "Looking for alternative containers..."
    MULTIQC_CONTAINER=$(find ${OUTPUT_DIR}/singularity -name "*multiqc*.img" | head -1)
    if [[ -z "${MULTIQC_CONTAINER}" ]]; then
        echo "ERROR: No MultiQC container found. Cannot proceed."
        exit 1
    fi
    echo "Using: ${MULTIQC_CONTAINER}"
fi

# Change to output directory (MultiQC will scan subdirectories)
cd "${OUTPUT_DIR}"

echo ""
echo "Running MultiQC aggregation..."
echo "This will scan all QC outputs and create summary report"
echo ""

# Run MultiQC with our fixed temp directory
# Note: Let MultiQC auto-detect modules instead of specifying them
#       This avoids version compatibility issues
singularity exec \
    --cleanenv \
    --bind "${PROJECT_BASE}:${PROJECT_BASE}" \
    --bind "${TMPDIR}:${TMPDIR}" \
    "${MULTIQC_CONTAINER}" \
    multiqc \
        --force \
        --title "Albopictus Diapause RNAseq - Complete Pipeline 3.19.0" \
        --filename "${MULTIQC_DIR}/multiqc_report.html" \
        --outdir "${MULTIQC_DIR}" \
        --dirs \
        --dirs-depth 2 \
        .

MULTIQC_EXIT=$?

echo ""
echo "============================================="

if [[ ${MULTIQC_EXIT} -eq 0 ]]; then
    echo "✅ MultiQC completed successfully!"
    echo ""
    echo "Report location: ${MULTIQC_DIR}/multiqc_report.html"
    echo ""
    echo "To view the report:"
    echo "  scp lcosme@cluster.hpcc.ucr.edu:${MULTIQC_DIR}/multiqc_report.html ."
    echo "  open multiqc_report.html"
    echo ""

    # List generated files
    echo "Generated files:"
    ls -lh "${MULTIQC_DIR}/"

    # Clean up temp directory
    echo ""
    echo "Cleaning up temp directory..."
    rm -rf "${TEMP_DIR}"
    echo "Temp directory removed: ${TEMP_DIR}"

else
    echo "❌ MultiQC failed with exit code ${MULTIQC_EXIT}"
    echo ""
    echo "Temp files preserved for debugging: ${TEMP_DIR}"
    echo "Check logs for details"
fi

echo "============================================="
echo "Finished at: $(date)"
echo "============================================="

exit ${MULTIQC_EXIT}
