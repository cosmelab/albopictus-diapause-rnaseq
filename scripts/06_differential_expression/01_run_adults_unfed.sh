#!/bin/bash
#SBATCH --job-name=adults_unfed_DE
#SBATCH --output=logs/adults_unfed_deseq2_%j.out
#SBATCH --error=logs/adults_unfed_deseq2_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc
#SBATCH --account=cosmelab

################################################################################
# Adult Female DESeq2 Analysis - Unfed Only
#
# Runs differential expression analysis on unfed adult females using DESeq2
# within Singularity container for complete reproducibility.
#
# Author: L. Cosme
# Date: October 22, 2025
################################################################################

set -euo pipefail

echo "========================================================================"
echo "Adult Unfed DESeq2 Analysis"
echo "========================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Project paths
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
CONTAINER="${PROJECT_BASE}/albopictus-diapause-rnaseq.sif"
SCRIPT="${PROJECT_BASE}/scripts/06_differential_expression/01_adults_deseq2_unfed.R"

# Create cache directories to avoid warnings
export TMPDIR="${PROJECT_BASE}/tmp"
export MPLCONFIGDIR="${PROJECT_BASE}/tmp/matplotlib"
export XDG_CACHE_HOME="${PROJECT_BASE}/tmp/cache"
export FONTCONFIG_PATH="${PROJECT_BASE}/tmp/fontconfig"

mkdir -p "${TMPDIR}"
mkdir -p "${MPLCONFIGDIR}"
mkdir -p "${XDG_CACHE_HOME}"
mkdir -p "${FONTCONFIG_PATH}"

echo "Project base: ${PROJECT_BASE}"
echo "Container: ${CONTAINER}"
echo "Script: ${SCRIPT}"
echo "Cache directory: ${TMPDIR}"
echo ""

# Load singularity
module load singularity-ce/3.9.3

# Verify container exists
if [ ! -f "${CONTAINER}" ]; then
    echo "ERROR: Container not found at ${CONTAINER}"
    exit 1
fi

# Verify script exists
if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: Script not found at ${SCRIPT}"
    exit 1
fi

echo "Starting DESeq2 analysis..."
echo "------------------------------------------------------------------------"

# Run analysis in container
singularity exec \
    --cleanenv \
    --no-home \
    --bind "${PROJECT_BASE}:${PROJECT_BASE}" \
    --bind "${TMPDIR}:${TMPDIR}" \
    "${CONTAINER}" \
    Rscript "${SCRIPT}"

EXIT_CODE=$?

echo "------------------------------------------------------------------------"
echo "Analysis completed with exit code: ${EXIT_CODE}"
echo "End time: $(date)"
echo "========================================================================"

exit ${EXIT_CODE}
