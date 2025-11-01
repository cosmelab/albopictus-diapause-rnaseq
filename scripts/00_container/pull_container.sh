#!/bin/bash
#SBATCH --job-name=pull_container
#SBATCH -o logs/pull_container_%j.o.txt
#SBATCH -e logs/pull_container_%j.e.txt
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc
#SBATCH --account=cosmelab

################################################################################
# Pull Container from GitHub Container Registry
#
# Pulls the latest albopictus-diapause-rnaseq container with:
#   - rdryad package for automated Dryad downloads
#   - All collaborator packages (apeglm, enhancedvolcano, pheatmap, etc.)
#
# IMPORTANT: Sets SINGULARITY_CACHEDIR to avoid filling home quota (20GB limit)
#
# Date: October 24, 2025
################################################################################

set -euo pipefail

echo "========================================================================"
echo "Container Pull from GitHub Container Registry"
echo "========================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo ""

# Project paths
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
CONTAINER_NAME="albopictus-diapause-rnaseq.sif"
REGISTRY_URL="docker://ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest"

# CRITICAL: Set cache directory to avoid filling home quota (20GB limit)
export SINGULARITY_CACHEDIR="${PROJECT_BASE}/output/singularity"
mkdir -p "${SINGULARITY_CACHEDIR}"

echo "Project base: ${PROJECT_BASE}"
echo "Container output: ${CONTAINER_NAME}"
echo "Registry: ${REGISTRY_URL}"
echo "Cache directory: ${SINGULARITY_CACHEDIR}"
echo ""

# Load singularity module
module load singularity

# Verify singularity loaded
if ! command -v singularity &> /dev/null; then
    echo "ERROR: Singularity module failed to load"
    exit 1
fi

echo "Singularity version: $(singularity --version)"
echo ""

# Change to project directory
cd "${PROJECT_BASE}"

echo "Starting container pull..."
echo "------------------------------------------------------------------------"

# Pull container
singularity pull "${CONTAINER_NAME}" "${REGISTRY_URL}"

EXIT_CODE=$?

echo "------------------------------------------------------------------------"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ Container pulled successfully"
    echo ""
    echo "Container details:"
    ls -lh "${CONTAINER_NAME}"
    echo ""
    echo "Cache directory usage:"
    du -sh "${SINGULARITY_CACHEDIR}"
else
    echo "❌ Container pull failed with exit code: ${EXIT_CODE}"
fi

echo ""
echo "End time: $(date)"
echo "========================================================================"

exit ${EXIT_CODE}
