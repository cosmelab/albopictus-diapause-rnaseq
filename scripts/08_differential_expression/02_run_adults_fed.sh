#!/bin/bash
#SBATCH --job-name=adults_fed_DE
#SBATCH --output=logs/adults_fed_deseq2_%j.out
#SBATCH --error=logs/adults_fed_deseq2_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc
#SBATCH --account=cosmelab

set -euo pipefail

echo "========================================================================"
echo "Adult Fed DESeq2 Analysis"
echo "========================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo ""

PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
CONTAINER="${PROJECT_BASE}/albopictus-diapause-rnaseq.sif"
SCRIPT="${PROJECT_BASE}/scripts/08_differential_expression/02_adults_deseq2_fed.R"

export TMPDIR="${PROJECT_BASE}/tmp"
export MPLCONFIGDIR="${PROJECT_BASE}/tmp/matplotlib"
export XDG_CACHE_HOME="${PROJECT_BASE}/tmp/cache"
export FONTCONFIG_PATH="${PROJECT_BASE}/tmp/fontconfig"

mkdir -p "${TMPDIR}" "${MPLCONFIGDIR}" "${XDG_CACHE_HOME}" "${FONTCONFIG_PATH}"

module load singularity-ce/3.9.3

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
