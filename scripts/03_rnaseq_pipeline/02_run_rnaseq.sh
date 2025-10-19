#!/bin/bash
#SBATCH --job-name=rnaseq_all
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=epyc
#SBATCH -o /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/rnaseq_main_%j.o.txt
#SBATCH -e /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/rnaseq_main_%j.e.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lcosme@ucr.edu

#############################################################################
# nf-core/rnaseq Pipeline - Production-Ready Launcher
#
# Purpose: Process all 44 samples with ONE Nextflow orchestrator
#          Nextflow submits individual SLURM jobs for each process
#
# Hardening Features:
#   - Version pinning for reproducibility
#   - Explicit work directory control
#   - Isolated container/temp directories
#   - Preflight validation mode
#   - Exit traps for audit trail
#   - Deterministic file permissions
#   - Bounded JVM and temp space
#
# Usage:
#   # Full run
#   sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
#
#   # Validation only (dry-run)
#   VALIDATE_ONLY=1 sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
#
# Date: October 18, 2025
#############################################################################

set -euo pipefail

#############################################################################
# Configuration
#############################################################################

# Project paths
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"

# Pipeline version pinning
PIPE_VER="3.19.0"
PIPE_VER_CLEAN="3_19_0"  # For run name (no periods allowed)
RUN_NAME="albopictus_diapause_rnaseq_${PIPE_VER_CLEAN}_$(date +%Y%m%d)"

# Validation mode (set VALIDATE_ONLY=1 to check params without running)
VALIDATE_ONLY="${VALIDATE_ONLY:-0}"

#############################################################################
# Environment Setup
#############################################################################

# Deterministic file permissions
umask 0022

# Load required modules
module load nextflow
module load singularity

# Nextflow JVM settings (bounded heap to avoid OOM)
export NXF_OPTS="-Xms1g -Xmx6g"

# Nextflow directories (isolated, project-specific)
export NXF_HOME="/bigdata/cosmelab/lcosme/.nextflow"
export NXF_WORK="${PROJECT_BASE}/output/nf-work"
export NXF_TEMP="${PROJECT_BASE}/output/nf-temp"

# Singularity settings (isolated cache and temp)
export SINGULARITY_CACHEDIR="${PROJECT_BASE}/output/singularity"
export SINGULARITY_TMPDIR="${PROJECT_BASE}/output/singularity_temp"
export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR}"

# Apptainer settings (future-proofing - Apptainer is Singularity's successor)
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR}"

# System temp directory override
export TMPDIR="${PROJECT_BASE}/output/tmp"

#############################################################################
# Create Required Directories
#############################################################################

mkdir -p "${PROJECT_BASE}/logs"
mkdir -p "${NXF_WORK}"
mkdir -p "${NXF_TEMP}"
mkdir -p "${SINGULARITY_CACHEDIR}"
mkdir -p "${SINGULARITY_TMPDIR}"
mkdir -p "${TMPDIR}"

#############################################################################
# Exit Trap for Provenance Logging
#############################################################################

finish() {
    EXIT_CODE=$?
    echo ""
    echo "==========================================="
    if [ ${EXIT_CODE} -eq 0 ]; then
        echo "✅ Pipeline completed successfully"
    else
        echo "❌ Pipeline failed with exit code ${EXIT_CODE}"
    fi
    echo "Finished at: $(date)"
    echo "Duration: ${SECONDS}s"
    echo "==========================================="
    exit ${EXIT_CODE}
}
trap finish EXIT

#############################################################################
# Provenance Logging
#############################################################################

cd "${PROJECT_BASE}"

echo "==========================================="
echo "nf-core/rnaseq Pipeline Launcher"
echo "==========================================="
echo "Run name: ${RUN_NAME}"
echo "Pipeline version: ${PIPE_VER}"
echo "Nextflow version: $(nextflow -version | head -n1)"
echo "Singularity version: $(singularity --version)"
echo "Started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: $(hostname)"
echo "User: $(whoami)"
echo ""
echo "Directories:"
echo "  Project: ${PROJECT_BASE}"
echo "  Work: ${NXF_WORK}"
echo "  Singularity cache: ${SINGULARITY_CACHEDIR}"
echo ""
echo "Input:"
echo "  Samplesheet: data/metadata/samplesheet.csv"
echo "  Samples: 44"
echo ""
echo "Validation mode: ${VALIDATE_ONLY}"
echo "==========================================="
echo ""

#############################################################################
# Build Nextflow Command
#############################################################################

NXF_CMD=(
    nextflow run nf-core/rnaseq
    -r "${PIPE_VER}"
    -name "${RUN_NAME}"
    -profile singularity
    -c scripts/03_rnaseq_pipeline/hpc_batch.conf
    -params-file scripts/03_rnaseq_pipeline/params.yaml
    -work-dir "${NXF_WORK}"
    -resume
    -with-report output/pipeline_report.html
    -with-timeline output/pipeline_timeline.html
    -with-trace output/pipeline_trace.txt
    -with-dag output/pipeline_dag.svg
    --multiqc_title "${RUN_NAME}"
)

# Note: All pipeline parameters are defined in params.yaml for reproducibility
# This ensures the exact run configuration is version-controlled and portable

#############################################################################
# Validation Mode
#############################################################################

if [ "${VALIDATE_ONLY}" = "1" ]; then
    echo "🔍 Running in VALIDATION mode (dry-run)"
    echo ""

    # Print full command
    echo "Command that would be executed:"
    echo "---"
    printf '%s \\\n' "${NXF_CMD[@]}"
    echo "---"
    echo ""

    # Validate samplesheet exists
    if [ ! -f data/metadata/samplesheet.csv ]; then
        echo "❌ ERROR: Samplesheet not found: data/metadata/samplesheet.csv"
        exit 1
    fi

    # Validate params file exists
    if [ ! -f scripts/03_rnaseq_pipeline/params.yaml ]; then
        echo "❌ ERROR: Parameters file not found: scripts/03_rnaseq_pipeline/params.yaml"
        exit 1
    fi

    # Validate references exist
    for REF in data/references/AalbF3_genome/AalbF3_genome.fa.gz \
               data/references/AalbF3_genome/AalbF3_annotation.gtf.gz \
               data/references/AalbF3_genome/rrna_db_manifest.txt; do
        if [ ! -f "${REF}" ]; then
            echo "❌ ERROR: Reference file not found: ${REF}"
            exit 1
        fi
    done

    echo "✅ All validation checks passed"
    echo ""
    echo "To run the pipeline, submit without VALIDATE_ONLY:"
    echo "  sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh"
    exit 0
fi

#############################################################################
# Execute Pipeline
#############################################################################

echo "🚀 Launching pipeline..."
echo ""

# Execute with full command array
"${NXF_CMD[@]}"

# Exit trap will handle completion message
