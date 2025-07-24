#!/bin/bash
# HPC Configuration Template
# Copy this file to hpc_config.sh and customize for your environment
# Source this file before running analysis scripts: source hpc_config.sh

# =============================================================================
# User-Specific Paths (CUSTOMIZE THESE)
# =============================================================================

# Base directory for user data (usually your home or allocation directory)
export USER_BASE_DIR="/bigdata/cosmelab/lcosme"

# Conda environment paths (for nf-core pipeline dependencies)
export CONDA_PKGS_DIRS="${USER_BASE_DIR}/conda/pkgs"
export CONDA_ENVS_PATH="${USER_BASE_DIR}/conda/envs"

# Nextflow home directory (for pipeline caching)
export NXF_HOME="${USER_BASE_DIR}/.nextflow"

# =============================================================================
# HPC Module Configuration (CUSTOMIZE FOR YOUR SYSTEM)
# =============================================================================

# Load required modules for your HPC system
# Uncomment and modify as needed for your cluster
module load nextflow
module load singularity
# module load java/21.0.7  # If nextflow requires specific Java version

# Alternative module commands for different HPC systems:
# module load nextflow/22.10.1
# module load singularity-ce/3.9.3
# module load apptainer/1.1.0

# =============================================================================
# Container Configuration
# =============================================================================

# Container registry preference (choose one)
export CONTAINER_REGISTRY="ghcr.io"  # GitHub Container Registry
# export CONTAINER_REGISTRY="docker.io"  # Docker Hub

# Container image
export CONTAINER_IMAGE="cosmelab/albopictus-diapause-rnaseq:latest"

# =============================================================================
# Resource Limits (ADJUST FOR YOUR HPC ALLOCATION)
# =============================================================================

# SLURM partition/queue
export DEFAULT_PARTITION="epyc"  # Change to your cluster's partition name

# Maximum resources for nf-core pipeline
export MAX_MEMORY="96.GB"
export MAX_CPUS="12"
export MAX_TIME="24:00:00"

# =============================================================================
# Project Paths (USUALLY NO NEED TO CHANGE)
# =============================================================================

# These are automatically set relative to the project directory
# export PROJECT_ROOT  # Set automatically by scripts
# export NXF_WORK       # Set automatically to ${PROJECT_ROOT}/work
# export NXF_TEMP       # Set automatically to ${PROJECT_ROOT}/temp

# =============================================================================
# Nextflow Configuration
# =============================================================================

export NXF_OPTS="-Xms512m -Xmx4g"

# Singularity cache directory (for downloaded containers)
# export NXF_SINGULARITY_CACHEDIR  # Set automatically to ${PROJECT_ROOT}/singularity

# =============================================================================
# Validation
# =============================================================================

echo "🔧 HPC Configuration Loaded:"
echo "   User Base: ${USER_BASE_DIR}"
echo "   Conda Packages: ${CONDA_PKGS_DIRS}"
echo "   Conda Environments: ${CONDA_ENVS_PATH}"
echo "   Nextflow Home: ${NXF_HOME}"
echo "   Container Registry: ${CONTAINER_REGISTRY}"
echo "   Default Partition: ${DEFAULT_PARTITION}"
echo "   Max Memory: ${MAX_MEMORY}"
echo "   Max CPUs: ${MAX_CPUS}"
echo ""
echo "💡 To customize, copy this template:"
echo "   cp hpc_config.template.sh hpc_config.sh"
echo "   # Edit hpc_config.sh with your specific paths"
echo "   source hpc_config.sh"