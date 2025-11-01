#!/usr/bin/bash
#SBATCH --job-name=deseq2_de
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deseq2_%j.log
#SBATCH --error=logs/deseq2_%j.e.txt

################################################################################
# Run DESeq2 Differential Expression Analysis
#
# Purpose: Perform DE analysis on ComBat-seq corrected counts
#
# Prerequisites: ComBat-seq batch correction must be complete
#   - results/count_matrices/salmon_gene_counts_combat_corrected.tsv must exist
#
# Author: RNA-seq Analysis Pipeline
# Date: October 20, 2025
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load Singularity module
module purge
module load singularity/3.9.3

# Project paths
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
CONTAINER="${PROJECT_BASE}/albopictus-diapause-rnaseq.sif"
R_SCRIPT="${PROJECT_BASE}/scripts/08_differential_expression/01_run_deseq2.R"

# Check prerequisites
COMBAT_COUNTS="${PROJECT_BASE}/results/count_matrices/salmon_gene_counts_combat_corrected.tsv"

echo "========================================================================"
echo "DESeq2 Differential Expression Analysis"
echo "========================================================================"
echo ""
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}M"
echo ""
echo "Container: ${CONTAINER}"
echo "R Script: ${R_SCRIPT}"
echo ""
echo "========================================================================"
echo ""

# Check if ComBat-corrected counts exist
if [[ ! -f "${COMBAT_COUNTS}" ]]; then
    echo "ERROR: ComBat-seq corrected counts not found!"
    echo "Expected file: ${COMBAT_COUNTS}"
    echo ""
    echo "Please ensure ComBat-seq batch correction has completed successfully."
    exit 1
fi

# Check if container exists
if [[ ! -f "${CONTAINER}" ]]; then
    echo "ERROR: Container not found: ${CONTAINER}"
    exit 1
fi

# Check if R script exists
if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${R_SCRIPT}"
    exit 1
fi

# Display input file info
echo "Input files:"
echo "  Count matrix: ${COMBAT_COUNTS}"
ls -lh "${COMBAT_COUNTS}"
echo ""

# Display container info
echo "Container information:"
singularity exec "${CONTAINER}" R --version | head -3
echo ""
echo "Key R packages:"
singularity exec "${CONTAINER}" Rscript -e "cat('  DESeq2:', as.character(packageVersion('DESeq2')), '\n')"
singularity exec "${CONTAINER}" Rscript -e "cat('  ggplot2:', as.character(packageVersion('ggplot2')), '\n')"
singularity exec "${CONTAINER}" Rscript -e "cat('  pheatmap:', as.character(packageVersion('pheatmap')), '\n')"
echo ""
echo "========================================================================"
echo ""

# Run DESeq2 analysis
echo "Running DESeq2 differential expression analysis..."
echo ""

singularity exec \
    --bind "${PROJECT_BASE}:${PROJECT_BASE}" \
    "${CONTAINER}" \
    Rscript "${R_SCRIPT}"

echo ""
echo "========================================================================"
echo "DESeq2 Analysis Complete!"
echo "========================================================================"
echo ""

# Display output files
echo "Output files generated:"
echo ""
echo "Results directory:"
ls -lh "${PROJECT_BASE}/results/differential_expression/" 2>/dev/null || echo "  (pending)"
echo ""
echo "Figures directory:"
ls -lh "${PROJECT_BASE}/results/de_figures/" 2>/dev/null || echo "  (pending)"
echo ""

echo "Next steps:"
echo "  1. Review DE results: results/differential_expression/deseq2_results_significant.tsv"
echo "  2. Check visualizations: results/de_figures/*.pdf"
echo "  3. Look up GWAS candidate genes in results"
echo "  4. Run GO enrichment analysis"
echo ""
echo "========================================================================"