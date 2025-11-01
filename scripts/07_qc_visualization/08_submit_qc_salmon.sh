#!/bin/bash
#SBATCH --job-name=qc_salmon
#SBATCH -o logs/qc_salmon_%j.o.txt
#SBATCH -e logs/qc_salmon_%j.e.txt
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=batch

# QC Visualization - Salmon
# Phase 2.5 - Quality control before differential expression
# Created: October 27, 2025

set -e
set -u

echo "======================================================================="
echo "QC VISUALIZATION - Salmon"
echo "======================================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo ""

# Load singularity
module load singularity

# Configuration
CONTAINER="albopictus-diapause-rnaseq.sif"
BIND_PATH="$PWD:/proj"
COUNT_TYPE="salmon"

# Create logs directory if needed
mkdir -p logs

# Run all QC scripts
echo "Running all QC scripts with Salmon data..."
echo ""

singularity exec --cleanenv --bind "$BIND_PATH" "$CONTAINER" \
  bash -c "cd /proj && bash scripts/07_qc_visualization/00_run_all_qc.sh salmon" <<< 'y'

echo ""
echo "======================================================================="
echo "QC VISUALIZATION COMPLETE"
echo "======================================================================="
echo "End time: $(date)"
echo ""
echo "Output directory: results/qc_visualization/salmon/"
echo ""
echo "Next steps:"
echo "  1. Review QC plots in results/qc_visualization/salmon/"
echo "  2. Check outlier report"
echo "  3. Document QC decisions"
echo "  4. Proceed to Phase 3 (Salmon validation)"
echo ""
