#!/bin/bash
#
# Master Script to Run All QC Visualization
#
# Purpose: Execute all QC scripts in correct order with container
#          This is Phase 2.5 - CRITICAL QC before differential expression
#
# Usage:
#   bash scripts/07_qc_visualization/RUN_ALL_QC.sh featurecounts
#   bash scripts/07_qc_visualization/RUN_ALL_QC.sh salmon
#
# Created: October 27, 2025
# Phase: 2.5 - QC Visualization

set -e  # Exit on error
set -u  # Exit on undefined variable

# =============================================================================
# Configuration
# =============================================================================

COUNT_TYPE=${1:-featurecounts}

if [[ "$COUNT_TYPE" != "featurecounts" && "$COUNT_TYPE" != "salmon" ]]; then
    echo "ERROR: Invalid count type. Use 'featurecounts' or 'salmon'"
    echo "Usage: bash $0 [featurecounts|salmon]"
    exit 1
fi

# Paths
SCRIPT_DIR="scripts/07_qc_visualization"
CONTAINER="albopictus-diapause-rnaseq.sif"
BIND_PATH="$PWD:/proj"

# =============================================================================
# Functions
# =============================================================================

run_r_script() {
    local script=$1
    local description=$2

    echo ""
    echo "======================================================================="
    echo "$description"
    echo "======================================================================="
    echo ""

    singularity exec --cleanenv --bind "$BIND_PATH" "$CONTAINER" \
        bash -c "cd /proj && Rscript $script $COUNT_TYPE"

    if [ $? -ne 0 ]; then
        echo "ERROR: Script failed: $script"
        exit 1
    fi
}

run_python_script() {
    local script=$1
    local description=$2

    echo ""
    echo "======================================================================="
    echo "$description"
    echo "======================================================================="
    echo ""

    singularity exec --cleanenv --bind "$BIND_PATH" "$CONTAINER" \
        bash -c "cd /proj && python $script"

    if [ $? -ne 0 ]; then
        echo "ERROR: Script failed: $script"
        exit 1
    fi
}

# =============================================================================
# Main
# =============================================================================

echo ""
echo "======================================================================="
echo "QC VISUALIZATION - ALL SCRIPTS"
echo "======================================================================="
echo ""
echo "Count type: $COUNT_TYPE"
echo "Container: $CONTAINER"
echo ""
echo "This will run 6 QC scripts in order:"
echo "  01. Global QC (all 45 samples)"
echo "  02. Adults QC (16 samples)"
echo "  03. Embryos QC (12 samples)"
echo "  04. Larvae QC (17 samples)"
echo "  05. MultiQC metrics extraction"
echo "  06. Outlier detection"
echo ""
echo "Estimated runtime: 15-20 minutes"
echo ""

read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted"
    exit 0
fi

START_TIME=$(date +%s)

# Run scripts in order
run_r_script "$SCRIPT_DIR/01_global_qc_plots.R" "01. Global QC - All 45 Samples"

run_r_script "$SCRIPT_DIR/02_adults_qc_plots.R" "02. Adults QC - 16 Samples"

run_r_script "$SCRIPT_DIR/03_embryos_qc_plots.R" "03. Embryos QC - 12 Samples"

run_r_script "$SCRIPT_DIR/04_larvae_qc_plots.R" "04. Larvae QC - 17 Samples"

# MultiQC integration (runs once, not per count type)
if [[ "$COUNT_TYPE" == "featurecounts" ]]; then
    run_python_script "$SCRIPT_DIR/05_extract_multiqc_metrics.py" "05. MultiQC Integration"
fi

run_r_script "$SCRIPT_DIR/06_outlier_detection.R" "06. Outlier Detection"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "======================================================================="
echo "QC VISUALIZATION COMPLETE"
echo "======================================================================="
echo ""
echo "Runtime: $((DURATION / 60)) minutes $((DURATION % 60)) seconds"
echo "Count type: $COUNT_TYPE"
echo ""
echo "Output directory: results/qc_visualization/$COUNT_TYPE"
echo ""
echo "NEXT STEPS:"
echo "  1. Review all PDF plots in results/qc_visualization/$COUNT_TYPE/"
echo "  2. Read outlier report: results/qc_visualization/$COUNT_TYPE/outlier_report.txt"
echo "  3. If outliers found, document decisions in qc_decisions.md"
echo "  4. Run for other count type if needed"
echo "  5. Proceed to Phase 3 (Salmon validation) only after QC review"
echo ""
