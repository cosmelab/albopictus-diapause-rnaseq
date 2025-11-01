#!/bin/bash
# Run QC analysis scripts using the Singularity container

# Load Singularity module
module load singularity/3.9.3

# Set project directory
PROJECT_DIR="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
CONTAINER="${PROJECT_DIR}/albopictus-diapause-rnaseq.sif"

# Check if container exists
if [ ! -f "$CONTAINER" ]; then
    echo "Error: Container not found at $CONTAINER"
    echo "Please build the container first"
    exit 1
fi

# Change to project directory
cd $PROJECT_DIR

echo "Starting QC analysis..."
echo "========================"

# Step 1: Extract QC metrics
echo "Step 1: Extracting QC metrics from all samples..."
singularity exec --cleanenv --bind $PROJECT_DIR:/proj $CONTAINER \
    bash -c "cd /proj && python scripts/07_qc_analysis/01_extract_qc_metrics.py"

# Check if extraction was successful
if [ -f "results/01_qc_analysis/metrics/all_samples_qc_metrics.csv" ]; then
    echo "✓ QC metrics extracted successfully"
    echo ""
    
    # Step 2: Generate figures
    echo "Step 2: Generating QC figures..."
    singularity exec --cleanenv --bind $PROJECT_DIR:/proj $CONTAINER \
        bash -c "cd /proj && python scripts/07_qc_analysis/02_generate_qc_figures.py"
    
    echo ""
    echo "QC analysis complete!"
    echo "===================="
    echo "Results saved to:"
    echo "  - results/01_qc_analysis/metrics/all_samples_qc_metrics.csv"
    echo "  - results/01_qc_analysis/metrics/per_sample/*.json"
    echo "  - results/01_qc_analysis/figures/*.pdf"
    echo "  - results/01_qc_analysis/tables/*.csv"
else
    echo "✗ Error: QC metrics extraction failed"
    exit 1
fi