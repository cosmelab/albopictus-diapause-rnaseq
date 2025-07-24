#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12        # Increased to 12 CPUs per task
#SBATCH --mem-per-cpu=8gb        # 8GB per CPU
#SBATCH --time=24:00:00          # Reduced time since we're processing one sample
#SBATCH --job-name=rnaseq
#SBATCH --array=1-44             # Array of 44 samples
#SBATCH -o logs/rnaseq_%A_%a.o.txt    # %A is job ID, %a is array task ID
#SBATCH -e logs/rnaseq_%A_%a.ERROR.txt

# Load required modules
module load nextflow
module load singularity

# Create logs directory if it doesn't exist
mkdir -p logs

# Set up environment - use variables for portability
# Conda environment paths (customizable via hpc_config.sh)
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/bigdata/cosmelab/lcosme/conda/pkgs}"
export CONDA_ENVS_PATH="${CONDA_ENVS_PATH:-/bigdata/cosmelab/lcosme/conda/envs}"

# Project-specific paths
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export NXF_HOME="${NXF_HOME:-/bigdata/cosmelab/lcosme/.nextflow}"
export NXF_WORK="${PROJECT_ROOT}/work"
export NXF_TEMP="${PROJECT_ROOT}/temp"
export NXF_OPTS="-Xms512m -Xmx4g"
export NXF_SINGULARITY_CACHEDIR="${PROJECT_ROOT}/singularity"

# Create and set permissions for work directories
mkdir -p $NXF_WORK
mkdir -p $NXF_TEMP
mkdir -p $NXF_SINGULARITY_CACHEDIR
chmod -R 755 $NXF_WORK
chmod -R 755 $NXF_TEMP
chmod -R 755 $NXF_SINGULARITY_CACHEDIR

# Change to project root directory
cd "${PROJECT_ROOT}"

# Get the sample name and FASTQ paths for this array task
# Add 1 to SLURM_ARRAY_TASK_ID to skip the header line
SAMPLE_LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" data/metadata/samplesheet.csv)
SAMPLE=$(echo "$SAMPLE_LINE" | cut -d',' -f1)
FASTQ1=$(echo "$SAMPLE_LINE" | cut -d',' -f2)
FASTQ2=$(echo "$SAMPLE_LINE" | cut -d',' -f3)

# Extract project ID and sample ID from sample name (assuming format: PRJNA158021_SRR458462)
PROJECT=$(echo "$SAMPLE" | cut -d'_' -f1)
SAMPLE_ID=$(echo "$SAMPLE" | cut -d'_' -f2)

# Generate unique run ID
RUN_ID="${SAMPLE_ID}_${SLURM_ARRAY_TASK_ID}_$(date +%Y%m%d_%H%M%S)"

echo "Processing sample: $SAMPLE"
echo "Project: $PROJECT"
echo "Sample ID: $SAMPLE_ID"
echo "Run ID: $RUN_ID"
echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"

# Create a temporary samplesheet with just this sample
echo "sample,fastq_1,fastq_2,strandedness" > samplesheet_${SLURM_ARRAY_TASK_ID}.csv
echo "$SAMPLE_LINE" >> samplesheet_${SLURM_ARRAY_TASK_ID}.csv

# Check if FASTQ files exist
if [ ! -f "$FASTQ1" ]; then
    echo "Error: FASTQ1 file not found: $FASTQ1"
    exit 1
fi

if [ ! -f "$FASTQ2" ]; then
    echo "Error: FASTQ2 file not found: $FASTQ2"
    exit 1
fi

# Check if reference files exist
FASTA="data/references/albo.fasta.gz"
GTF="data/references/genes_fixed.rsem.fixed.gtf"

if [ ! -f "$FASTA" ]; then
    echo "Error: FASTA file not found: $FASTA"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "Error: GTF file not found: $GTF"
    exit 1
fi

# Create project directory structure
mkdir -p output/${PROJECT}

# Run the pipeline for this sample
nextflow run nf-core/rnaseq \
    -profile singularity \
    -c nf-core/configs/hpc_batch.conf \
    --input samplesheet_${SLURM_ARRAY_TASK_ID}.csv \
    --outdir output/${PROJECT}/${SAMPLE_ID} \
    --fasta "$FASTA" \
    --gtf "$GTF" \
    --skip_markduplicates true \
    --skip_bigwig true \
    --skip_stringtie true \
    --skip_qualimap false \
    --skip_rseqc false \
    --skip_deseq2_qc true \
    --skip_dupradar true \
    --skip_preseq true \
    --skip_salmon false \
    --skip_featurecounts false \
    --save_trimmed true \
    --save_unaligned true \
    --save_intermediates true \
    --max_memory 96.GB \
    --max_cpus 12 \
    -with-trace \
    -with-timeline \
    -with-dag \
    -name "${RUN_ID}"

# Check if pipeline completed successfully
if [ $? -eq 0 ] && [ -d "output/${PROJECT}/${SAMPLE_ID}" ]; then
    echo "Pipeline completed successfully for sample ${SAMPLE_ID}"
    
    # Create logs directory if it doesn't exist
    mkdir -p logs/samplesheets
    
    # Move the samplesheet to logs directory with timestamp
    mv samplesheet_${SLURM_ARRAY_TASK_ID}.csv logs/samplesheets/${SAMPLE_ID}_${RUN_ID}.csv
else
    echo "Pipeline failed or output directory not found. Check the logs for details."
    # Move the samplesheet to logs directory even if pipeline fails
    mkdir -p logs/samplesheets
    mv samplesheet_${SLURM_ARRAY_TASK_ID}.csv logs/samplesheets/${SAMPLE_ID}_${RUN_ID}_FAILED.csv
    exit 1
fi