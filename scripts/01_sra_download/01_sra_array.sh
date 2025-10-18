#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=48:00:00
#SBATCH --job-name=sra_download
#SBATCH --array=1-3
#SBATCH -o logs/sra_download_%A_%a.o.txt
#SBATCH -e logs/sra_download_%A_%a.ERROR.txt

# Load singularity module (customize for your HPC)
module load singularity-ce/3.9.3

# Create logs directory if it doesn't exist
mkdir -p logs

# Set paths - customize these for your environment
PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
CONTAINER_PATH="${CONTAINER_PATH:-${PROJECT_DIR}/albopictus-diapause-rnaseq.sif}"
BIND_PATH="${PROJECT_DIR}:/proj"

# Define the datasets
DATASETS=("PRJNA268379" "PRJNA158021" "PRJNA187045")
DATASET=${DATASETS[$((SLURM_ARRAY_TASK_ID-1))]}

echo "Job started at: $(date)"
echo "Processing dataset: $DATASET"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Project directory: $PROJECT_DIR"
echo "Container: $CONTAINER_PATH"

# Download this dataset
singularity exec --cleanenv --bind $BIND_PATH $CONTAINER_PATH bash -c "cd /proj && python scripts/01_download_data/02_sra_download.py --dataset $DATASET --output-dir data --max-workers 6"

echo "Dataset $DATASET completed at: $(date)"