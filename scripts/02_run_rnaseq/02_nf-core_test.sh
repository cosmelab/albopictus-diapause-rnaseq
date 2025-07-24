#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=36:00:00
#SBATCH --partition=week
#SBATCH --job-name=rnaseq_test
#SBATCH -o rnaseq_test_%j.o.txt
#SBATCH -e rnaseq_test_%j.ERROR.txt

CONTAINER_PATH="/ycga-gpfs/project/caccone/lvc26/docker/diapause-rnaseq.sif"
BIND_PATH="/ycga-gpfs/project/caccone/lvc26/docker:/proj"

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

echo "Starting nf-core/rnaseq pipeline in container..."

# Run nf-core/rnaseq inside the container
apptainer exec --bind $BIND_PATH $CONTAINER_PATH bash -c "
    cd /proj
    mkdir -p logs
    rm -rf work/
    
    nextflow run nf-core/rnaseq \
        -c nf-core/configs/hpc_batch.conf \
        -work-dir work \
        -resume
"

echo "Pipeline finished at: $(date)"

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo "✅ Pipeline completed successfully!"
    echo "Results are in: output/"
    echo "MultiQC report: nf-core/results/multiqc/multiqc_report.html"
else
    echo "❌ Pipeline failed. Check logs for details."
    echo "Log file: logs/rnaseq_test_${SLURM_JOB_ID}.err"
fi
