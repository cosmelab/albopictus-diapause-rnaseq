#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=36:00:00
#SBATCH --job-name=rnaseq_test
#SBATCH -o rnaseq_test_%j.o.txt
#SBATCH -e rnaseq_test_%j.ERROR.txt

# Load required modules
module load singularity-ce/3.9.3

# Source HPC configuration
source scripts/hpc_config_local.sh

echo "Job started at: $(date)"

# Set environment variables
export BASE_DIR="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
export CONTAINER_PATH="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif"
export BIND_PATHS="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq:/proj"
export NXF_CONFIG="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/nf-core/configs/hpc_batch.conf"

# Run the pipeline with environment variables
singularity exec --cleanenv --bind $BIND_PATHS $CONTAINER_PATH bash -c "
    cd /proj
    nextflow run nf-core/rnaseq \
        -c nf-core/configs/hpc_batch.conf \
        -profile standard \
        --outdir /proj/output/test-results \
        -work-dir /tmp/nf-work \
        --base_dir $BASE_DIR \
        --container_path $CONTAINER_PATH \
        --bind_paths $BIND_PATHS \
        --input data/metadata/test_samplesheet.csv \
        --fasta data/references/albo.fasta.gz \
        --gff data/references/genes.gff.gz \
        --aligner star_salmon \
        --trimmer fastp \
        --extra_fastp_args '--qualified_quality_phred 15 --length_required 35' \
        --skip_markduplicates \
        --skip_stringtie \
        --without-suid
"

echo "Pipeline finished at: $(date)"