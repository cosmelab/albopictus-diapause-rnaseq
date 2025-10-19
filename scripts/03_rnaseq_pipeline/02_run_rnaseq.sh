#!/bin/bash
#SBATCH --job-name=rnaseq_all
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=epyc
#SBATCH -o /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/rnaseq_main_%j.o.txt
#SBATCH -e /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/rnaseq_main_%j.e.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lcosme@ucr.edu

#############################################################################
# nf-core/rnaseq Pipeline - Single Nextflow Orchestrator
#
# Purpose: Process all 44 samples with ONE Nextflow instance
#          Nextflow submits individual SLURM jobs for each process
#
# Usage: sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
#
# Date: October 18, 2025
#############################################################################

set -euo pipefail

# Load required modules
module load nextflow
module load singularity

# Set project base directory
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"

# Set up Nextflow environment
export NXF_OPTS="-Xms1g -Xmx4g"
export NXF_HOME="/bigdata/cosmelab/lcosme/.nextflow"
export NXF_WORK="${PROJECT_BASE}/output/nf-work"
export NXF_SINGULARITY_CACHEDIR="${PROJECT_BASE}/output/singularity"
export SINGULARITY_TMPDIR="${PROJECT_BASE}/output/singularity_temp"

# Create directories
mkdir -p ${PROJECT_BASE}/logs
mkdir -p $NXF_WORK
mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

# Change to project directory
cd "${PROJECT_BASE}"

echo "==========================================="
echo "nf-core/rnaseq Pipeline Started"
echo "==========================================="
echo "Date: $(date)"
echo "Samples: 44 (all datasets)"
echo "Samplesheet: data/metadata/samplesheet.csv"
echo "Output: output/"
echo "==========================================="
echo ""

# Run the pipeline with ALL 44 samples
# Nextflow will submit individual SLURM jobs for each process
nextflow run nf-core/rnaseq \
    -r 3.19.0 \
    -profile singularity \
    -c scripts/03_rnaseq_pipeline/hpc_batch.conf \
    -resume \
    -with-report output/pipeline_report.html \
    -with-timeline output/pipeline_timeline.html \
    -with-trace output/pipeline_trace.txt \
    --input data/metadata/samplesheet.csv \
    --outdir output \
    --fasta data/references/AalbF3_genome/AalbF3_genome.fa.gz \
    --gtf data/references/AalbF3_genome/AalbF3_annotation.gtf.gz \
    --remove_ribo_rna \
    --ribo_database_manifest data/references/AalbF3_genome/rrna_db_manifest.txt \
    --aligner star_salmon \
    --save_reference \
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
    --max_cpus 128 \
    --max_memory '500.GB' \
    --max_time '72.h'

echo ""
echo "==========================================="
echo "Pipeline completed at $(date)"
echo "==========================================="
