#!/bin/bash
#SBATCH --job-name=monitor_rnaseq
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/monitor_%j.o.txt
#SBATCH --error=/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/logs/monitor_%j.e.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lcosme@ucr.edu

#############################################################################
# RNA-seq Pipeline Monitor and Auto-Launcher
#
# Purpose: Monitor test job (array index 1), then automatically launch
#          remaining jobs (2-44) if successful, or email if failed
#
# Usage: sbatch scripts/03_rnaseq_pipeline/03_monitor_and_launch.sh <TEST_JOB_ID>
#
# Date: October 18, 2025
#############################################################################

set -euo pipefail

# Get test job ID from command line or use the current running test job
TEST_JOB_ID=${1:-20467416}
PROJECT_BASE="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
EMAIL="lcosme@ucr.edu"
CHECK_INTERVAL=300  # Check every 5 minutes

echo "=========================================="
echo "RNA-seq Pipeline Monitor Started"
echo "=========================================="
echo "Date: $(date)"
echo "Test Job ID: ${TEST_JOB_ID}"
echo "Email: ${EMAIL}"
echo "Project: ${PROJECT_BASE}"
echo "Check Interval: ${CHECK_INTERVAL} seconds"
echo "=========================================="
echo ""

cd ${PROJECT_BASE}

# Function to send email
send_email() {
    local subject="$1"
    local body="$2"
    echo -e "${body}" | mail -s "${subject}" ${EMAIL}
}

# Function to check if job is still running
is_job_running() {
    local job_id=$1
    if squeue -j ${job_id} 2>/dev/null | grep -q ${job_id}; then
        return 0  # Job is running
    else
        return 1  # Job is not running
    fi
}

# Function to check job exit status
get_job_status() {
    local job_id=$1
    # Get job status from sacct
    # Note: Job array format is JOBID_ARRAYINDEX
    local status=$(sacct -j ${job_id}_1 --format=State --noheader | head -1 | tr -d ' ')
    echo "${status}"
}

# Function to check if pipeline output looks successful
check_pipeline_success() {
    local sample_id="SRR458462"
    local project="PRJNA158021"
    local output_dir="${PROJECT_BASE}/output/${project}/${sample_id}"

    echo "Checking pipeline output directory: ${output_dir}"

    # Check if output directory exists
    if [ ! -d "${output_dir}" ]; then
        echo "ERROR: Output directory not found"
        return 1
    fi

    # Check for key output files
    local salmon_dir="${output_dir}/star_salmon"
    local multiqc_dir="${output_dir}/multiqc/star_salmon"

    # Check if Salmon quantification completed
    if [ ! -f "${salmon_dir}/salmon.merged.gene_counts.tsv" ]; then
        echo "ERROR: Salmon gene counts file not found"
        return 1
    fi

    # Check if gene counts are non-zero (critical check!)
    local gene_count=$(wc -l < "${salmon_dir}/salmon.merged.gene_counts.tsv")
    echo "Salmon gene count file has ${gene_count} lines"

    if [ ${gene_count} -lt 100 ]; then
        echo "ERROR: Salmon gene count file has too few entries (${gene_count})"
        return 1
    fi

    # Check if featureCounts completed (NEW - for validation)
    local featurecounts_file=$(find "${output_dir}/star_salmon" -name "*featureCounts.txt" 2>/dev/null | head -1)
    if [ -z "${featurecounts_file}" ]; then
        echo "ERROR: featureCounts output file not found"
        return 1
    fi

    local fc_count=$(wc -l < "${featurecounts_file}")
    echo "featureCounts file has ${fc_count} lines"

    if [ ${fc_count} -lt 100 ]; then
        echo "ERROR: featureCounts file has too few entries (${fc_count})"
        return 1
    fi

    # Check if MultiQC report exists
    if [ ! -f "${multiqc_dir}/multiqc_report.html" ]; then
        echo "WARNING: MultiQC report not found (not critical)"
    fi

    # Check log file for errors
    local log_file="${PROJECT_BASE}/logs/rnaseq_${TEST_JOB_ID}_1.o.txt"
    if grep -q "Pipeline completed with errors" "${log_file}"; then
        echo "ERROR: Pipeline log shows errors"
        return 1
    fi

    echo "Pipeline output looks good!"
    return 0
}

# Monitor the test job
echo "Monitoring test job ${TEST_JOB_ID}_1..."
echo "Waiting for job to complete..."
echo ""

while is_job_running ${TEST_JOB_ID}; do
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Test job still running... checking again in ${CHECK_INTERVAL}s"
    sleep ${CHECK_INTERVAL}
done

echo ""
echo "=========================================="
echo "Test job completed at $(date)"
echo "=========================================="
echo ""

# Wait a bit for files to be written
sleep 30

# Check job status
JOB_STATUS=$(get_job_status ${TEST_JOB_ID})
echo "Job status from SLURM: ${JOB_STATUS}"

# Check if job completed successfully
if [[ "${JOB_STATUS}" == "COMPLETED" ]]; then
    echo "SLURM reports job COMPLETED successfully"
    echo "Checking pipeline output files..."
    echo ""

    if check_pipeline_success; then
        echo ""
        echo "=========================================="
        echo "SUCCESS! Test job completed successfully!"
        echo "=========================================="
        echo ""
        echo "Launching remaining jobs (array indices 2-44)..."
        echo ""

        # Submit remaining jobs
        FULL_JOB_ID=$(sbatch --array=2-44 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh | awk '{print $NF}')

        echo "Submitted job array: ${FULL_JOB_ID}"
        echo "Array indices: 2-44 (43 samples)"
        echo ""
        echo "Monitor with: squeue -u $USER"
        echo "Check logs: tail -f logs/rnaseq_${FULL_JOB_ID}_*.o.txt"
        echo ""

        # Send success email
        EMAIL_SUBJECT="✅ RNA-seq Test PASSED - Full Array Launched"
        EMAIL_BODY="Success! Your RNA-seq test job (${TEST_JOB_ID}) completed successfully.

Test Sample: PRJNA158021_SRR458462
Status: COMPLETED
Gene counts detected: $(wc -l < output/PRJNA158021/SRR458462/star_salmon/salmon.merged.gene_counts.tsv) genes

The full array (samples 2-44) has been automatically launched:
Job ID: ${FULL_JOB_ID}
Array indices: 2-44
Total samples: 43

Monitor progress with:
  squeue -u $USER

View logs:
  ls -lht logs/rnaseq_${FULL_JOB_ID}_*.o.txt

Project directory:
  ${PROJECT_BASE}

This is an automated message from the RNA-seq pipeline monitor.
"
        send_email "${EMAIL_SUBJECT}" "${EMAIL_BODY}"

        echo "Success email sent to ${EMAIL}"

    else
        echo ""
        echo "=========================================="
        echo "FAILURE: Pipeline validation failed"
        echo "=========================================="
        echo ""

        # Send failure email with details
        EMAIL_SUBJECT="❌ RNA-seq Test FAILED - Pipeline Validation Error"
        EMAIL_BODY="ALERT: Your RNA-seq test job (${TEST_JOB_ID}) completed but failed validation checks.

Test Sample: PRJNA158021_SRR458462
SLURM Status: ${JOB_STATUS}
Issue: Pipeline output files missing or invalid

Please check:
1. Output directory: ${PROJECT_BASE}/output/PRJNA158021/SRR458462/
2. Log file: ${PROJECT_BASE}/logs/rnaseq_${TEST_JOB_ID}_1.o.txt
3. Error file: ${PROJECT_BASE}/logs/rnaseq_${TEST_JOB_ID}_1.ERROR.txt

The full array was NOT launched. Please investigate the issue before proceeding.

Check the monitor log for details:
  ${PROJECT_BASE}/logs/monitor_${SLURM_JOB_ID}.o.txt

This is an automated message from the RNA-seq pipeline monitor.
"
        send_email "${EMAIL_SUBJECT}" "${EMAIL_BODY}"

        echo "Failure email sent to ${EMAIL}"
        exit 1
    fi

else
    echo ""
    echo "=========================================="
    echo "FAILURE: Job did not complete successfully"
    echo "=========================================="
    echo "Job status: ${JOB_STATUS}"
    echo ""

    # Send failure email
    EMAIL_SUBJECT="❌ RNA-seq Test FAILED - Job Status: ${JOB_STATUS}"
    EMAIL_BODY="ALERT: Your RNA-seq test job (${TEST_JOB_ID}) failed.

Test Sample: PRJNA158021_SRR458462
SLURM Status: ${JOB_STATUS}

Please check the logs:
1. SLURM output: ${PROJECT_BASE}/logs/rnaseq_${TEST_JOB_ID}_1.o.txt
2. SLURM errors: ${PROJECT_BASE}/logs/rnaseq_${TEST_JOB_ID}_1.ERROR.txt
3. Nextflow log: ${PROJECT_BASE}/.nextflow.log

The full array was NOT launched. Please investigate the issue before proceeding.

Common issues:
- Check GTF file format
- Verify reference files exist
- Check for out-of-memory errors
- Review pipeline parameters

Project directory:
  ${PROJECT_BASE}

This is an automated message from the RNA-seq pipeline monitor.
"
    send_email "${EMAIL_SUBJECT}" "${EMAIL_BODY}"

    echo "Failure email sent to ${EMAIL}"
    exit 1
fi

echo ""
echo "=========================================="
echo "Monitor job completed at $(date)"
echo "=========================================="
