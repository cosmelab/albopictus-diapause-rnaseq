#!/bin/bash
# Session Management and Crash Recovery System
# Author: L. Cosme
# Date: October 2025
# Purpose: Track analysis sessions and enable easy recovery

BASE_DIR="/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq"
SESSION_DIR="${BASE_DIR}/logs/sessions"
CHECKPOINT_DIR="${BASE_DIR}/checkpoints"
TIMESTAMP=$(date +"%Y%m%d_%H%M")

# Function to start a new session
start_session() {
    SESSION_FILE="${SESSION_DIR}/${TIMESTAMP}_session_start.md"

    cat > ${SESSION_FILE} << EOF
# Analysis Session Started

**Date:** $(date +"%Y-%m-%d")
**Time:** $(date +"%H:%M:%S")
**User:** ${USER}
**Working Directory:** ${PWD}
**Container:** albopictus-diapause-rnaseq.sif

## Session Goals
$1

## Environment Check
- Singularity version: $(singularity --version 2>/dev/null || echo "Not found")
- R version: $(singularity exec ${BASE_DIR}/albopictus-diapause-rnaseq.sif R --version | head -1)
- Python version: $(singularity exec ${BASE_DIR}/albopictus-diapause-rnaseq.sif python --version)

## Active Checkpoints
$(ls -la ${CHECKPOINT_DIR}/*.checkpoint 2>/dev/null | tail -5 || echo "No checkpoints yet")

## Commands Log
Commands executed in this session will be logged to:
${SESSION_DIR}/${TIMESTAMP}_session_commands.sh

---
EOF

    echo "Session started: ${SESSION_FILE}"
    echo "# Commands executed on $(date)" > ${SESSION_DIR}/${TIMESTAMP}_session_commands.sh

    # Export session file for other scripts to use
    export CURRENT_SESSION="${TIMESTAMP}"
}

# Function to log commands
log_command() {
    if [ -z "${CURRENT_SESSION}" ]; then
        echo "No active session. Run: source session_manager.sh start 'Your goal here'"
        return 1
    fi

    echo "$@" >> ${SESSION_DIR}/${CURRENT_SESSION}_session_commands.sh
    echo "[$(date +%H:%M:%S)] Logged: $@"
}

# Function to create checkpoint
create_checkpoint() {
    CHECKPOINT_NAME=$1
    CHECKPOINT_FILE="${CHECKPOINT_DIR}/${TIMESTAMP}_${CHECKPOINT_NAME}.checkpoint"

    cat > ${CHECKPOINT_FILE} << EOF
Checkpoint: ${CHECKPOINT_NAME}
Date: $(date)
Session: ${CURRENT_SESSION}
Status: Complete

Files created:
$(find ${BASE_DIR}/results -type f -newermt "${TIMESTAMP:0:4}-${TIMESTAMP:4:2}-${TIMESTAMP:6:2}" 2>/dev/null | head -10)

Next steps:
$2
EOF

    echo "✓ Checkpoint created: ${CHECKPOINT_NAME}"
}

# Function to check analysis status
check_status() {
    echo "="
    echo "ANALYSIS STATUS REPORT"
    echo "="
    echo ""

    # Check data files
    echo "## Data Files"
    echo -n "Raw FASTQ files: "
    count=$(ls ${BASE_DIR}/data/raw/*/*.fastq.gz 2>/dev/null | wc -l)
    echo "${count} files"

    echo -n "Reference genome: "
    if [ -f "${BASE_DIR}/data/genome/AalbF3.fa" ]; then
        echo "✓ Present"
    else
        echo "✗ Missing"
    fi

    echo -n "GTF annotation: "
    if [ -f "${BASE_DIR}/data/annotation/AalbF3_final.gtf" ]; then
        echo "✓ Present"
    else
        echo "✗ Missing"
    fi

    echo ""
    echo "## Count Matrices"
    echo -n "Salmon counts: "
    if [ -f "${BASE_DIR}/results/count_matrices/salmon_gene_counts.tsv" ]; then
        echo "✓ Generated"
    else
        echo "✗ Not found"
    fi

    echo -n "FeatureCounts: "
    if [ -f "${BASE_DIR}/results/count_matrices/featurecounts_genelevel_counts.tsv" ]; then
        echo "✓ Generated"
    else
        echo "✗ Not found"
    fi

    echo -n "Stage-split counts: "
    for stage in adults embryos larvae; do
        if [ -f "${BASE_DIR}/results/count_matrices/${stage}_counts.tsv" ]; then
            echo -n "✓ ${stage} "
        else
            echo -n "✗ ${stage} "
        fi
    done
    echo ""

    echo ""
    echo "## Differential Expression"
    for stage in adults embryos larvae; do
        echo -n "${stage^} DE: "
        if [ -f "${BASE_DIR}/results/03_differential_expression/${stage}/deseq2_results_all.tsv" ]; then
            echo "✓ Complete"
        else
            echo "✗ Pending"
        fi
    done

    echo ""
    echo "## Recent Checkpoints"
    ls -lt ${CHECKPOINT_DIR}/*.checkpoint 2>/dev/null | head -5 | awk '{print $9}' | xargs -I {} basename {} .checkpoint

    echo ""
    echo "## Active Sessions"
    ls -lt ${SESSION_DIR}/*_session_start.md 2>/dev/null | head -3 | awk '{print $6, $7, $8, $9}'
}

# Function to recover from crash
recover_session() {
    echo "Analyzing recovery options..."

    # Find most recent session
    LAST_SESSION=$(ls -t ${SESSION_DIR}/*_session_start.md 2>/dev/null | head -1)

    if [ -z "${LAST_SESSION}" ]; then
        echo "No previous sessions found. Starting fresh."
        return 1
    fi

    echo "Last session: $(basename ${LAST_SESSION})"

    # Find last checkpoint
    LAST_CHECKPOINT=$(ls -t ${CHECKPOINT_DIR}/*.checkpoint 2>/dev/null | head -1)

    if [ -n "${LAST_CHECKPOINT}" ]; then
        echo "Last checkpoint: $(basename ${LAST_CHECKPOINT})"
        cat ${LAST_CHECKPOINT}
    fi

    # Show last commands executed
    SESSION_ID=$(basename ${LAST_SESSION} | cut -d_ -f1-2)
    COMMAND_FILE="${SESSION_DIR}/${SESSION_ID}_session_commands.sh"

    if [ -f "${COMMAND_FILE}" ]; then
        echo ""
        echo "Last 10 commands from previous session:"
        tail -10 ${COMMAND_FILE}
    fi

    echo ""
    echo "## Recovery Options:"
    echo "1. Continue from last checkpoint"
    echo "2. Start new session"
    echo "3. View full status report"

    echo ""
    echo "To continue, source this script with 'start' and your goal"
}

# Function to end session
end_session() {
    if [ -z "${CURRENT_SESSION}" ]; then
        echo "No active session"
        return 1
    fi

    SUMMARY_FILE="${SESSION_DIR}/${CURRENT_SESSION}_session_summary.md"

    cat > ${SUMMARY_FILE} << EOF
# Session Summary

**Session ID:** ${CURRENT_SESSION}
**End Time:** $(date)

## Completed Tasks
$(grep "create_checkpoint" ${SESSION_DIR}/${CURRENT_SESSION}_session_commands.sh 2>/dev/null | cut -d'"' -f2)

## Files Created
$(find ${BASE_DIR}/results -type f -newermt "${CURRENT_SESSION:0:4}-${CURRENT_SESSION:4:2}-${CURRENT_SESSION:6:2}" 2>/dev/null | wc -l) new files

## Commands Executed
$(wc -l ${SESSION_DIR}/${CURRENT_SESSION}_session_commands.sh | awk '{print $1}') commands

## Next Session Should
- Check status with: check_status
- Continue from latest checkpoint
- Review this summary: ${SUMMARY_FILE}
EOF

    echo "Session ended. Summary: ${SUMMARY_FILE}"
    unset CURRENT_SESSION
}

# Main logic
case "$1" in
    start)
        start_session "$2"
        ;;
    check)
        check_status
        ;;
    recover)
        recover_session
        ;;
    checkpoint)
        create_checkpoint "$2" "$3"
        ;;
    log)
        shift
        log_command "$@"
        ;;
    end)
        end_session
        ;;
    *)
        echo "Session Management System"
        echo ""
        echo "Usage:"
        echo "  source session_manager.sh start 'Goal description'  # Start new session"
        echo "  ./session_manager.sh check                          # Check status"
        echo "  ./session_manager.sh recover                        # Recover from crash"
        echo "  ./session_manager.sh checkpoint 'name' 'next steps' # Create checkpoint"
        echo "  ./session_manager.sh log 'command'                  # Log a command"
        echo "  source session_manager.sh end                       # End session"
        echo ""
        echo "Example:"
        echo "  source session_manager.sh start 'Run adult DE analysis'"
        echo "  ./session_manager.sh checkpoint 'adults_de_complete' 'Run embryos DE'"
        ;;
esac

# Export functions for use in current shell when sourced
if [ "${BASH_SOURCE[0]}" != "${0}" ]; then
    export -f log_command
    export -f create_checkpoint
fi