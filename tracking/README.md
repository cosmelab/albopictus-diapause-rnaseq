# Tracking System

This directory contains the project tracking system for managing analysis progress.

## Files in this directory

- `project_state.json` - Database tracking all phases and tasks
- `00_track.py` - CLI tool for checking status and marking tasks complete
- `visual_progress.py` - Generates visual progress dashboard
- `visual_progress.txt` - Generated progress dashboard (updated on demand)
- `project_summary.md` - Complete project information and experimental design
- `README.md` - This file

## How to use the tracking system

### Check current status
```bash
module load singularity
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py status"
```

### Get next task
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py next"
```

### Mark task complete
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py complete <task_id>"
```

### Generate visual progress dashboard
```bash
module load singularity
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/visual_progress.py"
```

This updates `visual_progress.txt` with current progress.

## What information is tracked

The tracking system monitors:
- **Phases**: Major project stages (data prep, QC, differential expression, etc.)
- **Tasks**: Individual steps within each phase
- **Scripts**: Which script files are used for each task
- **Output directories**: Where each task's results are saved
- **Status**: pending, in_progress, failed, or complete for each task
- **Reviewer requirements**: Critical deliverables for manuscript revision

## Files generated

All files in this directory are version-controlled except temporary outputs.
