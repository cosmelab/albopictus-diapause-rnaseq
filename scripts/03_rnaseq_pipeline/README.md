# 03 - RNA-seq Pipeline

**Purpose:** Run nf-core/rnaseq pipeline for quantification of all 44 samples.

**Status:** 🔄 Re-running (October 18, 2025) - Added featureCounts validation + fixed parallel execution conflicts

---

## Overview

This directory contains scripts for running the nf-core/rnaseq v3.19.0 pipeline on HPCC using SLURM array jobs. The pipeline processes 44 samples across 3 BioProjects in parallel.

## Pipeline Approach

**Strategy:** One sample per array job for maximum parallelization and resource efficiency.

**Benefits:**
- Each sample gets dedicated resources (12 CPUs, 96GB RAM)
- Failed samples don't affect others
- Easy to re-run individual samples
- Clear per-sample logs and outputs

### Two Execution Strategies

**Option 1: Test-Then-Launch (Conservative, Recommended)**
```bash
# Submit test job (sample 1)
sbatch --array=1 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh

# Use monitor script to validate and auto-launch remaining samples
sbatch scripts/03_rnaseq_pipeline/03_monitor_and_launch.sh <TEST_JOB_ID>
```
- **Use when:** First time running, testing new parameters, or validating pipeline changes
- **Benefit:** Catches errors early before wasting resources on all 44 samples
- **Time:** Adds ~2-4 hours for test + validation

**Option 2: Run All At Once (Fast, Time-Constrained)**
```bash
# Submit all 44 samples immediately
sbatch --array=1-44 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh
```
- **Use when:** Confident in configuration, time-constrained, or re-running with known-good parameters
- **Benefit:** Saves 2-4 hours by skipping test phase
- **Risk:** If configuration wrong, wastes resources on all samples
- **Current run (Oct 18):** Using this approach due to Christmas deadline

### Parallel Execution Fix (October 18, 2025)

**Problem Discovered:**
When running all 44 samples in parallel, Nextflow instances conflicted on shared directories:
- All jobs tried to lock same `.nextflow/cache/` directory → lock file errors
- All jobs wrote to same `.nextflow.log` in project root → file conflicts

**Solution Implemented:**
Each array task now gets completely isolated directories:

1. **Unique Nextflow home per task:**
   ```bash
   export NXF_HOME="${PROJECT_BASE}/output/nf-home-${SLURM_ARRAY_TASK_ID}"
   ```
   - Task 1 → `output/nf-home-1/.nextflow/`
   - Task 2 → `output/nf-home-2/.nextflow/`
   - Prevents lock file conflicts

2. **Unique working directory per task:**
   ```bash
   TASK_WORK_DIR="${PROJECT_BASE}/output/task-${SLURM_ARRAY_TASK_ID}"
   cd "${TASK_WORK_DIR}"
   ```
   - Task 1 cd's to `output/task-1/` → creates `.nextflow.log` there
   - Task 2 cd's to `output/task-2/` → creates `.nextflow.log` there
   - Prevents log file conflicts

3. **All paths made absolute:**
   - Since tasks no longer run from project root, all paths use `${PROJECT_BASE}/`

**Resume Still Works:**
- Each task's `SLURM_ARRAY_TASK_ID` is consistent across re-runs
- Task 5 always uses `output/task-5/` and `output/nf-home-5/`
- Resume finds the same cache directory automatically

---

## Scripts (Run in Order)

### 01_create_samplesheet.py

**Purpose:** Generate nf-core/rnaseq compatible samplesheet with all 44 samples.

**Input:**
- `data/metadata/sample_info.csv`
- FASTQ files in `data/raw_fastq/PRJNA*/`

**Output:**
- `data/metadata/samplesheet.csv`

**Samplesheet Format:**
```csv
sample,fastq_1,fastq_2,strandedness
PRJNA268379_SRR1719052,/path/to/SRR1719052_1.fastq.gz,/path/to/SRR1719052_2.fastq.gz,auto
```

**Usage:**
```bash
python scripts/03_rnaseq_pipeline/01_create_samplesheet.py
```

---

### 02_run_rnaseq_array.sh

**Purpose:** Main SLURM array job script to run nf-core/rnaseq pipeline for each sample.

**Key Features:**
- SLURM array job (1-44 samples)
- 12 CPUs per task, 8GB per CPU (96GB total RAM)
- 24-hour time limit per sample
- Automatic per-sample samplesheet creation
- Singularity containerization

**Pipeline Configuration:**

**Reference Files:**
```bash
FASTA="data/references/AalbF3_genome/AalbF3_genome.fa.gz"
GTF="data/references/AalbF3_genome/AalbF3_annotation.gtf.gz"
RRNA_MANIFEST="data/references/AalbF3_genome/rrna_db_manifest.txt"
```

**Key Parameters:**
- `--remove_ribo_rna` - Filter rRNA with SortMeRNA
- `--save_reference` - Save genome index for reuse
- `--skip_markduplicates true` - Skip duplicate marking (RNA-seq)
- `--skip_bigwig true` - Skip BigWig generation (large files)
- `--skip_stringtie true` - Skip transcript assembly
- `--skip_salmon false` - **Salmon quantification (primary method)**
- `--skip_featurecounts false` - **featureCounts enabled for validation**
- `--save_trimmed true` - Save trimmed FASTQ for QC
- `--save_unaligned true` - Save unaligned reads for troubleshooting
- `-resume` - Resume from cache if available (safe to include always)

**Important:** Configuration file `hpc_batch.conf` contains HPC-specific settings.

**Usage:**
```bash
# Test with one sample first
sbatch --array=1 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh

# Run all 44 samples
sbatch --array=1-44 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh
```

**Logs:**
- SLURM output: `logs/rnaseq_<JOBID>_<ARRAYID>.o.txt`
- SLURM errors: `logs/rnaseq_<JOBID>_<ARRAYID>.ERROR.txt`

**Runtime:** ~2 hours per sample (varies by sample size)

---

### 03_monitor_and_launch.sh

**Purpose:** Automated monitoring script - validates test job and auto-launches remaining samples.

**Workflow:**
1. Monitor test job (sample 1) until completion
2. Check SLURM job status (COMPLETED, FAILED, etc.)
3. Validate pipeline output:
   - Gene counts file exists
   - Gene count > 100 lines
   - No critical errors in log
4. If successful: Auto-launch jobs 2-44
5. If failed: Send email alert with error details
6. Send status email either way

**Key Validation Checks:**
```bash
# Check Salmon gene counts file
output/PRJNA158021/SRR458462/star_salmon/salmon.merged.gene_counts.tsv
# Expected: ~20,000-25,000 genes

# Check featureCounts output (NEW - for validation)
output/PRJNA158021/SRR458462/star_salmon/*featureCounts.txt
# Expected: ~20,000-25,000 genes

# Both methods should detect similar gene counts
# Previous runs with wrong GTF: 0 genes (FAILURE)
# Successful run: 22,176 genes (SUCCESS)
```

**Usage:**
```bash
# First, submit test job
TEST_JOB_ID=$(sbatch --array=1 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh | awk '{print $NF}')

# Then submit monitor job
sbatch scripts/03_rnaseq_pipeline/03_monitor_and_launch.sh $TEST_JOB_ID
```

**Email Notifications:**
- ✅ Success: "RNA-seq Test PASSED - Full Array Launched"
- ❌ Failure: "RNA-seq Test FAILED - [error details]"

**Actual Run (October 17, 2025):**
- 16:02 - Test job submitted (Job 20467416, sample 1)
- 17:32 - Monitor job submitted (Job 20468449)
- 17:57 - Test completed successfully (22,176 genes detected!)
- 17:58 - Monitor auto-launched remaining 43 samples (Job 20468495)
- 18:00 - All 44 samples processing in parallel

---

## HPC Configuration

### Critical Environment Variables

**Must set before running:**
```bash
export SINGULARITY_TMPDIR="${PROJECT_BASE}/output/singularity_temp"
mkdir -p $SINGULARITY_TMPDIR
```

**Why:** Singularity on HPC requires writable temp directory. Default `/tmp` often too small.

### Resource Allocation

**Per Sample:**
- CPUs: 12
- Memory: 96 GB (8GB per CPU)
- Time: 24 hours
- Partition: batch (default)

**Total for 44 samples:**
- CPUs: 528 (but run 10-15 at a time based on HPC availability)
- Wall time: ~2-4 hours for full array completion

---

## Output Structure

```
output/
├── PRJNA268379/
│   ├── SRR1719052/
│   │   ├── star_salmon/
│   │   │   ├── salmon.merged.gene_counts.tsv      # Gene counts
│   │   │   ├── salmon.merged.gene_tpm.tsv         # Gene TPM
│   │   │   └── quant.sf                            # Transcript quantification
│   │   ├── multiqc/
│   │   │   └── star_salmon/
│   │   │       ├── multiqc_report.html
│   │   │       └── multiqc_data/
│   │   ├── trimgalore/                             # Trimmed FASTQ
│   │   └── fastqc/                                 # QC reports
│   └── [... 15 more samples]
├── PRJNA158021/
│   └── [12 samples]
└── PRJNA187045/
    └── [16 samples]
```

**Key Output Files Per Sample:**
- `salmon.merged.gene_counts.tsv` - Gene-level counts (for DESeq2)
- `salmon.merged.gene_tpm.tsv` - Gene-level TPM (for visualization)
- `quant.sf` - Transcript-level quantification (for isoform analysis)
- `multiqc_report.html` - Comprehensive QC report

---

## Monitoring Progress

```bash
# Check running jobs
squeue -u $USER

# Check specific array job
squeue -j 20468495

# Watch logs in real-time
tail -f logs/rnaseq_20468495_2.o.txt

# Check completion status
sacct -j 20468495 --format=JobID,State,Elapsed,MaxRSS

# Count completed samples
ls -d output/PRJNA*/SRR*/star_salmon/ | wc -l
```

---

## Troubleshooting

### Common Issues

**1. Singularity temp directory errors**
```bash
export SINGULARITY_TMPDIR="${PROJECT_BASE}/output/singularity_temp"
```

**2. GTF validation failures**
- Ensure you ran all 3 steps in `02_annotation_prep/`
- Check GTF has gene_id and gene_biotype on all lines

**3. Zero genes detected**
- Previous issue with wrong annotation (fixed)
- Test job now shows 22,176 genes ✅

**4. Out of memory**
- Increase `--mem-per-cpu` in SLURM header
- Or reduce `--max_cpus` in nextflow params

**5. Sample fails but others succeed**
- Re-run specific sample: `sbatch --array=5 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh`
- Check sample-specific log for errors

---

## Next Step

After pipeline completion:
- **Proceed to:** `04_qc_analysis/` to extract and visualize QC metrics

---

## Success Metrics (October 17-18, 2025)

✅ **All 44 samples completed successfully**
- Test sample (SRR458462): 22,176 genes detected
- No pipeline failures
- All output files present
- MultiQC reports generated
- Total runtime: ~10 hours for all samples

**Key Achievement:** Fixed all 6 GTF annotation problems - zero gene detection issue resolved!

---

## Notes

- Container images cached in `output/singularity/` (3.8GB - reusable)
- Work directory `output/nf-work/` can be deleted after completion (saves ~2TB)
- Temporary samplesheets saved to `logs/samplesheets/` for reproducibility
- Pipeline version locked: nf-core/rnaseq v3.19.0
- All parameters documented in script for reproducibility
