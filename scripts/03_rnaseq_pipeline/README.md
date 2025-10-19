# 03 - RNA-seq Pipeline

**Purpose:** Run nf-core/rnaseq pipeline for quantification of all 44 samples.

**Status:** 🏃 RUNNING (October 18, 2025) - Using correct single Nextflow orchestrator approach (Job 20487053)

---

## Overview

This directory contains the script for running the nf-core/rnaseq v3.19.0 pipeline on HPCC with SLURM. The pipeline processes all 44 samples across 3 BioProjects.

## ⚠️ CRITICAL: The Correct Way to Run nf-core Pipelines

**WRONG (what we did for 10 months):**
```bash
# ❌ SLURM array job with one Nextflow instance per sample
sbatch --array=1-44 02_run_rnaseq_array.sh
# This runs 44 separate Nextflow instances → cache conflicts, wasted resources
```

**CORRECT (what we should have done):**
```bash
# ✅ ONE Nextflow orchestrator with all samples in samplesheet
sbatch 02_run_rnaseq.sh
# Nextflow submits its own SLURM jobs and manages parallelization
```

**Why the array approach was wrong:**
- Multiple Nextflow instances compete for same cache/work directories
- Lock file conflicts and task retries
- Rebuilds indices/containers repeatedly (wasted I/O)
- Defeats Nextflow's built-in parallelization

**How it should work:**
- Submit ONE orchestrator job (light: 4 CPUs, 16GB)
- Nextflow reads samplesheet with all 44 samples
- Nextflow submits individual SLURM jobs for each process (STAR, Salmon, featureCounts)
- Nextflow handles parallelization, caching, and resume automatically

## Current Approach (Corrected October 18, 2025)

**Single Nextflow Orchestrator:**
- ONE SLURM job runs Nextflow with all 44 samples in `data/metadata/samplesheet.csv`
- Nextflow's SLURM executor submits child jobs for each process
- Automatic parallelization across samples
- Proper caching and resume functionality

**Resources:**
- Orchestrator job: 4 CPUs, 16GB RAM, 72 hours
- Per-process jobs: Determined by nf-core pipeline labels
- Global caps: `--max_cpus 128 --max_memory 500.GB --max_time 72.h`

**Running:**
```bash
sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
```

---

## Scripts

### 01_create_samplesheet.py

**Purpose:** Generate nf-core/rnaseq compatible samplesheet with all 44 samples.

**Input:**
- `data/metadata/sample_info.csv`
- FASTQ files in `data/raw/PRJNA*/`

**Output:**
- `data/metadata/samplesheet.csv` (45 lines: 1 header + 44 samples)

**Samplesheet Format:**
```csv
sample,fastq_1,fastq_2,strandedness
PRJNA268379_SRR1719052,/path/to/SRR1719052_1.fastq.gz,/path/to/SRR1719052_2.fastq.gz,auto
```

**Usage:**
```bash
python scripts/03_rnaseq_pipeline/01_create_samplesheet.py
```

**Status:** ✅ Already run - samplesheet exists with all 44 samples

---

### 02_run_rnaseq.sh

**Purpose:** Single Nextflow orchestrator job that processes all 44 samples.

**Key Features:**
- ONE Nextflow instance with all samples in samplesheet
- Light orchestrator: 4 CPUs, 16GB RAM, 72 hour limit
- Nextflow submits child SLURM jobs for each process (STAR, Salmon, featureCounts)
- Automatic parallelization across samples
- Singularity containerization
- Proper caching and resume functionality

**Status:** ✅ Currently running (Job 20487198)

---

### 02_run_rnaseq_robust.sh

**Purpose:** Production-hardened launcher with all robustness and reproducibility improvements.

**Key Improvements Over Standard Version:**
- **Version Pinning:** Pipeline version locked as variable, run tagged with date
- **Profile Integration:** Uses `-profile slurm,singularity` for native SLURM executor
- **Reduced Orchestrator:** 2 CPUs, 8GB RAM (cheaper, child jobs get their own resources)
- **Explicit Paths:** `-work-dir` specified to prevent cache mixing across runs
- **Isolated Environments:** Dedicated temp/cache dirs prevent cross-user/run conflicts
- **Validation Mode:** `VALIDATE_ONLY=1` dry-run checks params before spending cluster time
- **Boolean Idiom:** Uses presence/absence instead of `=true/=false` (nf-core standard)
- **Exit Traps:** Captures success/failure status with timestamps for audit trails
- **Provenance Logging:** Records Nextflow version, paths, run metadata
- **Deterministic Permissions:** `umask 0022` ensures outputs readable across users
- **Bounded Resources:** JVM heap limits, controlled temp space prevent OOM/disk fills

**When to Use:**
- **New runs** on different HPCs (Yale, other collaborators)
- **Publication-quality** runs requiring full audit trail
- **Production deployments** where reproducibility is critical
- **Troubleshooting** with validation mode before full execution

**Usage:**
```bash
# Validation mode (check params, verify files exist)
VALIDATE_ONLY=1 sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq_robust.sh

# Full run
sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq_robust.sh
```

**Recommended For:**
- ✅ Yale HPC reproducibility test (tomorrow)
- ✅ Sharing with collaborators
- ✅ Future projects using nf-core pipelines

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
# Submit the orchestrator job (processes all 44 samples)
sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
```

**Logs:**
- Main orchestrator: `logs/rnaseq_main_<JOBID>.o.txt`
- Main errors: `logs/rnaseq_main_<JOBID>.e.txt`
- Nextflow log: `.nextflow.log`
- Per-process logs: `output/nf-work/`

**Runtime:** ~24-48 hours total for all 44 samples (parallel processing)

**Monitoring:**
```bash
# Watch orchestrator progress
tail -f logs/rnaseq_main_*.o.txt

# Check what SLURM jobs Nextflow has submitted
squeue -u $USER

# Check Nextflow work directory
ls -lht output/nf-work/
```

**⚠️ CRITICAL - Resume Functionality:**

Nextflow creates cache files needed for `-resume` to work:
- `.nextflow.log` - execution log (in project root)
- `.nextflow/` - cache directory (in project root)
- `output/nf-work/` - work directory with intermediate files

**DO NOT DELETE** these until analysis is completely done:
- If you delete them, `-resume` won't work
- Pipeline will start from scratch, wasting time/resources
- Keep them throughout the entire analysis

**After completely done:**
```bash
# Safe to move for archival (keeps resume capability)
mkdir -p logs/nextflow_cache
mv .nextflow.log logs/nextflow_cache/
mv .nextflow logs/nextflow_cache/

# Can delete work dir to save space (loses intermediate files but keeps results)
rm -rf output/nf-work/  # Saves ~3TB
```

**Lesson learned:** We deleted `.nextflow.log` after Run 3, breaking resume for Run 4. Had to delete 3.2TB cache and restart from scratch.

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
