# RNA-seq Pipeline - nf-core/rnaseq

**Status:** ✅ Production-Ready - Completed Successfully Oct 26, 2025 (Job 20677337)

---

## TO REPRODUCE THIS ANALYSIS - RUN THIS:

```bash
# Submit the pipeline job (processes all 45 samples)
sbatch scripts/05_rnaseq_pipeline/02_run_rnaseq.sh
```

**That's it.** One command runs everything end-to-end, including MultiQC.

---

## What the Script Does

### 02_run_rnaseq.sh - **THE COMPLETE PIPELINE**
- Submits ONE Nextflow orchestrator job to SLURM
- Nextflow reads all 45 samples from samplesheet
- Nextflow submits its own child jobs for STAR, Salmon, featureCounts, MultiQC
- Processes all samples in parallel automatically
- **Runtime:** ~24-48 hours for all 45 samples
- **Output:** BAM files, Salmon counts, featureCounts, QC reports, MultiQC summary
- **Production-hardened:** Handles disk space and temp directory correctly

---

## Configuration Files

### params.yaml
**What it controls:** WHAT to analyze
- Which samples (samplesheet path)
- Which reference genome
- Which GTF annotation
- Pipeline options (rRNA filtering, etc.)

### hpc_batch.conf
**What it controls:** HOW MUCH resources each step gets
- Memory limits per process
- CPU counts
- Time limits
- SLURM queue/partition

### custom.config (in project root)
**What it controls:** Nextflow executor settings
- How to submit jobs to SLURM
- Container settings
- Work directory location

---

## Monitoring Progress

```bash
# Check the main orchestrator job
tail -f logs/rnaseq_main_*.o.txt

# Check what SLURM jobs Nextflow created
squeue -u $USER

# Check completion
sacct -j <JOBID>
```

---

## Output Files

```
output/
├── star_salmon/                   # Main outputs
│   ├── *.bam                      # Alignments (45 files)
│   ├── */quant.sf                 # Salmon transcript counts (45 samples)
│   └── featurecounts/*.tsv        # Gene counts (45 files, 22,177 genes each)
└── multiqc/
    └── multiqc_report.html        # QC summary (16MB)
```

**Key files for next steps:**
- `output/star_salmon/featurecounts/*.tsv` → Use for DESeq2 (collaborator replication)
- Salmon counts → Will be combined in next step (06_count_matrix)

---

## Troubleshooting

### Pipeline fails immediately
- Check: `tail logs/rnaseq_main_*.e.txt`
- Common issue: Missing reference files or bad paths in params.yaml

### SortMeRNA fails (rRNA indexing)
- Check: `data/references/AalbF3/rrna_db_manifest.txt` has absolute path
- Fix: Change relative path to `/bigdata/cosmelab/lcosme/projects/.../combined_rRNA_sequences.fa`

### Out of memory
- Edit `hpc_batch.conf` and increase memory for the failing process

### Disk space issues
- **Solution already implemented:** `hpc_batch.conf` sets `scratch = false` to use bigdata TMPDIR
- This prevents "No space left on device" errors that can occur with node-local /scratch

---

## Important Notes

**DO NOT delete these until completely done:**
- `.nextflow.log` - Needed for `-resume` to work
- `.nextflow/` - Cache directory
- `output/nf-work/` - Intermediate files (can delete after to save space)

**If you delete them:** Pipeline will restart from scratch, wasting days of compute time.

**After analysis is complete:** Can archive to `logs/nextflow_cache/` or delete work dir to save ~3TB

---

## Success Criteria

✅ All 45 samples processed successfully (Job 20677337)
✅ 22,177 genes detected per sample
✅ MultiQC report generated (output/multiqc/star_salmon/multiqc_report.html)
✅ No pipeline failures
✅ Production-hardened and ready for replication

**Next step:** `scripts/06_count_matrix/` to combine counts

---

## Pipeline Improvements Made

**Oct 26, 2025 - Production Hardening:**
- Fixed disk space issue by setting `scratch = false` in hpc_batch.conf
- Removed `--containall` flag to allow proper TMPDIR access
- Pipeline now uses bigdata TMPDIR (4.8TB available) instead of node-local /scratch
- MultiQC now completes successfully as part of main pipeline run
- No manual intervention required - fully automated end-to-end
