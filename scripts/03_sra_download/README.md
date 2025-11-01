# 01 - SRA Download

**Purpose:** Download 44 RNA-seq samples from NCBI SRA across 3 BioProjects.

**Status:** ✅ Completed

---

## Overview

This directory contains scripts for downloading raw FASTQ files from the Sequence Read Archive (SRA) for three published RNA-seq studies on *Aedes albopictus* diapause.

## Datasets

### PRJNA268379 - Adult Females (Huang et al. 2015)
- **Samples:** 16 (4 conditions × 4 replicates)
- **Design:** 2×2 factorial (Photoperiod × Blood meal)
- **Platform:** Illumina HiSeq2000
- **Read Length:** Paired-end
- **Life Stage:** Adult females

### PRJNA158021 - Embryos (Poelchau et al. 2013a)
- **Samples:** 12 (4 timepoint-condition combinations × 3 replicates)
- **Design:** Time course (72-78h, 135-141h post-oviposition)
- **Platform:** Illumina GAIIx
- **Read Length:** Paired-end
- **Life Stage:** Embryos

### PRJNA187045 - Pharate Larvae (Poelchau et al. 2013b)
- **Samples:** 16 (6 timepoint-condition combinations × variable replicates)
- **Design:** Time course (11, 21, 40 days post-oviposition)
- **Platform:** Illumina HiSeq2000
- **Read Length:** Paired-end
- **Life Stage:** Pharate larvae (diapause vs quiescence)

**Total Samples:** 44

---

## Sample Metadata

Sample metadata is stored in:
- `data/metadata/samplesheet.csv` (main samplesheet for nf-core/rnaseq)
- `data/metadata/sample_info.csv` (detailed experimental metadata)

---

## Scripts (Run in Order)

### 01_sra_array.sh

**Purpose:** SLURM array job script to download FASTQ files in parallel from SRA.

**Requirements:**
- SRA Toolkit module loaded on HPC
- Sample list in `data/metadata/samplesheet.csv`

**Usage:**
```bash
# Download all 44 samples in parallel (uses SLURM array)
sbatch --array=1-44 scripts/01_sra_download/01_sra_array.sh

# Or download a single sample for testing
sbatch --array=1 scripts/01_sra_download/01_sra_array.sh
```

**Output:**
- Downloads FASTQ files to `data/raw_fastq/PRJNA*/`
- Logs saved to `logs/sra_download/`

**Runtime:** ~30 minutes per sample (depending on file size and network speed)

---

### 02_sra_download.py

**Purpose:** Python script called by the SLURM array job to handle actual download via SRA Toolkit.

**Called by:** `01_sra_array.sh` (not run directly)

**Functions:**
- Reads sample information from samplesheet
- Downloads FASTQ using `fasterq-dump`
- Validates download completion
- Compresses FASTQ files with `gzip`

---

### 03_info_sra_checker.py

**Purpose:** Verify that all SRA samples downloaded successfully and have non-zero file sizes.

**Usage:**
```bash
python scripts/01_sra_download/03_info_sra_checker.py
```

**Checks:**
- All 44 samples present
- R1 and R2 files exist for each sample
- Files are not empty (size > 0)
- Reports any missing or corrupted files

**Output:** Summary report to stdout

---

## Expected Output Structure

```
data/raw_fastq/
├── PRJNA268379/
│   ├── SRR1719052_1.fastq.gz
│   ├── SRR1719052_2.fastq.gz
│   ├── SRR1719053_1.fastq.gz
│   └── ...
├── PRJNA158021/
│   ├── SRR458462_1.fastq.gz
│   ├── SRR458462_2.fastq.gz
│   └── ...
└── PRJNA187045/
    ├── SRR629651_1.fastq.gz
    ├── SRR629651_2.fastq.gz
    └── ...
```

**Total Size:** ~150-200 GB (compressed FASTQ files)

---

## Troubleshooting

### Download Failures
If a sample fails to download:
1. Check network connectivity
2. Verify SRA accession is correct
3. Re-run specific array index: `sbatch --array=5 scripts/01_sra_download/01_sra_array.sh`

### Disk Space
Ensure sufficient disk space before downloading:
```bash
df -h $PWD
```

### SRA Toolkit Issues
If `fasterq-dump` fails, try:
```bash
module load sratoolkit
vdb-config --interactive  # Configure cache location
```

---

## Next Step

After successful download and verification:
- **Proceed to:** `02_annotation_prep/` to prepare GTF annotation files

---

## Notes

- FASTQ files are not tracked in git (in `.gitignore`)
- Files can be re-downloaded from SRA if needed (publicly available)
- Compressed FASTQ files are kept for pipeline input
- Download logs are saved for troubleshooting
