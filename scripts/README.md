# Scripts Directory - Analysis Pipeline

**This directory contains ALL scripts needed to replicate the analysis from raw data to final results.**

Scripts are numbered in execution order. Run them sequentially.

---

## 00. Container Management

**Directory:** `00_container/`

**Scripts:**
- `pull_container.sh` - Pull updated container from GitHub

**When:** When Dockerfile is updated

**Output:** `albopictus-diapause-rnaseq.sif` in project root

---

## 01. rRNA Download (Optional)

**Directory:** `01_download_rrna/`

**Purpose:** Download rRNA sequences for filtering (if needed)

**Status:** Optional - nf-core pipeline handles rRNA filtering

---

## 02. Reference Preparation

**Directory:** `02_reference_preparation/`

**⚠️ CRITICAL: MANUAL DRYAD DOWNLOAD REQUIRED**

**Dryad blocks ALL automated downloads - you MUST download manually:**

1. Open browser: https://datadryad.org/stash/dataset/doi:10.5061/dryad.mgqnk98z4
2. Click "Download Dataset" button
3. Extract to `data/references/AalbF3/`

**Required files:**
- `AalbF3_genome.fa.gz` - genome assembly
- `AalbF3_annotation.gff3.gz` - gene annotations

**Scripts in this directory:** Reference indexing and preparation

**Status:** COMPLETED (data already downloaded, references indexed)

---

## 03. SRA Data Download

**Directory:** `03_sra_download/`

**Scripts:**
- `01_sra_array.sh` - Download all 44 RNA-seq samples from SRA

**When:** Before running RNA-seq pipeline

**Output:** `data/sra/` - FASTQ files

**Status:** COMPLETED

---

## 04. Annotation Preprocessing

**Directory:** `04_annotation_preprocessing/`

**Scripts:**
- GFF3 to GTF conversion and transcript fixing

**When:** Before running RNA-seq pipeline

**Documentation:** See scripts in directory

**Status:** COMPLETED

---

## 05. RNA-seq Pipeline

**Directory:** `05_rnaseq_pipeline/`

**Scripts:**
- `02_run_rnaseq.sh` - Main pipeline launcher (nf-core/rnaseq)
- `03_rerun_multiqc.sh` - Re-run MultiQC if it fails

**When:** After references and data are ready

**Output:**
- `output/star_salmon/` - BAM files, quantifications
- `output/multiqc/` - Quality control report

**Status:** COMPLETED (Job 20644937)

---

## 06. Count Matrix Preparation

**Directory:** `06_count_matrix/`

**Scripts:**
- `01_combine_counts.py` - Combine Salmon/featureCounts into matrices

**When:** After RNA-seq pipeline completes

**Output:** `results/count_matrices/` - Combined count matrices

**Status:** PENDING

---

## 07. QC Analysis & Data Splitting

**Directory:** `07_qc_analysis/`

**Scripts:**
- Split counts by developmental stage
- Add sample metadata
- PCA visualization

**When:** After count matrices are combined

**Output:** Stage-specific count files

**Status:** PENDING

---

## 08. Differential Expression Analysis

**Directory:** `08_differential_expression/`

**Scripts:**

### Adults (Blood meal split required):
- `01_adults_deseq2_unfed.R` - Unfed females DESeq2
- `01_run_adults_unfed.sh` - SLURM submission script
- `02_adults_deseq2_fed.R` - Blood-fed females DESeq2
- `02_run_adults_fed.sh` - SLURM submission script

### Embryos (by timepoint):
- `03_embryos_72h_deseq2.R` - 72h embryos (TO CREATE)
- `04_embryos_135h_deseq2.R` - 135h embryos (TO CREATE)

### Larvae (by timepoint):
- `05_larvae_11d_deseq2.R` - 11 day larvae (TO CREATE)
- `06_larvae_21d_deseq2.R` - 21 day larvae (TO CREATE)
- `07_larvae_40d_deseq2.R` - 40 day larvae (TO CREATE)

**When:** After count matrices are split by stage

**Output:** `results/03_differential_expression/` - DESeq2 results per stage

**Status:** Scripts exist for adults, need to create embryo/larvae scripts

**Documentation:** See `06_differential_expression/README.md`

---

## 08. GWAS Candidate Analysis

**Scripts:** (TO CREATE)
- Extract GWAS candidate expression
- Fisher's exact test for enrichment
- Create candidate heatmaps

**When:** After DESeq2 completes for all stages

**Output:** GWAS validation results for reviewers

**Status:** PENDING - critical for Reviewer 2

---

## 09. Visualization & Tables

**Scripts:** (TO CREATE)
- Volcano plots with candidates highlighted
- Cross-stage heatmaps
- Supplementary tables

**When:** After DESeq2 and GWAS analysis complete

**Output:** Publication-ready figures and tables

**Status:** PENDING

---

## 10. Meta-Analysis

**Scripts:** (TO CREATE)
- Fisher's method for combining p-values
- Stouffer's weighted method
- Cross-stage summary

**When:** After all stages analyzed

**Output:** Combined results across stages

**Status:** PENDING

---

## 99. Archive Directories

**99_method_validation_archive/** - Archived method comparison scripts
**99_visualization_future/** - Future visualization scripts

**Purpose:** Archive of completed validation work and future development

---

## Critical Notes

1. **⚠️ ALWAYS USE THE CONTAINER - NEVER HPC TOOLS DIRECTLY:**
   ```bash
   module load singularity
   # For Python scripts:
   singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c "cd /proj && python script.py"
   # For R scripts:
   singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c "cd /proj && Rscript script.R"
   ```
   **NEVER use system Python, R, or any HPC-installed tools!**

2. **Track progress:**
   ```bash
   ./00_track.py status  # Before starting
   ./00_track.py complete <task_id>  # After finishing
   ```

3. **SLURM job logs:**
   - All logs go to `logs/` directory
   - Use `.o.txt` and `.e.txt` extensions (NOT .out/.err)

4. **Container-only policy:**
   - NEVER install software outside container
   - Update Dockerfile and rebuild if packages missing

---

## For Replication

To replicate this analysis from scratch:

```bash
# 1. Check where to start
./00_track.py status

# 2. Follow numbered script order
# 3. Mark each task complete
./00_track.py complete <task_id>

# 4. See detailed documentation in each directory's README.md
```

---

## Future: Snakemake Automation

**Note:** Once all individual scripts are tested and working, a Snakemake workflow will be created to automate the entire pipeline.

**Snakefile will:**
- Orchestrate all numbered scripts in proper order
- Handle dependencies between steps
- Enable parallel execution where possible
- Provide single-command reproduction: `snakemake --use-singularity`

**Location:** `Snakefile` (to be created in project root)

This ensures the workflow is:
1. **Documented** - Individual numbered scripts show each step
2. **Automated** - Snakemake runs everything in correct order
3. **Reproducible** - Single command replicates entire analysis

---

**Last Updated:** Oct 24, 2025
