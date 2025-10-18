# Albopictus Diapause RNA-seq Analysis Project

**Last Updated:** October 18, 2025 11:30
**Status:** ✅ ALL 44 SAMPLES COMPLETED! Ready for QC analysis!

---

## PROJECT OVERVIEW

**Goal:** Validate 34 GWAS candidate genes for diapause regulation in *Aedes albopictus*

**Datasets:** 3 independent RNA-seq experiments
- PRJNA268379: Adult females (16 samples) - Huang et al. 2015
- PRJNA158021: Embryos (12 samples) - Poelchau et al. 2013a
- PRJNA187045: Pharate larvae (16 samples) - Poelchau et al. 2013b
- **Total: 44 samples**

**Challenge:** Replicating collaborator's results from published manuscripts while addressing reviewer concerns

---

## CURRENT STATUS - OCTOBER 18, 2025 11:30

### ✅ PROJECT REORGANIZATION COMPLETE!
- [x] Cleaned up logs directory (removed old July runs)
- [x] Reorganized scripts into numbered workflow directories (00-07)
- [x] Created README.md files for each script directory
- [x] Updated all script cross-dependencies
- [x] Following bioinformatics best practices for reproducibility

### ✅ ALL 44 SAMPLES COMPLETED! ✓
- [x] Downloaded correct AalbF3 genome from Dryad (GCA_018104305.1)
- [x] Fixed chromosome naming (stripped "chr" prefix to match SNP chip format)
- [x] Downloaded 12 Ae. albopictus rRNA sequences for SortMeRNA
- [x] Created GTF with proper gene_id attributes for RSEM compatibility
- [x] Filtered 45 problematic transcripts (multi-chromosome/inconsistent strands)
- [x] **FIXED GTF newline bug** - regenerated GTF with proper formatting
- [x] Updated pipeline script to use new references and rRNA removal
- [x] Cleaned up 2TB of old cache files
- [x] **DOCUMENTED ALL ANNOTATION PROBLEMS FOR MANUSCRIPT METHODS**
- [x] **TEST JOB SUBMITTED AND RUNNING** (Job 20467416)

### Latest Fix (October 17, 2025 16:00)
**Problem**: GTF file had literal "\n" characters instead of actual newlines
**Root Cause**: Bug in `add_gene_biotype.py` line 74 - used `\\n` instead of `\n`
**Solution**:
1. Fixed the bug in add_gene_biotype.py
2. Created new script `convert_gff_to_gtf.py` for proper GFF→GTF conversion
3. Re-ran complete 3-step GTF generation workflow
4. GTF now validates correctly (9 columns, proper newlines)

**Test Result**: GTF_FILTER step now PASSES ✓ (previously failed)

### What Changed from Previous Runs
1. **Genome**: `albo.fasta.gz` (old) → `AalbF3_genome.fa.gz` (correct, chr prefix stripped)
2. **Annotation**: Multiple iterations → `AalbF3_annotation.gtf.gz` (final, properly formatted)
3. **Added rRNA removal**: SortMeRNA with 12 Ae. albopictus sequences
4. **Fixed all GTF formatting issues**: See "GTF Generation Workflow" section below

### Current Pipeline Status - COMPLETE! ✅
```bash
# Test job: COMPLETED ✓
# Job ID: 20467416 (sample 1: PRJNA158021_SRR458462)
# Result: Gene count = 22,176 genes detected! (was 0 before GTF fix)
# Status: GTF validation PASSED, all QC steps completed
# Runtime: 1h 55m

# Monitor job: COMPLETED ✓
# Job ID: 20468449
# Result: Validated output, auto-launched remaining samples
# Email sent: "✅ RNA-seq Test PASSED - Full Array Launched"

# Full array: COMPLETED ✅
# Job ID: 20468495 (samples 2-44)
# Start time: October 17, 2025 18:00
# End time: October 18, 2025 05:04 (last sample completed)
# Total runtime: ~11 hours for all 44 samples
# All samples completed successfully!

# Verify completion:
ls -d output/PRJNA*/SRR*/star_salmon/ | wc -l
# Expected: 44 (all samples)
# Actual: 44 ✓

# Check logs for any errors:
grep -i "error" logs/successful_runs/20241017_full_pipeline/rnaseq_*.o.txt
# No critical errors found ✓
```

### Automated Pipeline Launch (NEW!)
Created monitoring script that will:
1. ✅ Poll test job every 5 minutes until completion
2. ✅ Validate pipeline output (gene counts, output files)
3. ✅ If successful: Automatically launch jobs 2-44 (remaining 43 samples)
4. ✅ If failed: Email lcosme@ucr.edu with error details
5. ✅ Send status email either way

**Script location**: `scripts/02_run_rnaseq/monitor_and_launch.sh`

### Automation Results ✅
**Test succeeded! Here's what happened:**
- ✅ Monitor validated gene counts: 22,176 genes (expected ~15-20k)
- ✅ Monitor validated all output files exist
- ✅ Monitor automatically submitted: Job 20468495 (samples 2-44)
- ✅ Email sent: "✅ RNA-seq Test PASSED - Full Array Launched"
- ✅ All 44 samples now processing!

**Timeline:**
- 16:02 - Test job submitted (sample 1)
- 17:32 - Monitor job submitted
- 17:57 - Test completed successfully (runtime: 1h 55m)
- 17:58 - Monitor validated and auto-launched remaining 43 samples
- 18:00 - 10 samples running in parallel, 33 pending

---

## REFERENCE FILES (ZENODO UPLOAD LATER)

All large reference files are in `.gitignore`. Upload to Zenodo for publication.

### Current Reference Files (VERIFIED AND READY)
Location: `data/references/AalbF3_genome/`

| File | Size | Description | Status |
|------|------|-------------|--------|
| `AalbF3_genome.fa.gz` | 423 MB | Genome FASTA, chr prefix stripped | ✓ Ready |
| `AalbF3_annotation.gtf.gz` | 4.9 MB | Fixed GTF, 347,380 features, 32,889 transcripts | ✓ Ready |
| `AalbF3_annotation.gff.gz` | 6.3 MB | Original GFF3 (backup, not used) | ✓ Backup |
| `combined_rRNA_sequences.fa` | 37 KB | 12 Ae. albopictus rRNA sequences | ✓ Ready |
| `rrna_db_manifest.txt` | 57 bytes | SortMeRNA manifest (tracked in git) | ✓ Ready |

### Reference Source
- **Assembly**: Boyle et al. 2021, Insects 12:1-19 (PMID: 33669192)
  - Linkage-based genome assembly (AalbF3)
  - Based on 111,328 informative SNPs from RNAseq
  - 1.45 Gb assembly with 84.3% complete BUSCO genes
- **Annotation**: Palatini et al. 2020, BMC Biology
- **Accession**: GCA_018104305.1
- **Repository**: Dryad (doi:10.5061/dryad.mgqnk98z4)
- **Downloaded**: October 17, 2025

### Why This Reference?
- Same genome used by collaborators (srmarzec/albopictus_remapping)
- High-quality linkage-based chromosome-level assembly
- Complete GFF3 annotation with rRNA, tRNA, and all gene features
- Identifies chromosomal regions affecting photoperiodic diapause
- Public and citable (required for publication)

### Reference Preprocessing Steps
1. **Download from Dryad** (NCBI lacks complete annotation)
   - **Critical:** NCBI's GCA_018104305.1 does NOT include the complete GFF3 annotation
   - **Solution:** Use Dryad repository which has complete annotation with rRNA, tRNA, and all gene features
   - **URL:** https://datadryad.org/stash/dataset/doi:10.5061/dryad.mgqnk98z4
   - **Files needed:**
     - `AALBF3_assembly.fa.gz` (~423 MB)
     - `AALBF3_annotation.gff3.gz` (~6.3 MB)
   - **Download method:** Manual download from browser + `scp` to HPC (automated wget/curl may fail due to Dryad's download mechanism)
   - **Detailed instructions:** See `scripts/00_reference_preparation/README.md`

2. **Strip chr prefix** (match SNP chip naming convention)
   - Original: `>chr1.1`, `>chr2.206`
   - Fixed: `>1.1`, `>2.206`
   - Script: `scripts/01_download_data/strip_chr_prefix.sh`

3. **Convert GFF3 to GTF and fix for RSEM**
   - See "GFF3/GTF Annotation Problems" section below

4. **Download rRNA sequences**
   - 12 Ae. albopictus sequences from GenBank (Koh et al. 2023, eLife)
   - Script: `scripts/01_download_data/download_rrna_sequences.sh`
   - Headers cleaned for SortMeRNA compatibility

---

## COMPLETE TROUBLESHOOTING HISTORY

**Purpose:** Document ALL attempts, failures, and solutions so we never repeat the same mistakes.

### Attempt 1: Using NCBI Genome (FAILED)
**What we tried:** Download GCA_018104305.1 directly from NCBI
**Result:** FAILED - NCBI does NOT provide the complete GFF3 annotation file
**Lesson:** Must use Dryad repository for complete annotation

### Attempt 2: Automated Download from Dryad (FAILED)
**What we tried:** Use wget/curl to download files from Dryad
**Result:** FAILED - Dryad's download mechanism requires browser interaction
**Solution:** Manual download from browser + `scp` to HPC

### Attempt 3: Using GFF3 Directly with nf-core (FAILED)
**What we tried:** Use `AALBF3_annotation.gff3.gz` directly
**Result:** FAILED - nf-core requires `.gff.gz` extension, not `.gff3.gz`
**Solution:** Rename to `.gff.gz` OR convert to GTF

### Attempt 4: nf-core Internal GFF→GTF Conversion (FAILED)
**What we tried:** Let nf-core convert GFF to GTF using gffread
**Result:** FAILED - RSEM error "Cannot find gene_id!"
**Problem:** nf-core's gffread created 18,135 GTF lines without gene_id attributes
**Lesson:** Cannot rely on nf-core's internal conversion - must create proper GTF ourselves

### Attempt 5: Quick gene_id Fix (FAILED)
**What we tried:** Script `fix_gtf_gene_id.py` to add gene_id to missing lines
**Result:** FAILED - Created double semicolons `;;` (invalid GTF syntax)
**Lesson:** Need to strip trailing semicolons before adding attributes

### Attempt 6: Fixed gene_id Script (FAILED)
**What we tried:** Updated script to avoid double semicolons
**Result:** FAILED - RSEM error "transcript has exons on multiple chromosomes!"
**Problem:** NCBI annotation has 43 transcripts with exons spanning different scaffolds
**Lesson:** Must filter out problematic transcripts before using with RSEM

### Attempt 7: Custom GTF with Transcript Filtering (FAILED)
**What we tried:** Filter multi-chromosome transcripts, keep everything else
**Result:** FAILED - featureCounts error "failed to find gene_biotype"
**Problem:** Only 9.1% of lines had gene_biotype attribute
**Lesson:** featureCounts needs gene_biotype on ALL lines for categorization

### Attempt 8: Add gene_biotype from gbkey (FAILED)
**What we tried:** Script `add_gene_biotype.py` to map gbkey → gene_biotype
**Result:** FAILED - nf-core GTF_FILTER error "Expected 9 tab-separated columns"
**Problem:** Script had bug using `\\n` instead of `\n` - created literal "\n" text
**Root cause:** Line 74 of script wrote escaped newlines instead of actual newlines
**Lesson:** Always validate output file format, not just content

### Attempt 9: Complete 3-Step GTF Workflow (SUCCESS ✅)
**What we did:**
1. Created `convert_gff_to_gtf.py` - proper GFF3→GTF with all attributes
2. Created `gff_to_rsem_gtf.py` - filter 45 problematic transcripts
3. Fixed `add_gene_biotype.py` bug - changed `\\n` to `\n`
4. Re-ran all 3 steps to generate clean GTF

**Result:** SUCCESS!
- Test job: 22,176 genes detected (was 0 before)
- GTF validation: PASSED
- All 44 samples: COMPLETED

**Final GTF stats:**
- 344,588 features
- 53,619 valid transcripts (45 removed)
- All have gene_id ✓
- All have gene_biotype ✓
- Proper GTF9 format ✓
- No multi-chromosome transcripts ✓

---

## GFF3/GTF ANNOTATION PROBLEMS (FOR MANUSCRIPT METHODS)

**IMPORTANT:** Document these 6 problems in the manuscript methods section.

### Problem 1: GFF file extension validation
- **Error**: nf-core rejected `.gff3.gz` file extension
- **Cause**: Pipeline validates extensions and only accepts `.gff.gz` format
- **Solution**: Renamed `AalbF3_annotation.gff3.gz` → `AalbF3_annotation.gff.gz`
- **Impact**: Minor formatting issue, resolved immediately

### Problem 2: Missing gene_id attributes in GTF
- **Error**: RSEM failed with "Cannot find gene_id!" error
- **Cause**: nf-core's internal GFF→GTF conversion (gffread) created 18,135 GTF lines without gene_id attributes
- **Solution**: Created `scripts/utils/fix_gtf_gene_id.py` to extract gene_id from gene_name and add to all features
- **Impact**: Fixed 18,135 lines (5.2% of total features)
- **Key code**:
  ```python
  gene_match = re.search(r'gene[_ ]"([^"]+)"', line)
  gene_id = gene_match.group(1).replace('gene-', '')
  line = line.rstrip().rstrip(';') + f'; gene_id "{gene_id}";\n'
  ```

### Problem 3: Double semicolon syntax errors
- **Error**: RSEM failed with "Cannot locate the identifier from attribute ;!"
- **Cause**: Initial fix script added `; gene_id` to lines already ending with `;`, creating invalid `;;` syntax
- **Solution**: Modified script to strip trailing semicolons before adding gene_id attribute
- **Impact**: All 348,074 lines now have proper GTF9 format

### Problem 4: Transcripts on multiple chromosomes
- **Error**: RSEM validation failed with "transcript rna-XM_029864299.1 has exons on multiple chromosomes!"
- **Cause**: NCBI AalbF3 annotation contains 43 transcripts with exons spanning different scaffolds (annotation errors)
- **Solution**: Modified `scripts/utils/gff_to_rsem_gtf.py` to:
  - Track both strand consistency AND chromosome consistency per transcript
  - Filter out 45 problematic transcripts (43 multi-chromosome, 26 inconsistent strands, some with both)
  - Reduced from 348,074 to 347,380 features (0.2% loss)
- **Impact**: Minimal loss of features while ensuring RSEM compatibility
- **Examples of problematic transcripts**:
  - `rna-XM_029864299.1`: exons on scaffolds 2.198 and 2.205
  - `rna-XM_029866277.1`: inconsistent strands AND multiple chromosomes
  - `gene-LOC109431822`: spans scaffolds 2.175 and 3.83

### Problem 5: Missing gene_biotype attributes for featureCounts
- **Error**: featureCounts failed with "failed to find the gene identifier attribute 'gene_biotype'"
- **Cause**: Only 31,700 lines (9.1%) had `gene_biotype` attribute, but featureCounts needs it on all lines for gene categorization
- **Solution**: Created `scripts/utils/add_gene_biotype.py` to add `gene_biotype` to all lines based on `gbkey` mapping:
  - `gbkey "mRNA"` → `gene_biotype "protein_coding"`
  - `gbkey "rRNA"` → `gene_biotype "rRNA"`
  - `gbkey "tRNA"` → `gene_biotype "tRNA"`
  - `gbkey "ncRNA"` → `gene_biotype "ncRNA"`
  - `gbkey "misc_RNA"` → `gene_biotype "misc_RNA"`
- **Impact**: Added gene_biotype to all features
- **Why needed**: featureCounts uses gene_biotype to categorize reads (protein-coding vs rRNA vs other), essential for QC and comparison with published results

### Problem 6: Literal \n characters in GTF (Bug in add_gene_biotype.py)
- **Error**: nf-core GTF_FILTER failed with "Invalid GTF file: Expected 9 tab-separated columns"
- **Cause**: `add_gene_biotype.py` script had `\\n` (escaped backslash-n) instead of `\n` (newline), writing literal "\n" strings into the file and concatenating multiple lines together
- **Solution**:
  1. Fixed bug in `scripts/utils/add_gene_biotype.py` line 74: changed `\\n` to `\n`
  2. Created `scripts/utils/convert_gff_to_gtf.py` to properly convert GFF3 to GTF format
  3. Re-ran complete GTF generation workflow (see GTF Generation Workflow below)
- **Impact**: GTF now has proper newlines and exactly 9 tab-separated columns per line

### GTF Generation Workflow (October 17, 2025)

The complete workflow to generate a clean, RSEM-compatible GTF from the NCBI GFF3:

**Step 1: Convert GFF3 to GTF format**
```bash
python scripts/utils/convert_gff_to_gtf.py \
    data/references/AalbF3_genome/AalbF3_annotation.gff.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz
```
- Input: 376,867 lines (GFF3)
- Output: 356,246 features (GTF)
- Converts GFF3 attributes (key=value) to GTF format (key "value")
- Adds transcript_id and gene_id to all features
- Skips gene features (only transcript-level and below)

**Step 2: Filter problematic transcripts**
```bash
python scripts/utils/gff_to_rsem_gtf.py \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz
```
- Input: 356,246 features from 53,664 transcripts
- Output: 344,588 features (45 transcripts removed)
- Filters transcripts with inconsistent strands (26 transcripts)
- Filters transcripts on multiple chromosomes (43 transcripts)
- Some transcripts have both issues

**Step 3: Add gene_biotype attributes**
```bash
python scripts/utils/add_gene_biotype.py \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation.gtf.gz
```
- Input: 344,588 features
- Output: 344,588 features (all with gene_biotype added)
- Maps gbkey values to standard gene_biotype terms
- Required for featureCounts gene categorization

**Final GTF: 344,588 features from 53,619 valid transcripts (344,588 features / ~6.4 features per transcript)**

### Summary for Methods Section
The NCBI AalbF3 annotation (GCF_018104305.1) required preprocessing for compatibility with nf-core/rnaseq. We developed a three-step workflow:

1. **GFF3 to GTF conversion** - Converted NCBI GFF3 format to GTF, properly formatting attributes and adding transcript_id/gene_id to all features (356,246 features retained)

2. **Transcript filtering** - Removed 45 transcripts (0.08% of 53,664 total) with annotation errors:
   - 26 transcripts with inconsistent strand annotations
   - 43 transcripts with exons on multiple chromosomes (scaffolds)
   - Some transcripts had both issues

3. **Gene biotype annotation** - Added gene_biotype attributes to all features based on NCBI gbkey mapping for featureCounts compatibility

**Final annotation: 344,588 features from 53,619 valid transcripts, all with proper GTF formatting, transcript_id, gene_id, and gene_biotype attributes**

---

## NF-CORE PIPELINE CONFIGURATION

### Critical Configuration Notes

**ALWAYS SET SINGULARITY_TMPDIR** - Required for Singularity on HPC:
```bash
export SINGULARITY_TMPDIR="${PROJECT_BASE}/output/singularity_temp"
mkdir -p $SINGULARITY_TMPDIR
```
Configured in `scripts/02_run_rnaseq/01_run_rnaseq_array.sh` at line 26.

### Pipeline Script Location
`scripts/02_run_rnaseq/01_run_rnaseq_array.sh`

### Key Parameters
```bash
FASTA="data/references/AalbF3_genome/AalbF3_genome.fa.gz"
GTF="data/references/AalbF3_genome/AalbF3_annotation.gtf.gz"
RRNA_MANIFEST="data/references/AalbF3_genome/rrna_db_manifest.txt"

nextflow run nf-core/rnaseq \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --fasta "$FASTA" \
    --gtf "$GTF" \
    --remove_ribo_rna \
    --ribo_database_manifest "$RRNA_MANIFEST" \
    --save_reference
```

### Directory Management

**Cache directories (DELETE after each run to save space):**
- `output/nf-work/` - Nextflow work cache (can be 2TB!)
- `output/nf-temp/` - Temporary files
- `output/singularity_temp/` - Singularity runtime temp
- `output/temp_samplesheets/` - Temporary samplesheets

**Keep these directories:**
- `output/singularity/` - Container images (3.8GB, reusable)
- `output/previous_run_backup/` - Old results with wrong references

### Running the Pipeline
```bash
# Test with one sample first
sbatch --array=1 scripts/02_run_rnaseq/01_run_rnaseq_array.sh

# If test passes, run all 44 samples
sbatch --array=1-44 scripts/02_run_rnaseq/01_run_rnaseq_array.sh
```

### Check Logs
```bash
# SLURM output
tail -f logs/rnaseq_JOBID_ARRAYID.o.txt

# SLURM errors
tail -f logs/rnaseq_JOBID_ARRAYID.ERROR.txt
```

---

## SAMPLE INFORMATION

### Experimental Designs

**PRJNA268379 - Adult Females (Huang et al. 2015)**
- 2×2 factorial: Photoperiod (LD/SD) × Blood meal (NBF/BF)
- 4 replicates per condition = 16 samples
- Platform: HiSeq2000
- Life stage: Adult females

**PRJNA158021 - Embryos (Poelchau et al. 2013a)**
- Time course: 72-78h and 135-141h post-oviposition
- Conditions: Diapause vs Non-diapause
- 3 replicates per timepoint/condition = 12 samples
- Platform: GAIIx

**PRJNA187045 - Pharate Larvae (Poelchau et al. 2013b)**
- Time course: 11, 21, 40 days post-oviposition
- Conditions: Diapause vs Quiescence
- Variable replicates = 16 samples (not 17)
- Platform: HiSeq2000

### Batch Effects Concern
- Different sequencing platforms (GAIIx vs HiSeq2000)
- Must account for platform in statistical models
- Use ComBat-seq or include platform as covariate in DESeq2

---

## ANALYSIS WORKFLOW

### Phase 1: QC Assessment (IN PROGRESS)
1. Extract QC metrics from nf-core MultiQC outputs
2. Create manuscript-ready QC summary tables
3. Generate standard QC figures (depth, mapping rates, PCA)
4. Assess data quality and identify problematic samples

**Standard QC metrics:**
- Total reads
- Mapping rate (>80% expected)
- Uniquely mapped % vs multi-mapped %
- rRNA contamination %
- Insert size distribution
- Gene detection (# genes with >10 counts)

**Standard QC figures:**
- Sequencing depth bar chart
- Mapping rate stacked bar chart
- PCA plot (per-project and overall)
- Correlation heatmap
- Gene body coverage plot

### Phase 2: Count Matrix Preparation
1. Gather `salmon.merged.gene_counts.tsv` from all samples
2. Create combined count matrices by project
3. Import 34 GWAS candidate genes list
4. Verify gene IDs match annotation

### Phase 3: Differential Expression
**Analyze per-project first** (following collaborator's approach):

1. **PRJNA268379 (Adults)**: Factorial design
   ```r
   design = ~ photoperiod + blood_meal + photoperiod:blood_meal
   ```

2. **PRJNA158021 (Embryos)**: Time course
   ```r
   design = ~ timepoint + condition + timepoint:condition
   ```

3. **PRJNA187045 (Larvae)**: Time course
   ```r
   design = ~ timepoint + condition + timepoint:condition
   ```

**Statistical considerations:**
- Apply standard thresholds (FDR < 0.05, |log2FC| > 1)
- Document model choices for reviewers
- Compare with collaborator's published DEG lists

### Phase 4: Meta-analysis
1. Integrate results across life stages
2. Account for platform batch effects
3. Focus on 34 GWAS candidate genes
4. Prioritize genes consistent across datasets

### Phase 5: Publication
1. Generate volcano plots, MA plots, heatmaps
2. Gene set enrichment analysis
3. Prepare supplementary data files
4. Write methods section (include annotation preprocessing)

---

## CONTAINER USAGE

### Singularity Container
- **Image**: `albopictus-diapause-rnaseq.sif` (1.3GB)
- **Location**: `/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/`
- **Mount point**: Project directory → `/proj`

### Loading Container
```bash
module load singularity/3.9.3
cd /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq
```

### Interactive Shell
```bash
singularity shell --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif
cd /proj
# Now run Python/R scripts
```

### Execute Command
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
    bash -c "cd /proj && python scripts/03_qc_analysis/extract_qc_metrics.py"
```

### Container Contents
- Python 3.10
- R 4.3
- pandas, matplotlib, seaborn
- DESeq2, tximport, edgeR
- All dependencies for analysis scripts

---

## KEY FILE LOCATIONS

### Quantification Results
```
output/PRJNA*/SRR*/star_salmon/
├── quant.sf                              # Transcript-level
├── salmon.merged.gene_counts.tsv         # Gene-level counts
├── salmon.merged.gene_tpm.tsv            # Gene-level TPM
└── featurecounts/*.tsv                   # Alternative counts
```

### QC Reports
```
output/PRJNA*/SRR*/
├── multiqc/star_salmon/multiqc_report.html
├── multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt
├── fastqc/
└── star_salmon/qualimap/
```

### Logs
```
logs/
├── rnaseq_JOBID_*.o.txt                  # SLURM output
├── rnaseq_JOBID_*.ERROR.txt              # SLURM errors
└── nextflow_*.log                         # Nextflow logs
```

---

## COMPARISON WITH COLLABORATORS

### Their Approach (srmarzec/albopictus_remapping)
- Manual pipeline: Trimmomatic → STAR → HTSeq
- Project-specific DESeq2 scripts
- HEADCROP:15 trimming strategy
- Same AalbF3 reference genome
- Used older annotation (unknown version)

### Our Approach
- nf-core/rnaseq v3.19.0 (automated, containerized)
- Salmon quantification (faster, includes TPM)
- Comprehensive QC (MultiQC + Qualimap + RSeQC)
- Adaptive trimming
- rRNA removal with SortMeRNA
- Latest AalbF3 annotation (preprocessed for compatibility)

### Key Validation Steps

**1. Gene Count Comparison**
After pipeline completion, compare our gene-level counts with collaborators' HTSeq counts:
```bash
# Extract gene counts from our Salmon results
# Location: output/PRJNA*/SRR*/star_salmon/salmon.merged.gene_counts.tsv

# Compare with collaborator counts
# Location: data/collaborator_repos/albopictus_remapping/counts/
```

**Expected outcome**: High correlation (>0.95) for gene-level counts despite different methods

**2. Transcript-Level Resolution**
Our Salmon output provides transcript-level quantification that collaborators didn't have:
```bash
# Transcript counts: output/PRJNA*/SRR*/star_salmon/quant.sf
```

This allows us to:
- Identify isoform switching events during diapause
- Detect alternative splicing patterns
- Provide more granular gene expression analysis

**3. QC Metric Validation**
Compare key QC metrics to ensure data quality:
- Total reads processed (should match SRA)
- Mapping rates (expect >80%)
- rRNA contamination (compare with/without rRNA removal)
- Gene detection numbers (expect ~15,000-20,000 genes)

**4. Differential Expression Validation**
For the 34 GWAS candidate genes:
- Compare log2FC direction (up/down regulation)
- Compare significance (p-values)
- Identify any discrepancies and investigate causes

### Advantages of Our Approach
1. **Reproducibility**: Containerized with locked versions
2. **Speed**: Salmon is much faster than STAR+HTSeq
3. **QC**: More comprehensive metrics
4. **Automation**: Less manual intervention, fewer errors
5. **Best practices**: Following current RNA-seq standards
6. **Transcript resolution**: Isoform-level quantification
7. **Documented preprocessing**: All annotation fixes tracked

### Addressing Reviewers
When reviewers ask about differences:
1. Same samples and experimental designs
2. Modern pipeline improves accuracy and reproducibility
3. Salmon validated to produce comparable results to STAR+HTSeq (see [Patro et al. 2017](https://doi.org/10.1038/nmeth.4197))
4. Better QC filtering and improved algorithms
5. Any differences in DEGs likely due to improved methods
6. We will provide correlation analysis with original results
7. All preprocessing steps are documented and reproducible

### Planned Comparison Analysis
```r
# Script to compare with collaborator results
# Location: scripts/04_validation/compare_with_collaborators.R

# Steps:
# 1. Load our Salmon gene counts
# 2. Load collaborator HTSeq counts
# 3. Match sample IDs across datasets
# 4. Calculate Pearson/Spearman correlation
# 5. Create scatter plots (our counts vs. theirs)
# 6. Identify highly discrepant genes
# 7. Compare DE results for 34 candidate genes
# 8. Generate validation report
```

---

## REPRODUCIBLE WORKFLOW - SCRIPT ORGANIZATION

**Last Reorganized:** October 18, 2025
**Approach:** Numbered workflow directories following bioinformatics best practices
**Goal:** Prepare for Snakemake automation and manuscript publication

### Complete Workflow Structure

```
scripts/
├── 00_reference_preparation/      # ONE-TIME: Prepare genome and rRNA references
│   ├── README.md                  # Detailed instructions and Dryad links
│   ├── 01_strip_chr_prefix.sh     # Remove chr prefix from genome
│   ├── 02_download_rrna_sequences.sh  # Download 12 Ae. albopictus rRNA
│   └── 03_clean_rrna_headers.sh   # Simplify headers for SortMeRNA
│
├── 01_sra_download/               # Download 44 FASTQ files from SRA
│   ├── README.md                  # Sample metadata and download instructions
│   ├── 01_sra_array.sh            # SLURM array job for parallel download
│   ├── 02_sra_download.py         # Python script for fasterq-dump
│   └── 03_info_sra_checker.py     # Verify downloads completed
│
├── 02_annotation_prep/            # ONE-TIME: Fix GTF annotation for nf-core
│   ├── README.md                  # Documents all 6 GTF problems and solutions
│   ├── 01_convert_gff_to_gtf.py   # Step 1: GFF3 → GTF conversion
│   ├── 02_filter_problematic_transcripts.py  # Step 2: Remove 45 bad transcripts
│   └── 03_add_gene_biotype.py     # Step 3: Add gene_biotype for featureCounts
│
├── 03_rnaseq_pipeline/            # Run nf-core/rnaseq on all 44 samples
│   ├── README.md                  # Pipeline parameters and monitoring guide
│   ├── 01_create_samplesheet.py   # Generate nf-core samplesheet
│   ├── 02_run_rnaseq_array.sh     # Main SLURM array job (1-44)
│   ├── 03_monitor_and_launch.sh   # Auto-validate and launch remaining jobs
│   └── hpc_batch.conf             # HPC-specific Nextflow configuration
│
├── 04_qc_analysis/                # Extract and visualize QC metrics
│   ├── README.md                  # QC metrics documentation
│   ├── 01_extract_qc_metrics.py   # Parse MultiQC outputs
│   ├── 02_generate_qc_figures.py  # Create publication-ready QC figures
│   └── 03_run_qc_analysis.sh      # Wrapper script for Singularity container
│
├── 05_count_matrix/               # Combine count matrices
│   ├── README.md                  # Count matrix organization
│   ├── 01_combine_counts.py       # Merge Salmon gene counts
│   └── 02_create_summary_table.py # Generate summary statistics
│
├── 06_diff_expression/            # DESeq2 analysis per BioProject
│   ├── README.md                  # Statistical models and comparisons
│   ├── 01_de_analysis_PRJNA268379_adults.py      # Adults (factorial design)
│   ├── 02_de_analysis_PRJNA158021_embryos.py     # Embryos (time course)
│   ├── 03_de_analysis_PRJNA158021_comparison.py  # Embryo comparisons
│   ├── 04_embryo_detailed_analysis.py            # Detailed embryo analysis
│   └── 05_outlier_detection_PRJNA158021.py       # QC and outlier detection
│
├── 07_visualization/              # Publication figures
│   ├── README.md                  # Figure descriptions
│   └── 01_create_publication_figures.py  # Generate all manuscript figures
│
├── utils/                         # Helper scripts (not numbered)
│   ├── README.md                  # Utility documentation
│   ├── create_samplesheet.py      # ✅ Active: samplesheet generation
│   ├── update_sample_mapping.py   # ✅ Active: sample metadata updates
│   ├── organize_outputs.sh        # ⚠️ DEPRECATED: old laptop workflow
│   ├── fix_gtf_gene_id.py         # ⚠️ DEPRECATED: superseded by 02_annotation_prep
│   ├── convert_gff_to_gtf.py      # Reference copy (canonical in 02_annotation_prep)
│   ├── gff_to_rsem_gtf.py         # Reference copy (canonical in 02_annotation_prep)
│   └── add_gene_biotype.py        # Reference copy (canonical in 02_annotation_prep)
│
└── OLD_STRUCTURE/                 # Archived old organization (pre-Oct 18, 2025)
    ├── 01_download_data/
    ├── 02_run_rnaseq/
    ├── 03_qc_analysis/
    ├── 04_differential_expression/
    └── 05_visualization/
```

### Workflow Execution Order

**Phase 1: One-Time Setup - Reference Preparation**

**Step 00: Reference Genome Preparation**
1. Manually download from Dryad (NCBI lacks complete GFF3):
   - `AALBF3_assembly.fa.gz`
   - `AALBF3_annotation.gff3.gz`
   - Transfer to HPC via `scp`
2. Run `00_reference_preparation/01_strip_chr_prefix.sh` - Remove chr prefix
3. Run `00_reference_preparation/02_download_rrna_sequences.sh` - Get 12 rRNA sequences
4. Run `00_reference_preparation/03_clean_rrna_headers.sh` - Clean headers for SortMeRNA

**Step 02: GTF Annotation Preparation (CRITICAL - Multiple attempts required)**
**Problem:** Original GFF3 from Dryad had 6 critical issues that caused nf-core failures:
1. Wrong file extension (.gff3.gz → .gff.gz)
2. Missing gene_id attributes (RSEM error)
3. Double semicolon syntax errors
4. Transcripts on multiple chromosomes (45 transcripts)
5. Missing gene_biotype attributes (featureCounts error)
6. Literal `\n` characters instead of newlines

**Solution - 3-Step GTF Generation Workflow (October 17, 2025):**

```bash
# Step 1: Convert GFF3 to GTF format
python scripts/02_annotation_prep/01_convert_gff_to_gtf.py \
    data/references/AalbF3_genome/AalbF3_annotation.gff.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz
```
**Output:**
- Reading GFF3 file: data/references/AalbF3_genome/AalbF3_annotation.gff.gz
- Lines processed: 376,867
- Features written: 356,246

```bash
# Step 2: Filter 45 problematic transcripts
python scripts/02_annotation_prep/02_filter_problematic_transcripts.py \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz
```
**Output:**
- Total features processed: 356,246
- Total transcripts found: 53,664
- Transcripts with inconsistent strands: 26
- Transcripts on multiple chromosomes: 43
- Transcripts removed: 45 (some had both issues)
- Features written: 344,588

```bash
# Step 3: Add gene_biotype attributes
python scripts/02_annotation_prep/03_add_gene_biotype.py \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation.gtf.gz
```
**Output:**
- Lines with existing biotype: 31,700
- Lines added biotype: 312,888
- Total features: 344,588

**Final GTF Created:** `data/references/AalbF3_genome/AalbF3_annotation.gtf.gz`
- **Date created:** October 17, 2025 16:01
- **Size:** 4.5 MB (compressed)
- **Features:** 344,588
- **Valid transcripts:** 53,619 (45 removed from original 53,664)
- **Validation:** All lines have gene_id ✓, All lines have gene_biotype ✓, Proper GTF9 format ✓

**Test Result:** Used in nf-core pipeline test job → 22,176 genes detected (SUCCESS!)

**Scripts Used:**
- `scripts/02_annotation_prep/01_convert_gff_to_gtf.py` (ACTIVE VERSION)
- `scripts/02_annotation_prep/02_filter_problematic_transcripts.py` (ACTIVE VERSION)
- `scripts/02_annotation_prep/03_add_gene_biotype.py` (ACTIVE VERSION - bug fixed)

**Note:** Duplicate copies exist in `scripts/utils/` as reference/backup from development. The canonical working versions are in `scripts/02_annotation_prep/`.

**DO NOT re-run these scripts** - the GTF is complete and working. File is dated October 17, 2025 16:01.

**Step 01: SRA Download**
1. Download all 44 FASTQ files using SLURM array job:
   ```bash
   sbatch --array=1-44 scripts/01_sra_download/01_sra_array.sh
   ```
2. Verify downloads:
   ```bash
   python scripts/01_sra_download/03_info_sra_checker.py
   ```

**Phase 2: Pipeline Execution**

**Step 03: nf-core/rnaseq Pipeline**
1. Create samplesheet for nf-core (uses `scripts/utils/create_samplesheet.py` internally):
   ```bash
   python scripts/03_rnaseq_pipeline/01_create_samplesheet.py
   ```
2. Test with one sample first:
   ```bash
   sbatch --array=1 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh
   ```
3. Auto-monitor and launch remaining samples:
   ```bash
   sbatch scripts/03_rnaseq_pipeline/03_monitor_and_launch.sh <TEST_JOB_ID>
   ```

**COMPLETED:** All 44 samples processed successfully (Oct 17-18, 2025)
- Test job validated: 22,176 genes detected (was 0 before GTF fix!)
- Monitor script auto-launched remaining 43 samples
- Total runtime: ~11 hours

**Phase 3: Analysis (Steps 04-07)**
5. `04_qc_analysis/` - Extract QC metrics and generate figures (NEXT)
6. `05_count_matrix/` - Combine count matrices
7. `06_diff_expression/` - DESeq2 analysis
8. `07_visualization/` - Publication figures

### Script Numbering Convention

**Format:** `##_descriptive_name.ext`
- Leading zeros for proper sorting (01, 02, 03, ..., 10, 11)
- Descriptive names after number
- Grouped by logical workflow stages (directories)
- Each directory has README.md with detailed instructions

**Why This Approach:**
✅ Follows bioinformatics community best practices
✅ Natural sorting works correctly
✅ Easy to add scripts without renumbering everything
✅ Self-documenting workflow
✅ Prepares for Snakemake conversion
✅ Clear for manuscript reproducibility section

### Key Scripts - What Actually Works

**✅ ACTIVE SCRIPTS (Use these):**
- `scripts/02_annotation_prep/01_convert_gff_to_gtf.py` - GFF3→GTF conversion
- `scripts/02_annotation_prep/02_filter_problematic_transcripts.py` - Remove 45 bad transcripts
- `scripts/02_annotation_prep/03_add_gene_biotype.py` - Add gene_biotype (bug fixed)
- `scripts/utils/create_samplesheet.py` - Generate nf-core samplesheet (used by 03_rnaseq_pipeline)
- `scripts/01_sra_download/01_sra_array.sh` - Download FASTQ from SRA
- `scripts/01_sra_download/02_sra_download.py` - Called by array script
- `scripts/01_sra_download/03_info_sra_checker.py` - Verify downloads
- `scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh` - Main pipeline script
- `scripts/03_rnaseq_pipeline/03_monitor_and_launch.sh` - Auto-monitor and launch

**⚠️ DEPRECATED SCRIPTS (Do NOT use):**
- `scripts/utils/fix_gtf_gene_id.py` - Old approach, superseded by 3-step workflow
- `scripts/utils/organize_outputs.sh` - For old laptop workflow, not HPC

**📁 Reference Copies (Canonical versions in numbered directories):**
- `scripts/utils/convert_gff_to_gtf.py` - Reference copy (use version in 02_annotation_prep)
- `scripts/utils/gff_to_rsem_gtf.py` - Reference copy (use version in 02_annotation_prep)
- `scripts/utils/add_gene_biotype.py` - Reference copy (use version in 02_annotation_prep)

### Key Documentation

Each directory contains a `README.md` with:
- Purpose and overview
- Input/output files
- Script descriptions and usage
- Expected runtime
- Troubleshooting tips
- Next steps

**To reproduce full workflow:** Read `scripts/00_reference_preparation/README.md` and follow numbered directories in order.

---

## TODO LIST

### Immediate (October 2025)
- [ ] Test cleaned GTF with nf-core on single sample
- [ ] Verify gene detection no longer shows zeros
- [ ] If test passes, run full array for all 44 samples
- [ ] Extract QC metrics from new run

### QC Analysis
- [ ] Fix gene count detection in QC extraction script
- [ ] Create per-project QC summary tables
- [ ] Generate standard QC figures with project labels
- [ ] Assess batch effects between platforms

### Differential Expression
- [ ] Gather count matrices from all samples
- [ ] Import 34 GWAS candidate genes list
- [ ] Run DESeq2 per project
- [ ] Compare with collaborator's published results

### Validation with Collaborator Data
- [ ] Extract gene counts from Salmon outputs
- [ ] Load collaborator HTSeq count matrices
- [ ] Calculate correlation between our counts and theirs (expect >0.95)
- [ ] Create scatter plots comparing count methods
- [ ] Compare transcript counts (our advantage - isoform resolution)
- [ ] Validate 34 GWAS candidate gene expression patterns
- [ ] Document any discrepancies and investigate causes
- [ ] Create validation report for manuscript supplement

### Publication
- [ ] Write methods section including annotation preprocessing
- [ ] Create manuscript-ready figures
- [ ] Prepare supplementary data files
- [ ] Upload reference files to Zenodo

---

## KNOWN ISSUES AND SOLUTIONS

### Issue: Previous run showed 0 detected genes
**Cause**: Wrong reference annotation
**Solution**: Using correct AalbF3_annotation.gtf.gz with proper gene_id attributes

### Issue: rRNA showed 0% in previous run
**Cause**: No rRNA removal step
**Solution**: Added SortMeRNA with Ae. albopictus rRNA sequences

### Issue: RSEM failed with multi-chromosome transcripts
**Cause**: Annotation errors in NCBI AalbF3
**Solution**: Filtered 45 problematic transcripts (0.2% loss)

### Issue: Platform batch effects (GAIIx vs HiSeq2000)
**Cause**: Different experiments used different sequencers
**Solution**: Include platform as covariate in DESeq2 models or use ComBat-seq

---

## IMPORTANT PATHS

### Project Base
```
/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/
```

### Container
```
/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif
```

### References
```
/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/data/references/AalbF3_genome/
```

### Outputs
```
/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/output/PRJNA*/SRR*/
```

### Collaborator Code (for reference)
```
/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/data/collaborator_repos/albopictus_remapping/
```

---

## GITHUB AUTHENTICATION

**IMPORTANT:** Git authentication requires loading your saved token before push/pull operations.

```bash
# Load saved GitHub token (for private repos)
export GITHUB_TOKEN=$(cat ~/.github_token)

# Configure git to use the token
git config --global credential.helper store
git config --global url."https://${GITHUB_TOKEN}@github.com/".insteadOf "https://github.com/"

# Now you can push/pull
git push
```

**Token location:** `~/.github_token` (created October 16, 2025)
**Note:** Token is in .gitignore and will never be committed

---

## DATA MANAGEMENT STRATEGY

### Git Repository
- Track: Scripts, small config files, samplesheets, documentation
- Ignore: Large reference files, output data, containers, cache

### Zenodo Upload (Later)
- Reference genome: `AalbF3_genome.fa.gz` (423 MB)
- Reference annotation: `AalbF3_annotation.gtf.gz` (4.9 MB)
- rRNA sequences: `combined_rRNA_sequences.fa` (37 KB)
- Create download script to pull from Zenodo

### Reproducibility Package
Container will include:
1. All analysis scripts
2. Sample metadata
3. SRA download scripts
4. nf-core pipeline configuration
5. Reference to Zenodo for data files

---

## TIMELINE

### Completed (October 2025)
- [x] Downloaded all 44 SRA samples
- [x] Downloaded correct AalbF3 references
- [x] Fixed chromosome naming
- [x] Created rRNA database
- [x] Fixed GTF annotation for RSEM
- [x] Updated pipeline script
- [x] Documented all annotation problems

### In Progress
- [ ] Testing pipeline with fixed references (1 sample)
- [ ] Full re-run of all 44 samples

### Upcoming
- [ ] QC analysis and figures
- [ ] Differential expression per project
- [ ] Meta-analysis across life stages
- [ ] Manuscript preparation

---

## QUESTIONS FOR USER

### Pending Decisions
1. Where to get the 34 GWAS candidate genes list?
2. Which DESeq2 model for batch effect correction?
3. Threshold for excluding low-quality samples?

---

## NOTES

- All nf-core runs use consistent parameters (see pipeline script)
- Skipped some QC steps to speed up processing (bigwig, stringtie, dupradar, preseq)
- Saved trimmed reads and unaligned reads for potential re-analysis
- Container includes all dependencies for full reproducibility
- Every analytical decision documented for reviewer transparency
