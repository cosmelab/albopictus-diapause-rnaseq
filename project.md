# Albopictus Diapause RNA-seq Analysis Project

**Last Updated:** October 17, 2025
**Status:** Ready to re-run pipeline with fixed AalbF3 references

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

## CURRENT STATUS - OCTOBER 17, 2025

### Ready to Re-run Pipeline ✓
- [x] Downloaded correct AalbF3 genome from Dryad (GCA_018104305.1)
- [x] Fixed chromosome naming (stripped "chr" prefix to match SNP chip format)
- [x] Downloaded 12 Ae. albopictus rRNA sequences for SortMeRNA
- [x] Created GTF with proper gene_id attributes for RSEM compatibility
- [x] Filtered 45 problematic transcripts (multi-chromosome/inconsistent strands)
- [x] Updated pipeline script to use new references and rRNA removal
- [x] Cleaned up 2TB of old cache files
- [x] **DOCUMENTED ALL ANNOTATION PROBLEMS FOR MANUSCRIPT METHODS**

### What Changed from Previous Run
1. **Genome**: `albo.fasta.gz` (old) → `AalbF3_genome.fa.gz` (correct, chr prefix stripped)
2. **Annotation**: `genes_fixed.rsem.fixed.gtf` (old) → `AalbF3_annotation.gtf.gz` (fixed, RSEM-compatible)
3. **Added rRNA removal**: SortMeRNA with 12 Ae. albopictus sequences
4. **Fixed annotation errors**: See "GFF3/GTF Problems" section below

### Next Step
Test with one sample:
```bash
sbatch --array=1 scripts/02_run_rnaseq/01_run_rnaseq_array.sh
```

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
- **Assembly**: Matthews et al. 2018, G3 (Chromosome-level assembly)
- **Annotation**: Palatini et al. 2020, BMC Biology
- **Accession**: GCA_018104305.1
- **Repository**: Dryad (doi:10.5061/dryad.mgqnk98z4)
- **Downloaded**: October 17, 2025

### Why This Reference?
- Same genome used by collaborators (srmarzec/albopictus_remapping)
- Complete GFF3 annotation with rRNA, tRNA, and all gene features
- Chromosome-level assembly (better than fragmented assemblies)
- Public and citable (required for publication)

### Reference Preprocessing Steps
1. **Download from Dryad** (NCBI blocks automated downloads)
   - Manual download: https://datadryad.org/stash/dataset/doi:10.5061/dryad.mgqnk98z4
   - Files: `AALBF3_assembly.fa.gz`, `AALBF3_annotation.gff3.gz`

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

## GFF3/GTF ANNOTATION PROBLEMS (FOR MANUSCRIPT METHODS)

**IMPORTANT:** Document these 4 problems in the manuscript methods section.

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

### Summary for Methods Section
The NCBI AalbF3 annotation (GCF_018104305.1) required four preprocessing steps for RSEM compatibility:
1. File extension correction (`.gff3.gz` → `.gff.gz`)
2. Addition of 18,135 missing gene_id attributes
3. GTF syntax correction (removal of double semicolons)
4. Filtering of 45 transcripts with annotation errors (0.2% of transcripts)

**Final annotation: 347,380 features from 32,889 valid transcripts**

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

### Our Approach
- nf-core/rnaseq v3.19.0 (automated, containerized)
- Salmon quantification (faster, includes TPM)
- Comprehensive QC (MultiQC + Qualimap + RSeQC)
- Adaptive trimming
- rRNA removal with SortMeRNA

### Advantages of Our Approach
1. **Reproducibility**: Containerized with locked versions
2. **Speed**: Salmon is much faster than STAR+HTSeq
3. **QC**: More comprehensive metrics
4. **Automation**: Less manual intervention, fewer errors
5. **Best practices**: Following current RNA-seq standards

### Addressing Reviewers
When reviewers ask about differences:
1. Same samples and experimental designs
2. Modern pipeline improves accuracy and reproducibility
3. Salmon validated to produce comparable results to STAR+HTSeq
4. Better QC filtering and improved algorithms
5. Any differences in DEGs likely due to improved methods

---

## SCRIPTS AND UTILITIES

### Data Download Scripts
```
scripts/01_download_data/
├── download_rrna_sequences.sh        # Get 12 Ae. albopictus rRNA from GenBank
├── clean_rrna_headers.sh             # Simplify headers for SortMeRNA
└── strip_chr_prefix.sh               # Remove chr prefix from references
```

### Pipeline Scripts
```
scripts/02_run_rnaseq/
└── 01_run_rnaseq_array.sh            # Main SLURM array job script
```

### Utility Scripts
```
scripts/utils/
├── fix_gtf_gene_id.py                # Add missing gene_id attributes
└── gff_to_rsem_gtf.py                # Filter problematic transcripts
```

### QC Analysis Scripts
```
scripts/03_qc_analysis/
├── extract_qc_metrics.py             # Extract from MultiQC outputs
└── run_qc_analysis.sh                # Wrapper script
```

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

```bash
# Load saved GitHub token (for private repos)
export GITHUB_TOKEN=$(cat ~/.github_token)
# Token created: October 16, 2025
```

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
