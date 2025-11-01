# Albopictus Diapause RNA-seq Project - Complete Information

Last updated: October 31, 2025, 5:40 PM

## Today's Session (Oct 31, 2025)

**What we did:**
1. ✅ Consolidated tracking/ directory - removed 40+ duplicate .md files
2. ✅ Added correct biological information from the 3 manuscript PDFs
3. ✅ Moved visual_progress.py and visual_progress.txt to tracking/
4. ✅ Reorganized project_state.json - moved tracking setup to Phase 0
5. ✅ Updated visual_progress.py to show exact commands to run
6. ✅ Created tracking/README.md explaining how to use the system

**What we discovered:**
- results/ directory has MANY duplicate directories (qc_analysis, qc_figures, qc_metrics, count_matrices x4, validation x2)
- Scripts directory is well-organized and numbered correctly
- Need to clean up results/ to match analysis strategy

**What still needs to be done:**
1. Rerun failed QC jobs (priority)
2. Reorganize results/ directory structure
3. Update scripts to use new results/ paths
4. Continue with Phase 3 (Salmon validation)

## How to Use the Tracking System

Check current status:
```bash
module load singularity
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py status"
```

Get next task:
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py next"
```

Mark task complete:
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py complete <task_id>"
```

## Current Status - WHERE WE ARE RIGHT NOW (Oct 31, 2025 - 5:56 PM)

**AT: Phase 2 - Reorganize Results Directory Structure**
**Progress: 15% overall (honest assessment)**

**THE TRUTH:**
- We have a MESS in results/ directory with many duplicates
- We need to REORGANIZE the entire results/ structure
- We will need to RERUN scripts after reorganizing
- Previous "Phase 2 complete" was WRONG - resetting to Phase 2

**What needs to happen NEXT (in order):**
1. Delete or archive duplicate results directories
2. Create new clean results/ structure
3. Update all scripts to use new paths
4. Rerun count matrix organization scripts
5. Rerun ALL QC scripts to populate new structure

**AFTER reorganization:**
- Then run differential expression (Phase 4-6)
- Then validate methods (Phase 3)
- Then answer reviewers (Phase 10)
- Then GWAS candidates (Phase 11)

## What's Actually Complete (Honest Status)

### Phase 0: Setup - COMPLETE ✅
- Container pulled
- Reference genome downloaded from Dryad
- Tracking system setup

### Phase 1: Pipeline - COMPLETE ✅
- 45 samples downloaded from SRA
- nf-core/rnaseq pipeline executed (Job 20677337) - all 45 samples processed successfully
- Pipeline outputs in output/star_salmon/

### Phase 2: Reorganize Results - IN PROGRESS (0% done)
- NOT complete - we have a mess to clean up
- Need to reorganize entire results/ directory
- Need to rerun scripts after reorganization

### What exists:
- `output/star_salmon/featurecounts/*.tsv` - 22,177 genes x 45 samples
- `output/star_salmon/*/quant.sf` - Salmon quantification for all samples
- `output/count_matrices/` - Organized matrices by stage

## Project Overview

**Goal**: Re-analyze published RNA-seq data on new AalbF3 genome to validate GWAS candidate genes for diapause

**Data**: 45 samples from 3 independent studies
- Adults (16 samples) - PRJNA268379 - Huang et al. 2015
- Embryos (12 samples) - PRJNA158021 - Poelchau et al. 2013a
- Larvae (17 samples) - PRJNA187045 - Poelchau et al. 2013b

## The 45 Samples - Experimental Design

### Adults (16 samples) - 2x2 factorial
**Study**: Huang et al. 2015, PLOS Neglected Tropical Diseases
**Platform**: Illumina HiSeq 2000, Paired-end 101bp
**Factors**: Photoperiod (SD vs LD) x Blood meal (Non-blood-fed vs Blood-fed)
- SD_Unfed (n=4): SRR1663685, SRR1663687, SRR1663689, SRR1663695
- SD_Fed (n=4): SRR1663697, SRR1663700, SRR1663703, SRR1663707
- LD_Unfed (n=4): SRR1663709, SRR1663754, SRR1663769, SRR1663911
- LD_Fed (n=4): SRR1663913, SRR1663915, SRR1663917, SRR1663919

**Sample details**:
- Adult females (whole body, blood bolus removed from midgut)
- F3 laboratory generation from Manassas, VA population
- Reared at 21°C, 80% RH for 11 days under assigned photoperiods
- Blood-fed samples collected 26-28h post-blood meal
- Diapause incidence: SD = 81-95%, LD = 1-3%

**Key biology**: Photoperiod and blood feeding affect diapause induction in adults. Short-day photoperiods trigger early circadian and metabolic shifts. Blood feeding under short-day conditions amplifies oxidative and lipid metabolism for provisioning diapause eggs. SD+Blood fed produces diapause eggs, LD+Blood fed produces non-diapause eggs.

**Analysis**: Split by blood meal status (unfed vs fed analyzed separately)

### Embryos (12 samples) - 2x2 factorial time series
**Study**: Poelchau et al. 2013, Proceedings of the Royal Society B
**Platform**: Illumina HiSeq, Paired-end
**Factors**: Photoperiod (D vs ND) x Time (3d vs 6d post-oviposition)
- D_72h (n=3): SRR458464, SRR458465, SRR458466
- D_135h (n=3): SRR458467, SRR458468, SRR458469
- ND_72h (n=3): SRR458470, SRR458471, SRR458472
- ND_135h (n=3): SRR458473, SRR458474, SRR458475

**Sample details**:
- F9 laboratory strain
- Reared at 21°C, ~80% RH
- D (diapause, 8L:16D) photoperiod produces ~99% diapause eggs
- ND (non-diapause, 16L:8D) photoperiod produces non-diapause eggs
- 3d (72-78h): germ band retraction/dorsal closure stage
- 6d (135-141h): late embryogenesis (segmentation, head-thorax separation)

**Key biology**: Diapause preparation in embryos. At 3d, diapause embryos show early cell-cycle/ecdysone signaling shifts. At 6d, they reduce energy metabolism and fortify cuticular structures, preparing for metabolic suppression and desiccation resistance before dormancy.

**Analysis**: Split by timepoint (72h vs 135h analyzed separately)

### Larvae (17 samples) - 2x3 factorial (unbalanced)
**Study**: Poelchau et al. 2013b, Journal of Experimental Biology (PRJNA187045)
**Platform**: Illumina GA IIx, Paired-end 100bp
**Factors**: Photoperiod (SD vs LD) x Time (11d, 21d, 40d post-oviposition)
- SD_11d (n=3): SRR652090, SRR652091, SRR652092
- SD_21d (n=3): SRR652093, SRR652094, SRR652095
- SD_40d (n=3): SRR652096, SRR652097, SRR652098
- LD_11d (n=3): SRR652099, SRR652100, SRR652101
- LD_21d (n=3): SRR652116, SRR652117, SRR652118
- LD_40d (n=2): SRR652119, SRR652120

**Sample details**:
- Pharate larvae (fully developed larvae still in egg chorion)
- Comparing diapause (SD, 8L:16D) vs non-diapause/quiescence (LD, 16L:8D)
- Three timepoints: 11d, 21d, 40d post-oviposition
- Note: Unbalanced design at 40d (LD has n=2 instead of n=3)

**Key biology**: Study examines pharate larval stage comparing diapause vs quiescence at extended timepoints. Diapause larvae remain viable in egg chorion for extended periods. Non-diapause larvae would normally hatch around 3-4 days under favorable conditions.

**Analysis**: Each timepoint analyzed separately (11d, 21d, 40d)

## Critical Confounding

**Platform = Stage = Study** - Cannot compare across stages directly. Must analyze each stage independently.

## The 3-Part Analysis Plan

### Part 1: Replicate Collaborators (featureCounts)
- Use featureCounts gene counts (22,177 genes)
- Run collaborator's exact DESeq2 methods (see below)
- Goal: Validate our pipeline reproduces their results

### Part 1.5: Validate Salmon vs featureCounts (gate check)
- Must prove Salmon matches featureCounts
- Requirements: r > 0.95 correlation, >90% DEG overlap
- If fails: stop, do not proceed to Part 2

### Part 2: Answer Reviewers (Salmon - only after validation)
- Focus on genome-wide statistics
- Fisher's exact test for GWAS enrichment (4/5 genes DE)
- Sample information table
- Mapping statistics from MultiQC

### Part 3: GWAS Candidates (last)
- Extract candidate expression after Parts 1-2 complete

## Collaborator DESeq2 Methods (exact replication)

### Adults - Angela
- Split by blood meal: NB (unfed) and BM (fed) separately
- Design: `~ condition` where condition = SD vs LD
- Reference: "LD"
- LFC shrinkage: `lfcShrink(dds, coef="condition_SD_vs_LD", type="apeglm")`
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Significance: padj < 0.05 and abs(log2FC) > 0.58

### Embryos - Mackenzie
- Split by timepoint: 72h and 135h separately
- Design: `~ condition` where condition = DI vs NDI (not SD/LD!)
- Reference: "NDI"
- LFC shrinkage: `lfcShrink(dds, coef="condition_DI_vs_NDI", type="apeglm")`
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Significance: padj < 0.05 and abs(log2FC) > 0.58

### Larvae - Sarah
- Split by timepoint: 11d, 21d, 40d separately
- Design: `~ condition` where condition = D vs ND (not SD/LD!)
- Reference: "ND"
- LFC shrinkage: `lfcShrink(dds, coef="condition_D_vs_ND", type="apeglm")`
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Significance: padj < 0.05 and abs(log2FC) > 0.5 (note: 0.5, not 0.58)

## Container Usage (mandatory)

**Always start with**:
```bash
module load singularity
```

**Python scripts**:
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python script.py"
```

**R scripts**:
```bash
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && Rscript script.R"
```

Container location: `/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif`

## Directory Structure

```
output/
├── star_salmon/           # Pipeline outputs (complete)
│   ├── featurecounts/     # Gene counts for Part 1
│   └── */quant.sf         # Salmon for Part 2
├── count_matrices/        # Organized matrices (complete)
│   ├── featurecounts/
│   └── salmon/
└── results/
    └── qc_visualization/  # Empty - jobs failed
```

## Scripts Created (in execution order)

### Data preparation (complete):
1. `scripts/03_sra_download/01_sra_array.sh` - Download SRA data
2. `scripts/04_annotation_preprocessing/00_convert_gff3_to_gtf.py` - Convert GFF to GTF
3. `scripts/05_rnaseq_pipeline/01_create_samplesheet.py` - Create sample sheet
4. `scripts/05_rnaseq_pipeline/02_run_rnaseq.sh` - Run pipeline (Job 20677337)

### Count matrices (complete):
5. `scripts/06_count_matrix/01_combine_counts.py` - Combine counts
6. `scripts/06_count_matrix/03_sanity_check_counts.py` - Verify counts

### QC visualization (failed - need to rerun):
7. `scripts/07_qc_visualization/07_submit_qc_featurecounts.sh` - Failed job 20699269
8. `scripts/07_qc_visualization/08_submit_qc_salmon.sh` - Failed job 20699276

### Differential expression (not started):
9. `scripts/08_differential_expression/` - Multiple scripts created, not executed

## Reviewer Requirements

### Reviewer 2:
1. How many genes DE genome-wide at different thresholds
2. Fisher's exact test: Is 4/5 GWAS candidates being DE significant?

### Reviewer 3:
1. Sample information table (all 45 samples)
2. DESeq2 design explanation
3. Stage-specific expression patterns
4. Mapping statistics from MultiQC

## NEXT SESSION - START HERE (Tomorrow, Nov 1, 2025)

**PRIORITY: Reorganize results/ directory structure**

**Step 1: Review current tracking:**
```bash
cd /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq
python tracking/00_track.py status
cat tracking/visual_progress.txt
```

**Step 2: Decide on new results/ structure together**
- Look at what we have in results/
- Decide what to keep, what to archive, what to delete
- Design clean structure that matches analysis strategy

**Step 3: Execute reorganization**
- Move/rename/delete directories
- Update scripts with new paths
- Rerun scripts to populate new structure

**DO NOT just run QC jobs - we need to reorganize first!**

The old QC results are in a mess of duplicate directories. We need to clean up and start fresh with organized structure.