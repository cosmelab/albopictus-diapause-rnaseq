# Session Tracking - Current Work

**Last Updated:** Oct 23, 2025 18:45

---

## QUICK START GUIDE

**When resuming work, read these files in order:**

1. **THIS FILE (session_tracking.md)** - Current status, what's running, next steps
2. **assistant_rules.md** - Rules to follow (container-only, no caps in filenames, etc.)
3. **project.md** - Full project context, goals, samples, GWAS candidates
4. **stage_specific_analysis_roadmap.md** - Analysis strategy for each developmental stage
5. **technical.md** - HPC setup, container info, GitHub/Docker setup

**Current Status:**
- Pipeline Job 20644937 running (started 18:09, ~2-3 hours)
- Scripts cleaned up, committed (f758c74)
- Container rebuilding with rdryad

---

## THE ACTUAL ANALYSIS PLAN (STOP IGNORING THIS)

We have 3 experiments. We ran nf-core/rnaseq and got:
- **featureCounts** (HTSeq-compatible counts)
- **Salmon counts** (transcript-level quantification)

### Step 1: REPLICATE COLLABORATORS FIRST (Using HTSeq/featureCounts)

**Goal:** See if we can reproduce their findings using our new genome assembly

- ✅ Document their exact methods (DONE - see below)
- Use their exact DESeq2 approach:
  - Split by timepoint/treatment (NOT collapsed)
  - lfcShrink with apeglm
  - Their significance thresholds
- Generate same analyses they did
- *Optional:* Get their deposited counts if available, but maybe not necessary
- Compare patterns/results (accounting for genome assembly differences)

**Why?** Validate our pipeline produces sensible results before doing more complex models

### Step 2: PROPER STATISTICAL MODELS (Reviewers' Recommendations)

**Goal:** Robust analysis accounting for ALL experimental variables

**Use BOTH HTSeq and Salmon counts for:**

1. **Proper statistical framework**
   - Models that account for:
     - **Diapause** (SD vs LD / D vs ND / DI vs NDI)
     - **Time/Development** (embryo: 72h, 135h; larvae: 11d, 21d, 40d)
     - **Blood meal** (adults: unfed vs fed)
     - **Batch/Platform** (confounded with stage, need careful handling)

2. **Global gene expression analysis**
   - Patterns across all 3 life stages
   - Meta-analysis across experiments
   - Cross-stage comparisons

3. **Statistical power estimation**
   - For each comparison
   - Sample size considerations

4. **Method comparison**
   - HTSeq vs Salmon using proper models
   - Assess concordance

### Step 3: ANSWER ALL MOLECULAR ECOLOGY REVIEWER QUESTIONS

- Address each reviewer concern systematically
- Provide proper statistical justification
- Show reproducible workflow

---

## CURRENT STATUS

### Completed ✅
1. nf-core/rnaseq pipeline - got HTSeq and Salmon counts
2. Method validation: Salmon vs HTSeq (R²=0.98)

### COLLABORATOR METHODS - COMPLETE ✅

**Pipeline (all 3 experiments):**
1. Trimmomatic v0.39 (HEADCROP:15 - remove first 15bp)
2. STAR v2.7.1a (alignment)
3. HTSeq-count (gene counting)
4. DESeq2 (differential expression)

**Angela - Adults (data/collaborator_repos/scripts/Angela/Angela_DESeq.R):**
- **Split by blood meal:** Analyzed NB (no blood meal/unfed) and BM (blood meal/fed) SEPARATELY
- Design: `~ condition` (SD vs LD)
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Reference: "LD" (long day, non-diapause)
- **LFC shrinkage:** `lfcShrink(dds, coef="condition_SD_vs_LD", type="apeglm")`
- **Significance:** `padj < 0.05` AND `abs(log2FoldChange) > 0.58`
- **Volcano:** pCutoff = 10e-6, FCcutoff = 0.58

**Mackenzie - Embryos (data/collaborator_repos/scripts/Mackenzie/Mackenzie_DESeq.R):**
- **Split by timepoint:** Analyzed 72h and 135h SEPARATELY
- Design: `~ condition` (DI vs NDI)
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Reference: "NDI" (non-diapause-inducing)
- **LFC shrinkage:** `lfcShrink(dds, coef="condition_DI_vs_NDI", type="apeglm")`
- **Significance:** `padj < 0.05` AND `abs(log2FoldChange) > 0.58`
- **Volcano:** pCutoff = 10e-6, FCcutoff = 0.58

**Sarah - Larvae (data/collaborator_repos/scripts/Sarah/DESeq.R):**
- **Split by timepoint:** Analyzed 11d, 21d, 40d SEPARATELY
- Design: `~ condition` (D vs ND)
- Pre-filter: `rowSums(counts(dds)) >= 10`
- Reference: "ND" (non-diapause)
- **LFC shrinkage:** `lfcShrink(dds, coef="condition_D_vs_ND", type="apeglm")`
- **Significance:** `padj < 0.05` AND `abs(log2FoldChange) > 0.5`
- **Volcano:** pCutoff = 10e-6, FCcutoff = 0.5

**KEY INSIGHT:** All analyzed by timepoint/treatment SEPARATELY, NOT collapsed!

### WHAT WE HAVE NOW
- HTSeq counts from nf-core (featureCounts, compatible with HTSeq)
- Salmon counts from nf-core
- Both using new AalbF3 genome (collaborators used older assembly)

### OUR SCRIPTS STATUS

**Container:**
- ✅ Added all collaborator packages to Dockerfile (VERIFIED in conda channels):
  - **From bioconda:** bioconductor-apeglm, bioconductor-enhancedvolcano, bioconductor-mygene
  - **From conda-forge:** r-pheatmap, r-reshape2, r-viridis
  - **From CRAN:** GOplot (not in conda-forge/bioconda)
- ✅ Optimized installation order:
  1. r-base first
  2. R packages from conda (pheatmap, reshape2, viridis)
  3. Bioconductor packages
  4. CRAN packages (remaining)
- ✅ Pushed to GitHub (commit eda94db)
- 🔄 **Container build TRIGGERED** - check GitHub Actions
- ⏳ Once build completes: `singularity pull albopictus-diapause-rnaseq.sif docker://ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest`

**Adults (scripts/06_differential_expression/):**
- ✅ 01_adults_deseq2_unfed.R - UPDATED to match Angela's methods exactly
  - Added lfcShrink with apeglm
  - Added |log2FC| > 0.58 threshold for significance
  - Switched to pheatmap
  - Creates side-by-side MA plots (unshrunk vs shrunk)
- ✅ 02_adults_deseq2_fed.R - UPDATED to match Angela's methods exactly
  - Same updates as unfed script
- 🔄 Old jobs 20616368 & 20616372 may still be running - results won't match collaborators
- ⏳ READY to rerun once container is rebuilt

**Embryos:** Not created yet (need 72h and 135h separate analyses)
**Larvae:** Not created yet (need 11d, 21d, 40d separate analyses)

### PIPELINE RE-RUN STATUS (Oct 23, 2025 - 18:40)

**THE ACTUAL PROBLEM (VERIFIED):**
- Oct 18 run used: `-g gene_biotype` in featureCounts
- This groups ALL genes by biotype instead of counting individual genes
- Result: ~10 biotype categories instead of 22,176 individual gene counts
- **PROOF:** `/bigdata/.../output/nf-work/00/bf2f02386ad32b8f20d50a6737507b/.command.sh` line 3

**Job 20644937 RUNNING (started 18:09, running 29+ minutes)**
- **FIXED PARAMETER:** `featurecounts_group_type: gene_id` in params.yaml
- **References:** Using same files that worked Oct 18
  - FASTA: data/references/AalbF3/AalbF3_genome.fa.gz
  - GTF: data/references/AalbF3/AalbF3_annotation.gtf.gz
  - rRNA: data/references/AalbF3/combined_rRNA_sequences.fa
- **STATUS:** Pipeline running with -resume
  - ✅ Reusing trimmed reads (cached)
  - ✅ Reusing QC (cached)
  - 🔄 STAR genome generation in progress
  - Will re-run: featureCounts with correct parameter, then merge results
- **Expected completion:** ~2-3 hours (not 48-72 because -resume works)

### IMMEDIATE NEXT STEPS (when pipeline completes)

**Step 1: Verify Pipeline Outputs**
1. ⏳ Wait for Job 20644937 to complete (~2-3 hours from 18:09)
2. ⏳ Check featureCounts files have individual gene counts (not biotypes)
3. ⏳ Verify gene count: should be ~22,176 genes
4. ⏳ Check one featureCounts .command.sh to confirm `-g gene_id` was used

**Step 2: Combine Count Matrices**
5. ⏳ Run: `scripts/04b_count_matrix/01_combine_counts.py`
6. ⏳ Run: `scripts/04_qc_analysis/08_split_counts_by_stage.py`
7. ⏳ Run: `scripts/04_qc_analysis/09_split_counts_with_metadata.py`

**Step 3: Replicate Collaborators' DESeq2 Analyses**
8. ⏳ Adults unfed: `scripts/06_differential_expression/01_adults_deseq2_unfed.R`
9. ⏳ Adults fed: `scripts/06_differential_expression/02_adults_deseq2_fed.R`
10. ⏳ Create embryo scripts: 72h and 135h separate (Mackenzie's methods)
11. ⏳ Create larvae scripts: 11d, 21d, 40d separate (Sarah's methods)

**Step 4: GWAS Enrichment**
12. ⏳ Run GWAS enrichment analyses for all stages

### AFTER REPLICATION (Step 2: Reviewer-Recommended Models)

**Design proper models accounting for:**
- Diapause status + Time + Blood meal (adults) + Batch effects
- Use BOTH HTSeq and Salmon counts
- Statistical power estimation
- Global gene expression patterns across all 3 life stages
- Meta-analysis framework

### FINAL STEP (Step 3: Address Reviewers)

- Systematically answer each Molecular Ecology reviewer question
- Provide statistical justification for all analyses
- Document reproducible workflow

---

## SCRIPTS INVENTORY & CLEANUP PLAN

### Scripts Status by Directory:

**00_reference_preparation/** ✅ DONE (3 scripts)
- 01_strip_chr_prefix.sh, 02_download_rrna_sequences.sh, 03_clean_rrna_headers.sh

**01_data_acquisition/** ⏭️ OBSOLETE - FILES ALREADY EXIST
- download_dryad_genome.py - TESTED: Gets 403 Forbidden error
- process_dryad_genome.sh - Never run
- Files manually downloaded to: data/references/dryad_download/doi_10_5061_dryad_mgqnk98z4__v20210127/

**01_sra_download/** ✅ DONE (3 scripts)
- Downloaded all 44 SRA samples

**02_annotation_preprocessing/** ✅ DONE (2 scripts)
- Converted GFF3 to GTF, fixed transcript issues
- See: scripts/02_annotation_preprocessing/gtf_preprocessing_pipeline.md for details

**03_rnaseq_pipeline/** 🔄 RUNNING
- Job 20644937 running with corrected featureCounts parameter

**04_qc_analysis/** ⚠️ NEEDS CLEANUP - 14 scripts, only 3 needed
- KEEP: 08_split_counts_by_stage.py, 09_split_counts_with_metadata.py, 10_create_pca_plots.R
- DELETE: 01-07, 11-14 (obsolete/duplicates)

**05_batch_correction/** ❌ DELETE ENTIRE DIRECTORY
- Decision: NO batch correction (Platform = Stage = confounded)

**05_count_matrix/** ⏳ PENDING (2 scripts) - RENAME to 04b_count_matrix/

**06_differential_expression/** ⏳ MAIN - needs completion
- Have: adults unfed/fed scripts
- Need: embryos (72h, 135h), larvae (11d, 21d, 40d)

**06_diff_expression/** ❌ DELETE - superseded by R DESeq2 approach

**07_visualization/** ⏳ PENDING

**09_method_validation/** ✅ DONE

### Cleanup Status: ✅ DONE (Oct 23, 2025)
- Deleted: scripts/01_data_acquisition/, 02_annotation_prep/, 05_batch_correction/, 06_diff_expression/
- Renamed: 05_count_matrix/ → 04b_count_matrix/
- Cleaned 04_qc_analysis/: kept only scripts 08, 09, 10
- **Committed and pushed** (commit f758c74)
- **GitHub Actions building new container with rdryad**

---

## KEY INFO TO REMEMBER

### Analysis Strategy & Progress

**Adults (16 samples):** 2×2 (Photoperiod × BloodMeal)
- Strategy: SPLIT by blood meal (following Huang et al. 2015)
- ✅ Script 01_adults_deseq2_unfed.R created
- ✅ Script 01_run_adults_unfed.sh created
- 🔄 Job 20616368 RUNNING (submitted Oct 22, 16:35)
- ✅ Script 02_adults_deseq2_fed.R created
- ✅ Script 02_run_adults_fed.sh created
- 🔄 Job 20616372 RUNNING (submitted Oct 22, 16:35)
- Design: `~ Photoperiod` (SD vs LD, 4 vs 4 each)

**Embryos (12 samples):** 2×2 (Photoperiod × Time: 72h, 135h)
- Strategy: COLLAPSE timepoints → all LD vs all SD
- ⏳ Script 03 not created yet
- Design: `~ Photoperiod` (SD vs LD, 6 vs 6)

**Larvae (17 samples):** 2×3 (Photoperiod × Time: 11d, 21d, 40d)
- Strategy: COLLAPSE timepoints → all LD vs all SD
- ⏳ Script 04 not created yet
- Design: `~ Photoperiod` (SD vs LD, 9 vs 8)
- Note: Unbalanced design (LD_40d only n=2)

### Multiple Testing
- DESeq2 uses Benjamini-Hochberg FDR (padj column)
- Collaborators use: FDR < 0.05 AND |log2FC| > 0.5
- We currently use: FDR < 0.05 only

### Container Rule - CRITICAL
- NEVER install outside container
- **NEVER use HPC system Python directly** - home quota is only 20GB!
- **ALL scripts must run with:** singularity exec albopictus-diapause-rnaseq.sif python/Rscript
- pheatmap missing? Use base R or update Dockerfile

---

## FILES WE WORK WITH

**Main docs:**
- project.md - project details
- technical.md - HPC/container info
- assistant_rules.md - rules to follow
- session_tracking.md - THIS FILE

**Scripts created:**
- scripts/09_method_validation/01_compare_salmon_featurecounts.py ✅
- scripts/06_differential_expression/01_adults_deseq2_unfed.R (needs pheatmap fix)
- scripts/06_differential_expression/02_adults_deseq2_fed.R (needs pheatmap fix)
- scripts/06_differential_expression/run_adults_unfed_deseq2.sh ✅

**Results:**
- results/07_method_comparison/ - validation done ✅
- results/03_differential_expression/adults/ - empty, waiting for fixed scripts

---

## STOP GETTING LOST CHECKLIST

Before doing ANYTHING:
- [ ] Check this file for current status
- [ ] Check assistant_rules.md for rules
- [ ] Don't create new .md files
- [ ] Use existing files (project.md, technical.md)
- [ ] Update this file when done with tasks

---

## CURRENT SESSION SUMMARY (Oct 23, 2025)

### What Happened Today:

1. **Identified featureCounts problem** ✅
   - Oct 18 run used `-g gene_biotype` instead of `-g gene_id`
   - Result: Got biotype groups instead of individual gene counts
   - Verified in: `/output/nf-work/00/bf2f02386ad32b8f20d50a6737507b/.command.sh` line 3

2. **Fixed configuration** ✅
   - Updated `params.yaml`: `featurecounts_group_type: gene_id`
   - Using correct reference files (same ones that worked Oct 18)

3. **Cleaned up project** ✅
   - Deleted obsolete script directories
   - Cleaned 04_qc_analysis (kept only 08, 09, 10)
   - Updated .gitignore (ignore all data/references/)
   - Committed and pushed (commit f758c74)

4. **Added rdryad for automated downloads** ✅
   - Added rdryad R package to Dockerfile
   - GitHub Actions building new container
   - Test script: `scripts/01_data_acquisition/download_dryad_genome.py` gets 403
   - Solution: Use rdryad R package when container ready

5. **Pipeline re-running** 🔄
   - Job 20644937 started 18:09
   - Using -resume (reusing trimmed reads, QC)
   - Only re-running: STAR alignment, featureCounts with correct param
   - Expected completion: ~2-3 hours

### What's Next (when you return):

**Check pipeline completion:**
```bash
squeue -j 20644937
tail -50 logs/rnaseq_main_20644937.o.txt
```

**Verify featureCounts is correct:**
```bash
# Should see individual genes, not biotypes
head output/star_salmon/featurecounts/PRJNA268379_SRR1663709.featureCounts.tsv
wc -l output/star_salmon/salmon.merged.gene_counts.tsv  # Should be ~22,176
```

**Then start DESeq2 analyses** (see IMMEDIATE NEXT STEPS section above)
