# Session Tracking - Current Work

**Last Updated:** Oct 22, 2025 16:30

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

### PIPELINE RE-RUN STATUS (Oct 23, 2025)

**THE ACTUAL PROBLEM (VERIFIED):**
- Oct 18 run used: `-g gene_biotype` in featureCounts
- This groups ALL genes by biotype instead of counting individual genes
- Result: ~10 biotype categories instead of 22,176 individual gene counts
- **PROOF:** `/bigdata/.../output/nf-work/00/bf2f02386ad32b8f20d50a6737507b/.command.sh` line 3

**Job Running:** 20644924 (started 18:05)
- **FIXED PARAMETER:** `featurecounts_group_type: gene_id` in params.yaml
- **References:** Using same files that worked Oct 18
  - FASTA: data/references/AalbF3/AalbF3_genome.fa.gz
  - GTF: data/references/AalbF3/AalbF3_annotation.gtf.gz
  - rRNA: data/references/AalbF3/combined_rRNA_sequences.fa
- **STATUS:** Pipeline starting with -resume (will reuse trimmed reads and alignments)
- **Expected completion:** ~48-72 hours

### IMMEDIATE NEXT STEPS (Step 1: Replicate Collaborators)

1. ✅ Document collaborator methods (DONE - all 3 experiments)
2. ✅ Update Dockerfile with all collaborator packages (DONE)
3. ✅ Update adult scripts to match Angela's methods (DONE)
4. ⏳ **NEXT:** Rebuild container (user needs to update GitHub secrets)
5. Rerun adult analyses with updated container
6. Create embryo scripts (72h, 135h separately - matching Mackenzie's methods)
7. Create larvae scripts (11d, 21d, 40d separately - matching Sarah's methods)
8. *Optional:* Check if collaborators deposited their counts (GEO/SRA/Dryad)
9. Compare patterns (may not need exact counts comparison)

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

**02_annotation_prep/** ❌ DELETE - duplicate of 02_annotation_preprocessing/

**02_annotation_preprocessing/** ✅ DONE (2 scripts)
- Converted GFF3 to GTF, fixed transcript issues

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

### Cleanup Commands (run after pipeline completes):
```bash
rm -rf scripts/01_data_acquisition/
rm -rf scripts/02_annotation_prep/
rm -rf scripts/05_batch_correction/
rm -rf scripts/06_diff_expression/
mv scripts/05_count_matrix scripts/04b_count_matrix
# Clean 04_qc_analysis - keep only 08, 09, 10
```

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

## CRITICAL ISSUE FOUND (Oct 23, 2025)

### THE PROBLEM:
1. **Wrong GTF used**: Only protein-coding genes (15,542), missing:
   - 153 rRNA genes (can't measure contamination!)
   - 3,894 lncRNA genes (8 are GWAS candidates!)
   - 869 tRNA genes

2. **Wrong featureCounts parameter**: Used `-g gene_biotype` instead of `-g gene_id`
   - Result: Collapsed all genes into ~10 biotype categories
   - We need: Individual counts for 20,621 genes

### WHAT WAS FIXED:
1. ✅ Created complete GTF: `/data/references/AalbF3/AalbF3_all_gene_types.gtf`
2. ✅ Updated `scripts/03_rnaseq_pipeline/params.yaml`:
   - Points to complete GTF
   - Added `featurecounts_group_type: gene_id`
3. ✅ Cleaned up references (one directory: `/data/references/AalbF3/`)

### IMMEDIATE TODO (IN ORDER):

1. [ ] **Re-run nf-core pipeline with fixed params**
   ```bash
   cd /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq
   sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
   ```
   - Will take ~72 hours
   - Will generate proper counts for ALL 20,621 genes

2. [ ] **After pipeline completes, verify outputs:**
   - Check featureCounts has gene-level counts (not biotypes)
   - Check Salmon counts include all genes
   - Verify rRNA contamination metrics are calculated

3. [ ] **Create count matrices for DESeq2:**
   - Split by stage (adults, embryos, larvae)
   - Ensure integer counts
   - Match sample names to metadata

4. [ ] **Replicate collaborator analyses (Step 1)**
   - Adults: Split by blood meal status
   - Embryos: Split by timepoint (72h, 135h)
   - Larvae: Split by timepoint (11d, 21d, 40d)
   - Use their exact DESeq2 parameters

5. [ ] **Run proper statistical models (Step 2)**
   - Account for all variables
   - Compare HTSeq vs Salmon

6. [ ] **Address reviewer concerns (Step 3)**
