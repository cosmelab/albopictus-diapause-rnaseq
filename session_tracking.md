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
- ⏳ NEED TO REBUILD container before running analyses
- 📝 User needs to update GitHub secrets to trigger build

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

### Container Rule
- NEVER install outside container
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

## Current Todo List

1. [ ] Fix pheatmap in 01_adults_deseq2_unfed.R
2. [ ] Resubmit adults unfed job
3. [ ] Create batch script for adults fed
4. [ ] Submit adults fed job
5. [ ] Create embryos script
6. [ ] Create larvae script
7. [ ] Extract GWAS candidates
8. [ ] Fisher's enrichment test
