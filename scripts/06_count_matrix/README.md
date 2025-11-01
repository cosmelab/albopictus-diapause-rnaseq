# 06_count_matrix - Organize and Validate Count Matrices

**Phase:** Phase 2 - Organize Count Matrices
**Purpose:** Create organized, validated count matrices ready for differential expression analysis

---

## Overview

This directory contains scripts to process raw pipeline outputs into organized count matrices split by life stage and biological comparison. Both featureCounts (for collaborator replication) and Salmon (for reviewer requirements) are processed.

---

## Scripts in This Directory

### 01_combine_counts.py (Primary Script)
**Purpose:** Read pipeline outputs and create organized count matrices

**Input:**
- featureCounts: `output/star_salmon/featurecounts/*.featureCounts.tsv` (45 files)
- Salmon: `output/star_salmon/*/quant.sf` (45 directories)
- Metadata: `data/metadata/samples.csv`

**Output:**
- Organized matrices in `output/count_matrices/featurecounts/` and `output/count_matrices/salmon/`
- Split by stage (adults/embryos/larvae) and comparison (unfed/fed, timepoints)
- Metadata files for each comparison group

**Documentation:** See `README_01_combine_counts.md`

### 02_create_summary_table.py (Deprecated)
**Status:** Will be replaced with new stage-specific scripts
**Note:** Do not use this old script

### 03_sanity_check_counts.py (Validation Script)
**Purpose:** Comprehensive validation of count matrices

**Checks performed:**
1. File structure (all files exist, non-empty)
2. Sample counts (correct number per matrix)
3. Gene counts (~22,177 expected)
4. Missing data (< 1% NAs)
5. Replicate correlations (r > 0.9)
6. Metadata alignment (samples match)

**Exit code:**
- 0 = All checks passed, proceed to Phase 3
- 1 = Checks failed, fix issues before proceeding

**Documentation:** See `README_03_sanity_check.md`

---

## Workflow

### Step 1: Run combine counts script

```bash
module load singularity

singularity exec --cleanenv --bind $PWD:/proj \
  albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python scripts/06_count_matrix/01_combine_counts.py"
```

**Expected output:**
```
output/count_matrices/
├── featurecounts/
│   ├── adults/
│   │   ├── unfed_counts.tsv (8 samples)
│   │   ├── fed_counts.tsv (8 samples)
│   │   └── all_counts.tsv (16 samples)
│   ├── embryos/
│   │   ├── 72h_counts.tsv (6 samples)
│   │   ├── 135h_counts.tsv (6 samples)
│   │   ├── SD_timeseries_counts.tsv (6 samples)
│   │   ├── LD_timeseries_counts.tsv (6 samples)
│   │   └── all_counts.tsv (12 samples)
│   └── larvae/
│       ├── 11d_counts.tsv (6 samples)
│       ├── SD_timeseries_counts.tsv (9 samples)
│       └── all_counts.tsv (17 samples)
└── salmon/
    └── (same structure as featurecounts)
```

**Runtime:** 2-5 minutes

### Step 2: Run sanity checks

```bash
singularity exec --cleanenv --bind $PWD:/proj \
  albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python scripts/06_count_matrix/03_sanity_check_counts.py"
```

**Expected output:**
- All 6 checks pass
- Exit code 0

**Runtime:** 1-2 minutes

### Step 3: Update tracking system

```bash
# Mark tasks complete
singularity exec --cleanenv --bind $PWD:/proj \
  albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py complete combine_featurecounts"

singularity exec --cleanenv --bind $PWD:/proj \
  albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py complete combine_salmon"

singularity exec --cleanenv --bind $PWD:/proj \
  albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python tracking/00_track.py complete sanity_check_matrices"
```

### Step 4: Proceed to Phase 3

If all sanity checks passed, you're ready for:
**Phase 3: Validate Salmon vs featureCounts**

---

## Output Structure Details

### Count Matrix Files (_counts.tsv)
- Format: Tab-separated
- Rows: Genes (gene_id as index)
- Columns: Samples (sorted alphabetically)
- Values: Raw counts (integers)

Example:
```
gene_id	LD_Fed	LD_Fed	LD_Fed	LD_Fed	SD_Fed	SD_Fed	SD_Fed	SD_Fed
LOC109400916	99	105	98	102	85	89	91	87
LOC109400917	1234	1198	1256	1245	987	1002	995	1010
...
```

### Metadata Files (_metadata.tsv)
- Format: Tab-separated
- Same sample order as count matrix columns
- All columns from `data/metadata/samples.csv`

---

## Sample Group Definitions

### Adults (PRJNA268379) - 16 samples
- **unfed:** SD_Unfed (4), LD_Unfed (4) - No blood meal
- **fed:** SD_Fed (4), LD_Fed (4) - Blood fed (PRIMARY comparison)
- **all:** All 16 adult samples

### Embryos (PRJNA158021) - 12 samples
- **72h:** SD_72h (3), LD_72h (3) - Early development
- **135h:** SD_135h (3), LD_135h (3) - Diapause entry
- **SD_timeseries:** SD_72h (3), SD_135h (3) - Diapause programming
- **LD_timeseries:** LD_72h (3), LD_135h (3) - Normal development
- **all:** All 12 embryo samples

### Larvae (PRJNA187045) - 17 samples
- **11d:** SD_11d (3), LD_11d (3) - VALID comparison (before LD degradation)
- **SD_timeseries:** SD_11d (3), SD_21d (3), SD_40d (3) - Diapause depth
- **all:** All 17 larvae samples (includes degraded LD samples at 21d/40d)

---

## Troubleshooting

### Script fails with "Metadata file not found"
**Solution:** Run from project root where `data/metadata/samples.csv` exists

### Script fails with "No featureCounts files found"
**Solution:** Check pipeline completed: `ls output/star_salmon/featurecounts/*.tsv`

### Sanity check fails: low replicate correlation
**Possible causes:**
- Sample swap/mislabeling
- Technical failure for one replicate
- Check MultiQC report for that sample

**Action:**
1. Identify problematic sample
2. Review QC metrics
3. Consider excluding sample (document decision)

### Sanity check fails: wrong sample counts
**Possible causes:**
- Sample ID mismatch in metadata
- Missing samples in pipeline output

**Action:**
1. Check metadata: `cat data/metadata/samples.csv | grep <sample_id>`
2. Check pipeline output exists for that sample
3. Verify sample ID format matches

---

## Next Steps

After completing this phase:

1. **Phase 3:** Validate Salmon vs featureCounts
   - Calculate correlations (r > 0.95 expected)
   - Run same DESeq2 with both methods
   - Verify DEG overlap > 90%

2. **Phase 4-6:** DESeq2 analysis with featureCounts
   - Adults (split by blood meal)
   - Embryos (split by timepoint)
   - Larvae (11d comparison + SD time series)

3. **Phase 7-9:** Advanced analyses
   - Clustering
   - Time series
   - Meta-analysis

4. **Phase 10:** Reviewer requirements (Salmon - only after validation)

---

## References

- Analysis strategy: `tracking/analysis_strategy.md` (Phase 2.1)
- Experimental design: `tracking/experimental_design.md`
- Session tracking: `tracking/session_tracking.md`
- Quick start: `tracking/START_HERE.md`

---

**Created:** October 27, 2025
**Last Updated:** October 27, 2025
**Phase:** 2 - Organize Count Matrices
**Next Phase:** 3 - Validate Salmon vs featureCounts
