# 03_sanity_check_counts.py - Validate Count Matrices

**Script:** `03_sanity_check_counts.py`
**Phase:** Phase 2 - Organize Count Matrices (validation step)
**Purpose:** Comprehensive validation of count matrices before proceeding to differential expression analysis

---

## What This Script Does

This script performs 6 critical sanity checks on the count matrices created by `01_combine_counts.py`:

1. **File Structure:** All expected files exist and are non-empty
2. **Sample Counts:** Correct number of samples per matrix
3. **Gene Counts:** Gene counts in expected range (~22,177)
4. **Missing Data:** No excessive NAs or zeros
5. **Replicate Correlations:** Biological replicates correlate highly (r > 0.9)
6. **Metadata Alignment:** Samples in counts match metadata files

**If any check fails:** Script exits with error code and detailed message
**If all checks pass:** Script exits successfully, ready for Phase 3

---

## Why We Need This

Sanity checks catch errors early:
- Missing or corrupted files
- Sample mix-ups or duplicates
- Technical failures in count generation
- Metadata mismatches
- Poor quality samples (low correlations)

Finding these issues now (before DESeq2) saves hours of debugging later.

---

## Input Files

Expects the output structure from `01_combine_counts.py`:

```
output/count_matrices/
├── featurecounts/
│   ├── adults/
│   │   ├── unfed_counts.tsv
│   │   ├── unfed_metadata.tsv
│   │   ├── fed_counts.tsv
│   │   ├── fed_metadata.tsv
│   │   ├── all_counts.tsv
│   │   └── all_metadata.tsv
│   ├── embryos/
│   │   ├── 72h_counts.tsv
│   │   ├── 72h_metadata.tsv
│   │   └── ...
│   └── larvae/
│       ├── 11d_counts.tsv
│       ├── 11d_metadata.tsv
│       └── ...
└── salmon/
    └── (same structure)
```

---

## Output

### Console Output:

```
======================================================================
SANITY CHECK: COUNT MATRICES
======================================================================

This script validates the organized count matrices.
All checks must pass before proceeding to Phase 3 (Salmon validation).

----------------------------------------------------------------------
Check 1: File Structure
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/unfed_counts.tsv
        Size: 245.3 KB
  [PASS] ✓ adults/fed_counts.tsv
        Size: 243.1 KB
  ...

----------------------------------------------------------------------
Check 2: Sample Counts
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/unfed_counts.tsv
        8 samples
  [PASS] ✓ adults/fed_counts.tsv
        8 samples
  ...

----------------------------------------------------------------------
Check 3: Gene Counts
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/unfed_counts.tsv
        22177 genes
  ...

Gene count range: 22177 - 22177
Mean: 22177, Median: 22177

----------------------------------------------------------------------
Check 4: Missing Data
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/unfed_counts.tsv
        0.000% NAs
  ...

----------------------------------------------------------------------
Check 5: Replicate Correlations
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/SD_Unfed
        Mean r = 0.982, Min r = 0.975
  [PASS] ✓ adults/LD_Unfed
        Mean r = 0.985, Min r = 0.981
  ...

----------------------------------------------------------------------
Check 6: Metadata Alignment
----------------------------------------------------------------------

FEATURECOUNTS:
  [PASS] ✓ adults/unfed_counts.tsv
        8 samples match
  ...

======================================================================
SUMMARY
======================================================================
  [PASS] ✓ File Structure
  [PASS] ✓ Sample Counts
  [PASS] ✓ Gene Counts
  [PASS] ✓ Missing Data
  [PASS] ✓ Replicate Correlations
  [PASS] ✓ Metadata Alignment

All sanity checks passed!
Ready to proceed to Phase 3 (Validate Salmon vs featureCounts)
```

### Exit Codes:
- **0:** All checks passed
- **1:** One or more checks failed

---

## How To Run

### Using Container (Required):

```bash
# Load singularity module
module load singularity

# Run script with container
singularity exec --cleanenv --bind $PWD:/proj \
  /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python scripts/06_count_matrix/03_sanity_check_counts.py"
```

### Expected Runtime:
- 1-2 minutes for all checks

---

## Detailed Check Descriptions

### Check 1: File Structure
**What it checks:**
- All expected count matrix files exist
- All expected metadata files exist
- Files are non-empty (size > 0 bytes)

**Expected files per stage:**
- Adults: 3 matrices (unfed, fed, all)
- Embryos: 5 matrices (72h, 135h, SD_timeseries, LD_timeseries, all)
- Larvae: 3 matrices (11d, SD_timeseries, all)

**Total:** 11 matrices × 2 count types (featureCounts + Salmon) = 22 count files + 22 metadata files

**Failure causes:**
- `01_combine_counts.py` did not complete successfully
- Disk space issues during write
- File permissions issues

### Check 2: Sample Counts
**What it checks:**
- Each matrix has the expected number of samples

**Expected counts:**
- unfed: 8 (4 SD + 4 LD)
- fed: 8 (4 SD + 4 LD)
- 72h: 6 (3 SD + 3 LD)
- 135h: 6 (3 SD + 3 LD)
- SD_timeseries (embryos): 6 (3 @ 72h + 3 @ 135h)
- LD_timeseries: 6 (3 @ 72h + 3 @ 135h)
- SD_timeseries (larvae): 9 (3 @ 11d + 3 @ 21d + 3 @ 40d)
- 11d: 6 (3 SD + 3 LD)
- all (adults): 16
- all (embryos): 12
- all (larvae): 17

**Failure causes:**
- Samples missing from metadata
- Sample ID mismatches
- Logic error in `01_combine_counts.py`

### Check 3: Gene Counts
**What it checks:**
- Gene counts are in expected range (20,000 - 25,000)
- Reports min, max, mean, median across all matrices

**Expected:** ~22,177 genes (from featureCounts output)

**Why variation is OK:**
- Different filtering in different tools
- Salmon aggregates transcripts to genes (may differ slightly)
- Small variations (±1000 genes) are normal

**Failure causes:**
- Wrong gene ID field extracted
- File corruption
- Major pipeline failure

### Check 4: Missing Data
**What it checks:**
- Percentage of NA values
- Percentage of zero values

**Thresholds:**
- NAs: < 1% (stricter, should be ~0%)
- Zeros: < 50% (lenient, zeros are common in RNA-seq)

**Why zeros are normal:**
- Many genes not expressed in specific samples
- Low-abundance transcripts may have 0 counts
- 30-40% zeros is typical for RNA-seq

**Failure causes:**
- File read errors
- Corrupted data
- Missing samples

### Check 5: Replicate Correlations
**What it checks:**
- Pearson correlation between biological replicates
- Uses log2(count + 1) transformation

**Threshold:** r > 0.9 (minimum for any pair)

**What it tells us:**
- r > 0.95: Excellent replicates
- r = 0.9-0.95: Good replicates
- r < 0.9: Poor replicates (investigate)

**Compared groups:**
- SD_Unfed replicates (n=4)
- LD_Unfed replicates (n=4)
- SD_Fed replicates (n=4)
- LD_Fed replicates (n=4)
- SD_72h replicates (n=3)
- LD_72h replicates (n=3)
- SD_135h replicates (n=3)
- LD_135h replicates (n=3)
- SD_11d replicates (n=3)
- LD_11d replicates (n=3)

**Failure causes:**
- Sample swap or mislabeling
- Technical failure for one replicate
- Biological variation (rare in controlled experiments)

### Check 6: Metadata Alignment
**What it checks:**
- Sample IDs in count matrix match sample IDs in metadata
- No samples in counts missing from metadata
- No samples in metadata missing from counts

**Failure causes:**
- Metadata file incomplete
- Sample ID format mismatch
- Logic error in sample filtering

---

## What To Do If Checks Fail

### General Steps:
1. Read the error message carefully
2. Check which specific file/sample failed
3. Investigate the root cause
4. Fix the issue
5. Re-run `01_combine_counts.py` if needed
6. Re-run this sanity check script

### Specific Failures:

**File Structure Failed:**
- Check if `01_combine_counts.py` completed successfully
- Look for error messages in that script's output
- Verify disk space: `df -h`
- Check file permissions

**Sample Counts Failed:**
- Check metadata file: `data/metadata/samples.csv`
- Verify sample IDs are correct
- Check for typos in sample group definitions

**Gene Counts Failed:**
- If drastically different (< 15,000 or > 30,000): major issue, investigate pipeline
- If slightly different (20,000-25,000): likely OK, document and proceed

**Missing Data Failed:**
- Check for file corruption
- Verify pipeline completed for all samples
- Review MultiQC report for sample failures

**Replicate Correlations Failed:**
- Identify which specific replicates failed
- Check if one replicate is outlier
- Review that sample's QC metrics in MultiQC
- Consider excluding problematic sample

**Metadata Alignment Failed:**
- Check for sample ID typos
- Verify metadata file is up to date
- Check if sample was renamed

---

## What Comes Next

### If All Checks Pass:
1. Update `project_state.json`:
   ```bash
   singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
     bash -c "cd /proj && python tracking/00_track.py complete sanity_check_matrices"
   ```

2. Proceed to **Phase 3: Validate Salmon vs featureCounts**

### If Checks Fail:
1. Do NOT proceed to Phase 3
2. Investigate and fix issues
3. Re-run `01_combine_counts.py` if needed
4. Re-run this sanity check until all pass

---

## References

- Analysis strategy: `tracking/analysis_strategy.md` (Phase 2.1)
- Previous script: `01_combine_counts.py`
- Next phase: Phase 3 - Validate Salmon vs featureCounts

---

**Created:** October 27, 2025
**Last Updated:** October 27, 2025
**Next Phase:** Validate Salmon vs featureCounts (Phase 3)
