# 01_combine_counts.py - Organize Count Matrices by Stage and Comparison

**Script:** `01_combine_counts.py`
**Phase:** Phase 2 - Organize Count Matrices
**Purpose:** Read raw featureCounts and Salmon outputs and create organized count matrices split by life stage and biological comparison

---

## What This Script Does

This script processes the nf-core/rnaseq pipeline outputs and creates organized count matrices for differential expression analysis. It:

1. Reads featureCounts gene-level counts (45 files)
2. Reads Salmon transcript-level quantifications (45 directories)
3. Aggregates Salmon transcripts to gene level
4. Splits samples by life stage (adults, embryos, larvae)
5. Further splits by biological comparison (unfed/fed, timepoints)
6. Saves organized matrices with corresponding metadata

---

## Why We Need This

The nf-core/rnaseq pipeline outputs one file per sample. For DESeq2 analysis, we need:
- Combined count matrices (genes × samples)
- Matrices organized by comparison group
- Separate matrices for each analysis (adults unfed, adults fed, embryo 72h, etc.)
- Both featureCounts (for collaborator replication) and Salmon (for reviewers)

---

## Input Files

### featureCounts (Primary - for Part 1):
```
output/star_salmon/featurecounts/
├── PRJNA268379_SRR1663685.featureCounts.tsv  # One per sample
├── PRJNA268379_SRR1663687.featureCounts.tsv
├── ...
└── (45 files total)
```

**Format:**
- Tab-separated
- Columns: Geneid, Chr, Start, End, Strand, Length, counts
- ~22,177 genes per file

### Salmon (Secondary - for Part 2 after validation):
```
output/star_salmon/
├── PRJNA268379_SRR1663685/
│   └── quant.sf  # Transcript-level counts
├── PRJNA268379_SRR1663687/
│   └── quant.sf
├── ...
└── (45 directories total)
```

**Format:**
- Tab-separated quant.sf file
- Columns: Name, Length, EffectiveLength, TPM, NumReads
- Transcript-level (will be aggregated to gene level)

### Metadata:
```
data/metadata/samples.csv
```
**Key columns:**
- BioProject: PRJNA268379, PRJNA158021, or PRJNA187045
- Run: SRA run accession (e.g., SRR1663685)
- New abreviation: Sample ID (e.g., SD_Unfed, LD_72h, SD_11d)
- Sample collection: Life stage
- Photoperiod: 8L:16D (SD) or 16L:8D (LD)

---

## Output Files

### Directory Structure:
```
output/count_matrices/
├── featurecounts/
│   ├── adults/
│   │   ├── unfed_counts.tsv        # 8 samples (SD_Unfed + LD_Unfed)
│   │   ├── unfed_metadata.tsv
│   │   ├── fed_counts.tsv          # 8 samples (SD_Fed + LD_Fed)
│   │   ├── fed_metadata.tsv
│   │   ├── all_counts.tsv          # 16 samples (all adults)
│   │   └── all_metadata.tsv
│   ├── embryos/
│   │   ├── 72h_counts.tsv          # 6 samples
│   │   ├── 72h_metadata.tsv
│   │   ├── 135h_counts.tsv         # 6 samples
│   │   ├── 135h_metadata.tsv
│   │   ├── SD_timeseries_counts.tsv  # 6 SD samples
│   │   ├── SD_timeseries_metadata.tsv
│   │   ├── LD_timeseries_counts.tsv  # 6 LD samples
│   │   ├── LD_timeseries_metadata.tsv
│   │   ├── all_counts.tsv          # 12 samples
│   │   └── all_metadata.tsv
│   └── larvae/
│       ├── 11d_counts.tsv          # 6 samples (VALID comparison)
│       ├── 11d_metadata.tsv
│       ├── SD_timeseries_counts.tsv  # 9 SD samples
│       ├── SD_timeseries_metadata.tsv
│       ├── all_counts.tsv          # 17 samples
│       └── all_metadata.tsv
└── salmon/
    └── (same structure as featurecounts)
```

### File Format:

**Count matrices (_counts.tsv):**
- Tab-separated
- Rows: genes (gene_id as index)
- Columns: samples (sorted alphabetically)
- Values: raw counts

**Metadata files (_metadata.tsv):**
- Tab-separated
- Same sample order as count matrix columns
- All columns from samples.csv

---

## How To Run

### Using Container (Required):

```bash
# Load singularity module
module load singularity

# Run script with container
singularity exec --cleanenv --bind $PWD:/proj \
  /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && python scripts/06_count_matrix/01_combine_counts.py"
```

### Expected Runtime:
- 2-5 minutes for all 45 samples (both featureCounts and Salmon)

### Expected Output:
```
Loading metadata...
Loaded metadata: 45 samples

Processing featureCounts
Found 45 featureCounts files

Processing ADULTS
  Comparison: unfed
    Creating matrix for 8 samples: LD_Unfed, SD_Unfed
    Matrix shape: (22177, 8) (genes x samples)
    Saved: output/count_matrices/featurecounts/adults/unfed_counts.tsv
    Saved: output/count_matrices/featurecounts/adults/unfed_metadata.tsv

  Comparison: fed
    Creating matrix for 8 samples: LD_Fed, SD_Fed
    ...

(continues for all comparisons)

SUMMARY
featureCounts matrices created: 10
featureCounts metadata files: 10

Salmon matrices created: 10
Salmon metadata files: 10

Next step: Run 03_sanity_check_counts.py to validate outputs

Done!
```

---

## Sample Group Definitions

Defined in SAMPLE_GROUPS dictionary in script:

### Adults (PRJNA268379):
- **unfed:** SD_Unfed (n=4), LD_Unfed (n=4)
- **fed:** SD_Fed (n=4), LD_Fed (n=4)
- **all:** All 16 adult samples

### Embryos (PRJNA158021):
- **72h:** SD_72h (n=3), LD_72h (n=3)
- **135h:** SD_135h (n=3), LD_135h (n=3)
- **SD_timeseries:** SD_72h (n=3), SD_135h (n=3)
- **LD_timeseries:** LD_72h (n=3), LD_135h (n=3)
- **all:** All 12 embryo samples

### Larvae (PRJNA187045):
- **11d:** SD_11d (n=3), LD_11d (n=3) - VALID comparison only
- **SD_timeseries:** SD_11d (n=3), SD_21d (n=3), SD_40d (n=3)
- **all:** All 17 larvae samples

---

## Salmon Gene-Level Aggregation

Salmon outputs transcript-level counts. This script aggregates to gene level by:

1. Extracting gene_id from transcript name:
   - Pattern: `gene-LOC109400916-RA` → gene_id: `LOC109400916`
   - Uses regex: `r'(LOC\d+)'`

2. Summing NumReads across all transcripts per gene

3. Matching gene_ids with featureCounts gene_ids

**Note:** This is a simple aggregation. For more sophisticated methods (accounting for transcript length, etc.), use tximport in R. However, for validation purposes (Phase 3), this simple sum is sufficient.

---

## Sanity Checks Built-In

The script performs these checks automatically:

1. Metadata file exists and has required columns
2. featureCounts/Salmon files found (reports count)
3. Sample IDs in metadata match files
4. Warns if count files missing for any sample
5. Reports matrix dimensions (genes × samples)
6. Verifies output files created

---

## What Comes Next

After running this script:

1. **Run sanity checks:** `03_sanity_check_counts.py`
   - Verify sample counts correct
   - Check for missing data
   - Calculate replicate correlations
   - Verify gene counts ~22,177

2. **If sanity checks pass:** Proceed to Phase 3 (Validate Salmon vs featureCounts)

3. **If sanity checks fail:** Investigate and fix before proceeding

---

## Troubleshooting

### Error: "Metadata file not found"
**Solution:** Ensure you're running from project root directory where data/metadata/samples.csv exists

### Error: "No featureCounts files found"
**Solution:** Check that nf-core/rnaseq pipeline completed successfully. Files should be in output/star_salmon/featurecounts/

### Error: "No Salmon directories found"
**Solution:** Check output/star_salmon/ for directories like PRJNA268379_SRR1663685/

### Warning: "Count file not found for [sample]"
**Solution:** Check if that specific sample failed in the pipeline. Review pipeline logs.

### Matrix has fewer genes than expected
**Solution:** This is normal - gene counts may vary slightly between samples due to filtering. Sanity check script will verify this is acceptable.

---

## References

- Experimental design: `tracking/experimental_design.md`
- Analysis strategy: `tracking/analysis_strategy.md` (Phase 2.1)
- Sample metadata: `data/metadata/samples.csv`
- Pipeline outputs: `output/star_salmon/`

---

**Created:** October 27, 2025
**Last Updated:** October 27, 2025
**Next Script:** `03_sanity_check_counts.py`
