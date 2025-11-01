# 07_qc_visualization - Quality Control Before Differential Expression

**Phase:** Phase 2.5 - QC Visualization (CRITICAL GATE CHECK)
**Purpose:** Comprehensive visualization to verify data quality BEFORE proceeding to DE analysis

---

## Why This Phase is Essential

### Literature Consensus (DESeq2/DEGreport/HBC Training)
**PCA and correlation heatmaps are MANDATORY before differential expression analysis.**

**What QC visualization catches:**
- Sample outliers (technical failures)
- Sample swaps or mislabeling
- Batch effects
- Unexpected sources of variation
- Confir that experimental design is visible in the data

**Early detection saves time:** Finding issues now prevents hours of debugging DE results later.

### nf-core/rnaseq Best Practice
The pipeline already provides:
- Global PCA (all 45 samples) in MultiQC report
- FastQC, alignment stats, duplication analysis

**What we add:**
- **Stage-specific PCA** (within adults, embryos, larvae)
- **Biological context interpretation** (photoperiod, time, blood meal)
- **Outlier detection** specific to experimental groups
- **Integration** of MultiQC metrics with expression patterns

---

## Scripts in This Directory

### 01_global_qc_plots.R
**Purpose:** QC visualization across all 45 samples

**Input:**
- `output/count_matrices/featurecounts/adults/all_counts.tsv` (and embryos, larvae)
- `output/count_matrices/salmon/adults/all_counts.tsv` (and embryos, larvae)

**Analyses:**
1. **PCA (all 45 samples):** Colored by life stage
   - Expected: Strong clustering by stage (due to platform confounding)
   - Verify: No extreme outliers across all samples

2. **Correlation heatmap:** Hierarchical clustering
   - Expected: Blocks corresponding to life stages
   - Check: High within-stage correlations

3. **Library size boxplots:** Total counts per sample
   - Check: No extreme differences (>3 SD from mean)

4. **Detection rate:** Percent genes detected per sample
   - Check: Consistent detection rates within stages

**Output:**
```
results/qc_visualization/featurecounts/all_samples/
├── pca_all_45_samples.pdf
├── correlation_heatmap_all.pdf
├── sample_distance_heatmap.pdf
└── library_size_boxplot.pdf
```

**Runtime:** 2-3 minutes

---

### 02_adults_qc_plots.R
**Purpose:** QC for adults (16 samples) split by blood meal status

**Input:**
- `output/count_matrices/featurecounts/adults/unfed_counts.tsv` (8 samples)
- `output/count_matrices/featurecounts/adults/fed_counts.tsv` (8 samples)
- `output/count_matrices/featurecounts/adults/all_counts.tsv` (16 samples)

**Analyses:**

**1. Global adults PCA (all 16 samples):**
- Color by photoperiod (SD vs LD)
- Shape by blood meal status (Unfed vs Fed)
- Expected: Separation by both factors

**2. Unfed-specific PCA (8 samples):**
- SD_Unfed vs LD_Unfed
- Expected: Clear separation by photoperiod
- Check: Biological replicates cluster together

**3. Fed-specific PCA (8 samples):**
- SD_Fed vs LD_Fed (PRIMARY COMPARISON)
- Expected: Strong separation (diapause vs non-diapause oocytes)
- Check: Biological replicates cluster together

**4. Correlation heatmaps:**
- Within-group correlations (expect r > 0.9)
- Between-group differences

**5. Normalized count boxplots:**
- Distribution by treatment group
- Check for consistent normalization

**Expected results:**
- ✓ Replicates cluster together
- ✓ SD vs LD separation visible
- ✓ Blood meal effect visible
- ✓ No extreme outliers

**Output:**
```
results/qc_visualization/featurecounts/adults/
├── pca_all_adults.pdf
├── pca_by_photoperiod.pdf
├── pca_by_blood_meal.pdf
├── pca_unfed_only.pdf
├── pca_fed_only.pdf
├── correlation_heatmap_adults.pdf
├── correlation_heatmap_unfed.pdf
├── correlation_heatmap_fed.pdf
└── boxplot_by_group.pdf
```

**Runtime:** 3-4 minutes

---

### 03_embryos_qc_plots.R
**Purpose:** QC for embryos (12 samples) with time series analysis

**Input:**
- `output/count_matrices/featurecounts/embryos/72h_counts.tsv` (6 samples)
- `output/count_matrices/featurecounts/embryos/135h_counts.tsv` (6 samples)
- `output/count_matrices/featurecounts/embryos/all_counts.tsv` (12 samples)

**Analyses:**

**1. Global embryos PCA (all 12 samples):**
- Color by photoperiod (SD vs LD)
- Shape by timepoint (72h vs 135h)
- Expected: Separation by photoperiod, possible time effect

**2. 72h-specific PCA (6 samples):**
- SD_72h vs LD_72h
- Expected: Early diapause programming signature

**3. 135h-specific PCA (6 samples):**
- SD_135h vs LD_135h
- Expected: Diapause entry signature

**4. Interaction plot:**
- PC1 vs PC2, colored by photoperiod, shaped by time
- Check: Does photoperiod effect INCREASE over time?

**5. Time series trajectory:**
- Connect 72h → 135h within each photoperiod
- Visualize developmental progression

**6. Correlation heatmaps:**
- Within-timepoint correlations
- Check: Consistent replicates

**Expected results:**
- ✓ Photoperiod separation at both timepoints
- ✓ Possible interaction (photoperiod effect changes over time)
- ✓ Smooth developmental trajectory (not abrupt)

**Output:**
```
results/qc_visualization/featurecounts/embryos/
├── pca_all_embryos.pdf
├── pca_by_photoperiod.pdf
├── pca_by_timepoint.pdf
├── pca_interaction.pdf
├── pca_72h_only.pdf
├── pca_135h_only.pdf
├── time_series_trajectory.pdf
├── correlation_heatmap_embryos.pdf
├── correlation_heatmap_72h.pdf
└── correlation_heatmap_135h.pdf
```

**Runtime:** 3-4 minutes

---

### 04_larvae_qc_plots.R
**Purpose:** QC for larvae (17 samples) with degradation warning

**Input:**
- `output/count_matrices/featurecounts/larvae/11d_counts.tsv` (6 samples)
- `output/count_matrices/featurecounts/larvae/SD_timeseries_counts.tsv` (9 samples)
- `output/count_matrices/featurecounts/larvae/all_counts.tsv` (17 samples)

**Analyses:**

**1. Valid comparison PCA (11d only, 6 samples):**
- SD_11d vs LD_11d
- Expected: Clear separation
- This is the ONLY valid cross-photoperiod comparison

**2. SD time series PCA (9 samples):**
- SD_11d, SD_21d, SD_40d
- Expected: Gradual progression (diapause deepening)
- Check: No abrupt jumps

**3. LD degradation warning plot (all LD samples):**
- LD_11d, LD_21d, LD_40d
- Expected: LD_21d and LD_40d are outliers
- **WARNING ANNOTATION:** "LD samples at 21d/40d are degraded"

**4. Global larvae PCA (all 17 samples):**
- Visualize the degradation problem
- Color code: valid vs degraded

**5. Correlation heatmaps:**
- Check: LD_21d/40d have low correlation with LD_11d

**Expected results:**
- ✓ 11d comparison shows SD vs LD separation
- ✓ SD time series shows gradual changes
- ✗ LD 21d/40d are visible outliers (expected!)

**Output:**
```
results/qc_visualization/featurecounts/larvae/
├── pca_11d_valid_comparison.pdf
├── pca_SD_timeseries.pdf
├── pca_LD_degradation_warning.pdf
├── pca_all_larvae.pdf
├── correlation_heatmap_11d.pdf
├── correlation_heatmap_SD_timeseries.pdf
└── correlation_heatmap_larvae_with_warning.pdf
```

**Runtime:** 3-4 minutes

---

### 05_extract_multiqc_metrics.py
**Purpose:** Parse MultiQC report and integrate with expression QC

**Input:**
- `output/multiqc/star_salmon/multiqc_data/multiqc_data.json`

**Extracts:**
1. **Mapping statistics:** Total reads, uniquely mapped, multi-mapped
2. **Duplication rates:** From dupRadar
3. **RSeQC metrics:** Read distribution, junction annotation
4. **Library size:** Total counts per sample

**Combines with expression QC:**
- Correlate mapping rate with PC1/PC2
- Check if technical metrics explain variation
- Identify samples where technical issues drive outlier status

**Output:**
```
results/qc_visualization/multiqc_integration/
├── technical_metrics_summary.tsv
├── technical_vs_pca.pdf
├── mapping_rate_vs_expression.pdf
└── duplication_vs_counts.pdf
```

**Runtime:** 1-2 minutes

---

### 06_outlier_detection.R
**Purpose:** Automated flagging of problematic samples

**Criteria for flagging:**
1. **Low correlation:** r < 0.80 with other replicates
2. **Extreme library size:** >3 SD from group mean
3. **Low detection rate:** <60% of typical for that stage
4. **PCA outlier:** >3 SD from group centroid on PC1 or PC2
5. **Technical issues:** Very low mapping rate, extreme duplication

**For each flagged sample:**
- Sample ID
- Flagging reason(s)
- QC metrics
- Recommended action (review, possibly exclude)

**Output:**
```
results/qc_visualization/outlier_report.txt
results/qc_visualization/outlier_summary.tsv
```

**Example report:**
```
OUTLIER DETECTION REPORT
========================

Sample: PRJNA268379_SRR1663690 (Adults, SD_Fed replicate 3)
Flags:
  - Low correlation with other SD_Fed replicates (r = 0.78)
  - Extreme library size (6.2M reads, group mean = 23.4M)
  - Low mapping rate (65%, group mean = 88%)
Recommendation: EXCLUDE - Likely technical failure
```

**Runtime:** 2-3 minutes

---

## Workflow

### Step 1: Create all scripts
Run this command to check scripts exist:
```bash
ls -1 scripts/07_qc_visualization/*.{R,py}
```

Expected:
```
01_global_qc_plots.R
02_adults_qc_plots.R
03_embryos_qc_plots.R
04_larvae_qc_plots.R
05_extract_multiqc_metrics.py
06_outlier_detection.R
```

### Step 2: Run QC for featureCounts
```bash
module load singularity

# Run all QC scripts in order
for script in scripts/07_qc_visualization/0{1..6}*.{R,py}; do
  if [[ $script == *.R ]]; then
    echo "Running $script..."
    singularity exec --cleanenv --bind $PWD:/proj \
      albopictus-diapause-rnaseq.sif \
      bash -c "cd /proj && Rscript $script featurecounts"
  elif [[ $script == *.py ]]; then
    echo "Running $script..."
    singularity exec --cleanenv --bind $PWD:/proj \
      albopictus-diapause-rnaseq.sif \
      bash -c "cd /proj && python $script"
  fi
done
```

**Expected runtime:** 15-20 minutes total

### Step 3: Review QC results
**CRITICAL STEP - Manual review required**

Review each output directory:
```bash
# View all PDF plots
ls results/qc_visualization/featurecounts/*/*.pdf

# Read outlier report
cat results/qc_visualization/outlier_report.txt
```

**Questions to answer:**
1. Do biological replicates cluster together?
2. Is experimental design visible in PCA?
3. Are there any extreme outliers?
4. Do technical metrics explain variation?
5. Are LD larvae at 21d/40d outliers (expected)?

### Step 4: Document decisions
Create `results/qc_visualization/qc_decisions.md`:
```markdown
# QC Review Decisions

**Date:** 2025-10-27
**Reviewer:** [Your name]

## Overall Assessment
- [ ] All replicates cluster together
- [ ] Experimental factors visible in PCA
- [ ] No unexpected outliers

## Flagged Samples
| Sample | Issue | Decision | Rationale |
|--------|-------|----------|-----------|
| SRR... | Low r | EXCLUDE  | Technical failure |
| SRR... | ...   | KEEP     | Within acceptable range |

## Action Items
- [ ] Re-run count matrices excluding flagged samples (if any)
- [ ] Update metadata to mark excluded samples
- [ ] Proceed to Phase 3 (Salmon validation)
```

### Step 5: Run QC for Salmon (same process)
```bash
for script in scripts/07_qc_visualization/0{1..6}*.R; do
  echo "Running $script for Salmon..."
  singularity exec --cleanenv --bind $PWD:/proj \
    albopictus-diapause-rnaseq.sif \
    bash -c "cd /proj && Rscript $script salmon"
done
```

### Step 6: GATE CHECK
**DO NOT PROCEED to differential expression until:**
- [x] All QC plots reviewed
- [x] Outliers identified and documented
- [x] Biological separation confirmed
- [x] qc_decisions.md created
- [x] Any excluded samples documented

---

## Data Transformations Used

### For PCA
**Method:** Variance stabilizing transformation (VST) from DESeq2
```r
dds <- DESeqDataSetFromMatrix(counts, metadata, ~1)
vst_data <- vst(dds, blind=TRUE)
```

**Why VST:**
- Stabilizes variance across expression levels
- Recommended by DESeq2 for visualization
- Better than log2(count+1) for low-count genes

### For Correlation Heatmaps
**Method:** Pearson correlation on log2(count+1)
```r
log_counts <- log2(counts + 1)
cor_matrix <- cor(log_counts, method="pearson")
```

**Why log2:**
- Simple, interpretable
- Standard for correlation analysis
- Consistent with sanity check correlations

---

## Interpreting QC Plots

### "Good" QC Looks Like:
**PCA:**
- Biological replicates form tight clusters
- Experimental groups separate along PC1 or PC2
- PC1 explains majority of variance (>40%)

**Correlation heatmap:**
- Within-group correlations > 0.9
- Clear blocks in hierarchical clustering
- Dendogram groups by experimental design

**Library sizes:**
- Relatively consistent within stage
- No samples with <5M or >100M counts
- Coefficient of variation <50%

### "Bad" QC Looks Like:
**PCA:**
- Replicates scattered (not clustered)
- One replicate far from group
- Experimental design not visible

**Correlation heatmap:**
- Within-group correlations < 0.80
- No clear clustering pattern
- One sample doesn't group with replicates

**Library sizes:**
- Extreme outliers (>3 SD from mean)
- Bimodal distribution

---

## Troubleshooting

### Issue: Replicates don't cluster together
**Possible causes:**
- High biological variation
- Sample swap
- Technical failure for one replicate

**Actions:**
1. Check MultiQC for that sample
2. Review library size and mapping rate
3. Check sample metadata (correct labeling?)
4. Consider excluding outlier replicate

### Issue: Experimental design not visible
**Possible causes:**
- Treatment effect is subtle
- Other factors dominate (batch, time, etc.)
- Insufficient power (n too small)

**Actions:**
1. Check if effect is visible in stage-specific PCA
2. Review correlation heatmap (subtle separation?)
3. May still detect DEGs even if PCA overlap

### Issue: Many samples flagged as outliers
**Possible causes:**
- Thresholds too strict
- Systematic technical issue
- Natural biological variation

**Actions:**
1. Review flagging thresholds
2. Check MultiQC for batch effects
3. Consult with experimental biologist

---

## References

- **DESeq2 vignette:** http://bioconductor.org/packages/DESeq2
- **HBC Training:** https://hbctraining.github.io/DGE_workshop/
- **DEGreport package:** https://bioconductor.org/packages/DEGreport
- **nf-core/rnaseq QC:** https://nf-co.re/rnaseq/output

---

**Created:** October 27, 2025
**Phase:** 2.5 - QC Visualization
**Next Phase:** 3 - Validate Salmon vs featureCounts (only after QC gate check passed)
