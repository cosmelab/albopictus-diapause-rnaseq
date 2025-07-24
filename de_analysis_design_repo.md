# Differential Expression Analysis Design and Models

## Datasets

### Dataset 1: PRJNA268379 (Adult Females - Huang et al. 2015)
- **Design**: 2×2 Factorial (Photoperiod × Blood meal)
- **Conditions**: 
  - Photoperiod: Long-day (LD, 16L:8D) vs Short-day (SD, 8L:16D)
  - Blood meal: Non-blood-fed (NBF) vs Blood-fed (BF)
- **Replicates**: 4 biological replicates per condition (16 total libraries)
- **Platform**: HiSeq2000

### Dataset 2: PRJNA158021 (Embryos - Poelchau et al. 2013a)
- **Design**: Time course with diapause comparison
- **Timepoints**: 72-78h and 135-141h post-oviposition  
- **Conditions**: Diapause (D) vs Non-diapause (ND)
- **Focus**: Diapause preparation during embryonic development
- **Platform**: GAIIx

### Dataset 3: PRJNA187045 (Pharate Larvae - Poelchau et al. 2013b)
- **Design**: Time course comparison
- **Timepoints**: 11, 21, and 40 days post-oviposition
- **Conditions**: Diapause (D) vs Quiescence (Q)  
- **Focus**: Diapause maintenance vs general dormancy
- **Platform**: HiSeq2000

## Sample Metadata
| Dataset | Life_Stage | Condition | Timepoint | Photoperiod | Blood_meal | Comparison_Group |
|---------|------------|-----------|-----------|-------------|------------|------------------|
| PRJNA268379 | Adult | NBF_LD | Adult | Long_day | Non_fed | Non_diapause |
| PRJNA268379 | Adult | NBF_SD | Adult | Short_day | Non_fed | Diapause_induction |
| PRJNA268379 | Adult | BF_LD | Adult | Long_day | Blood_fed | Non_diapause |
| PRJNA268379 | Adult | BF_SD | Adult | Short_day | Blood_fed | Diapause_induction |
| PRJNA158021 | Embryo | D_72h | 72-78h_pov | Short_day | NA | Diapause_prep |
| PRJNA158021 | Embryo | ND_72h | 72-78h_pov | Long_day | NA | Non_diapause |
| PRJNA158021 | Embryo | D_141h | 135-141h_pov | Short_day | NA | Diapause_prep |
| PRJNA158021 | Embryo | ND_141h | 135-141h_pov | Long_day | NA | Non_diapause |
| PRJNA187045 | Pharate_larvae | D_11d | 11d_pov | Short_day | NA | Diapause_maintenance |
| PRJNA187045 | Pharate_larvae | Q_11d | 11d_pov | Long_day | NA | Quiescence |
| PRJNA187045 | Pharate_larvae | D_21d | 21d_pov | Short_day | NA | Diapause_maintenance |
| PRJNA187045 | Pharate_larvae | Q_21d | 21d_pov | Long_day | NA | Quiescence |
| PRJNA187045 | Pharate_larvae | D_40d | 40d_pov | Short_day | NA | Diapause_maintenance |
| PRJNA187045 | Pharate_larvae | Q_40d | 40d_pov | Long_day | NA | Quiescence |

## DESeq2 Analysis Models

### Primary Comparisons for GWAS Candidate Genes

#### 1. Adult Females (PRJNA268379)
```r
# Model 1: Diapause induction (main effect of photoperiod)
dds <- DESeqDataSetFromTximport(txi, 
                               colData = metadata, 
                               design = ~ blood_meal + photoperiod)

# Contrast: Short_day vs Long_day (diapause vs non-diapause)
results_adults <- results(dds, contrast = c("photoperiod", "Short_day", "Long_day"))

# Model 2: Blood meal interaction  
dds_interaction <- DESeqDataSetFromTximport(txi, 
                                          colData = metadata, 
                                          design = ~ photoperiod * blood_meal)

# Test interaction: does blood meal modify diapause response?
results_interaction <- results(dds_interaction, name = "photoperiodShort_day.blood_mealBlood_fed")
```

#### 2. Embryos (PRJNA158021)  
```r
# Model: Diapause preparation
dds_embryo <- DESeqDataSetFromTximport(txi, 
                                      colData = metadata, 
                                      design = ~ timepoint + condition)

# Contrast: Diapause vs Non-diapause
results_embryo <- results(dds_embryo, contrast = c("condition", "Diapause", "Non_diapause"))

# Time-specific contrasts
results_72h <- results(dds_embryo, contrast = list("condition_Diapause_vs_Non_diapause", 
                                                  "timepoint72h.conditionDiapause"))
```

#### 3. Pharate Larvae (PRJNA187045)
```r
# Model: Diapause maintenance vs quiescence
dds_pharate <- DESeqDataSetFromTximport(txi, 
                                       colData = metadata, 
                                       design = ~ timepoint + condition)

# Contrast: Diapause vs Quiescence  
results_pharate <- results(dds_pharate, contrast = c("condition", "Diapause", "Quiescence"))

# Time course analysis
dds_timecourse <- DESeqDataSetFromTximport(txi, 
                                          colData = metadata, 
                                          design = ~ condition * timepoint)
```

### Cross-Platform Comparison Strategy

#### Step 1: Individual Dataset Analysis
```r
# Analyze each dataset separately with appropriate models
# Extract results for the 34 candidate genes
candidate_genes <- c("LOC109397825", "LOC109405370", "LOC109398973", 
                    "LOC109397812", "LOC109405365", ...) # full candidate gene list

# Function to extract candidate gene results
extract_candidates <- function(results_obj, gene_list) {
  results_df <- as.data.frame(results_obj)
  results_df$gene_id <- rownames(results_df)
  return(results_df[results_df$gene_id %in% gene_list, ])
}
```

#### Step 2: Cross-Platform Normalization
```r
# Combine TPM data across datasets
# Account for platform differences using batch correction
library(sva)

# Combat batch correction for cross-platform comparison
tpm_combined <- cbind(tpm_adults, tpm_embryos, tpm_pharate)
metadata_combined <- rbind(meta_adults, meta_embryos, meta_pharate)

# Add batch variable for platform
metadata_combined$platform <- c(rep("HiSeq2000", nrow(meta_adults)),
                                rep("GAIIx", nrow(meta_embryos)), 
                                rep("HiSeq2000", nrow(meta_pharate)))

# Batch correction
tpm_corrected <- ComBat_seq(tpm_combined, 
                           batch = metadata_combined$platform,
                           group = metadata_combined$diapause_status)
```

#### Step 3: Meta-Analysis Approach
```r
# Combine effect sizes across life stages
library(metafor)

# For each gene, combine log2FC across studies
combine_effects <- function(gene_id) {
  # Extract log2FC and SE from each dataset
  adult_effect <- results_adults[gene_id, "log2FoldChange"]
  adult_se <- results_adults[gene_id, "lfcSE"]
  
  embryo_effect <- results_embryo[gene_id, "log2FoldChange"] 
  embryo_se <- results_embryo[gene_id, "lfcSE"]
  
  pharate_effect <- results_pharate[gene_id, "log2FoldChange"]
  pharate_se <- results_pharate[gene_id, "lfcSE"]
  
  # Meta-analysis
  meta_result <- rma(yi = c(adult_effect, embryo_effect, pharate_effect),
                     sei = c(adult_se, embryo_se, pharate_se),
                     method = "FE")
  
  return(meta_result)
}
```

## Key Comparisons for Analysis

### Primary Question: Are GWAS candidate genes differentially expressed in diapause?

**Main Contrasts:**
1. **Adults**: Short-day (diapause-inducing) vs Long-day (non-diapause)
2. **Embryos**: Diapause preparation vs Non-diapause development  
3. **Pharate larvae**: Diapause maintenance vs Quiescence

### Secondary Questions:
1. **Consistency across life stages**: Do genes show same direction of effect?
2. **Platform robustness**: Are results consistent after batch correction?
3. **Blood meal interaction**: Does nutrition modify diapause gene expression?
4. **Temporal dynamics**: How do expression patterns change over time?

## Expected Outputs

### For Each Gene:
- Log2 fold change in each life stage
- Statistical significance (adjusted p-value)
- Direction of effect consistency
- Meta-analysis combined effect size
- Platform-corrected expression levels

### Summary Tables:
1. **Candidate gene expression matrix** (TPM values across all samples)
2. **Differential expression results** (log2FC, p-values for each comparison)
3. **Meta-analysis summary** (combined effects across life stages)
4. **Gene prioritization ranking** (based on consistency + significance)

This approach systematically evaluates which of the 34 candidate genes show the strongest and most consistent evidence for involvement in diapause regulation across multiple life stages and platforms.