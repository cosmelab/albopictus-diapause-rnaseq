# WORKFLOW.md - Step-by-Step RNA-seq Analysis Pipeline

## 🎯 **Overview**

This document provides detailed step-by-step instructions for the complete RNA-seq analysis workflow, from data download to differential expression analysis.

## 📋 **Pipeline Summary**

1. **Data Acquisition** → 2. **HPC Setup** → 3. **nf-core/Salmon** → 4. **DESeq2 Analysis** → 5. **Meta-analysis**

---

## 🔽 **Step 1: Data Acquisition**

### 1.1 SRA Download Setup

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Set up SRA cache directory
mkdir -p ~/.ncbi
export NCBI_SETTINGS=~/.ncbi/user-settings.mkfg
```

### 1.2 Download Datasets

```bash
# Dataset 1: PRJNA268379 (Adult Females)
python scripts/sra_download.py --accession PRJNA268379 --output data/raw/PRJNA268379

# Dataset 2: PRJNA158021 (Embryos)  
python scripts/sra_download.py --accession PRJNA158021 --output data/raw/PRJNA158021

# Dataset 3: PRJNA187045 (Pharate Larvae)
python scripts/sra_download.py --accession PRJNA187045 --output data/raw/PRJNA187045
```

### 1.3 Verify Downloads

```bash
# Check file integrity and completeness
python scripts/info_sra_checker.py --input data/raw/ --output data/metadata/download_summary.csv
```

---

## 🖥️ **Step 2: HPC Environment Setup**

### 2.1 Create Conda Environment

```bash
# On HPC system
conda create -n nf-core-rnaseq python=3.10
conda activate nf-core-rnaseq

# Install nf-core
pip install nf-core
```

### 2.2 Download nf-core/rnaseq

```bash
# Download the pipeline
nf-core download nf-core/rnaseq

# Check available versions
nf-core list
```

### 2.3 Configure HPC Settings

```bash
# Copy HPC configuration
cp nf-core/configs/hpc_batch.conf ~/.nextflow/config

# Test configuration
nextflow run nf-core/rnaseq -profile test,hpc_batch
```

---

## 🐟 **Step 3: Salmon Quantification with nf-core**

### 3.1 Prepare Samplesheet

```bash
# Create samplesheet for nf-core
python scripts/create_samplesheet.py \
    --input data/raw/ \
    --output samplesheet.csv \
    --format paired
```

### 3.2 Download Reference Genome

```bash
# Download Aedes albopictus genome
wget -O data/references/Aalbopictus.fasta [GENOME_URL]
wget -O data/references/Aalbopictus.gtf [ANNOTATION_URL]

# Index for Salmon
salmon index -t data/references/Aalbopictus.fasta -i data/references/salmon_index
```

### 3.3 Run nf-core Pipeline

```bash
# Run with Salmon quantification
nextflow run nf-core/rnaseq \
    -profile hpc_batch \
    --input samplesheet.csv \
    --outdir results/nf-core \
    --aligner salmon \
    --fasta data/references/Aalbopictus.fasta \
    --gtf data/references/Aalbopictus.gtf \
    --salmon_index data/references/salmon_index \
    --skip_alignment \
    --skip_pseudo_alignment \
    --skip_trimming \
    --skip_fastqc \
    --skip_multiqc
```

### 3.4 Organize Outputs

```bash
# Organize results into standard structure
python scripts/organize_outputs.py \
    --input results/nf-core \
    --output results/organized
```

---

## 📊 **Step 4: Differential Expression Analysis**

### 4.1 Prepare R Environment

```r
# Load required libraries
library(DESeq2)
library(tximport)
library(sva)
library(metafor)
library(tidyverse)
```

### 4.2 Import Salmon Quantifications

```r
# Import transcript-level abundances
files <- file.path("results/organized/quant.sf")
names(files) <- sample_names

# Import with tximport
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Summarize to gene-level
txi_gene <- summarizeToGene(txi, tx2gene = tx2gene)
```

### 4.3 Dataset-Specific Analysis

#### Adults (PRJNA268379)

```r
# Create DESeq2 dataset
dds_adults <- DESeqDataSetFromTximport(txi_gene, 
                                      colData = metadata_adults,
                                      design = ~ blood_meal + photoperiod)

# Run DESeq2
dds_adults <- DESeq(dds_adults)

# Extract results
results_adults <- results(dds_adults, 
                         contrast = c("photoperiod", "Short_day", "Long_day"))
```

#### Embryos (PRJNA158021)

```r
# Create DESeq2 dataset
dds_embryo <- DESeqDataSetFromTximport(txi_gene,
                                      colData = metadata_embryo,
                                      design = ~ timepoint + condition)

# Run DESeq2
dds_embryo <- DESeq(dds_embryo)

# Extract results
results_embryo <- results(dds_embryo,
                         contrast = c("condition", "Diapause", "Non_diapause"))
```

#### Pharate Larvae (PRJNA187045)

```r
# Create DESeq2 dataset
dds_pharate <- DESeqDataSetFromTximport(txi_gene,
                                       colData = metadata_pharate,
                                       design = ~ timepoint + condition)

# Run DESeq2
dds_pharate <- DESeq(dds_pharate)

# Extract results
results_pharate <- results(dds_pharate,
                          contrast = c("condition", "Diapause", "Quiescence"))
```

---

## 🔄 **Step 5: Cross-Platform Integration**

### 5.1 Batch Correction

```r
# Combine TPM data
tpm_combined <- cbind(tpm_adults, tpm_embryos, tpm_pharate)
metadata_combined <- rbind(meta_adults, meta_embryos, meta_pharate)

# Add platform information
metadata_combined$platform <- c(rep("HiSeq2000", nrow(meta_adults)),
                                rep("GAIIx", nrow(meta_embryos)),
                                rep("HiSeq2000", nrow(meta_pharate)))

# Apply ComBat-seq correction
tpm_corrected <- ComBat_seq(tpm_combined,
                           batch = metadata_combined$platform,
                           group = metadata_combined$diapause_status)
```

### 5.2 Meta-Analysis

```r
# Function to combine effects across datasets
combine_effects <- function(gene_id) {
  # Extract effect sizes and standard errors
  effects <- c(results_adults[gene_id, "log2FoldChange"],
               results_embryo[gene_id, "log2FoldChange"],
               results_pharate[gene_id, "log2FoldChange"])
  
  ses <- c(results_adults[gene_id, "lfcSE"],
           results_embryo[gene_id, "lfcSE"],
           results_pharate[gene_id, "lfcSE"])
  
  # Fixed-effects meta-analysis
  meta_result <- rma(yi = effects, sei = ses, method = "FE")
  return(meta_result)
}

# Apply to candidate genes
candidate_genes <- c("LOC109397825", "LOC109405370", "LOC109398973", ...)
meta_results <- lapply(candidate_genes, combine_effects)
```

---

## 📈 **Step 6: Visualization and Reporting**

### 6.1 Quality Control Plots

```r
# PCA plots
pca_data <- plotPCA(vsd, intgroup = "condition")
ggsave("results/analysis/visualization/pca_plot.png", pca_data)

# Volcano plots
volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_point(data = candidate_results, color = "red")
ggsave("results/analysis/visualization/volcano_plot.png", volcano_plot)
```

### 6.2 Candidate Gene Analysis

```r
# Extract candidate gene results
candidate_results <- results_df[results_df$gene_id %in% candidate_genes, ]

# Create heatmap
candidate_expression <- tpm_corrected[candidate_genes, ]
pheatmap(candidate_expression, 
         annotation_col = metadata_combined,
         filename = "results/analysis/visualization/candidate_heatmap.png")
```

### 6.3 Generate Reports

```r
# Create summary tables
write.csv(candidate_results, "results/analysis/candidate_gene_results.csv")
write.csv(meta_summary, "results/analysis/meta_analysis_summary.csv")

# Generate publication figures
source("scripts/visualization/create_publication_figures.R")
```

---

## 🔍 **Step 7: Quality Control and Validation**

### 7.1 Check Alignment Quality

```bash
# Review MultiQC reports
open results/organized/multiqc/multiqc_report.html

# Check alignment statistics
cat results/organized/alignment_stats/summary.txt
```

### 7.2 Validate Results

```r
# Check for batch effects
plotPCA(vsd, intgroup = "platform")

# Verify candidate gene expression
plot_candidate_expression(candidate_genes, tpm_corrected, metadata_combined)
```

### 7.3 Reproducibility Checks

```bash
# Save session info
Rscript -e "sessionInfo()" > results/analysis/session_info.txt

# Document versions
conda list > results/analysis/conda_packages.txt
```

---

## 📋 **Expected Outputs**

### Files Generated

- `results/organized/quant.sf` - Salmon quantification files
- `results/analysis/differential_expression/` - DESeq2 results
- `results/analysis/visualization/` - Plots and figures
- `results/analysis/candidate_gene_results.csv` - Candidate gene analysis
- `results/analysis/meta_analysis_summary.csv` - Meta-analysis results

### Key Metrics

- **Alignment rates**: >80% for each dataset
- **Gene detection**: >10,000 genes per sample
- **Differential expression**: FDR < 0.05 for candidate genes
- **Batch correction**: PCA shows reduced platform separation

---

## ⚠️ **Troubleshooting**

### Common Issues

1. **SRA download failures**: Check network connection, retry with smaller batches
2. **HPC job failures**: Check resource limits, adjust memory/CPU requests
3. **DESeq2 convergence**: Remove low-count genes, check model specification
4. **Batch correction**: Verify platform information, check for residual effects

### Debugging Commands

```bash
# Check file integrity
md5sum data/raw/*/*.fastq.gz

# Monitor HPC jobs
squeue -u $USER

# Check pipeline logs
tail -f .nextflow.log
```

---

**Last Updated:** [Current Date]
**Pipeline Version:** 1.0
**Dependencies:** nf-core/rnaseq, Salmon, DESeq2, sva, metafor
