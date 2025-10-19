# 00 - Reference Preparation

**Purpose:** Prepare reference genome and annotation files for RNA-seq analysis.

**Status:** ✅ Completed (October 17, 2025)

---

## Overview

This directory contains scripts for preprocessing the *Aedes albopictus* AalbF3 reference genome and annotation files downloaded from Dryad.

## Data Source

**Repository:** Dryad Digital Repository
**DOI:** [10.5061/dryad.mgqnk98z4](https://doi.org/10.5061/dryad.mgqnk98z4)
**URL:** https://datadryad.org/stash/dataset/doi:10.5061/dryad.mgqnk98z4
**Citation:** Boyle et al. 2021, Insects 12:1-19 (PMID: 33669192); Palatini et al. 2020, BMC Biology (annotation lifted)
**Assembly:** GCA_018104305.1 (AalbF3 - linkage-based chromosome-level assembly)

### Files to Download from Dryad

1. **AALBF3_assembly.fa.gz** (genome FASTA) → rename to `AalbF3_genome_original.fa.gz`
2. **AALBF3_annotation.gff3.gz** (annotation GFF3) → rename to `AalbF3_annotation_original.gff3.gz`

**Download Location:** `data/references/AalbF3_genome/`

**Important:** NCBI's version of this assembly (GCA_018104305.1) does NOT include the complete GFF3 annotation file. The Dryad repository has the complete annotation with all gene features, rRNA, and tRNA annotations that are required for this analysis.

### Download Procedure

**Option 1: Manual Download (Recommended if automated download fails)**

1. Visit: https://datadryad.org/stash/dataset/doi:10.5061/dryad.mgqnk98z4
2. Click on each file to download:
   - `AALBF3_assembly.fa.gz` (~423 MB)
   - `AALBF3_annotation.gff3.gz` (~6.3 MB)
3. Transfer to HPC using `scp`:
   ```bash
   # From your local machine
   scp AALBF3_assembly.fa.gz username@hpcc:/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/data/references/AalbF3_genome/AalbF3_genome_original.fa.gz
   scp AALBF3_annotation.gff3.gz username@hpcc:/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/data/references/AalbF3_genome/AalbF3_annotation_original.gff3.gz
   ```

**Option 2: Automated Download (Try this first)**

Dryad may provide direct download URLs. Try:
```bash
# Create directory
mkdir -p data/references/AalbF3_genome/

# Attempt wget/curl download (you may need to get actual URLs from Dryad page)
# Visit the Dryad page, right-click "Download", and copy link address
wget -O data/references/AalbF3_genome/AalbF3_genome_original.fa.gz [DIRECT_URL_FROM_DRYAD]
wget -O data/references/AalbF3_genome/AalbF3_annotation_original.gff3.gz [DIRECT_URL_FROM_DRYAD]
```

**Why Manual Download May Be Required:**
- Dryad download links may require browser interaction
- Direct URLs may not be easily accessible via command line
- Large file sizes may cause timeout issues with automated tools
- Previous attempts at automated download failed, requiring manual intervention

---

## Scripts (Run in Order)

### 01_strip_chr_prefix.sh

**Purpose:** Remove "chr" prefix from chromosome names to match SNP chip naming convention.

**Input:**
- `data/references/AalbF3_genome/AalbF3_genome_original.fa.gz`

**Output:**
- `data/references/AalbF3_genome/AalbF3_genome.fa.gz`

**Changes:**
- `>chr1.1` → `>1.1`
- `>chr2.206` → `>2.206`
- etc.

**Why:** SNP chip data uses chromosome names without "chr" prefix. Consistency required for downstream GWAS integration.

**Usage:**
```bash
bash scripts/00_reference_preparation/01_strip_chr_prefix.sh
```

---

### 02_download_rrna_sequences.sh

**Purpose:** Download 12 *Aedes albopictus* rRNA sequences from GenBank for rRNA filtering with SortMeRNA.

**Source:** Koh et al. 2023, eLife (PMID: 37094422)

**Output:**
- `data/references/AalbF3_genome/combined_rRNA_sequences.fa`

**Sequences Downloaded:**
- 18S rRNA (5 sequences)
- 28S rRNA (4 sequences)
- 5.8S rRNA (2 sequences)
- 5S rRNA (1 sequence)

**Why:** Original run showed 0% rRNA contamination, indicating no rRNA removal was performed. Proper rRNA filtering improves accuracy of gene quantification.

**Usage:**
```bash
bash scripts/00_reference_preparation/02_download_rrna_sequences.sh
```

---

### 03_clean_rrna_headers.sh

**Purpose:** Simplify FASTA headers for SortMeRNA compatibility.

**Input:**
- `data/references/AalbF3_genome/combined_rRNA_sequences.fa`

**Output:**
- `data/references/AalbF3_genome/combined_rRNA_sequences.fa` (updated in place)
- `data/references/AalbF3_genome/rrna_db_manifest.txt` (SortMeRNA manifest file)

**Changes:**
- Remove GenBank metadata from headers
- Keep only accession and gene name
- Create manifest file for nf-core/rnaseq pipeline

**Usage:**
```bash
bash scripts/00_reference_preparation/03_clean_rrna_headers.sh
```

---

## Expected Output Files

After running all scripts, you should have:

```
data/references/AalbF3_genome/
├── AalbF3_genome.fa.gz                 # Genome FASTA (chr prefix stripped)
├── combined_rRNA_sequences.fa          # 12 rRNA sequences (cleaned headers)
└── rrna_db_manifest.txt               # SortMeRNA manifest
```

**Next Step:** Proceed to `01_sra_download/` to download raw sequencing data from SRA.

---

## Notes

- Original files with "chr" prefix are kept as `*_original.*` for reference
- rRNA sequences are specific to *Aedes albopictus* (not generic insect rRNA)
- All reference files are in `.gitignore` and will be uploaded to Zenodo for publication
- Total processing time: ~5 minutes
