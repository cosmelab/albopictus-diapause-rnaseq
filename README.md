# Aedes albopictus Diapause RNA-seq Analysis

<div align="center">

![Species](https://img.shields.io/badge/Species-Aedes%20albopictus-brightgreen?style=for-the-badge&logo=dna)
![Analysis](https://img.shields.io/badge/Analysis-RNA--seq-blue?style=for-the-badge&logo=chart-line)
![Container](https://img.shields.io/badge/Container-Ready-orange?style=for-the-badge&logo=docker)

**GWAS Candidate Gene Validation Through RNA-seq Analysis**

*Multi-stage differential expression analysis pipeline for validating diapause-associated genes*

[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Falbopictus--diapause--rnaseq-2496ED?style=flat&logo=docker&logoColor=white)](https://hub.docker.com/r/cosmelab/albopictus-diapause-rnaseq)
[![GHCR](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Falbopictus--diapause--rnaseq-181717?style=flat&logo=github&logoColor=white)](https://github.com/cosmelab/albopictus-diapause-rnaseq/pkgs/container/albopictus-diapause-rnaseq)

**📊 Analysis:** 34 Candidate Genes | **🐳 Container:** 1.3GB | **🏗️ Architecture:** AMD64

</div>

---

## Project Overview

RNA-seq analysis pipeline for validating GWAS-identified candidate genes associated with diapause regulation in *Aedes albopictus* across multiple life stages.

### Research Context

- **Species**: *Aedes albopictus* (Asian tiger mosquito)
- **Goal**: Validate 34 GWAS candidate genes through differential expression analysis
- **Question**: Are GWAS-identified genes differentially expressed between diapausing and non-diapausing conditions?
- **Approach**: Cross-platform integration and meta-analysis across three life stages

### Datasets

<div align="center">

| Dataset | Life Stage | Design | Platform |
|---------|------------|---------|----------|
| **PRJNA268379** | Adult Females | 2×2 Factorial (Photoperiod × Blood meal) | HiSeq2000 |
| **PRJNA158021** | Embryos | Time course (72-78h, 135-141h) | GAIIx |
| **PRJNA187045** | Pharate Larvae | Diapause vs Quiescence (11, 21, 40d) | HiSeq2000 |

</div>

## Project Structure

```
albopictus-diapause-rnaseq/
├── data/
│   ├── metadata/               # Sample information and candidate genes
│   ├── references/             # AalbF3 genome and annotation files
│   ├── collaborator_repos/     # Collaborator's analysis scripts (for comparison)
│   └── raw/                    # FASTQ files by dataset
├── scripts/                    # Analysis pipeline (numbered workflow)
│   ├── 00_reference_preparation/  # Download and prepare reference genome
│   ├── 01_sra_download/          # SRA data download
│   ├── 02_annotation_prep/       # Fix GTF annotation for nf-core
│   ├── 03_rnaseq_pipeline/       # nf-core RNA-seq pipeline
│   ├── 04_qc_analysis/           # Quality control extraction
│   ├── 05_count_matrix/          # Combine count matrices
│   ├── 06_diff_expression/       # DESeq2 differential expression
│   ├── 07_visualization/         # Publication figures
│   └── utils/                    # Utility scripts
├── logs/                      # SLURM job logs
├── output/                    # Pipeline outputs (44 samples)
└── albopictus-diapause-rnaseq.sif # Analysis container (1.3GB)
```

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/cosmelab/albopictus-diapause-rnaseq.git
cd albopictus-diapause-rnaseq
```

### 2. Build Container (Required First)
```bash
module load singularity-ce/3.9.3

# Choose either registry:
# From GitHub Container Registry:
singularity pull albopictus-diapause-rnaseq.sif docker://ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest

# OR from Docker Hub:
singularity pull albopictus-diapause-rnaseq.sif docker://cosmelab/albopictus-diapause-rnaseq:latest

# Test the container works:
singularity shell --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif
```

### 3. Clone Collaborator Reference Repository
```bash
# Clone the collaborator's analysis repository for comparison
mkdir -p data/collaborator_repos
cd data/collaborator_repos
git clone https://github.com/srmarzec/albopictus_remapping.git
cd ../..

# This contains their original analysis scripts using STAR + HTSeq
```

### 4. Run Analysis Pipeline
```bash
# Download SRA data (44 samples in parallel)
sbatch --array=1-44 scripts/01_sra_download/01_sra_array.sh

# Run nf-core RNA-seq pipeline (completed Oct 17-18, 2025)
sbatch --array=1-44 scripts/03_rnaseq_pipeline/02_run_rnaseq_array.sh
```

## Analysis Workflow

**Current Status:** Pipeline complete (Oct 17-18, 2025) - All 44 samples processed successfully!

The complete workflow consists of 8 numbered stages:

1. **00_reference_preparation**: Download AalbF3 genome from Dryad, prepare rRNA database
2. **01_sra_download**: Download 44 FASTQ samples from SRA (3 datasets)
3. **02_annotation_prep**: Fix GTF annotation (remove 45 problematic transcripts, add gene_biotype)
4. **03_rnaseq_pipeline**: nf-core/rnaseq with STAR alignment + Salmon quantification
5. **04_qc_analysis**: Extract QC metrics, run featureCounts for validation (IN PROGRESS)
6. **05_count_matrix**: Combine count matrices across samples
7. **06_diff_expression**: DESeq2 analysis per dataset + meta-analysis
8. **07_visualization**: Generate publication figures

**Next Step:** Run featureCounts on existing BAM files to validate Salmon results

## Key Tools

**RNA-seq Processing (nf-core containers):**
- **nf-core/rnaseq**: Standardized RNA-seq pipeline
- **Salmon**: Transcript quantification

**Data Analysis (custom container):**
- Analysis tools available for future development

## Documentation

- **[project.md](project.md)**: Complete project documentation, status, and analysis strategy
- **[technical.md](technical.md)**: HPC setup, git hooks, GitHub/Docker Hub integration
- **[assistant_rules.md](assistant_rules.md)**: AI assistant behavior rules for this project

**See also:**
- Each `scripts/XX_*/README.md` contains detailed instructions for that workflow step
- `scripts/00_reference_preparation/README.md` - Reference genome download and preparation
- `scripts/02_annotation_prep/README.md` - GTF annotation fixes and validation

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Project:** Albopictus Diapause RNA-seq Analysis  
**Purpose:** GWAS candidate gene validation through multi-stage RNA-seq analysis