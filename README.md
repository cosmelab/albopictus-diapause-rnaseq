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
│   ├── raw/                    # FASTQ files by dataset
│   ├── metadata/               # Sample information and candidate genes
│   ├── references/             # Genome and annotation files
│   └── sra/                    # SRA cache files
├── scripts/                    # Analysis pipeline
│   ├── 01_download_data/       # SRA data download
│   ├── 02_run_rnaseq/         # nf-core RNA-seq pipeline
│   ├── 03_qc_analysis/        # Quality control
│   ├── 04_differential_expression/ # DESeq2 analysis
│   └── 05_visualization/      # Publication figures
├── logs/                      # Job logs
├── output/                    # Pipeline outputs
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
# Download SRA data (3 datasets in parallel) - requires the .sif container
sbatch scripts/01_download_data/01_sra_array.sh

# Run nf-core RNA-seq pipeline - uses nf-core containers
sbatch scripts/02_run_rnaseq/01_run_rnaseq_array.sh
```

## Analysis Workflow

The current pipeline consists of two main steps:

1. **Data Download**: Retrieve SRA datasets using parallel array jobs
2. **RNA-seq Processing**: nf-core/rnaseq pipeline with Salmon quantification

## Key Tools

**RNA-seq Processing (nf-core containers):**
- **nf-core/rnaseq**: Standardized RNA-seq pipeline
- **Salmon**: Transcript quantification

**Data Analysis (custom container):**
- Analysis tools available for future development

## Documentation

- **[workflow.md](workflow.md)**: Detailed analysis pipeline
- **[project_rules.md](project_rules.md)**: Analysis strategy and statistical methods
- **[hpcc_nfcore_install.md](hpcc_nfcore_install.md)**: HPC deployment guide

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Project:** Albopictus Diapause RNA-seq Analysis  
**Purpose:** GWAS candidate gene validation through multi-stage RNA-seq analysis