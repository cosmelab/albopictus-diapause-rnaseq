# 🧬 Aedes Diapause RNA-seq Analysis

<div align="center">

<!-- Primary Badges -->
![Species](https://img.shields.io/badge/Species-Aedes%20albopictus-brightgreen?style=for-the-badge&logo=dna)
![Analysis](https://img.shields.io/badge/Analysis-RNA--seq-blue?style=for-the-badge&logo=chart-line)
![Container](https://img.shields.io/badge/Container-Ready-orange?style=for-the-badge&logo=docker)

<!-- Secondary Badges -->
![Python](https://img.shields.io/badge/Python-3.10-blue?style=flat-square&logo=python)
![R](https://img.shields.io/badge/R-4.3.2-blue?style=flat-square&logo=r)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Platform](https://img.shields.io/badge/Platform-Linux%20HPC-lightgrey?style=flat-square)

**🧬 GWAS Candidate Gene Validation Through RNA-seq Analysis**

*Multi-stage differential expression analysis pipeline for validating diapause-associated genes*

<!-- Container Registry Badges -->
[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-cosmelab%2Falbopictus--diapause--rnaseq-2496ED?style=flat&logo=docker&logoColor=white)](https://hub.docker.com/r/cosmelab/albopictus-diapause-rnaseq)
[![GHCR](https://img.shields.io/badge/GHCR-ghcr.io%2Fcosmelab%2Falbopictus--diapause--rnaseq-181717?style=flat&logo=github&logoColor=white)](https://github.com/cosmelab/albopictus-diapause-rnaseq/pkgs/container/albopictus-diapause-rnaseq)

<!-- Quick Stats -->
**📊 Analysis Power:** 34 Candidate Genes | **🐳 Container Size:** ~4GB | **🏗️ Architecture:** AMD64 (Optimized)

</div>

---

## 📋 **Table of Contents**

- [🎯 Project Overview](#-project-overview)
- [🚀 Quick Start](#-quick-start)
- [🏗️ Project Structure](#️-project-structure)
- [🐳 Container Usage](#-container-usage)
- [🛠️ Analysis Pipeline](#️-analysis-pipeline)
- [🔗 Remote Development Setup](#-remote-development-setup)
- [📚 Documentation](#-documentation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)

---

## 🎯 **Project Overview**

This repository contains the RNA-seq analysis pipeline for validating GWAS-identified candidate genes associated with diapause regulation in *Aedes albopictus* across multiple life stages and experimental conditions.

### Research Context

- **Species**: *Aedes albopictus* (Asian tiger mosquito)
- **Research Goal**: Validate 34 GWAS candidate genes through differential expression analysis
- **Main Question**: Are genes identified in GWAS for diapause regulation differentially expressed between diapausing and non-diapausing conditions?
- **Analysis Focus**: Cross-platform integration and meta-analysis across three life stages
- **Methodology**: DESeq2-based differential expression with batch correction and meta-analysis

<div align="center">

### 📊 **Project Statistics**

<table>
<tr>
<td align="center" width="16.66%">

**🧬 Datasets**
<br>
`3 Independent Studies`

</td>
<td align="center" width="16.66%">

**📈 Life Stages**
<br>
`Adult, Embryo, Pharate`

</td>
<td align="center" width="16.66%">

**🔬 Candidate Genes**
<br>
`34 GWAS Genes`

</td>
<td align="center" width="16.66%">

**🐳 Container Size**
<br>
`~4GB`

</td>
<td align="center" width="16.66%">

**🏗️ Architecture**
<br>
`AMD64 (Optimized)`

</td>
<td align="center" width="16.66%">

**📚 Documentation**
<br>
`Complete`

</td>
</tr>
</table>

</div>

### 📊 **Experimental Design**

<div align="center">

<table>
<tr>
<th width="25%">📁 Dataset</th>
<th width="25%">🧬 Life Stage</th>
<th width="25%">🔬 Design</th>
<th width="25%">📈 Platform</th>
</tr>
<tr>
<td><b>PRJNA268379</b></td>
<td>Adult Females</td>
<td>2×2 Factorial (Photoperiod × Blood meal)</td>
<td>HiSeq2000</td>
</tr>
<tr>
<td><b>PRJNA158021</b></td>
<td>Embryos</td>
<td>Time course (72-78h, 135-141h)</td>
<td>GAIIx</td>
</tr>
<tr>
<td><b>PRJNA187045</b></td>
<td>Pharate Larvae</td>
<td>Diapause vs Quiescence (11, 21, 40d)</td>
<td>HiSeq2000</td>
</tr>
</table>

</div>

## 👥 **For New Collaborators**

<div align="center">

### 🚀 **Quick Start (4 Steps)**

**New to the project?** Get up and running in minutes:

</div>

<details>
<summary><b>📋 Step-by-Step Setup Guide</b> (Click to expand)</summary>

<table>
<tr>
<td align="center" width="60">

### 1️⃣

</td>
<td>

**🔧 Clone Repository & Setup**
```bash
git clone https://github.com/cosmelab/albopictus-diapause-rnaseq.git
cd albopictus-diapause-rnaseq
./setup.sh  # Creates directories + environment setup
```
*Creates project structure and prepares analysis environment*

</td>
</tr>
<tr>
<td align="center" width="60">

### 2️⃣

</td>
<td>

**🐳 Build Container (HPC Required)**
```bash
# Load Singularity module (adjust for your HPC)
module load singularity-ce/3.9.3

# Set up cache and build container
mkdir singularity_temp_cache
export SINGULARITY_CACHEDIR=$PWD/singularity_temp_cache
singularity build albopictus-diapause-rnaseq.sif docker://ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest

# Clean up cache
rm -rf singularity_temp_cache
```
*Container is needed for SRA download and analysis scripts*

</td>
</tr>
<tr>
<td align="center" width="60">

### 3️⃣

</td>
<td>

**⚙️ Configure HPC Settings**
```bash
# Copy and customize HPC configuration
cp hpc_config.template.sh hpc_config.sh
# Edit hpc_config.sh with your specific paths

# Source configuration (before running scripts)
source hpc_config.sh
```
*Customizes paths and modules for your HPC environment*

</td>
</tr>
<tr>
<td align="center" width="60">

### 4️⃣

</td>
<td>

**🧪 Test Everything Works**
```bash
# Test the container
singularity shell --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif
# Inside container:
cd /proj && python --version && R --version
```
*Verify container and tools are working*

</td>
</tr>
</table>

</details>

<div align="center">

**🎉 That's it!** You now have:

| Status | Component |
|--------|-----------|
| ✅ | Project directory structure |
| ✅ | All RNA-seq analysis tools in container |
| ✅ | Scripts and configurations |
| ✅ | Ready to start analysis |

</div>

### 📚 **Next Steps for Collaborators**

<table>
<tr>
<td width="33%">

**🖥️ HPC Users**
- See [hpcc_nfcore_install.md](hpcc_nfcore_install.md)
- SLURM job templates
- Resource recommendations

</td>
<td width="33%">

**📊 Analysis Guide**
- Check [WORKFLOW.md](WORKFLOW.md)
- Analysis workflows
- Step-by-step pipeline

</td>
<td width="33%">

**📖 Project Rules**
- [PROJECT_RULES.md](PROJECT_RULES.md)
- Analysis strategy
- Statistical models

</td>
</tr>
</table>

---

## 🏗️ **Project Structure**

```
albopictus-diapause-rnaseq/
├── data/
│   ├── raw/                    # FASTQ files by dataset (PRJNA*)
│   ├── metadata/               # Sample information and candidate genes
│   ├── references/             # Genome and annotation files (from Zenodo)
│   └── sra/                    # SRA cache files
├── scripts/                    # 🔢 Numbered Workflow Scripts
│   ├── 01_download_data/       # SRA data download
│   ├── 02_run_rnaseq/         # nf-core RNA-seq pipeline
│   ├── 03_qc_analysis/        # Quality control and plots
│   ├── 04_differential_expression/ # DESeq2 analysis
│   ├── 05_visualization/      # Publication figures
│   └── utils/                 # Utility scripts
├── nf-core/
│   └── configs/               # HPC configuration files
├── logs/                      # SLURM job logs
├── work/                      # Nextflow work directory
├── temp/                      # Temporary files
├── singularity/               # Singularity container cache
├── results/                   # Analysis results
├── output/                    # Final outputs and figures
├── Dockerfile                 # Container definition
└── hpc_config.template.sh     # HPC configuration template
```

## 🐳 **Container Usage**

This project provides a comprehensive RNA-seq analysis environment through Docker and Singularity containers with optimized performance.

### Available Containers

The analysis environment is available from two registries (choose either):

#### GitHub Container Registry (GHCR)
```bash
# Pull from GHCR (GitHub integration)
docker pull ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest
```

#### Docker Hub
```bash
# Pull from Docker Hub (widely supported)
docker pull cosmelab/albopictus-diapause-rnaseq:latest
```

**Both registries work great!** Choose based on your preference or institutional requirements.

### Singularity on HPC

For HPC systems without Docker, use Singularity with either registry:

```bash
# Load Singularity module (if required)
module load singularity-ce/3.9.3

# Choose either registry:
# From GHCR:
singularity pull albopictus-diapause-rnaseq.sif docker://ghcr.io/cosmelab/albopictus-diapause-rnaseq:latest
# OR from Docker Hub:
singularity pull albopictus-diapause-rnaseq.sif docker://cosmelab/albopictus-diapause-rnaseq:latest
```

### Running the Container

#### Interactive Shell
```bash
# Start interactive shell with project directory mounted
singularity shell \
    --cleanenv \
    --bind /path/to/your/project:/proj \
    albopictus-diapause-rnaseq.sif

# Once inside the container, navigate to the project directory
cd /proj
```

#### Execute Commands
```bash
# Run DESeq2 analysis
singularity exec \
    --cleanenv \
    --bind /path/to/your/project:/proj \
    albopictus-diapause-rnaseq.sif \
    Rscript /proj/scripts/differential_expression/PRJNA268379_adult_females_Angela_DESeq.R
```

### SLURM Batch Script Examples

#### Example 1: nf-core RNA-seq Pipeline
```bash
#!/bin/bash
#SBATCH --job-name=rnaseq-pipeline
#SBATCH --output=rnaseq_%j.out
#SBATCH --error=rnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=epyc

# Load required modules
module load singularity-ce/3.9.3

# Set paths
PROJECT_DIR="/path/to/your/project"
CONTAINER="albopictus-diapause-rnaseq.sif"

# Run nf-core RNA-seq pipeline
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    nextflow run nf-core/rnaseq \
    -profile singularity \
    --input /proj/samplesheet.csv \
    --outdir /proj/results/nf-core \
    --genome custom \
    --fasta /proj/data/references/genome.fa \
    --gtf /proj/data/references/annotation.gtf
```

#### Example 2: Differential Expression Analysis
```bash
#!/bin/bash
#SBATCH --job-name=deseq2-analysis
#SBATCH --output=deseq2_%j.out
#SBATCH --error=deseq2_%j.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=epyc

module load singularity-ce/3.9.3

PROJECT_DIR="/path/to/your/project"
CONTAINER="albopictus-diapause-rnaseq.sif"

# Run differential expression analysis
singularity exec \
    --cleanenv \
    --bind ${PROJECT_DIR}:/proj \
    ${CONTAINER} \
    bash -c "cd /proj && Rscript scripts/differential_expression/PRJNA268379_adult_females_Angela_DESeq.R"
```

## 🛠️ **Analysis Pipeline**

<div align="center">

### 🔢 **Numbered Workflow Execution**

*Step-by-step pipeline with organized scripts*

</div>

### 🚀 **Running the Complete Workflow**

```bash
# Step 1: Download SRA data (3 datasets in parallel)
sbatch scripts/01_download_data/01_sra_array.sh

# Step 2: Run nf-core RNA-seq pipeline (after Step 1 completes)
sbatch scripts/02_run_rnaseq/01_run_rnaseq_array.sh

# Step 3: Quality control analysis
python scripts/03_qc_analysis/01_create_qc_plots.py
python scripts/03_qc_analysis/02_create_qc_plots_advanced.py
python scripts/03_qc_analysis/04_create_qc_table.py

# Step 4: Differential expression analysis
Rscript scripts/04_differential_expression/PRJNA268379_adult_females_Angela_DESeq.R
python scripts/04_differential_expression/de_analysis_PRJNA268379.py

# Step 5: Create publication figures
python scripts/05_visualization/01_combine_counts.py
python scripts/05_visualization/02_create_summary_table.py
python scripts/05_visualization/03_create_publication_figures.py
```

<div align="center">

### 🔬 **Complete RNA-seq Analysis Suite**

*Comprehensive pipeline from raw reads to differential expression and meta-analysis*

</div>

<details>
<summary><b>🧬 Core RNA-seq Tools</b> (Click to expand)</summary>

<table>
<tr>
<th width="20%">🔧 Tool</th>
<th width="15%">📦 Version</th>
<th width="65%">🎯 Purpose</th>
</tr>
<tr>
<td><b>🐟 Salmon</b></td>
<td><code>1.10.1</code></td>
<td>Ultra-fast transcript quantification</td>
</tr>
<tr>
<td><b>📊 DESeq2</b></td>
<td><code>1.42.0</code></td>
<td>Differential expression analysis</td>
</tr>
<tr>
<td><b>📥 tximport</b></td>
<td><code>1.30.0</code></td>
<td>Import transcript-level estimates for gene-level analysis</td>
</tr>
<tr>
<td><b>🔧 sva</b></td>
<td><code>3.50.0</code></td>
<td>Batch effect correction (ComBat-seq)</td>
</tr>
<tr>
<td><b>📈 metafor</b></td>
<td><code>4.4-0</code></td>
<td>Meta-analysis across datasets</td>
</tr>
<tr>
<td><b>⚡ FastQC</b></td>
<td><code>0.12.1</code></td>
<td>Quality control of sequencing data</td>
</tr>
<tr>
<td><b>📋 MultiQC</b></td>
<td><code>1.15</code></td>
<td>Aggregate quality control reports</td>
</tr>
</table>

</details>

<details>
<summary><b>🐍 Python Environment</b> (Click to expand)</summary>

**Python 3.10** with comprehensive data science packages:

<table>
<tr>
<td width="50%">

**📊 Data Analysis**
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scipy` - Scientific computing
- `scikit-learn` - Machine learning

</td>
<td width="50%">

**📈 Visualization**
- `matplotlib` - Plotting
- `seaborn` - Statistical visualization
- `plotly` - Interactive plots
- `bokeh` - Web-based visualization

</td>
</tr>
<tr>
<td>

**🧬 Bioinformatics**
- `biopython` - Sequence analysis
- `pysam` - SAM/BAM processing
- `requests` - Data download
- `beautifulsoup4` - Web scraping

</td>
<td>

**🔧 Utilities**
- `xmltodict` - XML processing
- `lxml` - XML/HTML parsing
- `argparse` - Command line interfaces
- `logging` - Progress tracking

</td>
</tr>
</table>

</details>

<details>
<summary><b>📊 R Environment</b> (Click to expand)</summary>

**R 4.3.2** with specialized RNA-seq analysis packages:

<table>
<tr>
<td width="50%">

**🧬 RNA-seq Analysis**
- `DESeq2` - Differential expression
- `tximport` - Transcript import
- `sva` - Batch correction
- `metafor` - Meta-analysis
- `edgeR` - Alternative DE analysis
- `limma` - Linear modeling

</td>
<td width="50%">

**📈 Visualization & Stats**
- `ggplot2` - Data visualization
- `qqman` - Manhattan plots
- `qqplotr` - Q-Q plots
- `pheatmap` - Heatmaps
- `RColorBrewer` - Color palettes
- `ggrepel` - Text positioning

</td>
</tr>
<tr>
<td>

**🔧 Data Manipulation**
- `tidyverse` - Data science toolkit
- `data.table` - Fast data processing
- `readxl` / `writexl` - Excel files
- `here` - Path management

</td>
<td>

**🧮 Bioconductor**
- `AnnotationDbi` - Annotation databases
- `biomaRt` - Gene annotation
- `GenomicRanges` - Genomic intervals
- `S4Vectors` - Vector operations

</td>
</tr>
</table>

</details>

## 🔬 **Analysis Workflow**

<div align="center">

### 🔄 **Complete Analysis Pipeline**

</div>

<details>
<summary><b>📊 Step-by-Step Workflow</b> (Click to expand)</summary>

<table>
<tr>
<td align="center" width="25%">

### 📥 **Data Acquisition**

**🎯 Purpose**
<br>
Download SRA datasets for three life stages

**🔧 Tools**
<br>
SRA-tools, Python scripts

**📈 Output**
<br>
FASTQ files organized by dataset

</td>
<td align="center" width="25%">

### 🧬 **Quantification**

**🎯 Purpose**
<br>
Transcript-level abundance estimation

**🔧 Tools**
<br>
Salmon, nf-core/rnaseq pipeline

**📈 Output**
<br>
Gene and transcript count matrices

</td>
<td align="center" width="25%">

### 📊 **Differential Expression**

**🎯 Purpose**
<br>
Identify differentially expressed genes

**🔧 Tools**
<br>
DESeq2, tximport

**📈 Output**
<br>
Log2FC, p-values for each comparison

</td>
<td align="center" width="25%">

### 🔗 **Meta-Analysis**

**🎯 Purpose**
<br>
Combine results across life stages

**🔧 Tools**
<br>
metafor, sva (batch correction)

**📈 Output**
<br>
Combined effect sizes and rankings

</td>
</tr>
</table>

</details>

<details>
<summary><b>⚙️ Detailed Analysis Methods</b> (Click to expand)</summary>

#### 🧬 **Per-Dataset Analysis**
- **Adults (PRJNA268379)**: 2×2 factorial design (photoperiod × blood meal)
- **Embryos (PRJNA158021)**: Time course with diapause comparison  
- **Pharate Larvae (PRJNA187045)**: Diapause vs quiescence comparison

#### 🔗 **Cross-Platform Integration**
- **Batch Correction**: ComBat-seq for platform effects (HiSeq2000 vs GAIIx)
- **Meta-Analysis**: Fixed-effects combination of log2FC across datasets
- **Gene Prioritization**: Ranking based on consistency and significance

#### 📊 **Statistical Models**
- **DESeq2**: Negative binomial GLM with appropriate design matrices
- **Batch Correction**: sva package ComBat-seq method
- **Meta-Analysis**: metafor package fixed-effects models

</details>

## 📊 **Expected Results**

<div align="center">

### 🎯 **Key Deliverables**

</div>

<details>
<summary><b>📈 Primary Outputs</b> (Click to expand)</summary>

<table>
<tr>
<td width="50%">

**🧬 Gene Expression Data**
<br>
- Normalized count matrices (TPM values)
- Platform-corrected expression levels
- Quality control metrics and plots

**📊 Differential Expression**
<br>
- Log2 fold changes for each life stage
- Statistical significance (adjusted p-values)
- Effect size consistency across datasets

</td>
<td width="50%">

**🔗 Meta-Analysis Results**
<br>
- Combined effect sizes across life stages
- Gene prioritization rankings
- Cross-platform validation results

**📈 Visualization**
<br>
- PCA plots for sample relationships
- Volcano plots for differential expression
- Heatmaps for candidate gene patterns

</td>
</tr>
</table>

</details>

<details>
<summary><b>🎯 Candidate Gene Focus</b> (Click to expand)</summary>

**Primary Analysis Target**: 34 GWAS-identified genes with LOC identifiers

#### Expected Validation Outcomes:
- **Consistent Direction**: Genes showing same direction of effect across life stages
- **Statistical Significance**: FDR < 0.05 in at least one life stage  
- **Effect Size**: Meaningful log2FC values (|log2FC| > 0.5)
- **Cross-Platform Validation**: Reproducible results after batch correction

#### Gene Prioritization Criteria:
1. **Consistency Score**: Same direction across datasets
2. **Statistical Strength**: Combined p-value from meta-analysis  
3. **Effect Magnitude**: Average absolute log2FC across stages
4. **Platform Independence**: Robust to batch correction

</details>

## 🔗 **Remote Development Setup**

### SSH Configuration for HPC Access

To access your HPC system remotely via Cursor, VS Code, or other SSH clients:

#### 1. SSH Config Setup
Add this to your `~/.ssh/config` file:

```bash
# HPC SSH Configuration
Host ucr-hpc
  HostName cluster.hpcc.ucr.edu
  User your-username
  IdentityFile ~/.ssh/id_rsa
  AddKeysToAgent yes
  UseKeychain yes
  ControlMaster auto
  ControlPath ~/.ssh/cm-%r@%h:%p
  ControlPersist 600
```

#### 2. SSH Key Setup
```bash
# Generate SSH key (if you don't have one)
ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa

# Copy public key to HPC
ssh-copy-id -i ~/.ssh/id_rsa.pub your-username@cluster.hpcc.ucr.edu

# Test connection
ssh ucr-hpc
```

### SSH Tunneling for Jupyter

```bash
# Create SSH tunnel for Jupyter
ssh -L 8888:localhost:8888 ucr-hpc

# On HPC, start Jupyter
singularity exec --cleanenv --bind /path/to/project:/proj albopictus-diapause-rnaseq.sif \
    bash -c "cd /proj && jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root"

# Open in browser: http://localhost:8888
```

## 📚 **Documentation**

- **[workflow.md](workflow.md)**: Step-by-step analysis pipeline and methodology
- **[project_rules.md](project_rules.md)**: Analysis strategy and statistical models  
- **[hpcc_nfcore_install.md](hpcc_nfcore_install.md)**: HPC deployment and Nextflow installation
- **[github_dockerhub_setup.md](github_dockerhub_setup.md)**: GitHub and Docker Hub integration setup
- **[de_analysis_design_repo.md](de_analysis_design_repo.md)**: Differential expression analysis design
- **[assistant_rules.md](assistant_rules.md)**: AI assistant guidelines and rules

## 🤝 **Contributing**

This project is designed for publication purposes. For questions or issues, please open a GitHub issue or contact the [Cosme Lab](https://github.com/cosmelab).

## 📄 **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 **Acknowledgments**

- **Cosme Lab team** for project development
- **Bioconductor** for R package ecosystem  
- **nf-core** for standardized RNA-seq pipelines
- **Conda-forge** for Python package management
- **Singularity** for HPC containerization

---

**Last Updated:** 2025-01-18  
**Project:** Albopictus Diapause RNA-seq Analysis  
**Purpose:** GWAS candidate gene validation through RNA-seq analysis