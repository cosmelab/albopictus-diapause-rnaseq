# Nextflow + nf‑core/rnaseq Installation on UCR HPCC

This guide documents the steps to install Nextflow and run the nf‑core/rnaseq pipeline on UCR’s HPCC without using containers.

## 1. Load Java

```bash
module load java/21.0.7
```

## 2. Install Nextflow

1. Create a personal `bin` directory:
   ```bash
   mkdir -p /path/to/your/allocation/bin
   ```
2. Download the Nextflow launcher:
   ```bash
   curl -fsSL https://get.nextflow.io -o /path/to/your/allocation/bin/nextflow
   ```
3. Make it executable:
   ```bash
   chmod +x /path/to/your/allocation/bin/nextflow
   ```
4. Add it to your `PATH` (e.g. in `~/.zshrc`):
   ```bash
   export PATH=/path/to/your/allocation/bin:$PATH
   ```

## 3. Redirect Conda Storage to Bigdata

```bash
export CONDA_ENVS_PATH=/path/to/your/allocation/conda/envs
export CONDA_PKGS_DIRS=/path/to/your/allocation/conda/pkgs
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"
```

## 4. Create & Activate Nextflow Conda Environment

```bash
conda create -y -n nfnextflow -c bioconda -c conda-forge nextflow
conda activate nfnextflow
```

## 5. Run nf‑core/rnaseq

```bash
nextflow run nf-core/rnaseq \
  -profile conda \
  --reads '/path/to/your/allocation/data/*_R{1,2}.fastq.gz' \
  --genome GRCh38 \
  --outdir /path/to/your/allocation/results
```

## 6. Verify Installation

- Check Nextflow version:
  ```bash
  nextflow -v
  ```
- Confirm pipeline launch:
  ```bash
  nextflow run nf-core/rnaseq -resume
  ```

