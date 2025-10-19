# Technical Setup Documentation

This file consolidates all technical setup documentation for the albopictus-diapause-rnaseq project.

---

## Table of Contents

1. [HPC Setup - nf-core Installation](#hpc-setup---nf-core-installation)
2. [Git Hook Setup - AI Attribution Blocker](#git-hook-setup---ai-attribution-blocker)
3. [GitHub and Docker Hub Setup](#github-and-docker-hub-setup)

---

## HPC SETUP - NF-CORE INSTALLATION

Nextflow + nf-core/rnaseq Installation on UCR HPCC

This section documents the steps to install Nextflow and run the nf-core/rnaseq pipeline on UCR's HPCC.

### 1. Load Java

```bash
module load java/21.0.7
```

### 2. Install Nextflow

1. Create a personal `bin` directory:
   ```bash
   mkdir -p /bigdata/cosmelab/lcosme/bin
   ```

2. Download the Nextflow launcher:
   ```bash
   curl -fsSL https://get.nextflow.io -o /bigdata/cosmelab/lcosme/bin/nextflow
   ```

3. Make it executable:
   ```bash
   chmod +x /bigdata/cosmelab/lcosme/bin/nextflow
   ```

4. Add it to your `PATH` (e.g. in `~/.zshrc`):
   ```bash
   export PATH=/bigdata/cosmelab/lcosme/bin:$PATH
   ```

### 3. Redirect Conda Storage to Bigdata

```bash
export CONDA_ENVS_PATH=/bigdata/cosmelab/lcosme/conda/envs
export CONDA_PKGS_DIRS=/bigdata/cosmelab/lcosme/conda/pkgs
mkdir -p "$CONDA_ENVS_PATH" "$CONDA_PKGS_DIRS"
```

### 4. Create & Activate Nextflow Conda Environment

```bash
conda create -y -n nfnextflow -c bioconda -c conda-forge nextflow
conda activate nfnextflow
```

### 5. Run nf-core/rnaseq

```bash
nextflow run nf-core/rnaseq \
  -profile conda \
  --reads '/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/data/*_R{1,2}.fastq.gz' \
  --genome GRCh38 \
  --outdir /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/results
```

### 6. Verify Installation

- Check Nextflow version:
  ```bash
  nextflow -v
  ```

- Confirm pipeline launch:
  ```bash
  nextflow run nf-core/rnaseq -resume
  ```

---

## GIT HOOK SETUP - AI ATTRIBUTION BLOCKER

### What Was Done

Created a global git commit-msg hook that automatically rejects any commit messages containing AI attribution (Claude, Anthropic, etc.).

### Problem Solved

Prevents accidental commits that include:
- "Claude"
- "Anthropic"
- "Co-Authored-By: Claude"
- "Generated with Claude"
- AI assistant references
- 🤖 emoji

### Implementation Steps

#### 1. Created Global Template Directory

```bash
mkdir -p ~/.git-templates/hooks
```

#### 2. Created commit-msg Hook

```bash
cat > ~/.git-templates/hooks/commit-msg << 'EOF'
#!/bin/bash
# Reject commits with AI attribution

commit_msg_file="$1"
commit_msg=$(cat "$commit_msg_file")

# Check for forbidden terms
if echo "$commit_msg" | grep -iE "(Claude|Anthropic|Co-Authored-By: Claude|Generated with.*Claude|🤖|AI assistant)" > /dev/null; then
    echo "ERROR: Commit message contains AI attribution!"
    echo "This violates repository rules."
    echo "Please remove any references to Claude, Anthropic, or AI assistance."
    exit 1
fi

exit 0
EOF
```

#### 3. Made Hook Executable

```bash
chmod +x ~/.git-templates/hooks/commit-msg
```

#### 4. Configured Git Globally

```bash
git config --global init.templateDir ~/.git-templates
```

#### 5. Applied to Current Repository

```bash
git init
```

### How It Works

- **New repositories**: Automatically get the hook when created or cloned
- **Existing repositories**: Run `git init` in the repo directory to apply the template
- **Current repository**: Already protected

### Testing the Hook

Try to commit with AI attribution:
```bash
git commit -m "Test commit

Co-Authored-By: Claude <noreply@anthropic.com>"
```

Expected result:
```
ERROR: Commit message contains AI attribution!
This violates repository rules.
Please remove any references to Claude, Anthropic, or AI assistance.
```

### Applying to Other Existing Repositories

Navigate to any existing repository and run:
```bash
cd /path/to/your/repo
git init
```

This will copy the hook from the template directory without affecting your existing git history.

### Location of Files

- **Global template**: `~/.git-templates/hooks/commit-msg`
- **Repository hook**: `.git/hooks/commit-msg` (in each repo)
- **Git config**: Run `git config --global --get init.templateDir` to verify

### Modifying the Hook

To update the hook for all future repositories:
1. Edit `~/.git-templates/hooks/commit-msg`
2. Run `git init` in existing repositories to update them

### Reference

This hook enforces rules from `assistant_rules.md`:
- Rule 2: Never include AI assistant names, "Anthropic", or any AI attribution in commits
- Rule 3: Never use commit messages that mention AI assistance

---

## GITHUB AND DOCKER HUB SETUP

### Overview

This section explains how to set up the RNA-seq analysis project on GitHub and connect it to Docker Hub for automated image builds.

### Step 1: Create GitHub Repository

#### Option A: Web Browser (Recommended)

1. Go to [github.com](https://github.com) and log in
2. Click the "+" icon in the top right corner
3. Select "New repository"
4. Fill in the details:
   - **Repository name**: `albopictus-diapause-rnaseq`
   - **Description**: `RNA-seq analysis for validating GWAS candidate genes in Aedes albopictus diapause`
   - **Visibility**: Choose Public or Private
   - **DO NOT** check "Add a README file", "Add .gitignore", or "Choose a license"
5. Click "Create repository"

#### Option B: GitHub CLI

```bash
# Install GitHub CLI
brew install gh

# Authenticate
gh auth login

# Create repository
gh repo create albopictus-diapause-rnaseq --public --description "RNA-seq analysis for diapause GWAS validation" --source=. --remote=origin --push
```

### Step 2: Push Your Code to GitHub

After creating the repository, run these commands in your terminal:

```bash
# Load GitHub token
export GITHUB_TOKEN=$(cat ~/.github_token)

# Configure git to use the token
git config --global credential.helper store
git config --global url."https://${GITHUB_TOKEN}@github.com/".insteadOf "https://github.com/"

# Add all files to git
git add .

# Make your first commit
git commit -m "Initial commit: RNA-seq analysis environment"

# Add the GitHub repository as remote
git remote add origin https://github.com/cosmelab/albopictus-diapause-rnaseq.git

# Push to GitHub
git push -u origin main
```

### Step 3: Set Up Docker Hub

#### Create Docker Hub Account

1. Go to [hub.docker.com](https://hub.docker.com)
2. Click "Sign Up" and create an account
3. Verify your email address

#### Create Docker Hub Repository

1. Log in to Docker Hub
2. Click "Create Repository"
3. Fill in the details:
   - **Repository name**: `albopictus-diapause-rnaseq`
   - **Description**: `RNA-seq analysis environment for Aedes albopictus diapause studies`
   - **Visibility**: Choose Public or Private
4. Click "Create"

### Step 4: Connect GitHub to Docker Hub

#### Method A: GitHub Actions (Recommended)

1. **Create GitHub Actions Workflow**
   Create the file `.github/workflows/docker-build.yml`:

```yaml
name: Build and Push Docker Image

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  build:
    runs-on: ubuntu-latest

    steps:
    - name: Checkout code
      uses: actions/checkout@v3

    - name: Set up Docker Buildx
      uses: docker/setup-buildx-action@v2

    - name: Login to Docker Hub
      uses: docker/login-action@v2
      with:
        username: ${{ secrets.DOCKERHUB_USERNAME }}
        password: ${{ secrets.DOCKERHUB_TOKEN }}

    - name: Build and push
      uses: docker/build-push-action@v4
      with:
        context: .
        push: true
        tags: ${{ secrets.DOCKERHUB_USERNAME }}/albopictus-diapause-rnaseq:latest
        cache-from: type=gha
        cache-to: type=gha,mode=max
```

2. **Add Docker Hub Secrets to GitHub**
   - Go to your GitHub repository
   - Click "Settings" → "Secrets and variables" → "Actions"
   - Click "New repository secret"
   - Add these secrets:
     - `DOCKERHUB_USERNAME`: Your Docker Hub username
     - `DOCKERHUB_TOKEN`: Your Docker Hub access token

#### Method B: Docker Hub Automated Builds (Legacy)

1. In Docker Hub, go to your repository
2. Click "Builds" → "Configure Automated Builds"
3. Connect your GitHub account
4. Select your repository
5. Configure build settings:
   - **Source Type**: Branch
   - **Source**: main
   - **Docker Tag**: latest
   - **Dockerfile location**: /Dockerfile
6. Click "Create"

### Step 5: Create Docker Hub Access Token

1. Log in to Docker Hub
2. Go to "Account Settings" → "Security"
3. Click "New Access Token"
4. Give it a name (e.g., "GitHub Actions")
5. Copy the token (you won't see it again)
6. Add it to GitHub secrets as `DOCKERHUB_TOKEN`

### Step 6: Test the Setup

#### Test GitHub Repository

```bash
# Make a small change
echo "# Test" >> README.md

# Commit and push
git add README.md
git commit -m "Test commit"
git push
```

#### Test Docker Build

1. Check GitHub Actions tab to see if the build runs
2. Check Docker Hub to see if the image is built
3. Test pulling the image:
```bash
docker pull your-dockerhub-username/albopictus-diapause-rnaseq:latest
```

### Troubleshooting

#### Common Issues

1. **Authentication Errors**
   - Verify Docker Hub username and token in GitHub secrets
   - Ensure token has write permissions

2. **Build Failures**
   - Check Dockerfile syntax
   - Verify all required files are in the repository
   - Check GitHub Actions logs for specific errors

3. **Push Failures**
   - Ensure you have write access to the Docker Hub repository
   - Verify the repository name matches exactly

#### Useful Commands

```bash
# Check git remote
git remote -v

# Check Docker Hub login
docker login

# Test local build
docker build -t albopictus-diapause-rnaseq .

# Test local run
docker run -it albopictus-diapause-rnaseq
```

### Next Steps

1. **Set up branch protection** in GitHub
2. **Configure automated testing** in GitHub Actions
3. **Set up version tagging** for releases
4. **Configure multi-platform builds** if needed

### Security Notes

- Never commit secrets directly to your repository
- Use GitHub secrets for sensitive information
- Regularly rotate Docker Hub access tokens
- Consider using Docker Hub's vulnerability scanning
- GitHub token stored in `~/.github_token` (in .gitignore)

### Support

- [GitHub Documentation](https://docs.github.com/)
- [Docker Hub Documentation](https://docs.docker.com/docker-hub/)
- [GitHub Actions Documentation](https://docs.github.com/en/actions)

---

## PROJECT-SPECIFIC NOTES

- **HPC Path**: `/bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq/`
- **Container**: `albopictus-diapause-rnaseq.sif` (Singularity, 1.3GB)
- **GitHub Token**: Stored in `~/.github_token` (never committed)
- **Nextflow Config**: `scripts/03_rnaseq_pipeline/hpc_batch.conf`
