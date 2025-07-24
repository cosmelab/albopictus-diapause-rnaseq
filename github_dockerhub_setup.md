# GitHub and Docker Hub Setup Guide

## Overview
This guide explains how to set up your RNA-seq analysis project on GitHub and connect it to Docker Hub for automated image builds.

## Step 1: Create GitHub Repository

### Option A: Web Browser (Recommended)
1. Go to [github.com](https://github.com) and log in
2. Click the "+" icon in the top right corner
3. Select "New repository"
4. Fill in the details:
   - **Repository name**: `DiapauseRNASeq`
   - **Description**: `RNA-seq analysis environment for diapause studies`
   - **Visibility**: Choose Public or Private
   - **DO NOT** check "Add a README file", "Add .gitignore", or "Choose a license"
5. Click "Create repository"

### Option B: GitHub CLI
```bash
# Install GitHub CLI
brew install gh

# Authenticate
gh auth login

# Create repository
gh repo create DiapauseRNASeq --public --description "RNA-seq analysis environment" --source=. --remote=origin --push
```

## Step 2: Push Your Code to GitHub

After creating the repository, run these commands in your terminal:

```bash
# Add all files to git
git add .

# Make your first commit
git commit -m "Initial commit: RNA-seq analysis environment"

# Add the GitHub repository as remote
git remote add origin https://github.com/cosmelab/albopictus-diapause-rnaseq.git

# Push to GitHub
git push -u origin main
```

## Step 3: Set Up Docker Hub

### Create Docker Hub Account
1. Go to [hub.docker.com](https://hub.docker.com)
2. Click "Sign Up" and create an account
3. Verify your email address

### Create Docker Hub Repository
1. Log in to Docker Hub
2. Click "Create Repository"
3. Fill in the details:
   - **Repository name**: `diapause-rnaseq`
   - **Description**: `RNA-seq analysis environment for diapause studies`
   - **Visibility**: Choose Public or Private
4. Click "Create"

## Step 4: Connect GitHub to Docker Hub

### Method A: GitHub Actions (Recommended)

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
        tags: ${{ secrets.DOCKERHUB_USERNAME }}/diapause-rnaseq:latest
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

### Method B: Docker Hub Automated Builds (Legacy)
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

## Step 5: Create Docker Hub Access Token

1. Log in to Docker Hub
2. Go to "Account Settings" → "Security"
3. Click "New Access Token"
4. Give it a name (e.g., "GitHub Actions")
5. Copy the token (you won't see it again)
6. Add it to GitHub secrets as `DOCKERHUB_TOKEN`

## Step 6: Test the Setup

### Test GitHub Repository
```bash
# Make a small change
echo "# Test" >> README.md

# Commit and push
git add README.md
git commit -m "Test commit"
git push
```

### Test Docker Build
1. Check GitHub Actions tab to see if the build runs
2. Check Docker Hub to see if the image is built
3. Test pulling the image:
```bash
docker pull your-dockerhub-username/diapause-rnaseq:latest
```

## Troubleshooting

### Common Issues

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

### Useful Commands

```bash
# Check git remote
git remote -v

# Check Docker Hub login
docker login

# Test local build
docker build -t diapause-rnaseq .

# Test local run
docker run -it diapause-rnaseq
```

## Next Steps

1. **Set up branch protection** in GitHub
2. **Configure automated testing** in GitHub Actions
3. **Set up version tagging** for releases
4. **Configure multi-platform builds** if needed

## Security Notes

- Never commit secrets directly to your repository
- Use GitHub secrets for sensitive information
- Regularly rotate Docker Hub access tokens
- Consider using Docker Hub's vulnerability scanning

## Support

- [GitHub Documentation](https://docs.github.com/)
- [Docker Hub Documentation](https://docs.docker.com/docker-hub/)
- [GitHub Actions Documentation](https://docs.github.com/en/actions) 