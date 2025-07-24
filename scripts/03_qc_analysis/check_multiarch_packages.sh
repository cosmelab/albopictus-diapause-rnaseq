#!/bin/bash

# Multi-architecture package compatibility checker
# Quick check for RNA-seq analysis packages

echo "🔍 Multi-architecture Package Compatibility Checker"
echo "=================================================="

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to check conda package
check_conda_package() {
    local package=$1
    local channel=${2:-"conda-forge"}
    
    echo -n "  Checking ${channel}::${package}... "
    
    # Check x86_64
    x86_result=$(conda search "${channel}::${package}" --platform linux-64 --json 2>/dev/null | jq -r ".[\"${package}\"] | length" 2>/dev/null || echo "0")
    
    # Check ARM64
    arm_result=$(conda search "${channel}::${package}" --platform linux-aarch64 --json 2>/dev/null | jq -r ".[\"${package}\"] | length" 2>/dev/null || echo "0")
    
    if [ "$x86_result" -gt 0 ] && [ "$arm_result" -gt 0 ]; then
        echo -e "${GREEN}✅ Both architectures${NC}"
        return 0
    elif [ "$x86_result" -gt 0 ]; then
        echo -e "${YELLOW}⚠️  x86_64 only${NC}"
        return 1
    elif [ "$arm_result" -gt 0 ]; then
        echo -e "${YELLOW}⚠️  ARM64 only${NC}"
        return 2
    else
        echo -e "${RED}❌ Not found${NC}"
        return 3
    fi
}

# Initialize counters
total=0
compatible=0
x86_only=0
arm_only=0
not_found=0

echo -e "\n📦 Checking Conda packages (conda-forge)..."
conda_packages=(
    "python"
    "jupyter"
    "jupyterlab"
    "notebook"
    "cmake"
    "make"
    "gcc"
    "gxx"
    "eza"
    "starship"
    "openjdk"
    "datamash"
    "r-base"
)

for package in "${conda_packages[@]}"; do
    check_conda_package "$package" "conda-forge"
    case $? in
        0) ((compatible++));;
        1) ((x86_only++));;
        2) ((arm_only++));;
        3) ((not_found++));;
    esac
    ((total++))
done

echo -e "\n🧬 Checking Bioconductor packages..."
bioc_packages=(
    "bioconductor-deseq2"
    "bioconductor-tximport"
    "bioconductor-annotationdbi"
    "bioconductor-biomart"
    "bioconductor-sva"
    "sra-tools"
    "entrez-direct"
    "bcftools"
    "samtools"
    "bedtools"
    "tabix"
    "vcftools"
    "snakemake"
)

for package in "${bioc_packages[@]}"; do
    check_conda_package "$package" "bioconda"
    case $? in
        0) ((compatible++));;
        1) ((x86_only++));;
        2) ((arm_only++));;
        3) ((not_found++));;
    esac
    ((total++))
done

echo -e "\n📊 Checking R packages (CRAN via conda-forge)..."
r_packages=(
    "r-here"
    "r-data.table"
    "r-metafor"
    "r-tidyverse"
    "r-ggplot2"
    "r-qqman"
    "r-qqplotr"
    "r-reticulate"
    "r-broom"
    "r-readxl"
    "r-writexl"
    "r-knitr"
    "r-rmarkdown"
)

for package in "${r_packages[@]}"; do
    check_conda_package "$package" "conda-forge"
    case $? in
        0) ((compatible++));;
        1) ((x86_only++));;
        2) ((arm_only++));;
        3) ((not_found++));;
    esac
    ((total++))
done

echo -e "\n🐍 Checking Python packages (PyPI)..."
python_packages=(
    "pandas"
    "numpy"
    "scipy"
    "matplotlib"
    "seaborn"
    "upsetplot"
    "radian"
    "pygments"
    "prompt-toolkit"
    "biopython"
    "scikit-learn"
    "requests"
    "beautifulsoup4"
    "xmltodict"
    "lxml"
)

for package in "${python_packages[@]}"; do
    check_conda_package "$package" "conda-forge"
    case $? in
        0) ((compatible++));;
        1) ((x86_only++));;
        2) ((arm_only++));;
        3) ((not_found++));;
    esac
    ((total++))
done

echo -e "\n=================================================="
echo -e "📋 SUMMARY"
echo -e "=================================================="
echo -e "Total packages: ${total}"
echo -e "Multi-arch compatible: ${compatible} ($(echo "scale=1; ${compatible}*100/${total}" | bc)%)"
echo -e "x86_64 only: ${x86_only}"
echo -e "ARM64 only: ${arm_only}"
echo -e "Not found: ${not_found}"

if [ $compatible -eq $total ]; then
    echo -e "\n${GREEN}🎉 All packages are multi-architecture compatible!${NC}"
else
    echo -e "\n${YELLOW}⚠️  $((total - compatible)) packages may have architecture issues.${NC}"
    echo -e "   Consider checking these packages manually or finding alternatives."
fi

echo -e "\n💡 Quick test commands for your Docker image:"
echo -e "   # Test x86_64:"
echo -e "   docker run --platform linux/amd64 cosmelab/albopictus-diapause-rnaseq:latest python -c \"import pandas; print('x86_64 works')\""
echo -e "   # Test ARM64:"
echo -e "   docker run --platform linux/arm64 cosmelab/albopictus-diapause-rnaseq:latest python -c \"import pandas; print('ARM64 works')\"" 