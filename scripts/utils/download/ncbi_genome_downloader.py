#!/usr/bin/env python3
"""
Download Aedes albopictus reference genome and annotation from NCBI
and prepare for nf-core/rnaseq pipeline
"""
import os
import subprocess
import sys
import shutil
from pathlib import Path
import argparse
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_command(cmd, check=True):
    """Run shell command with error handling"""
    logger.info(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=check, 
                              capture_output=True, text=True)
        if result.stdout:
            logger.info(f"Output: {result.stdout.strip()}")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e}")
        logger.error(f"STDERR: {e.stderr}")
        if check:
            sys.exit(1)
        return e

def check_required_tools():
    """Check if required tools are available"""
    required_tools = ['datasets', 'gffread']
    missing_tools = []
    
    for tool in required_tools:
        try:
            result = subprocess.run([tool, '--version'], 
                                  capture_output=True, text=True)
            logger.info(f"✓ Found {tool}: {result.stdout.strip()}")
        except FileNotFoundError:
            logger.warning(f"✗ {tool} not found")
            missing_tools.append(tool)
    
    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}")
        logger.error("Please install them with:")
        logger.error("micromamba install -c conda-forge -c bioconda ncbi-datasets-cli gffread")
        return False
    
    return True

def download_genome_data(accession, output_dir):
    """Download genome data using NCBI datasets"""
    logger.info(f"Downloading genome data for {accession}")
    
    # Create output directory
    ref_dir = Path(output_dir) / "references"
    ref_dir.mkdir(parents=True, exist_ok=True)
    
    # Change to references directory
    original_dir = os.getcwd()
    os.chdir(ref_dir)
    
    try:
        # Download with datasets command
        cmd = f"datasets download genome accession {accession} --include gff3,genome,seq-report"
        run_command(cmd)
        
        # Unzip the download
        if Path("ncbi_dataset.zip").exists():
            run_command("unzip -o ncbi_dataset.zip")
            logger.info("Successfully unzipped genome data")
        else:
            logger.error("Download zip file not found!")
            return False
            
        return True
        
    finally:
        os.chdir(original_dir)

def organize_reference_files(accession, output_dir):
    """Organize downloaded files for nf-core/rnaseq"""
    ref_dir = Path(output_dir) / "references"
    ncbi_data_dir = ref_dir / "ncbi_dataset" / "data" / accession
    
    if not ncbi_data_dir.exists():
        logger.error(f"Expected NCBI data directory not found: {ncbi_data_dir}")
        return False
    
    logger.info("Organizing reference files...")
    
    # Find genome FASTA file
    genome_files = list(ncbi_data_dir.glob("*_genomic.fna"))
    if genome_files:
        genome_file = genome_files[0]
        target_genome = ref_dir / "genome.fa"
        shutil.copy2(genome_file, target_genome)
        logger.info(f"Copied genome: {genome_file.name} -> genome.fa")
    else:
        logger.error("Genome FASTA file not found!")
        return False
    
    # Find GFF annotation file
    gff_files = list(ncbi_data_dir.glob("genomic.gff"))
    if gff_files:
        gff_file = gff_files[0]
        target_gff = ref_dir / "annotation.gff3"
        shutil.copy2(gff_file, target_gff)
        logger.info(f"Copied annotation: {gff_file.name} -> annotation.gff3")
    else:
        logger.error("GFF annotation file not found!")
        return False
    
    # Copy sequence report if available
    seq_report_files = list(ncbi_data_dir.glob("sequence_report.jsonl"))
    if seq_report_files:
        seq_report = seq_report_files[0]
        target_report = ref_dir / "sequence_report.jsonl"
        shutil.copy2(seq_report, target_report)
        logger.info(f"Copied sequence report: {seq_report.name}")
    
    return True

def create_reference_info(accession, output_dir):
    """Create reference information file"""
    ref_dir = Path(output_dir) / "references"
    info_file = ref_dir / "reference_info.txt"
    
    with open(info_file, 'w') as f:
        f.write(f"# Aedes albopictus Reference Genome Information\n")
        f.write(f"# Downloaded from NCBI on: {subprocess.check_output(['date']).decode().strip()}\n")
        f.write(f"\n")
        f.write(f"Accession: {accession}\n")
        f.write(f"Species: Aedes albopictus\n")
        f.write(f"Source: NCBI GenBank\n")
        f.write(f"\n")
        f.write(f"Files:\n")
        f.write(f"- genome.fa: Reference genome FASTA\n")
        f.write(f"- annotation.gff3: Gene annotations (GFF3 format)\n")
        f.write(f"- sequence_report.jsonl: Assembly metadata\n")
        f.write(f"\n")
        f.write(f"For nf-core/rnaseq pipeline:\n")
        f.write(f"--fasta {ref_dir}/genome.fa --gff {ref_dir}/annotation.gff3\n")
    
    logger.info(f"Created reference info file: {info_file}")

def validate_files(output_dir):
    """Validate downloaded and organized files"""
    ref_dir = Path(output_dir) / "references"
    
    required_files = ["genome.fa", "annotation.gff3"]
    optional_files = ["sequence_report.jsonl"]
    
    logger.info("Validating downloaded files...")
    
    all_good = True
    for file_name in required_files:
        file_path = ref_dir / file_name
        if file_path.exists() and file_path.stat().st_size > 0:
            logger.info(f"✓ {file_name} - {file_path.stat().st_size:,} bytes")
        else:
            logger.error(f"✗ {file_name} - missing or empty")
            all_good = False
    
    for file_name in optional_files:
        file_path = ref_dir / file_name
        if file_path.exists():
            logger.info(f"✓ {file_name} - {file_path.stat().st_size:,} bytes")
        else:
            logger.info(f"? {file_name} - not available")
    
    return all_good

def main():
    parser = argparse.ArgumentParser(
        description="Download Aedes albopictus reference genome from NCBI"
    )
    parser.add_argument(
        "--accession", 
        default="GCA_018104305.1",
        help="NCBI genome accession (default: GCA_018104305.1)"
    )
    parser.add_argument(
        "--output-dir",
        default="data",
        help="Output directory (default: data)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files exist"
    )
    
    args = parser.parse_args()
    
    logger.info("🧬 NCBI Genome Reference Downloader")
    logger.info("=" * 50)
    logger.info(f"Accession: {args.accession}")
    logger.info(f"Output directory: {args.output_dir}")
    
    # Check required tools
    if not check_required_tools():
        logger.error("Please install missing tools first")
        sys.exit(1)
    
    # Check if files already exist
    ref_dir = Path(args.output_dir) / "references"
    genome_file = ref_dir / "genome.fa"
    annotation_file = ref_dir / "annotation.gff3"
    
    if genome_file.exists() and annotation_file.exists() and not args.force:
        logger.info("Reference files already exist. Use --force to re-download")
        if validate_files(args.output_dir):
            logger.info("✅ All reference files are ready!")
            return
    
    # Download genome data
    if not download_genome_data(args.accession, args.output_dir):
        logger.error("Failed to download genome data")
        sys.exit(1)
    
    # Organize files
    if not organize_reference_files(args.accession, args.output_dir):
        logger.error("Failed to organize reference files")
        sys.exit(1)
    
    # Create info file
    create_reference_info(args.accession, args.output_dir)
    
    # Validate final files
    if validate_files(args.output_dir):
        logger.info("✅ Reference genome download completed successfully!")
        logger.info(f"Files are ready in: {ref_dir}")
        
        # Print nf-core command suggestion
        logger.info("\n📋 For nf-core/rnaseq pipeline:")
        logger.info(f"--fasta {ref_dir.absolute()}/genome.fa --gff {ref_dir.absolute()}/annotation.gff3")
    else:
        logger.error("❌ Some files are missing or invalid")
        sys.exit(1)

if __name__ == "__main__":
    main()