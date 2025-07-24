#!/usr/bin/env python3
"""
SRA Data Download Script for Diapause RNA-seq Analysis
Downloads fastq files for three datasets: adults, embryos, and pharate larvae
Uses correct accessions from NCBI
"""

import os
import sys
import subprocess
import json
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import List, Dict, Optional
import argparse
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sra_download.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class SRAAccession:
    """Class to store SRA accession information"""
    accession: str
    dataset: str
    condition: str
    life_stage: str
    description: str

# Dataset definitions with CORRECT accessions from NCBI
DATASETS = {
    "PRJNA268379": {
        "name": "Adult_Females",
        "description": "Adult females - photoperiod and blood meal experiment",
        "accessions": [
            # Long day, no blood meal
            SRAAccession("SRR1664192", "PRJNA268379", "LD_NB", "Adult", "Ae_albo_ld_nb_4"),
            SRAAccession("SRR1664190", "PRJNA268379", "LD_NB", "Adult", "Ae_albo_ld_nb_3"),
            SRAAccession("SRR1663916", "PRJNA268379", "LD_NB", "Adult", "Ae_albo_ld_nb_2"),
            SRAAccession("SRR1663913", "PRJNA268379", "LD_NB", "Adult", "Ae_albo_ld_nb_1"),
            # Long day, blood meal
            SRAAccession("SRR1663911", "PRJNA268379", "LD_BM", "Adult", "Ae_albo_ld_bm_4"),
            SRAAccession("SRR1663843", "PRJNA268379", "LD_BM", "Adult", "Ae_albo_ld_bm_3"),
            SRAAccession("SRR1663769", "PRJNA268379", "LD_BM", "Adult", "Ae_albo_ld_bm_2"),
            SRAAccession("SRR1663754", "PRJNA268379", "LD_BM", "Adult", "Ae_albo_ld_bm_1"),
            # Short day, no blood meal
            SRAAccession("SRR1663709", "PRJNA268379", "SD_NB", "Adult", "Ae_albo_sd_nb_4"),
            SRAAccession("SRR1663707", "PRJNA268379", "SD_NB", "Adult", "Ae_albo_sd_nb_3"),
            SRAAccession("SRR1663703", "PRJNA268379", "SD_NB", "Adult", "Ae_albo_sd_nb_2"),
            SRAAccession("SRR1663700", "PRJNA268379", "SD_NB", "Adult", "Ae_albo_sd_nb_1"),
            # Short day, blood meal
            SRAAccession("SRR1663697", "PRJNA268379", "SD_BM", "Adult", "Ae_albo_sd_bm_4"),
            SRAAccession("SRR1663689", "PRJNA268379", "SD_BM", "Adult", "Ae_albo_sd_bm_3"),
            SRAAccession("SRR1663687", "PRJNA268379", "SD_BM", "Adult", "Ae_albo_sd_bm_2"),
            SRAAccession("SRR1663685", "PRJNA268379", "SD_BM", "Adult", "Ae_albo_sd_bm_1"),
        ]
    },
    "PRJNA158021": {
        "name": "Embryos",
        "description": "Embryonic diapause preparation",
        "accessions": [
            SRAAccession("SRR458462", "PRJNA158021", "DI_72h", "Embryo", "Diapause inducing, 72-78h post-oviposition, replicate 1"),
            SRAAccession("SRR458463", "PRJNA158021", "DI_135h", "Embryo", "Diapause inducing, 135-141h post-oviposition, replicate 1"),
            SRAAccession("SRR458464", "PRJNA158021", "NDI_72h", "Embryo", "Non-diapause inducing, 72-78h post-oviposition, replicate 1"),
            SRAAccession("SRR458465", "PRJNA158021", "NDI_135h", "Embryo", "Non-diapause inducing, 135-141h post-oviposition, replicate 1"),
            SRAAccession("SRR458466", "PRJNA158021", "DI_72h", "Embryo", "Diapause inducing, 72-78h post-oviposition, replicate 2"),
            SRAAccession("SRR458467", "PRJNA158021", "DI_135h", "Embryo", "Diapause inducing, 135-141h post-oviposition, replicate 2"),
            SRAAccession("SRR458468", "PRJNA158021", "NDI_72h", "Embryo", "Non-diapause inducing, 72-78h post-oviposition, replicate 2"),
            SRAAccession("SRR458469", "PRJNA158021", "NDI_135h", "Embryo", "Non-diapause inducing, 135-141h post-oviposition, replicate 2"),
            SRAAccession("SRR458470", "PRJNA158021", "DI_72h", "Embryo", "Diapause inducing, 72-78h post-oviposition, replicate 3"),
            SRAAccession("SRR458471", "PRJNA158021", "DI_135h", "Embryo", "Diapause inducing, 135-141h post-oviposition, replicate 3"),
            SRAAccession("SRR458472", "PRJNA158021", "NDI_72h", "Embryo", "Non-diapause inducing, 72-78h post-oviposition, replicate 3"),
            SRAAccession("SRR458473", "PRJNA158021", "NDI_135h", "Embryo", "Non-diapause inducing, 135-141h post-oviposition, replicate 3"),
        ]
    },
    "PRJNA187045": {
        "name": "Pharate_Larvae",
        "description": "Pharate larvae - diapause maintenance",
        "accessions": [
            # Diapause-inducing conditions
            SRAAccession("SRR652072", "PRJNA187045", "D_11d", "Pharate_larvae", "Ae. albopictus D 11d pharate larvae, rep 1b"),
            SRAAccession("SRR652088", "PRJNA187045", "D_11d", "Pharate_larvae", "Ae. albopictus D 11d pharate larvae, rep 4"),
            SRAAccession("SRR652090", "PRJNA187045", "D_11d", "Pharate_larvae", "Ae. albopictus D 11d pharate larvae, rep 2"),
            SRAAccession("SRR652091", "PRJNA187045", "D_21d", "Pharate_larvae", "Ae. albopictus D 21d pharate larvae, rep 1"),
            SRAAccession("SRR652092", "PRJNA187045", "D_21d", "Pharate_larvae", "Ae. albopictus D 21d pharate larvae, rep 2"),
            SRAAccession("SRR652097", "PRJNA187045", "D_21d", "Pharate_larvae", "Ae. albopictus D 21d pharate larvae, rep 4"),
            SRAAccession("SRR652098", "PRJNA187045", "D_40d", "Pharate_larvae", "Ae. albopictus D 40d pharate larvae, rep 2"),
            SRAAccession("SRR652099", "PRJNA187045", "D_40d", "Pharate_larvae", "Ae. albopictus D 40d pharate larvae, rep 3"),
            SRAAccession("SRR652100", "PRJNA187045", "D_40d", "Pharate_larvae", "Ae. albopictus D 40d pharate larvae, rep 4"),
            # Non-diapause inducing conditions
            SRAAccession("SRR652101", "PRJNA187045", "ND_11d", "Pharate_larvae", "Ae. albopictus ND 11d pharate larvae, rep 1"),
            SRAAccession("SRR652102", "PRJNA187045", "ND_11d", "Pharate_larvae", "Ae. albopictus ND 11d pharate larvae, rep 2"),
            SRAAccession("SRR652103", "PRJNA187045", "ND_11d", "Pharate_larvae", "Ae. albopictus ND 11d pharate larvae, rep 4"),
            SRAAccession("SRR652116", "PRJNA187045", "ND_21d", "Pharate_larvae", "Ae. albopictus ND 21d pharate larvae, rep 1"),
            SRAAccession("SRR652118", "PRJNA187045", "ND_21d", "Pharate_larvae", "Ae. albopictus ND 21d pharate larvae, rep 2"),
            SRAAccession("SRR652119", "PRJNA187045", "ND_21d", "Pharate_larvae", "Ae. albopictus ND 21d pharate larvae, rep 3"),
            SRAAccession("SRR652120", "PRJNA187045", "ND_40d", "Pharate_larvae", "Ae. albopictus ND 40d pharate larvae, rep 3"),
            SRAAccession("SRR672801", "PRJNA187045", "ND_40d", "Pharate_larvae", "Ae. albopictus ND 40d pharate larvae, rep 4"),
        ]
    }
}

class SRADownloader:
    """Main class for downloading SRA data"""
    
    def __init__(self, base_dir: str = "data", max_workers: int = 4):
        self.base_dir = Path(base_dir)
        self.max_workers = max_workers
        self.failed_downloads = []
        self.successful_downloads = []
        
        # Configure SRA cache to use data/sra directory
        self.configure_sra_cache()
        
    def configure_sra_cache(self):
        """Configure SRA toolkit to use data/sra for cache"""
        sra_cache_dir = self.base_dir / "sra"
        sra_cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Skip vdb-config if running as root to avoid warnings
        if os.getuid() == 0:
            logger.info(f"Running as root - skipping vdb-config, using --output-directory instead")
            logger.info(f"SRA cache will be directed to: {sra_cache_dir}")
            return
        
        try:
            # Set SRA cache directory for non-root users
            result = subprocess.run([
                "vdb-config", "--set", 
                f"/repository/user/main/public/root={sra_cache_dir.absolute()}"
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"SRA cache configured to use: {sra_cache_dir}")
            else:
                logger.warning(f"Could not configure SRA cache: {result.stderr}")
        except Exception as e:
            logger.warning(f"SRA cache configuration failed: {e}")

        
    def setup_directories(self):
        """Use existing directory structure from setup.sh"""
        # Ensure raw data directories exist for each dataset
        for dataset_id in DATASETS.keys():
            raw_dir = self.base_dir / "raw" / dataset_id
            raw_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"FASTQ files will go to: {raw_dir}")
        
        # Ensure SRA cache directory exists
        sra_dir = self.base_dir / "sra"
        sra_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"SRA cache files will go to: {sra_dir}")
        
    def check_sra_tools(self) -> bool:
        """Check if SRA toolkit is available"""
        try:
            result = subprocess.run(["fastq-dump", "--version"], 
                                  capture_output=True, text=True, timeout=10)
            logger.info(f"SRA toolkit version: {result.stdout.strip()}")
            return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            logger.error("SRA toolkit not found. Please install sra-tools.")
            return False
    
    def download_single_accession(self, accession: SRAAccession, fastq_output_dir: Path, test_mode: bool = False) -> bool:
        """Download a single SRA accession"""
        logger.info(f"Starting download: {accession.accession} ({accession.description})")
        
        # Ensure FASTQ output directory exists
        fastq_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if FASTQ files already exist - handle mixed single/paired-end properly
        single_end_file = fastq_output_dir / f"{accession.accession}.fastq.gz"
        file_1 = fastq_output_dir / f"{accession.accession}_1.fastq.gz"
        file_2 = fastq_output_dir / f"{accession.accession}_2.fastq.gz"
        
        # Skip if ANY valid FASTQ files exist for this accession
        if (single_end_file.exists() and single_end_file.stat().st_size > 0) or \
           (file_1.exists() and file_1.stat().st_size > 0):
            logger.info(f"FASTQ files already exist for {accession.accession}, skipping...")
            return True
        
        try:
            # Step 1: Prefetch SRA file (force it to go to cache directory)
            logger.info(f"Prefetching SRA file for {accession.accession}")
            sra_cache_dir = self.base_dir / "sra"
            
            prefetch_cmd = [
                "prefetch", 
                "--output-directory", str(sra_cache_dir),
                accession.accession
            ]
            
            result = subprocess.run(prefetch_cmd, capture_output=True, text=True, timeout=7200)
            if result.returncode != 0:
                logger.error(f"Prefetch failed for {accession.accession}: {result.stderr}")
                return False
            
            # Step 2: Convert SRA to FASTQ (output goes to data/raw/)
            logger.info(f"Converting SRA to FASTQ for {accession.accession}")
            
            # Point to the SRA file in cache (handle both .sra and .sralite)
            sra_cache_path = sra_cache_dir / accession.accession
            sra_file_path = None
            
            # Look for .sra file first, then .sralite
            if (sra_cache_path / f"{accession.accession}.sra").exists():
                sra_file_path = sra_cache_path / f"{accession.accession}.sra"
            elif (sra_cache_path / f"{accession.accession}.sralite").exists():
                sra_file_path = sra_cache_path / f"{accession.accession}.sralite"
                logger.info(f"Found .sralite file, will try fastq-dump anyway")
            
            if not sra_file_path:
                logger.error(f"No SRA file found for {accession.accession}")
                return False
            
            fastq_cmd = [
                "fastq-dump", 
                "--split-files",  # Create separate files for forward/reverse reads
                "--gzip",         # Compress output files
                "--outdir", str(fastq_output_dir),  # Output FASTQ files to data/raw/
            ]
            
            # Add test mode option for small files (only first 100 reads)
            if test_mode:
                fastq_cmd.extend(["--maxSpotId", "100"])
                logger.info(f"TEST MODE: Only converting first 100 reads")
            
            # Add the SRA file path
            fastq_cmd.append(str(sra_file_path))
            
            result = subprocess.run(fastq_cmd, capture_output=True, text=True, timeout=7200)
            
            if result.returncode != 0:
                logger.error(f"fastq-dump failed for {accession.accession}: {result.stderr}")
                return False
            
            # Step 3: Check what FASTQ files were actually created
            actual_files = list(fastq_output_dir.glob(f"{accession.accession}*.fastq.gz"))
            
            if not actual_files:
                logger.error(f"No FASTQ files created for {accession.accession}")
                return False
            
            # Determine if single-end or paired-end based on what was created
            single_end_file = fastq_output_dir / f"{accession.accession}.fastq.gz"
            file_1 = fastq_output_dir / f"{accession.accession}_1.fastq.gz"
            file_2 = fastq_output_dir / f"{accession.accession}_2.fastq.gz"
            
            if single_end_file.exists() and single_end_file.stat().st_size > 0:
                # Single-end sequencing
                logger.info(f"Single-end sequencing detected for {accession.accession}")
                file_info = f"{single_end_file.name}: {single_end_file.stat().st_size / 1e6:.1f} MB"
                logger.info(f"Successfully created FASTQ file for {accession.accession}: {file_info}")
                
            elif file_1.exists() and file_1.stat().st_size > 0:
                # Paired-end sequencing
                if file_2.exists() and file_2.stat().st_size > 0:
                    logger.info(f"Paired-end sequencing detected for {accession.accession}")
                    file_info = []
                    for f in [file_1, file_2]:
                        size_mb = f.stat().st_size / 1e6
                        file_info.append(f"{f.name}: {size_mb:.1f} MB")
                    logger.info(f"Successfully created FASTQ files for {accession.accession}: {', '.join(file_info)}")
                else:
                    logger.info(f"Single-end sequencing detected for {accession.accession} (only _1 file created)")
                    file_info = f"{file_1.name}: {file_1.stat().st_size / 1e6:.1f} MB"
                    logger.info(f"Successfully created FASTQ file for {accession.accession}: {file_info}")
            else:
                logger.error(f"No valid FASTQ files found for {accession.accession}")
                return False
            
            logger.info(f"Files location: {fastq_output_dir}")
            return True
                
        except subprocess.TimeoutExpired:
            logger.error(f"Timeout downloading {accession.accession}")
            return False
        except Exception as e:
            logger.error(f"Error downloading {accession.accession}: {str(e)}")
            return False
    
    def download_dataset(self, dataset_id: str, test_mode: bool = False):
        """Download all accessions for a dataset"""
        if dataset_id not in DATASETS:
            logger.error(f"Unknown dataset: {dataset_id}")
            return
        
        dataset_info = DATASETS[dataset_id]
        # FASTQ files go to data/raw/DATASET_ID/
        fastq_output_dir = self.base_dir / "raw" / dataset_id
        accessions = dataset_info["accessions"]
        
        if test_mode:
            accessions = accessions[:1]  # Only download first 1 for testing
            logger.info(f"TEST MODE: Only downloading first 1 accession with 100 reads")
        
        logger.info(f"Downloading dataset {dataset_id}")
        logger.info(f"SRA cache directory: {self.base_dir / 'sra'}")
        logger.info(f"FASTQ output directory: {fastq_output_dir}")
        logger.info(f"Number of accessions: {len(accessions)}")
        
        # Sequential downloads to avoid overwhelming the connection
        for acc in accessions:
            try:
                success = self.download_single_accession(acc, fastq_output_dir, test_mode)
                if success:
                    self.successful_downloads.append(acc.accession)
                else:
                    self.failed_downloads.append(acc.accession)
            except Exception as e:
                logger.error(f"Exception downloading {acc.accession}: {str(e)}")
                self.failed_downloads.append(acc.accession)
        
        # Generate metadata for this specific dataset
        self.generate_metadata(dataset_id)
    
    def generate_metadata(self, dataset_id: str = None):
        """Generate metadata file for downloaded samples"""
        if dataset_id:
            # Generate metadata for specific dataset only
            datasets_to_process = {dataset_id: DATASETS[dataset_id]}
            metadata_file = self.base_dir / f"sample_metadata_{dataset_id}.json"
            csv_file = self.base_dir / f"sample_mapping_{dataset_id}.csv"
        else:
            # Generate metadata for all datasets
            datasets_to_process = DATASETS
            metadata_file = self.base_dir / "sample_metadata_all.json"
            csv_file = self.base_dir / "sample_mapping_all.csv"
        
        metadata = []
        csv_lines = ["accession,dataset,condition,life_stage,description,fastq_files,sequencing_type"]
        
        for ds_id, dataset_info in datasets_to_process.items():
            fastq_dir = self.base_dir / "raw" / ds_id
            
            for acc in dataset_info["accessions"]:
                # Check what files actually exist for this accession
                single_end_file = fastq_dir / f"{acc.accession}.fastq.gz"
                file_1 = fastq_dir / f"{acc.accession}_1.fastq.gz"
                file_2 = fastq_dir / f"{acc.accession}_2.fastq.gz"
                
                entry = {
                    "accession": acc.accession,
                    "dataset": acc.dataset,
                    "condition": acc.condition, 
                    "life_stage": acc.life_stage,
                    "description": acc.description,
                }
                
                # Add file paths based on what actually exists
                if single_end_file.exists():
                    entry["fastq"] = f"{dataset_info['name']}/{acc.accession}.fastq.gz"
                    entry["sequencing_type"] = "single-end"
                    fastq_files = f"{acc.accession}.fastq.gz"
                elif file_1.exists():
                    entry["fastq_1"] = f"{dataset_info['name']}/{acc.accession}_1.fastq.gz"
                    entry["sequencing_type"] = "single-end" if not file_2.exists() else "paired-end"
                    if file_2.exists():
                        entry["fastq_2"] = f"{dataset_info['name']}/{acc.accession}_2.fastq.gz"
                        fastq_files = f"{acc.accession}_1.fastq.gz;{acc.accession}_2.fastq.gz"
                    else:
                        fastq_files = f"{acc.accession}_1.fastq.gz"
                else:
                    entry["sequencing_type"] = "not_downloaded"
                    fastq_files = "not_downloaded"
                
                metadata.append(entry)
                
                # Add to CSV
                csv_lines.append(f"{acc.accession},{acc.dataset},{acc.condition},{acc.life_stage},{acc.description},{fastq_files},{entry['sequencing_type']}")
        
        # Write JSON metadata
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Write CSV mapping
        with open(csv_file, 'w') as f:
            f.write('\n'.join(csv_lines))
        
        logger.info(f"Metadata saved to {metadata_file}")
        logger.info(f"CSV mapping saved to {csv_file}")
    
    def print_summary(self):
        """Print download summary"""
        total = len(self.successful_downloads) + len(self.failed_downloads)
        logger.info("\n" + "="*50)
        logger.info("DOWNLOAD SUMMARY")
        logger.info("="*50)
        logger.info(f"Total accessions processed: {total}")
        logger.info(f"Successful downloads: {len(self.successful_downloads)}")
        logger.info(f"Failed downloads: {len(self.failed_downloads)}")
        
        if self.failed_downloads:
            logger.info("\nFailed accessions:")
            for acc in self.failed_downloads:
                logger.info(f"  - {acc}")
        
        logger.info("="*50)

def main():
    parser = argparse.ArgumentParser(description="Download SRA data for diapause RNA-seq analysis")
    parser.add_argument("--dataset", choices=list(DATASETS.keys()) + ["all"], 
                       default="all", help="Dataset to download")
    parser.add_argument("--output-dir", default="data", 
                       help="Output directory for downloads")
    parser.add_argument("--max-workers", type=int, default=4,
                       help="Maximum number of parallel downloads")
    parser.add_argument("--test", action="store_true",
                       help="Test mode - only download first 2 accessions")
    
    args = parser.parse_args()
    
    # Initialize downloader
    downloader = SRADownloader(base_dir=args.output_dir, max_workers=args.max_workers)
    
    # Check SRA tools
    if not downloader.check_sra_tools():
        sys.exit(1)
    
    # Setup directories
    downloader.setup_directories()
    
    # Download data
    if args.dataset == "all":
        for dataset_id in DATASETS.keys():
            downloader.download_dataset(dataset_id, test_mode=args.test)
        # Also generate combined metadata
        downloader.generate_metadata()
    else:
        downloader.download_dataset(args.dataset, test_mode=args.test)
    
    # Print summary
    downloader.print_summary()

if __name__ == "__main__":
    main()