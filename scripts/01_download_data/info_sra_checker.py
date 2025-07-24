#!/usr/bin/env python3
"""
Simplified SRA Checker - Get all info from single API call
"""

import requests
import json
import csv
import time
import xml.etree.ElementTree as ET
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_project_info(project_id):
    """Get all SRA info for a project in one call"""
    logger.info(f"Getting info for {project_id}")
    
    # Search for project
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        'db': 'sra',
        'term': f'{project_id}[BioProject]',
        'retmode': 'json',
        'retmax': '1000'
    }
    
    time.sleep(1)
    response = requests.get(search_url, params=search_params, timeout=30)
    if response.status_code != 200:
        logger.error(f"Search failed: {response.status_code}")
        return []
    
    search_data = response.json()
    ids = search_data.get('esearchresult', {}).get('idlist', [])
    logger.info(f"Found {len(ids)} records")
    
    if not ids:
        return []
    
    # Get summary for all records
    time.sleep(2)
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    summary_params = {
        'db': 'sra',
        'id': ','.join(ids),
        'retmode': 'json'
    }
    
    response = requests.get(summary_url, params=summary_params, timeout=120)
    if response.status_code != 200:
        logger.error(f"Summary failed: {response.status_code}")
        return []
    
    summary_data = response.json()
    
    records = []
    for uid in ids:
        if uid in summary_data.get('result', {}):
            record = summary_data['result'][uid]
            
            # Parse the XML data that contains everything we need
            expxml = record.get('expxml', '')
            runs_xml = record.get('runs', '')
            
            # Extract experiment info from expxml
            sample_title = "Unknown"
            platform = "Unknown"
            layout = "Unknown"
            
            try:
                if expxml:
                    exp_root = ET.fromstring(f"<root>{expxml}</root>")
                    
                    # Get sample title from experiment title
                    title_elem = exp_root.find('.//Title')
                    if title_elem is not None:
                        sample_title = title_elem.text or "Unknown"
                    
                    # Get platform
                    platform_elem = exp_root.find('.//Platform')
                    if platform_elem is not None:
                        platform = platform_elem.text or "Unknown"
                    
                    # Get layout
                    paired_elem = exp_root.find('.//PAIRED')
                    single_elem = exp_root.find('.//SINGLE')
                    if paired_elem is not None:
                        layout = "PAIRED"
                    elif single_elem is not None:
                        layout = "SINGLE"
                    
            except Exception as e:
                logger.warning(f"Could not parse expxml for {uid}: {e}")
            
            # Extract run info from runs XML
            try:
                if runs_xml:
                    runs_root = ET.fromstring(f"<root>{runs_xml}</root>")
                    for run_elem in runs_root.findall('Run'):
                        acc = run_elem.get('acc', '')
                        spots = int(run_elem.get('total_spots', 0))
                        bases = int(run_elem.get('total_bases', 0))
                        
                        if acc and (acc.startswith('SRR') or acc.startswith('ERR') or acc.startswith('DRR')):
                            # Estimate size (bases * 0.25 bytes per base, roughly)
                            size_mb = (bases * 0.3) / (1024 * 1024)  # Rough estimate
                            avg_length = bases / spots if spots > 0 else 0
                            
                            # Parse condition information from sample title
                            condition_info = parse_condition_info(sample_title, project_id)
                            
                            record = {
                                'accession': acc,
                                'dataset': project_id,
                                'spots': spots,
                                'bases': bases,
                                'size_mb': f"{size_mb:.1f}",
                                'avg_length': f"{avg_length:.1f}",
                                'platform': platform,
                                'layout': layout,
                                'sample_title': sample_title,
                                'bioproject': project_id
                            }
                            
                            # Add condition-specific fields
                            record.update(condition_info)
                            
                            records.append(record)
                            
                            logger.info(f"Found: {acc} - {spots:,} spots, {bases:,} bases, ~{size_mb:.1f}MB")
                            
            except Exception as e:
                logger.warning(f"Could not parse runs XML for {uid}: {e}")
    
    return records

def parse_condition_info(sample_title, project_id):
    """Parse meaningful condition information from sample titles"""
    info = {
        'condition': 'unknown',
        'photoperiod': 'unknown', 
        'blood_meal': 'unknown',
        'replicate': 0,
        'comparison_group': 'unknown',
        'life_stage': 'unknown',
        'time_point': 'unknown'
    }
    
    title_lower = sample_title.lower()
    
    if project_id == "PRJNA268379":  # Adults
        info['life_stage'] = 'adult'
        
        # Parse photoperiod
        if 'ld' in title_lower or 'long' in title_lower:
            info['photoperiod'] = 'long_day'
            info['comparison_group'] = 'non_diapause_inducing'  # Long day = non-diapause inducing
        elif 'sd' in title_lower or 'short' in title_lower:
            info['photoperiod'] = 'short_day'
            info['comparison_group'] = 'diapause_inducing'  # Short day = diapause inducing
            
        # Parse blood meal
        if 'nb' in title_lower or 'no_blood' in title_lower:
            info['blood_meal'] = 'no_blood'
        elif 'bm' in title_lower or 'blood' in title_lower:
            info['blood_meal'] = 'blood_meal'
            
        # Parse replicate
        for i in range(1, 10):
            if f'_{i}' in title_lower or f'rep{i}' in title_lower or f'replicate_{i}' in title_lower:
                info['replicate'] = i
                break
                
        # Create condition code
        photoperiod_code = 'LD' if info['photoperiod'] == 'long_day' else 'SD'
        blood_code = 'BM' if info['blood_meal'] == 'blood_meal' else 'NB'
        info['condition'] = f"{photoperiod_code}_{blood_code}"
        
    elif project_id == "PRJNA158021":  # Embryos
        info['life_stage'] = 'embryo'
        
        # Parse diapause condition
        if 'non' in title_lower and 'diapause' in title_lower:
            info['condition'] = 'non_diapause'
            info['comparison_group'] = 'non_diapause_inducing'
        elif 'diapause' in title_lower:
            info['condition'] = 'diapause'
            info['comparison_group'] = 'diapause_inducing'
            
        # Parse time points
        if '72' in title_lower or '78' in title_lower:
            info['time_point'] = '72-78h'
        elif '135' in title_lower or '141' in title_lower:
            info['time_point'] = '135-141h'
            
        # Parse replicate
        for i in range(1, 10):
            if f'_{i}' in title_lower or f'rep{i}' in title_lower:
                info['replicate'] = i
                break
                
    elif project_id == "PRJNA187045":  # Pharate larvae
        info['life_stage'] = 'pharate_larvae'
        
        # Parse diapause condition - looking for "D " or "ND " in title
        if ' nd ' in title_lower or title_lower.startswith('nd '):
            info['condition'] = 'non_diapause'
            info['comparison_group'] = 'non_diapause_inducing'
            info['photoperiod'] = 'non_diapause_inducing'
        elif ' d ' in title_lower or (title_lower.startswith('d ') and not title_lower.startswith('nd')):
            info['condition'] = 'diapause'  
            info['comparison_group'] = 'diapause_inducing'
            info['photoperiod'] = 'diapause_inducing'
            
        # Parse time points - looking for "11d", "21d", "40d" in title
        if ' 11d ' in title_lower:
            info['time_point'] = '11_days_post_oviposition'
        elif ' 21d ' in title_lower:
            info['time_point'] = '21_days_post_oviposition'
        elif ' 40d ' in title_lower:
            info['time_point'] = '40_days_post_oviposition'
            
        # Parse replicate - looking for "rep 1", "rep 2", "rep 1b", etc.
        import re
        rep_match = re.search(r'rep (\d+[a-z]?)', title_lower)
        if rep_match:
            rep_str = rep_match.group(1)
            # Extract just the number part
            rep_num = re.search(r'\d+', rep_str)
            if rep_num:
                info['replicate'] = int(rep_num.group())
                
        # Create a more specific condition code
        condition_code = 'D' if info['condition'] == 'diapause' else 'ND'
        time_code = info['time_point'].split('_')[0] if info['time_point'] != 'unknown' else 'unk'
        info['condition'] = f"{condition_code}_{time_code}d"
    
    return info
    
    return records

def main():
    projects = ["PRJNA268379", "PRJNA158021", "PRJNA187045"]
    all_records = []
    
    for project in projects:
        records = get_project_info(project)
        all_records.extend(records)
        time.sleep(3)  # Be nice to NCBI
    
    # Write to CSV
    output_file = Path("data/sra_metadata.csv")
    output_file.parent.mkdir(exist_ok=True)
    
    if all_records:
        fieldnames = ['accession', 'dataset', 'spots', 'bases', 'size_mb', 'avg_length', 'platform', 'layout', 
                     'sample_title', 'bioproject', 'condition', 'photoperiod', 'blood_meal', 'replicate', 
                     'comparison_group', 'life_stage', 'time_point']
        
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_records)
        
        logger.info(f"Wrote {len(all_records)} records to {output_file}")
        
        # Print summary
        total_size = sum(float(r['size_mb']) for r in all_records)
        total_spots = sum(r['spots'] for r in all_records)
        
        print(f"\nSUMMARY:")
        print(f"Total accessions: {len(all_records)}")
        print(f"Estimated total size: {total_size/1024:.1f} GB")
        print(f"Total spots: {total_spots:,}")
        
        for project in projects:
            project_records = [r for r in all_records if r['dataset'] == project]
            project_size = sum(float(r['size_mb']) for r in project_records)
            print(f"{project}: {len(project_records)} files, ~{project_size/1024:.1f} GB")
    
    else:
        logger.error("No records found")

if __name__ == "__main__":
    main()