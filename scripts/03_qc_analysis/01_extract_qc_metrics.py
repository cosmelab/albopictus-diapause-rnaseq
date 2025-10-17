#!/usr/bin/env python3
"""
Extract QC metrics from nf-core RNA-seq outputs
Works directly with the output directory structure
"""

import pandas as pd
import json
import os
from pathlib import Path
import glob
import re

def extract_multiqc_stats(sample_path):
    """Extract metrics from MultiQC general stats"""
    stats_file = Path(sample_path) / "multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt"
    
    if not stats_file.exists():
        return None
        
    df = pd.read_csv(stats_file, sep='\t')
    
    # Get the main sample row (not _1 or _2)
    main_row = df[~df['Sample'].str.contains('_1|_2')].iloc[0]
    
    # Get read 1 and 2 rows for PE-specific metrics
    read1_rows = df[df['Sample'].str.endswith('_1')]
    read2_rows = df[df['Sample'].str.endswith('_2')]
    
    metrics = {
        # Basic info
        'sample_name': main_row['Sample'],
        
        # Sequencing metrics
        'total_reads_millions': float(main_row['star-total_reads']),
        'reads_mapped_percent': float(main_row['star-mapped_percent']),
        'uniquely_mapped_percent': float(main_row['star-uniquely_mapped_percent']),
        'multi_mapped_percent': float(main_row['star-mapped_percent'] - main_row['star-uniquely_mapped_percent']),
        
        # RNA-seq specific
        'exonic_percent': float(main_row['qualimap_rnaseq-reads_aligned_exonic']),
        'intronic_percent': float(main_row['qualimap_rnaseq-reads_aligned_intronic']),
        'intergenic_percent': float(main_row['qualimap_rnaseq-reads_aligned_intergenic']),
        'rrna_percent': float(main_row['custom_content_biotype_counts-percent_rRNA']) if 'custom_content_biotype_counts-percent_rRNA' in main_row else 0.0,
        'five_three_bias': float(main_row['qualimap_rnaseq-5_3_bias']),
        
        # Quality metrics
        'error_rate_percent': float(main_row['samtools_stats-error_rate']),
        'properly_paired_percent': float(main_row['samtools_stats-reads_properly_paired_percent']),
        'insert_size_mean': float(main_row['samtools_stats-insert_size_average']),
    }
    
    # Add PE-specific metrics if available
    if not read1_rows.empty and not read2_rows.empty:
        read1 = read1_rows.iloc[0]
        read2 = read2_rows.iloc[0]
        
        metrics.update({
            'read_length': float(read1['fastqc_raw-avg_sequence_length']),
            'gc_percent': (float(read1['fastqc_raw-percent_gc']) + float(read2['fastqc_raw-percent_gc'])) / 2,
            'duplication_percent': (float(read1['fastqc_raw-percent_duplicates']) + float(read2['fastqc_raw-percent_duplicates'])) / 2,
            'adapter_percent': (float(read1['cutadapt-percent_trimmed']) + float(read2['cutadapt-percent_trimmed'])) / 2,
        })
    
    return metrics

def extract_qualimap_metrics(sample_path):
    """Extract additional metrics from Qualimap RNA-seq QC"""
    qc_file = list(Path(sample_path).glob("star_salmon/qualimap/*/rnaseq_qc_results.txt"))
    
    if not qc_file:
        return {}
        
    metrics = {}
    with open(qc_file[0], 'r') as f:
        for line in f:
            if 'reads aligned' in line and '=' in line:
                value = line.split('=')[1].strip()
                if '(' in value:
                    number = value.split('(')[0].strip().replace(',', '')
                    metrics['total_aligned_reads'] = int(number)
            elif 'total alignments' in line and '=' in line:
                value = line.split('=')[1].strip().replace(',', '')
                metrics['total_alignments'] = int(value)
    
    return metrics

def extract_gene_counts_stats(sample_path):
    """Extract gene detection statistics from salmon counts"""
    counts_file = Path(sample_path) / "star_salmon/salmon.merged.gene_counts.tsv"
    
    if not counts_file.exists():
        return {}
        
    df = pd.read_csv(counts_file, sep='\t', index_col=0)
    
    # Get the sample column (should be only one)
    if df.shape[1] > 0:
        counts = pd.to_numeric(df.iloc[:, 0], errors='coerce').fillna(0)
        
        metrics = {
            'total_genes': int(len(counts)),
            'genes_detected_gt0': int((counts > 0).sum()),
            'genes_detected_gt10': int((counts > 10).sum()),
            'genes_detected_gt100': int((counts > 100).sum()),
            'genes_detected_gt1000': int((counts > 1000).sum()),
        }
    else:
        metrics = {}
    
    return metrics

def extract_strand_specificity(sample_path):
    """Extract strand specificity from RSeQC infer_experiment"""
    infer_file = list(Path(sample_path).glob("star_salmon/rseqc/infer_experiment/*.infer_experiment.txt"))
    
    if not infer_file:
        return {}
        
    metrics = {}
    with open(infer_file[0], 'r') as f:
        for line in f:
            if 'Fraction of reads explained by' in line:
                if '++,--' in line:
                    metrics['strand_forward_percent'] = float(line.split(':')[1].strip())
                elif '+-,-+' in line:
                    metrics['strand_reverse_percent'] = float(line.split(':')[1].strip())
    
    return metrics

def process_sample(sample_path, project_id, sample_id):
    """Process a single sample and extract all QC metrics"""
    
    print(f"Processing {project_id}/{sample_id}...")
    
    # Initialize metrics dictionary
    all_metrics = {
        'project_id': project_id,
        'sample_id': sample_id,
        'srr_id': sample_id,  # Assuming sample_id is the SRR ID
    }
    
    # Extract metrics from different sources
    multiqc_metrics = extract_multiqc_stats(sample_path)
    if multiqc_metrics:
        all_metrics.update(multiqc_metrics)
    
    qualimap_metrics = extract_qualimap_metrics(sample_path)
    all_metrics.update(qualimap_metrics)
    
    gene_stats = extract_gene_counts_stats(sample_path)
    all_metrics.update(gene_stats)
    
    strand_metrics = extract_strand_specificity(sample_path)
    all_metrics.update(strand_metrics)
    
    return all_metrics

def main():
    """Main function to process all samples"""
    
    # Use simplified directory structure
    output_dir = Path("results/01_qc_analysis/metrics")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "per_sample").mkdir(exist_ok=True)
    
    all_samples = []
    
    # Process each project
    for project_path in sorted(Path("output").glob("PRJNA*")):
        project_id = project_path.name
        
        # Process each sample in the project
        for sample_path in sorted(project_path.glob("SRR*")):
            sample_id = sample_path.name
            
            metrics = process_sample(sample_path, project_id, sample_id)
            
            # Save individual sample JSON
            output_file = output_dir / "per_sample" / f"{project_id}_{sample_id}_qc.json"
            with open(output_file, 'w') as f:
                json.dump(metrics, f, indent=2)
            
            all_samples.append(metrics)
    
    # Create combined dataframe
    df = pd.DataFrame(all_samples)
    
    # Save combined CSV
    df.to_csv(output_dir / "all_samples_qc_metrics.csv", index=False)
    
    print(f"\nProcessed {len(all_samples)} samples")
    print(f"Results saved to {output_dir}")
    
    # Print summary statistics
    print("\nSummary by project:")
    summary = df.groupby('project_id').agg({
        'total_reads_millions': ['mean', 'std', 'min', 'max'],
        'uniquely_mapped_percent': ['mean', 'std', 'min', 'max'],
        'exonic_percent': ['mean', 'std', 'min', 'max']
    }).round(2)
    print(summary)

if __name__ == "__main__":
    main()