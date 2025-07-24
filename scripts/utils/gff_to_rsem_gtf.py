#!/usr/bin/env python3

import sys
from collections import defaultdict

def fix_rsem_gtf(input_file, output_file):
    """Fix RSEM GTF file by ensuring all exons of a transcript are on the same strand."""
    print(f"Reading input file: {input_file}")
    
    # Track transcripts and their features
    transcript_features = defaultdict(list)
    transcript_strands = defaultdict(set)
    total_features = 0
    
    # Read input file
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            total_features += 1
            feature_type = fields[2]
            strand = fields[6]
            attributes = fields[8]
            
            # Parse attributes to get transcript_id
            transcript_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    transcript_id = attr.split('"')[1]
                    break
            
            if not transcript_id:
                continue
            
            # Track strand and feature
            transcript_strands[transcript_id].add(strand)
            transcript_features[transcript_id].append(line)
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"Total features processed: {total_features}")
    print(f"Total transcripts found: {len(transcript_features)}")
    print(f"Transcripts with inconsistent strands: {sum(1 for s in transcript_strands.values() if len(s) > 1)}")
    
    # Write output file
    print(f"\nWriting output file: {output_file}")
    written_features = 0
    skipped_transcripts = 0
    
    with open(output_file, 'w') as f:
        for transcript_id, strands in transcript_strands.items():
            if len(strands) > 1:
                print(f"Skipping transcript {transcript_id} with inconsistent strands: {strands}")
                skipped_transcripts += 1
                continue
                
            for line in transcript_features[transcript_id]:
                f.write(line)
                written_features += 1
    
    print(f"\nOutput Summary:")
    print(f"Total features written: {written_features}")
    print(f"Total transcripts skipped: {skipped_transcripts}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gff_to_rsem_gtf.py input.gtf output.gtf")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    fix_rsem_gtf(input_file, output_file) 