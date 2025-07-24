#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR="data/raw"

# Check if directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Directory $FASTQ_DIR does not exist"
    exit 1
fi

# Initialize counters
total_size=0
file_count=0

echo "Checking FASTQ files in $FASTQ_DIR..."
echo "----------------------------------------"

# Find all FASTQ files and calculate sizes
while IFS= read -r file; do
    size=$(du -h "$file" | cut -f1)
    size_bytes=$(du -b "$file" | cut -f1)
    total_size=$((total_size + size_bytes))
    file_count=$((file_count + 1))
    echo "File: $file"
    echo "Size: $size"
    echo "----------------------------------------"
done < <(find "$FASTQ_DIR" -type f -name "*.fastq.gz" -o -name "*.fq.gz")

# Convert total size to human readable format
total_size_human=$(numfmt --to=iec-i --suffix=B --format="%.2f" $total_size)

echo "Summary:"
echo "Total number of FASTQ files: $file_count"
echo "Total size: $total_size_human"

# Estimate runtime
# Assuming:
# - 1GB of data takes ~5 minutes to process on a modern CPU
# - Using 8 cores for parallel processing
# - Including overhead for file I/O and other operations

# Convert to GB for calculation (1GB = 1024*1024*1024 bytes)
total_size_gb=$((total_size / 1024 / 1024 / 1024))
# Add 7% to convert from GiB to GB (1024^3 vs 1000^3)
total_size_gb=$((total_size_gb * 107 / 100))
# Calculate minutes: (total_size_gb * 5 minutes per GB) / 8 cores
estimated_minutes=$((total_size_gb * 5 / 8))
# Convert to hours
estimated_hours=$((estimated_minutes / 60))
# Calculate remaining minutes
remaining_minutes=$((estimated_minutes % 60))

echo "----------------------------------------"
echo "Runtime Estimation:"
echo "Total size in GB: $total_size_gb GB"
echo "Estimated processing time: $estimated_hours hours and $remaining_minutes minutes"
echo "Note: This is a rough estimate based on:"
echo "- 1GB of data takes ~5 minutes to process"
echo "- Using 8 cores for parallel processing"
echo "- Including overhead for file I/O and other operations"
echo "Actual runtime may vary based on:"
echo "- System specifications"
echo "- Available memory"
echo "- Disk I/O speed"
echo "- Network speed (if data is remote)"
echo "- Pipeline configuration" 