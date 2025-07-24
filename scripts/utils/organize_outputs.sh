#!/bin/bash

# Enhanced organization script combining both approaches
# Adds cross-project aggregation and manuscript-ready outputs

# Function to create directory if it doesn't exist
create_dir() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        echo "Created directory: $1"
    fi
}

# Function to copy files with verbose output
copy_file() {
    if [ -f "$1" ]; then
        cp -v "$1" "$2"
    else
        echo "Warning: Source file not found: $1"
    fi
}

# Function to aggregate MultiQC reports
aggregate_multiqc() {
    local project_dir="$1"
    local zenodo_dir="$2"
    local project_id=$(basename "$project_dir")
    
    echo "Aggregating QC reports for $project_id..."
    
    # Create MultiQC directory directly under project
    local multiqc_dir="$zenodo_dir/$project_id/multiqc"
    create_dir "$multiqc_dir"
    
    # Process each sample directory for MultiQC
    for sample_dir in "$project_dir"/SRR*; do
        if [ -d "$sample_dir" ]; then
            sample_id=$(basename "$sample_dir")
            echo "Processing QC files for sample: $sample_id"
            
            # Copy MultiQC files with prefix
            if [ -d "$sample_dir/multiqc/star_salmon/multiqc_report_data" ]; then
                for file in "$sample_dir/multiqc/star_salmon/multiqc_report_data/"*; do
                    filename=$(basename "$file")
                    cp "$file" "$multiqc_dir/${project_id}_${sample_id}_${filename}"
                done
            fi
        fi
    done
}

# Main script
echo "Starting enhanced organization of output files..."

# Create main zenodo directory
zenodo_dir="results/zenodo"
create_dir "$zenodo_dir"

# Process each project
for project_dir in results/PRJNA*; do
    if [ -d "$project_dir" ]; then
        project_id=$(basename "$project_dir")
        echo "Processing project: $project_id"
        
        # Create project directory structure
        create_dir "$zenodo_dir/$project_id/gene_counts"
        create_dir "$zenodo_dir/$project_id/gene_tpm"
        create_dir "$zenodo_dir/$project_id/transcript_counts"
        create_dir "$zenodo_dir/$project_id/transcript_tpm"
        create_dir "$zenodo_dir/$project_id/metadata"
        
        # Process each sample directory for quantification files
        for sample_dir in "$project_dir"/SRR*; do
            if [ -d "$sample_dir" ]; then
                sample_id=$(basename "$sample_dir")
                echo "Processing quantification files for sample: $sample_id"
                
                # Copy gene and transcript files with sample-specific names
                if [ -f "$sample_dir/star_salmon/salmon.merged.gene_counts.tsv" ]; then
                    cp "$sample_dir/star_salmon/salmon.merged.gene_counts.tsv" \
                       "$zenodo_dir/$project_id/gene_counts/${project_id}_${sample_id}_gene_counts.tsv"
                fi
                if [ -f "$sample_dir/star_salmon/salmon.merged.gene_tpm.tsv" ]; then
                    cp "$sample_dir/star_salmon/salmon.merged.gene_tpm.tsv" \
                       "$zenodo_dir/$project_id/gene_tpm/${project_id}_${sample_id}_gene_tpm.tsv"
                fi
                if [ -f "$sample_dir/star_salmon/salmon.merged.transcript_counts.tsv" ]; then
                    cp "$sample_dir/star_salmon/salmon.merged.transcript_counts.tsv" \
                       "$zenodo_dir/$project_id/transcript_counts/${project_id}_${sample_id}_transcript_counts.tsv"
                fi
                if [ -f "$sample_dir/star_salmon/salmon.merged.transcript_tpm.tsv" ]; then
                    cp "$sample_dir/star_salmon/salmon.merged.transcript_tpm.tsv" \
                       "$zenodo_dir/$project_id/transcript_tpm/${project_id}_${sample_id}_transcript_tpm.tsv"
                fi
            fi
        done
        
        # Copy metadata
        echo "Copying metadata files..."
        find . -name "sample_mapping_${project_id}.csv" -exec cp -v {} "$zenodo_dir/$project_id/metadata/" \;
        
        # Aggregate QC reports
        echo "Aggregating QC reports..."
        aggregate_multiqc "$project_dir" "$zenodo_dir"
        
        echo "Completed processing $project_id"
        echo "----------------------------------------"
    fi
done

# Display summary
echo ""
echo "=========================================="
echo "ORGANIZATION COMPLETE!"
echo "=========================================="
echo ""
echo "Data organized in:"
echo "  - $zenodo_dir/"
echo ""
echo "Next steps:"
echo "1. Review the organized data"
echo "2. Download to local machine for analysis"
echo "3. Run R scripts for manuscript figures"
echo "4. Compress and upload to Zenodo"
echo ""