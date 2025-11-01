#!/usr/bin/env python3
"""
Create Manuscript Tables for Aedes albopictus Diapause RNA-seq Reanalysis
"""

import pandas as pd
import os
from pathlib import Path
import sys

# Configuration
output_dir = "output/tables"
projects = ["PRJNA268379", "PRJNA158021", "PRJNA187045"]  # Ordered biologically: adult females → embryos → pharate larvae
samples_file = "data/metadata/samples.csv"

# Ensure output directory exists
Path(output_dir).mkdir(parents=True, exist_ok=True)

def validate_samples_df(df):
    """Validate the samples dataframe has required columns."""
    required_columns = ['BioProject', 'Run', 'Sample collection', 'Photoperiod', 'Description', 'Publication']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"❌ Missing required columns in samples file: {', '.join(missing_columns)}")
        sys.exit(1)
    return True

print("🧬 Creating manuscript tables...")

# 1. LOAD SAMPLE METADATA ============================================
print("📂 Loading sample metadata...")

if os.path.exists(samples_file):
    try:
        samples_df = pd.read_csv(samples_file)
        validate_samples_df(samples_df)
        print(f"  ✅ Loaded {len(samples_df)} samples from {samples_file}")
    except Exception as e:
        print(f"  ❌ Error reading {samples_file}: {str(e)}")
        sys.exit(1)
else:
    print(f"  ❌ {samples_file} not found!")
    sys.exit(1)

# 2. CREATE TABLE 1: EXPERIMENTAL DESIGN OVERVIEW ==================
print("📋 Creating Table 1: Experimental Design Overview...")

try:
    # Create a summary for each project
    table1_data = []
    
    for project in projects:
        project_df = samples_df[samples_df['BioProject'] == project]
        
        # Get unique conditions and their counts with abbreviations
        conditions = project_df.groupby(['Photoperiod', 'Description', 'New abreviation']).size().reset_index(name='Count')
        
        # Sort conditions by time points if they exist
        if '72h' in conditions['New abreviation'].iloc[0] or '11d' in conditions['New abreviation'].iloc[0]:
            # Extract time points for sorting
            conditions['time_point'] = conditions['New abreviation'].str.extract(r'(\d+[hd])')
            conditions['time_value'] = conditions['time_point'].str.extract(r'(\d+)').astype(int)
            conditions['time_unit'] = conditions['time_point'].str.extract(r'([hd])')
            # Convert hours to days for consistent sorting
            conditions['time_days'] = conditions.apply(
                lambda x: x['time_value'] if x['time_unit'] == 'd' else x['time_value']/24, 
                axis=1
            )
            conditions = conditions.sort_values(['Photoperiod', 'time_days'])
        
        conditions_str = ', '.join([
            f"{row['New abreviation']} ({row['Count']})" 
            for _, row in conditions.iterrows()
        ])
        
        # Get unique life stages
        life_stages = ', '.join(sorted(project_df['Sample collection'].unique()))
        
        # Get publication
        publication = project_df['Publication'].iloc[0]
        
        # Add to table data
        table1_data.append({
            'BioProject': project,
            'Samples from': life_stages,
            'Samples per Treatment': conditions_str,
            'Number of samples': len(project_df),
            'Publication': publication
        })
    
    # Create DataFrame and save
    table1 = pd.DataFrame(table1_data)
    table1.to_csv(f"{output_dir}/Table_1_Experimental_Design_Overview.csv", index=False)
    print(f"  ✅ Created Table 1 with {len(table1)} projects")
except Exception as e:
    print(f"  ❌ Error creating Table 1: {str(e)}")
    sys.exit(1)

# 3. CREATE TABLE S1: DETAILED SAMPLE INFORMATION ==================
print("📝 Creating Table S1: Detailed Sample Information...")

try:
    # Select and reorder columns, removing those with identical values
    columns_to_include = [
        'BioProject', 'Run', 'New abreviation', 'Sample collection',
        'Photoperiod', 'Description',
        'Replicate'
    ]
    
    # Create table
    table_s1 = samples_df[columns_to_include].copy()
    
    # Rename columns
    column_mapping = {
        'New abreviation': 'Sample ID',
        'Sample collection': 'Samples from',
        'Description': 'Treatment',
        'Photoperiod': 'Light Cycle'
    }
    table_s1 = table_s1.rename(columns=column_mapping)
    
    # Replace photoperiod values with SD/LD
    table_s1['Light Cycle'] = table_s1['Light Cycle'].replace({
        '8L:16D': 'SD',
        '16L:8D': 'LD'
    })
    
    # Sort by project (maintaining biological order) and photoperiod
    project_order = {project: i for i, project in enumerate(projects)}
    table_s1['project_order'] = table_s1['BioProject'].map(project_order)
    table_s1 = table_s1.sort_values(['project_order', 'Light Cycle', 'Run'])
    table_s1 = table_s1.drop('project_order', axis=1)
    
    # Save
    table_s1.to_csv(f"{output_dir}/Table_S1_Detailed_Sample_Information.csv", index=False)
    print(f"  ✅ Created Table S1 with {len(table_s1)} samples")
except Exception as e:
    print(f"  ❌ Error creating Table S1: {str(e)}")
    sys.exit(1)

# Create markdown file with table descriptions
with open(f"{output_dir}/table_descriptions.md", 'w') as f:
    f.write("""# Table Descriptions

## Table 1: Experimental Design Overview
Experimental Design Overview of the three RNA-seq studies re-analyzed using nf-core/rnaseq pipeline. The studies are ordered biologically from adult females to embryos to pharate larvae. Photoperiod conditions are indicated as SD (short day, 8L:16D) or LD (long day, 16L:8D). For adult females, samples were collected from fed and unfed conditions. Embryo samples were collected at 72h and 135h post-oviposition, while pharate larvae were sampled at 11, 21, and 40 days post-oviposition. The number of biological replicates is indicated in parentheses for each condition.

## Table S1: Detailed Sample Information
The aim of these experiments was to capture gene expression patterns at key moments of diapause: its onset in adult females, early embryonic development, and maintenance in pharate larvae. Adult females were maintained under short (SD) or long day (LD) photoperiod at 21°C since pupal stage, with samples collected after blood feeding during oocyte production. Embryonic samples were collected during early development (72h and 135h post-oviposition), while pharate larvae were sampled at different developmental stages (11, 21, and 40 days). Each sample is identified by its NCBI BioProject and SRA Run accession numbers.
""")

# 4. SUMMARY ==========================================================
print("\n" + "="*60)
print("📊 MANUSCRIPT TABLES SUMMARY")
print("="*60)
print(f"📁 Tables saved to: {output_dir}")
print("\n📄 Files created:")

try:
    for file in sorted(os.listdir(output_dir)):
        if file.endswith('.csv'):
            file_path = os.path.join(output_dir, file)
            df = pd.read_csv(file_path)
            print(f"  ✅ {file}")
            print(f"     └─ {len(df)} rows × {len(df.columns)} columns")

    print(f"\n📈 Dataset overview:")
    print(f"  • Total projects: {samples_df['BioProject'].nunique()}")
    print(f"  • Total samples: {len(samples_df)}")
    print(f"  • Photoperiod conditions: {', '.join(sorted(samples_df['Photoperiod'].unique()))}")
    print(f"  • Sampled life stages: {', '.join(sorted(samples_df['Sample collection'].unique()))}")

    print("\n🎉 All manuscript tables created successfully!")
    print("💡 Tables are ready for Excel formatting and manuscript inclusion")
except Exception as e:
    print(f"❌ Error generating summary: {e}")
    sys.exit(1)