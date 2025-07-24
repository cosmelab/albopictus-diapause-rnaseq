import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style and create output directory
plt.style.use('default')
sns.set_palette("husl")
Path("output/figures").mkdir(parents=True, exist_ok=True)

# Read the data
df = pd.read_csv("output/tables/Table_S2_QC_Metrics.csv")

# Define sample order
sample_order = [
    # PRJNA268379 (Fed/Unfed)
    'SD_Unfed', 'LD_Unfed', 'SD_Fed', 'LD_Fed',
    # PRJNA158021 (72h/135h)
    'SD_72h', 'LD_72h', 'SD_135h', 'LD_135h',
    # PRJNA187045 (11d/21d/40d)
    'SD_11d', 'LD_11d', 'SD_21d', 'LD_21d', 'SD_40d', 'LD_40d'
]

# Create a mapping for sorting
def get_sample_order(sample_id):
    for i, pattern in enumerate(sample_order):
        if pattern in sample_id:
            return i
    return len(sample_order)  # Put unmatched samples at the end

# Sort the dataframe
df['sort_order'] = df['Sample ID'].apply(get_sample_order)
df = df.sort_values('sort_order')

# 1. RNA-seq Feature Distribution
plt.figure(figsize=(15, 8))
metrics = ['Exonic_Reads (%)', 'Intronic_Reads (%)', 'Intergenic_Reads (%)']
colors = ['#2ecc71', '#f1c40f', '#e74c3c']  # Green, Yellow, Red

# Plot stacked bars
bottom = np.zeros(len(df))
for metric, color in zip(metrics, colors):
    plt.bar(range(len(df)), df[metric], bottom=bottom, label=metric.replace('_', ' '), color=color)
    bottom += df[metric]

plt.title('RNA-seq Feature Distribution Across Samples', pad=20, fontsize=14)
plt.xlabel('Samples', labelpad=10, fontsize=12)
plt.ylabel('Percentage (%)', labelpad=10, fontsize=12)
plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=3)
plt.grid(True, axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('output/figures/QC_Feature_Distribution.pdf', format='pdf', bbox_inches='tight')
plt.savefig('output/figures/QC_Feature_Distribution.eps', format='eps', bbox_inches='tight')
plt.savefig('output/figures/QC_Feature_Distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Mapping Statistics
plt.figure(figsize=(15, 8))
mapping_metrics = ['Uniquely_Mapped (%)', 'Multi_Mapped (%)']
mapping_colors = ['#3498db', '#e67e22']  # Blue, Orange

bottom = np.zeros(len(df))
for metric, color in zip(mapping_metrics, mapping_colors):
    plt.bar(range(len(df)), df[metric], bottom=bottom, label=metric.replace('_', ' '), color=color)
    bottom += df[metric]

plt.title('Read Mapping Distribution Across Samples', pad=20, fontsize=14)
plt.xlabel('Samples', labelpad=10, fontsize=12)
plt.ylabel('Percentage (%)', labelpad=10, fontsize=12)
plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=2)
plt.grid(True, axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('output/figures/QC_Mapping_Distribution.pdf', format='pdf', bbox_inches='tight')
plt.savefig('output/figures/QC_Mapping_Distribution.eps', format='eps', bbox_inches='tight')
plt.savefig('output/figures/QC_Mapping_Distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. Quality Metrics
plt.figure(figsize=(15, 8))
quality_metrics = ['Duplication_Rate (%)', 'Adapter_Contamination (%)', 'Error_Rate (%)']
quality_colors = ['#9b59b6', '#34495e', '#e74c3c']  # Purple, Dark Blue, Red

for metric, color in zip(quality_metrics, quality_colors):
    plt.plot(range(len(df)), df[metric], marker='o', label=metric.replace('_', ' '), color=color)

plt.title('Quality Metrics Across Samples', pad=20, fontsize=14)
plt.xlabel('Samples', labelpad=10, fontsize=12)
plt.ylabel('Percentage (%)', labelpad=10, fontsize=12)
plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.legend(bbox_to_anchor=(0.5, 1.15), loc='upper center', ncol=3)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('output/figures/QC_Quality_Metrics.pdf', format='pdf', bbox_inches='tight')
plt.savefig('output/figures/QC_Quality_Metrics.eps', format='eps', bbox_inches='tight')
plt.savefig('output/figures/QC_Quality_Metrics.png', dpi=300, bbox_inches='tight')
plt.close()

# 4. Sequencing Depth and Insert Size
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12), sharex=True)

# Total Reads
ax1.bar(range(len(df)), df['Total_Reads (M)'], color='#2ecc71')
ax1.set_title('Sequencing Depth Across Samples', pad=20, fontsize=14)
ax1.set_ylabel('Total Reads (M)', labelpad=10, fontsize=12)
ax1.grid(True, axis='y', alpha=0.3)

# Insert Size
ax2.bar(range(len(df)), df['Insert_Size (bp)'], color='#3498db')
ax2.set_title('Insert Size Distribution Across Samples', pad=20, fontsize=14)
ax2.set_xlabel('Samples', labelpad=10, fontsize=12)
ax2.set_ylabel('Insert Size (bp)', labelpad=10, fontsize=12)
ax2.grid(True, axis='y', alpha=0.3)

plt.xticks(range(len(df)), df['Sample ID'], rotation=45, ha='right')
plt.tight_layout()
plt.savefig('output/figures/QC_Depth_and_InsertSize.pdf', format='pdf', bbox_inches='tight')
plt.savefig('output/figures/QC_Depth_and_InsertSize.eps', format='eps', bbox_inches='tight')
plt.savefig('output/figures/QC_Depth_and_InsertSize.png', dpi=300, bbox_inches='tight')
plt.close()