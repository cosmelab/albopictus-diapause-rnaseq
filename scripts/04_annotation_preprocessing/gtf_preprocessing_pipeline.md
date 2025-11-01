# GTF Preprocessing Pipeline Documentation

## Overview

The AalbF3 annotation from Dryad required extensive preprocessing to be compatible with nf-core/rnaseq pipeline. This document details the issues encountered and solutions implemented.

---

## Issues Identified and Fixed

### 1. File Format Conversion (GFF3 → GTF)

**Problem:** Dryad provides GFF3 format, but nf-core/rnaseq requires GTF format.

**Solution:**
```bash
gffread AalbF3.gff3 -T -o AalbF3.gtf
```

### 2. Missing gene_id Attributes (18,135 lines)

**Problem:** Some features lacked required `gene_id` attribute, causing pipeline failures.

**Solution:** Used parent gene_id from transcript features:
```python
# In convert_gff_to_gtf.py
if 'gene_id' not in attributes and 'Parent' in attributes:
    attributes['gene_id'] = attributes['Parent']
```

### 3. Double Semicolon Syntax Errors

**Problem:** Some lines had `;;` causing parsing errors.

**Solution:** Simple string replacement:
```python
line = line.replace(';;', ';')
```

### 4. Transcripts Spanning Multiple Chromosomes (45 problematic transcripts)

**Problem:** 45 transcripts had exons on different chromosomes - biologically impossible and breaks aligners.

**Identified transcripts:**
```
rna-XR_003913530.1  (exons on NC_050624.1 and NC_050629.1)
rna-XR_003911803.1  (exons on NC_050622.1 and NC_050628.1)
... [43 more]
```

**Solution:** Created `filter_problematic_transcripts.py` to remove these entirely:
```python
# Identify problematic transcripts
for transcript_id, chroms in transcript_chromosomes.items():
    if len(chroms) > 1:
        problematic_transcripts.add(transcript_id)

# Remove all features for problematic transcripts
if transcript_id not in problematic_transcripts:
    output_file.write(line)
```

### 5. Missing gene_biotype Attributes

**Problem:** featureCounts requires `gene_biotype` for proper quantification.

**Solution:** Added biotype based on gene features:
```python
# In add_gene_biotype.py
biotype_mapping = {
    'mRNA': 'protein_coding',
    'tRNA': 'tRNA',
    'rRNA': 'rRNA',
    'ncRNA': 'ncRNA',
    'misc_RNA': 'misc_RNA'
}
```

### 6. Literal \n Characters in Output

**Problem:** Early script versions wrote literal `\n` instead of newlines.

**Solution:** Fixed string formatting to use actual newlines.

---

## Complete Processing Pipeline

### Step 1: Download from Dryad
```bash
python scripts/01_data_acquisition/download_dryad_genome.py \
    --output-dir data/genome \
    --use-api
```

### Step 2: Convert GFF3 to GTF
```bash
singularity exec albopictus-diapause-rnaseq.sif \
    gffread data/genome/AalbF3.gff3 \
    -T -o data/annotation/AalbF3_raw.gtf
```

### Step 3: Fix GTF Issues
```bash
# Fix missing gene_id and syntax errors
singularity exec albopictus-diapause-rnaseq.sif \
    python scripts/02_annotation_preprocessing/01_fix_gtf_gene_ids.py \
    --input data/annotation/AalbF3_raw.gtf \
    --output data/annotation/AalbF3_fixed.gtf

# Remove problematic transcripts
singularity exec albopictus-diapause-rnaseq.sif \
    python scripts/02_annotation_preprocessing/02_filter_problematic_transcripts.py \
    --input data/annotation/AalbF3_fixed.gtf \
    --output data/annotation/AalbF3_filtered.gtf

# Add gene_biotype for featureCounts
singularity exec albopictus-diapause-rnaseq.sif \
    python scripts/02_annotation_preprocessing/03_add_gene_biotype.py \
    --input data/annotation/AalbF3_filtered.gtf \
    --output data/annotation/AalbF3_final.gtf
```

### Step 4: Validate Final GTF
```bash
# Check for required attributes
grep -c "gene_id" AalbF3_final.gtf
grep -c "gene_biotype" AalbF3_final.gtf

# Verify no problematic transcripts remain
for transcript in $(cat problematic_transcripts.txt); do
    grep -c "$transcript" AalbF3_final.gtf
done
```

---

## Final Statistics

| Metric | Original GFF3 | Final GTF |
|--------|--------------|-----------|
| Total features | 476,234 | 475,189 |
| Genes | 22,175 | 22,175 |
| Transcripts | 32,891 | 32,846 |
| Removed transcripts | - | 45 |
| Lines with gene_id | 458,099 | 475,189 |
| Lines with gene_biotype | 0 | 475,189 |

---

## Scripts Created

1. **`01_fix_gtf_gene_ids.py`**
   - Fixes missing gene_id attributes
   - Removes double semicolons
   - Ensures proper GTF format

2. **`02_filter_problematic_transcripts.py`**
   - Identifies transcripts spanning multiple chromosomes
   - Removes all features for problematic transcripts
   - Creates report of removed transcripts

3. **`03_add_gene_biotype.py`**
   - Adds gene_biotype attribute required by featureCounts
   - Maps transcript types to standard biotypes
   - Preserves all other attributes

---

## Validation with nf-core/rnaseq

The preprocessed GTF was successfully used with:
- nf-core/rnaseq v3.19.0
- STAR aligner (index generation and alignment)
- Salmon quantification
- featureCounts gene-level quantification

All 45 samples processed without GTF-related errors.

---

## Reproducibility Notes

1. **Always use the container** for GTF processing to ensure consistent tool versions
2. **Keep original files** - never modify AalbF3.gff3 directly
3. **Document removed features** - the 45 transcripts are listed in `problematic_transcripts.txt`
4. **Version control scripts** - all preprocessing scripts in GitHub

---

## Common Issues and Solutions

### Issue: "Error: no gene_id attribute"
**Solution:** Run `01_fix_gtf_gene_ids.py`

### Issue: "Transcript has exons on multiple chromosomes"
**Solution:** Run `02_filter_problematic_transcripts.py`

### Issue: "featureCounts only outputs 5 categories"
**Solution:** Ensure `gene_biotype` attribute exists and use `-g gene_id` not `-g gene_biotype`

### Issue: "GTF parsing error at line X"
**Solution:** Check for double semicolons or missing quotes around attribute values

---

## References

- gffread: Part of Cufflinks suite
- GTF format specification: https://www.ensembl.org/info/website/upload/gff.html
- featureCounts manual: Subread documentation

---

*Last updated: October 21, 2025*

---

## CRITICAL UPDATE: October 23, 2025

### MAJOR ISSUE DISCOVERED

The GTF created on Oct 21 was **MISSING 5,079 GENES**:
- **Only had**: 15,542 protein-coding genes
- **Missing**:
  - 153 rRNA genes (couldn't measure contamination!)
  - 3,894 lncRNA genes (8 are GWAS candidates!)
  - 869 tRNA genes
  - Other ncRNAs

### WHY IT HAPPENED

The filtering scripts removed ALL non-protein-coding genes when they shouldn't have.

### WHAT WAS FIXED (Oct 23)

Created `AalbF3_all_gene_types.gtf` from the original GFF3 that:
1. Includes ALL 20,621 genes (all biotypes)
2. Removes "gene-" prefix from gene IDs (LOC109400916, not gene-LOC109400916)
3. Removes "rna-" prefix from transcript IDs
4. Maintains all required GTF attributes

### FILES NOW IN `/data/references/AalbF3/`
- `AalbF3_all_gene_types.gtf` - **USE THIS ONE** (20,621 genes, all types)
- `AalbF3_only_protein_coding.gtf` - The broken one (15,542 genes only)
- `AalbF3_protein_coding_with_prefixes.gtf` - Also broken (has prefixes)

### PIPELINE FIXES APPLIED
1. Updated `scripts/03_rnaseq_pipeline/params.yaml`:
   - Points to `AalbF3_all_gene_types.gtf`
   - Added `featurecounts_group_type: gene_id` (was using gene_biotype!)

### NEXT STEPS
```bash
cd /bigdata/cosmelab/lcosme/projects/albopictus-diapause-rnaseq
sbatch scripts/03_rnaseq_pipeline/02_run_rnaseq.sh
```

The `-resume` flag will detect the GTF change and re-run what's needed (~72 hours).

*Last updated: October 23, 2025*