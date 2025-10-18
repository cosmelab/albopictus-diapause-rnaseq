# 02 - Annotation Preparation

**Purpose:** Convert and fix NCBI GFF3 annotation for compatibility with RSEM and featureCounts.

**Status:** ✅ Completed (October 17, 2025)

---

## Overview

The NCBI AalbF3 annotation (GCA_018104305.1) required extensive preprocessing for compatibility with nf-core/rnaseq pipeline. This directory contains scripts for a three-step GTF generation workflow that addresses multiple annotation problems.

## Background: Why Annotation Preprocessing Was Needed

The original NCBI GFF3 annotation had **6 critical issues** that caused pipeline failures:

1. **Wrong file extension** - nf-core requires `.gff.gz` not `.gff3.gz`
2. **Missing gene_id attributes** - RSEM requires `gene_id` on all GTF lines
3. **Double semicolon syntax errors** - Invalid GTF9 format (`;;` instead of `;`)
4. **Transcripts spanning multiple chromosomes** - RSEM validation failed
5. **Missing gene_biotype attributes** - featureCounts requires gene categorization
6. **Literal `\n` characters** - GTF had text "\\n" instead of actual newlines

All problems are documented in detail in `project.md` (lines 152-268).

---

## Three-Step GTF Generation Workflow

### Step 1: Convert GFF3 to GTF Format

**Script:** `01_convert_gff_to_gtf.py`

**Purpose:** Convert NCBI GFF3 format to GTF format with proper attribute formatting.

**Input:**
- `data/references/AalbF3_genome/AalbF3_annotation.gff.gz` (original NCBI GFF3)

**Output:**
- `data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz`

**Processing:**
- Converts GFF3 attributes (`key=value`) to GTF format (`key "value"`)
- Adds `transcript_id` and `gene_id` to all features
- Skips gene features (only keeps transcript-level and below)
- 376,867 GFF3 lines → 356,246 GTF features

**Usage:**
```bash
python scripts/02_annotation_prep/01_convert_gff_to_gtf.py \
    data/references/AalbF3_genome/AalbF3_annotation.gff.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz
```

---

### Step 2: Filter Problematic Transcripts

**Script:** `02_filter_problematic_transcripts.py`

**Purpose:** Remove transcripts with annotation errors that cause RSEM validation failures.

**Input:**
- `data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz`

**Output:**
- `data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz`

**Filters Applied:**
1. **Strand consistency** - Remove transcripts with exons on different strands (26 transcripts)
2. **Chromosome consistency** - Remove transcripts with exons on multiple chromosomes (43 transcripts)

**Examples of Problematic Transcripts:**
- `rna-XM_029864299.1` - exons on scaffolds 2.198 AND 2.205
- `rna-XM_029866277.1` - inconsistent strands AND multiple chromosomes
- `gene-LOC109431822` - spans scaffolds 2.175 AND 3.83

**Processing:**
- 356,246 features from 53,664 transcripts
- → 344,588 features from 53,619 valid transcripts
- **Loss:** 45 transcripts (0.08% of total)

**Usage:**
```bash
python scripts/02_annotation_prep/02_filter_problematic_transcripts.py \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz
```

---

### Step 3: Add gene_biotype Attributes

**Script:** `03_add_gene_biotype.py`

**Purpose:** Add `gene_biotype` attribute to all features for featureCounts gene categorization.

**Input:**
- `data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz`

**Output:**
- `data/references/AalbF3_genome/AalbF3_annotation.gtf.gz` (FINAL)

**Mapping:**
```
gbkey "mRNA"     → gene_biotype "protein_coding"
gbkey "rRNA"     → gene_biotype "rRNA"
gbkey "tRNA"     → gene_biotype "tRNA"
gbkey "ncRNA"    → gene_biotype "ncRNA"
gbkey "misc_RNA" → gene_biotype "misc_RNA"
```

**Why Needed:** featureCounts uses `gene_biotype` to categorize reads (protein-coding vs rRNA vs other), essential for QC comparisons with published results.

**Processing:**
- 344,588 features → 344,588 features (all with gene_biotype added)
- **Bug fixed:** Original script had `\\n` (literal backslash-n) instead of `\n` (newline)

**Usage:**
```bash
python scripts/02_annotation_prep/03_add_gene_biotype.py \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation.gtf.gz
```

---

## Complete Workflow Example

Run all three steps in sequence:

```bash
# Step 1: GFF3 to GTF conversion
python scripts/02_annotation_prep/01_convert_gff_to_gtf.py \
    data/references/AalbF3_genome/AalbF3_annotation.gff.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz

# Step 2: Filter problematic transcripts
python scripts/02_annotation_prep/02_filter_problematic_transcripts.py \
    data/references/AalbF3_genome/AalbF3_annotation_step1.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz

# Step 3: Add gene_biotype attributes
python scripts/02_annotation_prep/03_add_gene_biotype.py \
    data/references/AalbF3_genome/AalbF3_annotation_step2.gtf.gz \
    data/references/AalbF3_genome/AalbF3_annotation.gtf.gz
```

**Total Runtime:** ~2-3 minutes

---

## Final Output Statistics

**Final GTF:** `data/references/AalbF3_genome/AalbF3_annotation.gtf.gz`

- **Features:** 344,588
- **Transcripts:** 53,619 (45 problematic transcripts removed)
- **Genes:** ~32,889
- **Size:** 4.9 MB (compressed)

**Validation:**
- ✅ All lines have proper gene_id attributes
- ✅ All lines have proper gene_biotype attributes
- ✅ No transcripts span multiple chromosomes
- ✅ No transcripts have inconsistent strands
- ✅ Proper GTF9 format (9 tab-separated columns)
- ✅ Actual newlines (not literal `\n` characters)

---

## For Manuscript Methods Section

When writing the methods, include:

> The NCBI AalbF3 annotation (GCA_018104305.1) required preprocessing for compatibility with nf-core/rnaseq. We developed a three-step workflow:
>
> 1. **GFF3 to GTF conversion** - Converted NCBI GFF3 format to GTF, properly formatting attributes and adding transcript_id/gene_id to all features (356,246 features retained)
>
> 2. **Transcript filtering** - Removed 45 transcripts (0.08% of 53,664 total) with annotation errors: 26 transcripts with inconsistent strand annotations, and 43 transcripts with exons on multiple chromosomes/scaffolds (some had both issues)
>
> 3. **Gene biotype annotation** - Added gene_biotype attributes to all features based on NCBI gbkey mapping for featureCounts compatibility
>
> Final annotation: 344,588 features from 53,619 valid transcripts, all with proper GTF formatting, transcript_id, gene_id, and gene_biotype attributes.

---

## Next Step

After successful annotation preparation:
- **Proceed to:** `03_rnaseq_pipeline/` to run nf-core/rnaseq pipeline

---

## Notes

- Intermediate files (`*_step1.gtf.gz`, `*_step2.gtf.gz`) are kept for reproducibility
- All annotation fixes are tracked in git history
- Original GFF3 is preserved as backup
- This preprocessing ensures zero gene detection issues (previous runs showed 0 detected genes)
