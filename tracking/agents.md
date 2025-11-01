# Agent-Assisted Workflow Implementation Guide

**Project:** Albopictus Diapause RNA-seq Analysis
**Last Updated:** 2025-11-01
**Environment:** HPC with Singularity container

## Overview

This document outlines strategic implementation points for using specialized AI agents (`/agents` command) throughout the analysis pipeline. Each agent type serves specific phases where autonomous task completion, code generation, or architectural planning provides maximum value.

## Container Strategy Context

**Current workflow:**
- Primary analysis container: `albopictus-diapause-rnaseq.sif` (1.5GB)
- Contains all tools: Python, R, DESeq2, tidyverse, plotting libraries
- Used for all computational steps via Singularity on HPC

**Future deliverable (Phase 13):**
- Manuscript-only container: minimal image with final analysis scripts
- Purpose: Enable reproduction of published figures/tables without full pipeline
- Contents: Final R/Python scripts, processed data inputs, output generation only
- Goal: Reduce barrier to replication for reviewers/readers

## Agent Types Available

### Software Architect
- Design system architecture and directory structures
- Plan workflow organization and dependencies
- Create implementation strategies for complex features

### Code Writer
- Generate scripts from specifications
- Implement algorithms and statistical tests
- Create automation and pipeline code

### Code Reviewer
- Validate code correctness and adherence to specifications
- Identify bugs, edge cases, and optimization opportunities
- Verify reproducibility and documentation quality

---

## Recommended Implementation Points

### Phase 2: Reorganize Results Directory Structure (CURRENT - PRIORITY)

#### Software Architect

**Task:** Design clean, hierarchical results/ structure

**Rationale:**
- Current structure has duplicates (qc_analysis, qc_figures, qc_metrics, count_matrices x4)
- Need organization that matches 3-part analysis strategy
- Must support stage-specific analyses and reviewer requirements
- Should facilitate manuscript container extraction later

**Implementation:**
```bash
/agents Software Architect
```

**Prompt:**
> Design a comprehensive results/ directory structure for RNA-seq differential expression analysis with the following requirements:
>
> **Analysis Strategy (3 parts):**
> 1. Part 1: Replicate collaborators with featureCounts
> 2. Part 1.5: Validate Salmon vs featureCounts (gate check)
> 3. Part 2: Answer reviewers with Salmon (only if validation passes)
>
> **Data Organization:**
> - 3 life stages analyzed independently (adults, embryos, larvae)
> - Stage-specific subdirectories within each analysis type
> - QC metrics, count matrices, differential expression results
> - Method validation outputs (correlation plots, DEG overlap)
> - GWAS candidate expression tables
> - Reviewer response tables and figures
>
> **Constraints:**
> - All paths referenced in R/Python scripts must update
> - Structure should enable easy extraction for manuscript-only container
> - Must avoid current duplication issues
> - Support both featureCounts and Salmon workflows in parallel
>
> Provide: Directory tree, file naming conventions, and migration plan from current structure.

**Output Expected:**
- Detailed directory tree structure
- File organization conventions
- Migration strategy from current mess
- Script path update checklist

---

### Phase 3: Validate Salmon vs featureCounts

#### Code Writer

**Task:** Generate correlation and validation analysis scripts

**Rationale:**
- Critical gate check: must prove r > 0.95 correlation, >90% DEG overlap
- Complex multi-part validation (per-sample correlation, DEG overlap, effect size)
- Standardized statistical approach needed across all 45 samples
- Must generate publication-quality diagnostic plots

**Implementation:**
```bash
/agents Code Writer
```

**Prompt:**
> Write R script to validate Salmon quantification against featureCounts gene counts:
>
> **Inputs:**
> - featureCounts: `output/count_matrices/featurecounts/*_counts.tsv`
> - Salmon: `output/count_matrices/salmon/*_counts.tsv`
> - Metadata: `data/metadata/samples.csv`
>
> **Analysis Requirements:**
> 1. Per-sample Pearson correlation (all 45 samples)
> 2. Global correlation plot with r² annotations
> 3. Run identical DESeq2 on adults_fed dataset (SD vs LD, n=4 vs 4)
> 4. Compare DEG lists: overlap percentage, shared/unique genes
> 5. Correlate log2FC values between methods
>
> **Validation Thresholds (CRITICAL):**
> - Per-sample correlation: r > 0.95 for all samples
> - DEG overlap: >90% agreement
> - Effect size correlation: r > 0.9
>
> **Output:**
> - Validation report: `results/validation/salmon_vs_featurecounts_report.txt`
> - Correlation plots: `results/validation/figures/`
> - DEG comparison tables
> - PASS/FAIL determination with clear messaging
>
> **Technical Requirements:**
> - Use DESeq2 from BiocManager
> - Handle gene ID matching between methods
> - Log-transform counts for correlation (log2(count + 1))
> - Generate publication-quality ggplot2 figures
> - Must run in Singularity container: `singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif bash -c "cd /proj && Rscript script.R"`

**Output Expected:**
- Complete R script with documentation
- SLURM submission script
- Expected runtime and memory requirements

---

### Phase 4-6: Differential Expression Analysis

#### Code Reviewer

**Task:** Validate DESeq2 scripts match collaborator specifications exactly

**Rationale:**
- Part 1 requires exact replication of collaborator methods
- Multiple stage-specific specifications (adults, embryos, larvae)
- Different thresholds per collaborator (0.58 vs 0.5 log2FC)
- Critical for validating pipeline on new genome

**Implementation:**
```bash
/agents Code Reviewer
```

**Prompt:**
> Review the following DESeq2 analysis scripts for exact adherence to collaborator specifications:
>
> **Scripts to review:**
> - `scripts/08_differential_expression/01_adults_unfed_deseq2.R`
> - `scripts/08_differential_expression/02_adults_fed_deseq2.R`
> - `scripts/08_differential_expression/03_embryos_72h_deseq2.R`
> - `scripts/08_differential_expression/04_embryos_135h_deseq2.R`
> - `scripts/08_differential_expression/06_larvae_11d_deseq2.R`
>
> **Exact specifications per collaborator:**
>
> Adults (Angela):
> - Split by blood meal status (unfed/fed analyzed separately)
> - Design: `~ condition` (SD vs LD)
> - Reference level: "LD"
> - LFC shrinkage: `lfcShrink(dds, coef="condition_SD_vs_LD", type="apeglm")`
> - Pre-filter: `rowSums(counts(dds)) >= 10`
> - Significance: `padj < 0.05 AND abs(log2FC) > 0.58`
>
> Embryos (Mackenzie):
> - Split by timepoint (72h/135h separately)
> - Design: `~ condition` (DI vs NDI) ← NOTE: NOT SD/LD!
> - Reference level: "NDI"
> - LFC shrinkage: `lfcShrink(dds, coef="condition_DI_vs_NDI", type="apeglm")`
> - Pre-filter: `rowSums(counts(dds)) >= 10`
> - Significance: `padj < 0.05 AND abs(log2FC) > 0.58`
>
> Larvae (Sarah):
> - 11d timepoint only (LD valid at this stage)
> - Design: `~ condition` (D vs ND) ← NOTE: NOT SD/LD!
> - Reference level: "ND"
> - LFC shrinkage: `lfcShrink(dds, coef="condition_D_vs_ND", type="apeglm")`
> - Pre-filter: `rowSums(counts(dds)) >= 10`
> - Significance: `padj < 0.05 AND abs(log2FC) > 0.5` ← NOTE: 0.5, not 0.58!
>
> **Review Checklist:**
> - [ ] Correct design formula per stage
> - [ ] Correct reference level (LD/NDI/ND)
> - [ ] Correct coefficient for lfcShrink
> - [ ] apeglm shrinkage method specified
> - [ ] Pre-filtering threshold correct
> - [ ] Significance thresholds match (especially larvae 0.5 vs others 0.58)
> - [ ] Condition labels match specifications (DI/NDI for embryos, D/ND for larvae)
> - [ ] Output paths follow new directory structure
> - [ ] Metadata subsetting correct per stage/timepoint
>
> **Critical:** Flag any deviations from specifications. These must match exactly to replicate collaborator results.

**Output Expected:**
- Detailed review report with pass/fail per script
- Specific line-by-line corrections needed
- Verification of statistical approach
- Documentation improvements

---

### Phase 10: Reviewer Requirements

#### Code Writer

**Task:** Generate statistical tests and tables for reviewer responses

**Rationale:**
- Reviewer 2 requires Fisher's exact test for GWAS enrichment
- Reviewer 3 requires genome-wide statistics at multiple thresholds
- Standard formatting needed across all response tables
- Must be publication-ready with clear statistical reporting

**Implementation 1: Fisher's Exact Test**
```bash
/agents Code Writer
```

**Prompt:**
> Write R script to perform Fisher's exact test for GWAS candidate enrichment in differential expression results:
>
> **Scientific Question:**
> Is the observation that 4 out of 5 GWAS candidate genes are differentially expressed significantly enriched compared to genome-wide expectation?
>
> **Inputs:**
> - GWAS candidates: `data/metadata/gwas_candidates_final_34.txt` (gene IDs)
> - DESeq2 results: `results/differential_expression/*/deseq2_results.tsv`
> - Total genes tested: 22,177 (from featureCounts)
>
> **Analysis per stage:**
> - Adults unfed: SD_Unfed vs LD_Unfed
> - Adults fed: SD_Fed vs LD_Fed
> - Embryos 72h: DI_72h vs NDI_72h
> - Embryos 135h: DI_135h vs NDI_135h
> - Larvae 11d: D_11d vs ND_11d
>
> **Fisher's Exact Test Setup:**
>
> |                    | GWAS Candidates | Non-GWAS Genes |
> |--------------------|-----------------|----------------|
> | DE (padj < 0.05)   | a               | b              |
> | Not DE             | c               | d              |
>
> Test: `fisher.test(matrix(c(a,b,c,d), nrow=2))`
>
> **Output Requirements:**
> 1. Summary table with columns:
>    - Stage/comparison
>    - GWAS candidates DE / total GWAS candidates
>    - Genome-wide DE / total genes tested
>    - Odds ratio
>    - Fisher's exact test p-value
>    - Interpretation (enriched/depleted/not significant)
>
> 2. Meta-analysis across stages:
>    - Combined evidence using Fisher's method (combine p-values)
>    - Overall enrichment conclusion
>
> 3. Sensitivity analysis:
>    - Test at multiple thresholds: padj < 0.01, 0.05, 0.1
>    - Test with/without log2FC threshold
>
> **Output Files:**
> - `results/reviewer_response/gwas_enrichment_test.tsv`
> - `results/reviewer_response/gwas_enrichment_summary.txt`
>
> **Technical:**
> - Handle missing gene IDs gracefully
> - Report exact p-values (not "< 0.001")
> - Include 95% confidence intervals for odds ratios
> - Clear interpretation for non-statisticians

**Implementation 2: Genome-wide Statistics**
```bash
/agents Code Writer
```

**Prompt:**
> Write Python script to generate genome-wide differential expression statistics table for reviewer response:
>
> **Inputs:**
> - DESeq2 results from all stages: `results/differential_expression/*/deseq2_results.tsv`
> - Stage metadata to organize table
>
> **Table Requirements:**
>
> Count DEGs at multiple threshold combinations:
>
> | Stage       | Comparison | Total Genes | padj<0.05 | padj<0.05 & \|log2FC\|>0.5 | padj<0.05 & \|log2FC\|>1 | padj<0.01 | padj<0.01 & \|log2FC\|>1 |
> |-------------|------------|-------------|-----------|---------------------------|-------------------------|-----------|-------------------------|
> | Adults      | Unfed      | 22,177      | ...       | ...                       | ...                     | ...       | ...                     |
> | Adults      | Fed        | 22,177      | ...       | ...                       | ...                     | ...       | ...                     |
> | Embryos     | 72h        | 22,177      | ...       | ...                       | ...                     | ...       | ...                     |
> | Embryos     | 135h       | 22,177      | ...       | ...                       | ...                     | ...       | ...                     |
> | Larvae      | 11d        | 22,177      | ...       | ...                       | ...                     | ...       | ...                     |
>
> **Additional columns:**
> - Up-regulated count (log2FC > 0)
> - Down-regulated count (log2FC < 0)
> - Percent of genome differentially expressed
>
> **Output:**
> - Formatted markdown table: `results/reviewer_response/genome_wide_stats.md`
> - CSV for supplementary: `results/reviewer_response/genome_wide_stats.csv`
> - LaTeX table: `results/reviewer_response/genome_wide_stats.tex`
>
> **Technical:**
> - Use pandas for table manipulation
> - Handle NA/NaN values in padj column
> - Verify gene count consistency across stages
> - Generate publication-ready formatting

**Output Expected:**
- Complete scripts with error handling
- Output tables in multiple formats
- Clear statistical reporting
- README documenting how to run

---

### Phase 11: GWAS Candidate Analysis

#### Software Architect + Code Writer

**Task:** Design and implement candidate gene extraction and visualization system

**Rationale:**
- Must extract expression for GWAS candidates across all stages
- Create integrated cross-stage visualization
- Support iterative updates when user re-does GWAS
- Prepare for manuscript figures

**Implementation (2-step):**

**Step 1: Software Architect**
```bash
/agents Software Architect
```

**Prompt:**
> Design a flexible system for extracting and visualizing GWAS candidate gene expression across RNA-seq stages:
>
> **Requirements:**
> - Input: List of gene IDs (updated periodically as GWAS analysis refines)
> - Extract from 5 comparisons: adults_unfed, adults_fed, embryos_72h, embryos_135h, larvae_11d
> - Generate cross-stage heatmap showing log2FC and significance
> - Create volcano plots with candidates highlighted
> - Support both featureCounts and Salmon results
> - Prepare for manuscript figure generation
>
> **Design Considerations:**
> - Easy to update when candidate list changes
> - Automated regeneration of all plots
> - Consistent color schemes and formatting across figures
> - Export publication-quality PDFs (300 DPI, vector graphics)
> - Document which candidates are DE in which stages
>
> **Integration:**
> - Should plug into existing results/differential_expression/ structure
> - Output to results/gwas_candidates/ for manuscript container later
> - Must work within Singularity container environment
>
> Provide: System architecture, file organization, workflow diagram, and implementation checklist.

**Step 2: Code Writer** (after architecture approved)
```bash
/agents Code Writer
```

Use approved architecture to generate implementation scripts.

---

### Phase 13: Manuscript-Only Container (NEW PHASE)

#### Software Architect

**Task:** Design minimal container for manuscript figure/table reproduction

**Rationale:**
- Full analysis container is 1.5GB with many dependencies
- Reviewers/readers only need to reproduce final figures and tables
- Reduce replication barrier by providing minimal environment
- Extract only essential scripts and processed inputs

**Implementation:**
```bash
/agents Software Architect
```

**Prompt:**
> Design a minimal Docker container for reproducing manuscript figures and tables only:
>
> **Full Pipeline Container (current):**
> - `albopictus-diapause-rnaseq.sif` (1.5GB)
> - Contents: Python, R, DESeq2, STAR, Salmon, MultiQC, tidyverse, ggplot2, etc.
> - Purpose: Run entire pipeline from raw reads to final results
>
> **Manuscript Container (goal):**
> - Minimal image with only visualization/table generation dependencies
> - Input: Processed DESeq2 results, count matrices (already computed)
> - Output: All manuscript figures and supplementary tables
> - Target size: <500MB
>
> **Required Components:**
> - R with ggplot2, pheatmap, ComplexHeatmap for figures
> - Python with pandas, matplotlib, seaborn for tables
> - Scripts: Only final visualization/table generation (10-15 scripts)
> - Data: DESeq2 results, metadata, GWAS candidates, QC metrics (not raw reads)
>
> **Deliverables:**
> 1. Dockerfile for manuscript-only container
> 2. Directory structure for input data
> 3. Master script to regenerate all figures/tables
> 4. README for users: download container → run script → get outputs
> 5. Zenodo upload plan: container image + input data archive
>
> **Constraints:**
> - Must work independently of HPC environment
> - No Singularity requirement (use Docker directly)
> - Clear documentation for non-expert users
> - All outputs match exactly what's in manuscript
>
> Provide: Dockerfile, directory structure, automation script, and deployment plan.

**Output Expected:**
- Minimal Dockerfile
- Input data organization
- One-command reproduction script
- Zenodo/DockerHub deployment plan

#### Code Writer (after architecture approved)

**Task:** Implement manuscript container automation

Generate master script that:
1. Validates input data presence
2. Runs all figure generation scripts
3. Runs all table generation scripts
4. Validates outputs match manuscript
5. Generates reproduction report

---

## Workflow Integration

### Suggested Agent Usage Pattern

1. **Planning Phase** → Software Architect
   - Directory structure design
   - Workflow architecture
   - System integration planning

2. **Implementation Phase** → Code Writer
   - Script generation from specifications
   - Test automation
   - Pipeline scripting

3. **Validation Phase** → Code Reviewer
   - Verify correctness against specifications
   - Code quality and reproducibility checks
   - Documentation review

### Container Execution Context

All generated scripts must work within container:
```bash
module load singularity
singularity exec --cleanenv --bind $PWD:/proj albopictus-diapause-rnaseq.sif \
  bash -c "cd /proj && <command>"
```

Agents should be instructed to generate container-compatible code.

---

## Priority Recommendations

### Immediate (Phase 2):
1. **Software Architect**: Design results/ directory structure
2. **Code Writer**: Generate migration scripts for reorganization

### Near-term (Phases 3-6):
3. **Code Writer**: Salmon validation scripts
4. **Code Reviewer**: Validate DESeq2 replication scripts

### Mid-term (Phase 10):
5. **Code Writer**: Reviewer response statistical tests
6. **Code Writer**: Genome-wide statistics tables

### Long-term (Phases 11-13):
7. **Software Architect**: GWAS candidate visualization system
8. **Software Architect**: Manuscript-only container design
9. **Code Writer**: Container automation and deployment

---

## Documentation Standards

When using agents, maintain documentation:
- Save agent prompts to `tracking/agent_prompts/`
- Document agent outputs in `tracking/agent_outputs/`
- Update `project_state.json` when agent completes tasks
- Record agent usage in session notes

This ensures reproducibility of the agent-assisted workflow itself.

---

**Next Steps:**
1. Review this implementation guide
2. Execute Software Architect for Phase 2 directory reorganization
3. Update `project_state.json` with Phase 13 (manuscript container)
4. Begin systematic agent integration per phase

