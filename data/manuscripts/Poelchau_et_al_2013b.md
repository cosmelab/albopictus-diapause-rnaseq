# Technical Summary — Poelchau et al. 2013 (Journal of Experimental Biology)

**Title:** RNA-Seq reveals early distinctions and late convergence of gene expression between diapause and quiescent embryos in the mosquito *Aedes albopictus*  
**Authors:** Monica F. Poelchau, Xin Huang, David A. Schmidt, Kristie L. Harshman, Peter A. Armbruster  
**Journal:** Journal of Experimental Biology, 216(21): 4082–4090 (2013)  
**DOI:** [10.1242/jeb.089508](https://doi.org/10.1242/jeb.089508)  
**BioProject:** PRJNA201231  
**Data:** SRA accession SRP026208  

---

## 1. Biological Question and Objective
**Objective:** Examine the transcriptional basis of **embryonic diapause** in *Aedes albopictus* by comparing gene expression of diapause-destined (short-day) and non-diapause (long-day) embryos during development.  
**Main question:** How early do gene expression differences appear between diapause and non-diapause embryos, and what biological processes characterize diapause entry?

---

## 2. Biological Samples
- **Organism:** *Aedes albopictus*  
- **Stage:** Embryos during early development (0–5 days post-oviposition)  
- **Colony origin:** Field-derived from Manassas, Virginia  
- **Treatment:** Two photoperiod regimes to induce diapause vs. non-diapause development  
  - **Short-day (SD; 8L:16D)** → Diapause-inducing  
  - **Long-day (LD; 16L:8D)** → Non-diapause control  
- **Temperature:** 21 °C constant  
- **Samples per condition:** 4 biological replicates per stage and condition  

---

## 3. Experimental Treatments and Sampling Design

| Treatment | Photoperiod | Developmental Stage | Description | Replicates |
|------------|-------------|--------------------|--------------|-------------|
| SD | 8L:16D | 0–5 days post-oviposition | Diapause-inducing embryos | 4 |
| LD | 16L:8D | 0–5 days post-oviposition | Non-diapause embryos | 4 |

**Sampling timepoints:**  
1. **Day 0** – freshly laid eggs (pre-blastoderm)  
2. **Day 1** – early embryogenesis (blastoderm formation)  
3. **Day 2** – gastrulation  
4. **Day 3** – organogenesis  
5. **Day 4** – late organogenesis  
6. **Day 5** – pre-hatching stage  

Eggs were collected daily for each photoperiod and pooled by stage.

---

## 4. Sample Preparation
1. **Collection:** Embryos separated by age (0–5 days post-oviposition).  
2. **Preservation:** RNAlater, stored at −80 °C.  
3. **RNA extraction:** TRI Reagent (Sigma) followed by Turbo DNase treatment.  
4. **Quality control:** RNA integrity checked via Agilent 2100 Bioanalyzer.  
5. **Library preparation:** Illumina TruSeq RNA Sample Prep Kit v2.  
6. **Replicates:** 4 per photoperiod × 6 timepoints → 48 total libraries.

---

## 5. Data Collection and Sequencing Platform

| Parameter | Details |
|------------|----------|
| **Platform** | Illumina HiSeq 2000 |
| **Read type** | Paired-end 101 bp |
| **Insert size** | ~200 bp |
| **Total reads generated** | 1.0 × 10⁹ paired reads across 48 libraries |
| **Reads per library (avg)** | 20–25 million |
| **Data accession** | SRA SRP026208 (BioProject PRJNA201231) |

---

## 6. Computational Workflow and Analysis

| Step | Tool | Parameters |
|------|------|-------------|
| Adapter trimming | ssaha2 vs UniVec | ≥95% identity, score ≥18 |
| Quality filtering | SolexaQA v2.2 | phred ≥30, min length 50 bp |
| Digital normalization | khmer | 20-mer, coverage cutoff = 20 |
| Assembly | Trinity (r2012-10-05) | Default parameters |
| Redundancy reduction | CD-HIT-EST | 99% identity |
| Reference merging | CAP3 + EXONERATE | ≥70% identity |
| Quantification | RSEM v1.2.4 | Expectation–maximization |
| Normalization | edgeR (TMM) | logCPM ≥ 1 in ≥ 4 libraries |
| Differential expression | edgeR | |log₂FC| > 0.5, FDR < 0.05 |
| Pathway analysis | GOseq + KEGG (*Aedes aegypti*) | FDR < 0.05; ≥ 5 DE genes per set |

---

## 7. Statistical Overview
- **DE threshold:** |log₂FC| > 0.5 and FDR < 0.05 (Benjamini–Hochberg).  
- **Normalization:** TMM via edgeR.  
- **Replicates:** 4 biological replicates per stage per condition.  
- **Validation:** qPCR on representative transcripts confirmed RNA-seq patterns.

---

## 8. Main Findings

### A. Early-Stage Differences (Days 0–2)
- **Diapause embryos:** Strong upregulation of genes linked to *signal transduction*, *hormonal pathways*, and *cell cycle regulation*.  
- **Non-diapause embryos:** Enrichment in *cell proliferation* and *metabolic activity*.  
- Indicates that diapause-destined embryos enter a **distinct transcriptional program early** in development.

### B. Mid-to-Late Stages (Days 3–5)
- Expression profiles of SD and LD embryos **converge**, reflecting morphological completion.  
- Most DE genes at early stages lose differential expression as development proceeds.

### C. Enriched Pathways
- Diapause embryos show enrichment in:
  - **MAPK signaling**
  - **Insulin signaling**
  - **Cytoskeletal remodeling**
  - **Energy metabolism and lipid storage**
- Suggests **metabolic reallocation and developmental arrest** as hallmarks of diapause entry.

---

## 9. Biological Meaning
The results demonstrate that **diapause and quiescent embryos differ transcriptionally early** in development, with diapause induction involving major metabolic and signaling reprogramming.  
By day 5, gene expression between the two treatments converges, consistent with both reaching a quiescent state.  
Thus, diapause is **not merely a delayed developmental arrest**, but a **regulated alternative developmental trajectory** initiated soon after oviposition.

---

## 10. Reproducibility Checklist

| Component | Description |
|------------|-------------|
| Species | *Aedes albopictus* |
| Stage | Embryos (0–5 days post-oviposition) |
| Photoperiods | 8L:16D (SD) vs 16L:8D (LD) |
| Temperature | 21 °C constant |
| Replicates | 4 biological per stage and treatment |
| RNA extraction | TRI Reagent + Turbo DNase |
| Library prep | Illumina TruSeq RNA v2 |
| Sequencing | HiSeq 2000, PE 101 bp |
| QC | phred ≥ 30, length ≥ 50 bp |
| Assembly | Trinity + CD-HIT-EST + CAP3 merges |
| Quantification | RSEM v1.2.4 |
| DE analysis | edgeR, |log₂FC| > 0.5, FDR < 0.05 |
| Enrichment | GOseq + KEGG (*Ae. aegypti*) |
| Accession | PRJNA201231 |

---

## 11. Summary
Photoperiodic cues induce early transcriptional reprogramming in *Aedes albopictus* embryos.  
**Diapause-destined embryos** activate signaling and metabolic control genes immediately after oviposition, establishing a distinct developmental pathway from non-diapause embryos.  
Later, both states converge transcriptionally as embryos reach quiescence, highlighting diapause as a **programmed developmental alternative** rather than a passive dormancy.

---
