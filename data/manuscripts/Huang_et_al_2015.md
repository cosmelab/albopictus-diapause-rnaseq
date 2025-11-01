# Technical Summary — Huang et al. 2015 (PLOS Neglected Tropical Diseases)

**Title:** Global Transcriptional Dynamics of Diapause Induction in Non-Blood-Fed and Blood-Fed *Aedes albopictus*  
**Authors:** Xin Huang, Monica F. Poelchau, Peter A. Armbruster  
**Journal:** PLOS Neglected Tropical Diseases 9(4): e0003724 (2015)  
**DOI:** [10.1371/journal.pntd.0003724](https://doi.org/10.1371/journal.pntd.0003724)  
**BioProject:** PRJNA268379  

---

## 1. Biological Question and Objective
The study investigated **how photoperiod and blood-feeding influence global transcriptional changes** in *Aedes albopictus* adult females during diapause induction.  
**Main question:** Which molecular pathways and genes are regulated during the induction of diapause, and how do photoperiod and nutritional state (blood-fed vs. non-blood-fed) affect these changes?

---

## 2. Biological Samples
- **Organism:** *Aedes albopictus* (Asian tiger mosquito).  
- **Stage:** Adult females, whole-body samples.  
- **Population origin:** >200 larvae collected from used tires, Manassas, VA, USA.  
- **Generations:** Laboratory F3 generation reared at 21 °C and 80% relative humidity.  
- **Treatment duration before sampling:** 11 days under assigned photoperiods.

---

## 3. Experimental Treatments and Design

| Treatment | Photoperiod | Feeding State | Description | Biological Replicates |
|------------|--------------|---------------|--------------|------------------------|
| SD-NB | 8L:16D (short day) | Non-blood-fed | Diapause-inducing, non-fed | 4 |
| SD-BM | 8L:16D (short day) | Blood-fed | Diapause-inducing, post-feeding | 4 |
| LD-NB | 16L:8D (long day) | Non-blood-fed | Non-diapause, non-fed | 4 |
| LD-BM | 16L:8D (long day) | Blood-fed | Non-diapause, post-feeding | 4 |

**Total libraries:** 16 (4 per treatment).  
**Blood-feeding protocol:** Human host feeding at Zeitgeber time 3–4 h. Females sampled **26–28 h post-blood meal** (ZT 6–8 h).  
**Phenotypic validation:** Diapause incidence (SD = 81–95%; LD = 1–3%).

---

## 4. Sample Preparation
1. **Collection:** Blood bolus dissected from midguts. Whole-body females stored in RNAlater at 4 °C for 24 h, then frozen at −80 °C.  
2. **RNA extraction:** TRI Reagent (Sigma) followed by Turbo DNase (Ambion) digestion.  
3. **Quality control:** RNA quality verified using Agilent 2100 Bioanalyzer.  
4. **Library preparation:** Illumina TruSeq RNA v2 kit; unique barcodes per sample.  
5. **Pooling:** Barcoded libraries multiplexed across sequencing lanes.

---

## 5. Data Collection and Sequencing Platform

| Parameter | Details |
|------------|----------|
| **Sequencing Platform** | Illumina HiSeq 2000 |
| **Read Type** | Paired-end, 101 bp reads |
| **Insert Size** | ~203 bp |
| **Lanes** | Two |
| **Total Raw Reads** | 472,006,080 read pairs |
| **Reads Retained after QC** | 182,036,362 paired reads (+108,891,552 single reads) |
| **Mapping Rate** | 89.4% mapped to reference transcriptome |
| **Accession** | SRA PRJNA268379 |

---

## 6. Computational Workflow and Analysis

| Step | Tool / Method | Key Parameters |
|------|----------------|----------------|
| Adapter & rRNA removal | ssaha2 vs UniVec and *Ae. albopictus* rRNA | ≥95% identity; alignment score ≥18 |
| Quality filtering | SolexaQA v2.2 | phred ≥30; read length ≥50 bp |
| Digital normalization | khmer | 20-mer; coverage cutoff = 20 |
| Assembly | Trinity (Feb 2013) | Default; min k-mer coverage = 2 |
| Redundancy reduction | CD-HIT-EST | 99% identity |
| Reference-guided merging | CAP3 + EXONERATE | ≥70% identity; BLASTX E ≤ 1e−6 |
| Annotation reference | OrthoDB v7 | Diptera orthologs |
| Read mapping / quantification | RSEM v1.2.4 | Gene-level expression estimates |
| Normalization | edgeR (TMM) | Filter: logCPM ≥ 1 in ≥ 4 libraries |
| Differential expression | edgeR | |log₂FC| > 0.5; FDR < 0.05 (Benjamini–Hochberg) |
| Functional enrichment | GOseq + KEGG (*Ae. aegypti*) | FDR < 0.05; ≥ 5 DE genes per pathway |
| Visualization | limma-voom, MDS, heatmaps | Z-score standardized |

---

## 7. Statistical Overview
- **Total samples analyzed:** 16  
- **Thresholds:**  
  - |log₂FC| > 0.5  
  - FDR < 0.05  
  - logCPM ≥ 1 in ≥ 4 libraries retained  
- **Validation:** qRT-PCR of selected genes correlated r = 0.92 with RNA-seq.  

---

## 8. Key Findings

### A. Non-blood-fed (SD vs LD)
- Upregulated: *timeless*, *cryptochrome 1*, *pepck*, *Δ9-desaturases*, *JH-inducible proteins*.  
- Downregulated: DNA replication and cell-cycle genes (*pcna*, *gadd45*).  
- Enriched pathways: amino acid biosynthesis, alanine/aspartate/glutamate metabolism.  
- **Interpretation:** Early metabolic reprogramming and circadian response initiate diapause induction.

### B. Blood-fed (SD vs LD)
- Upregulated: *fatty acid synthase*, desaturases, *CYP302A1*, *CYP314A1*.  
- Enriched: oxidative phosphorylation, fatty acid synthesis, β-alanine metabolism.  
- **Interpretation:** Photoperiodic cue plus blood meal enhance energy metabolism and lipid provisioning for diapause eggs.

### C. Shared responses
- Overlapping SD-upregulated genes in NB + BM include circadian, oxidative, and metabolic regulators.  
- Downregulated categories center on cell cycle and replication control.

---

## 9. Biological Meaning
- **Photoperiodic diapause induction** begins before egg development, signaled by circadian and metabolic shifts in adults.  
- **Energy allocation** and **lipid metabolism** dominate transcriptional programs preparing for diapause.  
- Blood feeding reinforces these adjustments, indicating maternal provisioning for diapause eggs.  
- Provides the first transcriptome-wide view of diapause induction in mosquitoes.

---

## 10. Reproducibility Checklist
| Component | Description |
|------------|-------------|
| Species | *Aedes albopictus* |
| Stage | Adult female (whole body, blood bolus removed) |
| Photoperiods | 8L:16D (SD) vs 16L:8D (LD) |
| Feeding | NB (non-blood-fed) vs BM (blood-fed) |
| Replicates | 4 per group (16 total libraries) |
| RNA extraction | TRI Reagent + Turbo DNase |
| Library kit | Illumina TruSeq RNA v2 |
| Sequencing | HiSeq 2000, paired-end 101 bp |
| QC | phred ≥ 30, length ≥ 50 bp |
| Assembly | Trinity → CD-HIT-EST → CAP3 merges |
| Mapping | RSEM v1.2.4 |
| Normalization | edgeR TMM |
| DE thresholds | |log₂FC| > 0.5, FDR < 0.05 |
| Enrichment | GOseq + KEGG, FDR < 0.05, ≥ 5 genes |
| Data accession | PRJNA268379 |

---

## 11. Summary
Photoperiod and nutritional status shape distinct transcriptional programs in adult females during diapause induction.  
**Short-day photoperiods** trigger early circadian and metabolic shifts, while **blood feeding under short-day conditions** amplifies oxidative and lipid metabolism for diapause egg provisioning.  
These processes reveal the integrated metabolic and endocrine reprogramming underlying mosquito diapause.

---
