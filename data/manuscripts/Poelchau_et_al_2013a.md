# Technical Summary — Poelchau et al. 2013 (Proceedings of the Royal Society B)

**Title:** Deep sequencing reveals complex mechanisms of diapause preparation in the invasive mosquito, *Aedes albopictus*  
**Authors:** Monica F. Poelchau, Julie A. Reynolds, Christine G. Elsik, David L. Denlinger, Peter A. Armbruster  
**Journal:** Proceedings of the Royal Society B, 280: 20130143 (2013)  
**DOI:** https://doi.org/10.1098/rspb.2013.0143  
**Open access:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3619507/  
**Data:** SRA submissions SRA044835 and SRA051478; assembly at http://AlbopictusExpression.org

---

## 1) Biological Question & Objective
- **Objective:** Identify transcriptional mechanisms of **pre-diapause (diapause preparation)** in embryos destined for diapause versus non-diapause.  
- **Central question:** Which pathways and genes differ between **short-day (D)** and **long-day (ND)** embryos at key developmental landmarks, and do results support a conserved “diapause genetic toolkit”?

---

## 2) Biological Samples
- **Species:** *Aedes albopictus*.  
- **Stage:** Embryos collected at two ages post-oviposition (pov): **3 days (72–78 h)** and **6 days (135–141 h)**.  
- **Colony:** Laboratory strain (F9).  
- **Environment:** 21 °C, ~80% RH.  
- **Photoperiod during induction:**  
  - **D (8L:16D)** → produces diapause-destined eggs.  
  - **ND (16L:8D)** → produces non-diapause eggs.  
- **Diapause incidence check:** D treatment produced **~99% diapause** eggs.

---

## 3) Experimental Design & Treatments
- Rear six cohorts at 16L:8D to pupation.  
- Upon pupation, split to **three biological replicates under D** and **three under ND** (≈500 individuals per cohort).  
- Collect eggs over a **6 h oviposition window**, then age and harvest embryos at:
  - **3d pov** (germ band retraction/dorsal closure landmark).  
  - **6d pov** (late embryogenesis: segmentation and head–thorax separation).  
- **Design matrix (RNA-seq):** 2 photoperiods (D/ND) × 2 timepoints (3d/6d) × 3 replicates = **12 libraries**.

---

## 4) Sample Preparation
1. **Staging:** Bleach/clear subsets to verify morphology; confirm D and ND are morphologically matched at each timepoint.  
2. **Diapause incidence:** Reserve egg subsets per replicate; score diapause per established protocol.  
3. **RNA:** Extract total RNA from each embryo sample; build strand-agnostic mRNA-seq libraries.

---

## 5) Data Collection — Sequencing
- **Platform:** Illumina **HiSeq**.  
- **Layout:** Paired-end mRNA-Seq libraries.  
- **Lanes:** Two flow-cell lanes (all **3d** libraries on lane 1; all **6d** libraries on lane 2).  
- **Accessions:** **SRA044835**, **SRA051478**.  
- **Reference resource:** Transcriptome at **AlbopictusExpression.org**.

---

## 6) Computational Workflow (Replicable Steps)
| Step | Tool / Resource | Key settings / notes |
|---|---|---|
| Read QC / screening | Adapter/rRNA filtering; standard RNA-seq QC | Ensure uniform handling across 12 libs |
| Quantification | **RSEM** | Gene (unigene) abundance estimates |
| Normalization | **edgeR TMM** | Library size & composition normalization |
| DE testing | **limma** (Bioconductor) | Linear models with empirical Bayes |
| DE thresholds | — | **FDR (BH) < 0.05** and **\|log2FC\| > 0.5** |
| Enrichment (primary) | **DAVID** | **KEGG** + **GO FAT**, **FDR < 0.05**; annotation clustering to reduce redundancy |
| Enrichment (custom) | **GENEMERGE v1.2** | Curated sets: **insulin signaling**, **heat-shock/stress**, **ecdysone signaling** |
| Orthology | **BioMart**, OrthoDB | 1:1 and apparent 1:1 orthologs (Diptera) |
| Developmental context | **D. melanogaster** RNA-seq time series | Compare *Ae. albopictus* DE gene trends with **8–12 h** and **18–22 h** Dmel stages using z-score heatmaps |

---

## 7) Statistical Analysis & Thresholds
- **Units:** Gene-level counts (RSEM).  
- **Normalization:** **TMM** (edgeR).  
- **DE definition:** **FDR < 0.05** and **\|log2FC\| > 0.5** (limma).  
- **Enrichment significance:** **FDR < 0.05** (DAVID annotation clusters; representative, non-redundant term per cluster).  
- **Replicates:** **n = 3** per condition/timepoint; **12 total libraries**.

---

## 8) Main Findings (Quantitative & Pathways)
- **More DE later:** **6d** shows **4,337** DE genes vs **2,506** at **3d** (D vs ND).  
- **3d (early pre-diapause):**  
  - Enriched: **cell cycle/mitosis**, **DNA replication/repair**, **DNA binding**; strong signal for **ecdysone signaling** among DE genes.  
  - Pattern: Many **positive cell-cycle regulators** show **higher expression in D** at 3d.  
- **6d (late pre-diapause):**  
  - Enriched: **mitochondrion**, **oxidation–reduction**, **energy production**, **transport** → numerous transcripts **down** in D (consistent with depressed metabolism before diapause).  
  - Enriched and **up** in D: **cuticle structure** and **polysaccharide/chitin** categories → consistent with **serosal/cuticular provisioning** and desiccation resistance.  
- **Candidate “toolkit” genes:**  
  - **Pepck** (gluconeogenesis) **up in D** at both 3d and 6d → cross-taxon marker for diapause energetics.  
  - **PCNA** (cell-cycle progression) shows stage-dependent regulation consistent with diapause-linked cell-cycle modulation.  
- **Not enriched at transcript level in pre-diapause:** curated **insulin signaling** and **stress (HSP)** sets.  
- **Developmental comparison with Dmel:** DE gene trends suggest a **timing delay** of certain processes in D embryos (physiological delay on a **fixed morphological** trajectory).

---

## 9) Biological Meaning
- Embryos destined for diapause **reprogram transcription early** (3d) through **cell-cycle/ecdysone** axes, then **reduce energy metabolism** and **fortify cuticular structures** by 6d, aligning with **metabolic suppression** and **desiccation resistance** before dormancy.  
- Recurrent involvement of **Pepck** and regulated **cell-cycle genes** supports elements of a **conserved diapause toolkit** spanning insect taxa.

---

## 10) Replicability Checklist (Copy-paste)
- **Design:** D (8L:16D) vs ND (16L:8D); **3d** and **6d** pov; **n = 3** replicates per cell → **12 libraries**.  
- **Rearing:** 21 °C, ~80% RH; F9 lab strain.  
- **Sequencing:** Illumina **HiSeq**, paired-end; **two lanes** (3d on lane 1; 6d on lane 2).  
- **Accessions:** **SRA044835**, **SRA051478**; transcriptome at **AlbopictusExpression.org**.  
- **Quantification:** **RSEM**; **TMM** normalization (edgeR).  
- **DE calling:** **limma**; **FDR < 0.05**, **\|log2FC\| > 0.5**.  
- **Enrichment:** **DAVID** (KEGG/GO FAT; FDR < 0.05; representative term/cluster).  
- **Custom gene-set tests:** **GENEMERGE** for insulin, HSP/stress, ecdysone signaling.  
- **Developmental benchmarking:** 1:1 orthologs vs **D. melanogaster** RNA-seq time series; z-score heatmaps.

---

## 11) One-paragraph Wrap
Pre-diapause embryos display **early endocrine/cell-cycle shifts** and **late metabolic suppression + cuticle provisioning**, producing a transcriptomic trajectory consistent with **programmed diapause preparation** rather than simple delay. **Pepck** and **PCNA** emerge as cross-study candidates for a conserved diapause toolkit, while insulin/stress signatures do not dominate at this preparatory stage.

---
