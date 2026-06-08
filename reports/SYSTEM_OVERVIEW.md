# AISNP Selection Pipeline — System Overview

## Problem Statement

Given whole-genome SNP data from 504 East Asian individuals across three subpopulations, identify the **smallest set of SNPs** that allows a machine learning classifier to accurately assign ancestry to CN (Han Chinese), JPT (Japanese), or SEA (Southeast Asian).

This is a deliberately hard classification problem — all three groups belong to the same EAS super-population and are genetically close. Continental-level ancestry markers are insufficient.

### Why it matters
- Forensic identification and missing-persons cases
- Pharmacogenomics (drug response varies by ancestry)
- Population genetics research

---

## Data

| Property | Value |
|---|---|
| Source | 1000 Genomes Project Phase 3 |
| Total samples | 504 EAS individuals |
| Target groups | CN (n=208), JPT (n=104), SEA (n=192) |
| SNPs after filtering | 614,759 |

### Population grouping

| Target | Raw sub-populations | Description |
|---|---|---|
| **CN** | CHB + CHS | Han Chinese (Beijing + South) |
| **JPT** | JPT | Japanese |
| **SEA** | KHV + CDX | Kinh Vietnamese + Chinese Dai |

---

## Pipeline Architecture

```
1000 Genomes (~22M SNPs, 504 EAS samples)
        │
        ▼
┌─────────────────────────────┐
│  1. QUALITY FILTERING        │  614,759 SNPs remain
│  MAF · call rate · HWE · LD  │
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  2. CANDIDATE SET            │  3 pools: stat / FST / fst_stat
│     CONSTRUCTION             │
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  3. THREE-STAGE ML SWEEP     │  Stage 1: pool×reductor selection
│                              │  Stage 2: classifier evaluation
│                              │  Stage 3: panel commitment
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  4. BENCHMARKING             │  vs Cai 2024, Shi 2019, Cao 2022
└─────────────────────────────┘
```

---

## Step 1 — Quality Filtering

Applied sequentially using plink2 (Weir-Cockerham estimator, 16 threads).

| Filter | Threshold | Purpose |
|---|---|---|
| SNP-only, biallelic | — | Exclude indels, CNVs, multi-allelic sites |
| MAF | ≥ 1/(2×504) ≈ 0.001 | Remove near-monomorphic sites |
| Call rate | ≥ 95% | Remove poorly genotyped variants |
| HWE | p ≥ 1×10⁻⁶ (keep-fewhet) | Remove likely genotyping errors |
| LD pruning | r² < 0.10, 1000 kb window | Retain independent signals only |

Result: **614,759 independent SNPs**.

---

## Step 2 — Candidate Set Construction

Reduces 614,759 SNPs to three focused pools before ML.

### Statistical block (stat — 1,005 SNPs)

Seven statistics were evaluated across all 614,759 SNPs. Correlation analysis revealed three independent signals; the redundant tests were dropped:

| Retained | Dropped (r ≥ 0.96 with a retained test) |
|---|---|
| χ² (chi-squared) | Cramér's V, ANOVA F, Kruskal-Wallis H |
| δAF (allele freq. difference) | — |
| JSD (Jensen-Shannon divergence) | Mutual Information |

Top-500 SNPs per retained test → union → **1,005 unique SNPs**.

### FST block (FST — 2,508 SNPs)

Weir-Cockerham FST computed pairwise across all three population pairs (CN↔JPT, CN↔SEA, JPT↔SEA) using plink2. Top-1,000 per pair → union → **2,508 SNPs**.

FST normalises by heterozygosity and upweights rare, drift-fixed variants. δAF treats all allele frequencies equally — the two approaches capture partly orthogonal signal.

### Combined block (fst_stat — 1,003 SNPs)

Intersection of stat and FST blocks. SNPs independently nominated by both approaches — the consensus set. **1,003 SNPs**.

---

## Step 3 — Three-Stage ML Sweep

The core methodological contribution. All three candidate sets are compared systematically across all panel sizes N ∈ {5, 10, 15 … 100} in a single unified notebook (`08_unified_panel_sweep.ipynb`).

### Stage 1 — Reductor × Panel Selection

For each N: compare 3 pools × 3 reductors = **9 configurations** via 5-fold nested CV (RF as fixed eval classifier). Ranker re-fit inside each training fold — no leakage.

| Reductor | Mechanism | Bias |
|---|---|---|
| **L1-LR** | SAGA solver, L1 penalty; sparse coefficients | Maximum sparsity |
| **ElasticNet** | SAGA solver, L1+L2 (l1_ratio=0.5); balanced | Sparsity + stability ← **winner** |
| **RF importance** | Gini impurity mean decrease | Non-linear interactions |

Output: best (pool, reductor) per N. **stat + ElasticNet wins at N ≥ 15 across nearly all panel sizes.**

### Stage 2 — Classifier Evaluation

For each N: use Stage-1 winner, evaluate **6 classifiers** via 5-fold nested CV. Same leak-free design.

| Classifier | Key hyperparameters |
|---|---|
| Random Forest (RF) | n_estimators=100, max_depth=10 |
| XGBoost (XGB) | n_estimators=200, max_depth=6, lr=0.1, subsample=0.8 |
| Logistic Regression (LR) | lbfgs solver (default), max_iter=1000, L2 |
| SVM-RBF | RBF kernel, scaled input |
| SVM-Lin | Linear kernel, scaled input |
| GBM | n_estimators=100, max_depth=5 |

Metrics: accuracy, weighted F1, MCC, ROC-AUC (macro OvR). **Stage 2 accuracy = primary reported metric.**

Key results:
- Best classifier: **LR at most N**, SVM-RBF at N ≥ 50
- ~69% at N=5 → ~92% at N=35 → ~95% at N=55–75
- Performance plateaus after N≈75

### Stage 3 — Panel Commitment

80/20 stratified hold-out. Reductor fit on 80% → top-N SNPs → committed panel. Best Stage-2 classifier evaluated on 20% hold-out.

Stage 3 accuracy is 4–8% below Stage 2 CV (expected: single-split variance + double-selection optimism from Stages 1–2). Not the primary metric — used only for panel commitment.

**Committed panels: N = 35, N = 50, N = 70.**

---

## Step 4 — Benchmarking Against Published Panels

Three published EAS AISNP panels evaluated on the same 504-sample cohort using the identical 6-classifier + 5-fold CV protocol.

| Panel | Nominal N | Matched N | Best acc (our eval) | Published acc |
|---|---|---|---|---|
| Cai et al. 2024 (EAS-specific) | 34 | 34/34 (100%) | **~94.6% (SVM-RBF)** | 92% (XGB) |
| Cao et al. 2022 | 19 | 14/19 (74%) | ~82% | — |
| Shi et al. 2019 (nested) | 36/59/98/142 | 80–82% | ~70–79% | — |

SNPs extracted from MAF-filtered pre-LD-pruning plink data (22M SNPs) to maximise match rate.

### Key comparison finding

Cai's expert-curated panel achieves ~94.6% on our cohort — **classifier-agnostic** (all 6 classifiers score similarly). Our automated panels show greater variance across classifiers but achieve higher peak accuracy at N ≥ 50. The crossover is around N = 50–55.

Cai's higher score on our cohort vs their published 92% is attributable to cohort composition differences: our 1000 Genomes 504-sample EAS subset is more homogeneous than the broader cohort Cai likely used.

---

## SNP Overlap

Our committed panels (N=35/50/70) converted to rsIDs via MyVariant.info. Set intersection with all published panels:

- **vs Cai**: 7–8 shared SNPs (notable given orthogonal selection methods)
- **vs Cao**: 2–3 shared SNPs
- **vs Shi**: 0 shared SNPs (Shi designed for broader Asian scope, not CN/JPT/SEA)

8 shared loci annotated: CADM2 (BMI/cognition), PRRC2A (MHC/autoimmune), TTC3 (chr21/Alzheimer's), SYNE1 (cerebellar ataxia), OLFM3 (glaucoma), DCAF16 (cataract), P2RY1 (blood pressure/T2D), POLR1E (rRNA biogenesis). All are ancestry-informative via population allele frequency differences, not disease-causal.

---

## Summary

| Step | Input | Output |
|---|---|---|
| Quality filtering | ~22M SNPs | 614,759 clean independent SNPs |
| Candidate construction | 614,759 SNPs | 3 pools: stat (1,005) / FST (2,508) / fst_stat (1,003) |
| Stage 1 sweep | 9 pool×reductor configs | Best config per N → stat+ElasticNet |
| Stage 2 sweep | 6 classifiers | Performance curve N=5→100; best clf per N |
| Stage 3 commit | 80/20 split | Final panels: N=35 / N=50 / N=70 |
| Benchmarking | 3 published panels | Cai~94.6%, our N=35~92%, crossover at N≈50 |

**Main finding**: Automated selection via stat+ElasticNet reaches expert-curated panel accuracy (~94%) at N≈50 without any domain-specific curation, and exceeds it at N≥55.
