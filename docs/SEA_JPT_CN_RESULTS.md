# SEA-JPT-CN Population Analysis Results

Comprehensive analysis results for distinguishing Chinese (CN), Japanese (JPT), and Southeast Asian (SEA) populations using Ancestry-Informative SNPs.

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Dataset Overview](#dataset-overview)
3. [Quality Control Results](#quality-control-results)
4. [Statistical SNP Selection](#statistical-snp-selection)
5. [Machine Learning Performance](#machine-learning-performance)
6. [Optimal SNP Panel Analysis](#optimal-snp-panel-analysis)
7. [Key Findings](#key-findings)

---

## Executive Summary

### Objective

Identify a minimal set of Ancestry-Informative SNPs (AISNPs) capable of distinguishing between three population groups: Chinese (CN), Japanese (JPT), and Southeast Asian (SEA).

### Key Results

| Metric | Value |
|--------|-------|
| **Final Consensus SNPs** | 31 |
| **Best Model** | SVM (RBF kernel) |
| **Best Accuracy (50 SNPs)** | 93.07% (Test), 96.28% (CV) |
| **Best Accuracy (31 SNPs)** | 84.72% ± 2.94% |
| **Best Accuracy (45 SNPs)** | 96.27% (TabPFN CV) |
| **Optimal SNP Count** | 35-50 SNPs |
| **Target Populations** | CN, JPT, SEA |

### Pipeline Summary

```
Raw VCF → Hard Filtering → Situational Filtering → Statistical Selection → ML Validation
           (SNPs-only,      (HWE, LD pruning)       (χ², MI, IG, KL)        (10+ models)
            MAF, Call rate)
```

---

## Dataset Overview

### Target Populations

| Population Code | Description | Sample Size |
|-----------------|-------------|-------------|
| **CN** | Chinese | - |
| **JPT** | Japanese | - |
| **SEA** | Southeast Asian | - |

### Population Genetic Similarity (F<sub>ST</sub>)

Pairwise Hudson F<sub>ST</sub> values indicate close genetic relationships:

| Population Pair | Hudson F<sub>ST</sub> | Interpretation |
|-----------------|----------------------|----------------|
| CN ↔ JPT | 0.0069 | Very low differentiation |
| CN ↔ SEA | 0.0044 | Very low differentiation |
| JPT ↔ SEA | 0.0134 | Low differentiation |

> **Note**: F<sub>ST</sub> values between 0.004-0.013 indicate minimal genetic differentiation, making subpopulation classification challenging. The selected AISNPs must capture subtle frequency differences.

---

## Quality Control Results

### Stage 1: Hard Filtering

**Filters Applied**:

| Filter | Parameters | Purpose |
|--------|------------|---------|
| SNP-only | `--snps-only` | Focus on single nucleotide variants |
| Biallelic | `--max-alleles 2` | Ensure computational simplicity |
| MAF | `--maf 0.0016` | Remove singletons/rare variants |
| Call rate | `--geno 0.05` | ≥95% genotyping completeness |

### Stage 2: Situational Filtering

| Filter | Parameters | Purpose |
|--------|------------|---------|
| Hardy-Weinberg Equilibrium | p < 1e-6, keep-fewhet | Remove genotyping errors |
| Unique IDs | Deduplicate | Ensure variant uniqueness |
| LD Pruning | window=1000kb, step=1, r²=0.1 | Remove redundant variants |

---

## Statistical SNP Selection

### Four-Test Consensus Approach

SNPs were required to pass **ALL 4 statistical tests** to be included in the final consensus set:

1. **Pearson Chi-Squared Test (χ²)**: FDR-corrected q-value < 0.05
2. **Mutual Information (MI)**: Top 500 SNPs
3. **Information Gain (IG)**: Top 500 SNPs
4. **Kullback-Leibler Divergence (KL)**: Top 500 SNPs

### Consensus SNPs

**Final count**: **31 SNPs** passed all four statistical tests.

#### Top 10 Consensus SNPs by Composite Rank

| Rank | SNP ID | χ² p-value | MI | KL Divergence |
|------|--------|------------|-----|---------------|
| 1 | 4:17813761[b37]G,A | 9.58e-39 | 0.192 | 3.88 |
| 2 | 1:12387655[b37]G,A | 1.84e-35 | 0.174 | 3.87 |
| 3 | 5:41181491[b37]G,T | 1.44e-31 | 0.155 | 3.28 |
| 4 | 14:96938945[b37]A,T | 7.55e-31 | 0.161 | 2.86 |
| 5 | 11:112053732[b37]T,A | 5.88e-23 | 0.119 | 2.15 |
| 6 | 1:168205652[b37]G,A | 1.10e-21 | 0.105 | 2.41 |
| 7 | 7:156743124[b37]G,T | 3.69e-17 | 0.087 | 3.23 |
| 8 | 16:46417894[b37]G,A | 7.71e-19 | 0.101 | 1.53 |
| 9 | 20:5547557[b37]A,T | 2.32e-17 | 0.085 | 1.72 |
| 10 | 1:79095386[b37]A,T | 2.58e-17 | 0.089 | 1.63 |

#### Full Ranked Consensus SNP List (31 SNPs)

| Rank | SNP ID |
|------|--------|
| 1 | 3:130239945[b37]G,A |
| 2 | 16:46417894[b37]G,A |
| 3 | 1:12387655[b37]G,A |
| 4 | 4:17813761[b37]G,A |
| 5 | 10:3043329[b37]A,G |
| 6 | 5:41181491[b37]G,T |
| 7 | 16:88769665[b37]T,G |
| 8 | 6:37497412[b37]C,T |
| 9 | 10:21179885[b37]C,T |
| 10 | 8:142363071[b37]T,C |
| 11 | 2:123492479[b37]A,G |
| 12 | 14:96938945[b37]A,T |
| 13 | 8:95173897[b37]T,C |
| 14 | 1:143543213[b37]G,A |
| 15 | 1:213582667[b37]C,A |
| 16 | 14:104193238[b37]G,C |
| 17 | 9:127101899[b37]C,A |
| 18 | 16:56678700[b37]G,C |
| 19 | 20:29647723[b37]G,A |
| 20 | 18:60094992[b37]T,C |
| 21 | 11:112053732[b37]T,A |
| 22 | 7:156743124[b37]G,T |
| 23 | 1:6133700[b37]C,G |
| 24 | 4:173803959[b37]A,G |
| 25 | 14:106237499[b37]G,A |
| 26 | 1:168205652[b37]G,A |
| 27 | 16:55844911[b37]G,A |
| 28 | 1:79095386[b37]A,T |
| 29 | 2:238277264[b37]T,G |
| 30 | 20:5547557[b37]A,T |
| 31 | 9:127990148[b37]A,T |

---

## Machine Learning Performance

### Cross-Validation Results (31 Consensus SNPs)

10 classifiers evaluated using 5-fold stratified cross-validation:

| Rank | Model | Accuracy | F1-Score |
|------|-------|----------|----------|
| 1 | XGBoost | 84.52% | 0.845 |
| 2 | Logistic Regression | 84.71% | 0.847 |
| 3 | SVM (RBF) | 82.93% | 0.829 |
| 4 | Random Forest | 82.73% | 0.827 |
| 5 | SVM (Linear) | 82.14% | 0.821 |
| 6 | MLP Neural Network | 80.95% | 0.809 |
| 7 | Gradient Boosting | 80.35% | 0.803 |
| 8 | AdaBoost | 78.57% | 0.786 |
| 9 | Naive Bayes | 70.82% | 0.708 |
| 10 | K-Nearest Neighbors | 70.03% | 0.700 |

---

## Optimal SNP Panel Analysis

### Performance vs. Number of SNPs (25-50 SNPs)

#### Accuracy by SNP Count (Test Set)

| N_SNPs | Random Forest | XGBoost | Logistic Regression | SVM (RBF) | Gradient Boosting | TabPFN |
|--------|---------------|---------|---------------------|-----------|-------------------|--------|
| 25 | 82.18% | 82.18% | 82.18% | 84.16% | 80.20% | 81.19% |
| 30 | 89.11% | 87.13% | 82.18% | 90.10% | 86.14% | 88.12% |
| 35 | 91.09% | 92.08% | 90.10% | 92.08% | 86.14% | 89.11% |
| 40 | 90.10% | 91.09% | 91.09% | 93.07% | 90.10% | 93.07% |
| 45 | 92.08% | 91.09% | 92.08% | 92.08% | 92.08% | 91.09% |
| 50 | 92.08% | 92.08% | 91.09% | **93.07%** | 93.07% | 91.09% |

#### Cross-Validation Accuracy by SNP Count

| N_SNPs | Random Forest | XGBoost | Logistic Regression | SVM (RBF) | Gradient Boosting | TabPFN |
|--------|---------------|---------|---------------------|-----------|-------------------|--------|
| 25 | 88.10% | 88.35% | 89.08% | 90.07% | 84.88% | 90.83% |
| 30 | 89.33% | 90.08% | 89.83% | 91.31% | 88.84% | **92.81%** |
| 35 | 91.82% | 90.57% | 92.31% | 93.55% | 90.33% | 92.07% |
| 40 | 92.31% | 91.32% | 93.56% | 94.05% | 90.83% | **94.79%** |
| 45 | 92.29% | 91.56% | 94.54% | 95.04% | 90.59% | **96.27%** |
| 50 | 93.55% | 91.56% | 94.80% | **96.28%** | 91.82% | 96.03% |

### Best Configuration Analysis

| Metric | Best Model | N_SNPs | Score |
|--------|------------|--------|-------|
| **Best Test Accuracy** | SVM (RBF) / Gradient Boosting | 50 | 93.07% |
| **Best CV Accuracy** | SVM (RBF) | 50 | 96.28% |
| **Best CV (45 SNPs)** | TabPFN | 45 | 96.27% |
| **Best CV (40 SNPs)** | TabPFN | 40 | 94.79% |

### Model Performance Summary (Averaged Across All SNP Counts)

| Model | Mean Accuracy | Std | Max Accuracy | Best SNP Count |
|-------|---------------|-----|--------------|----------------|
| SVM (RBF) | 91.32% | 3.65% | 93.07% | 50 |
| Logistic Regression | 89.80% | 4.23% | 92.08% | 45 |
| TabPFN | 89.63% | 4.31% | 93.07% | 40 |
| XGBoost | 89.28% | 3.95% | 92.08% | 35, 50 |
| Random Forest | 89.44% | 3.93% | 92.08% | 45, 50 |
| Gradient Boosting | 87.96% | 4.71% | 93.07% | 50 |

---

## Key Findings

### 1. Population Differentiation Challenge

- F<sub>ST</sub> values range from 0.0044 to 0.0134, indicating **minimal genetic differentiation**
- CN-SEA show the lowest differentiation (F<sub>ST</sub> = 0.0044)
- JPT-SEA show the highest differentiation (F<sub>ST</sub> = 0.0134)

### 2. Optimal Panel Size

- **31 consensus SNPs** achieve ~85% accuracy
- **35 SNPs** show strong improvement (>90% accuracy)
- **40-50 SNPs** reach peak performance (~93-96%)
- **Recommended panel size**: 40-50 SNPs for optimal accuracy

### 3. Best Performing Models

1. **SVM (RBF)**: Most consistent performer, best at 50 SNPs (93.07% test, 96.28% CV)
2. **TabPFN**: Excellent CV performance, especially at 40-45 SNPs
3. **Logistic Regression**: Stable and interpretable
4. **XGBoost**: Strong overall performance

### 4. Accuracy per SNP Efficiency

| SNP Count | Best Accuracy | Accuracy per SNP |
|-----------|---------------|------------------|
| 31 | 84.72% | 2.73% |
| 35 | 93.55% | 2.67% |
| 40 | 94.79% | 2.37% |
| 45 | 96.27% | 2.14% |
| 50 | 96.28% | 1.93% |

### 5. Conclusions

1. The four-test statistical consensus approach effectively identifies population-discriminating variants despite low F<sub>ST</sub> values
2. **40-50 SNPs** provide optimal accuracy (>93%) for CN/JPT/SEA classification
3. SVM (RBF) and TabPFN are the recommended classifiers
4. The panel offers a cost-effective solution for East Asian population substructure analysis

---

## Output Files

| File | Description |
|------|-------------|
| `statistical_all4_snps_02b.csv` | Full statistics for consensus SNPs |
| `consensus_snps_ranked.txt` | Ranked consensus SNP list |
| `consensus_snps_cv_results.csv` | Cross-validation results |
| `varying_snp_results.csv` | Performance by SNP count |
| `varying_snp_performance.png` | Performance visualization |

---

*Analysis conducted using the AISNP-Research pipeline. Data from 1000 Genomes Project Phase 3.*
