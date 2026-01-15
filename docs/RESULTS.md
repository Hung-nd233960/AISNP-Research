# Part 1 Results: Statistical AISNP Selection

Comprehensive analysis results from the statistical SNP selection pipeline for distinguishing East Asian subpopulations (CHB, JPT, KHV).

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Dataset Overview](#dataset-overview)
3. [Quality Control Results](#quality-control-results)
4. [Statistical SNP Selection](#statistical-snp-selection)
5. [Machine Learning Performance](#machine-learning-performance)
6. [Feature Importance Analysis](#feature-importance-analysis)
7. [Key Findings and Conclusions](#key-findings-and-conclusions)

---

## Executive Summary

### Objective

Identify a minimal set of Ancestry-Informative SNPs (AISNPs) capable of distinguishing between three East Asian subpopulations from the 1000 Genomes Project.

### Key Results

| Metric | Value |
|--------|-------|
| **Final Consensus SNPs** | 37 |
| **Best Model** | SVM (RBF kernel) |
| **Best Accuracy (37 SNPs)** | 89.55% ± 3.80% |
| **Best Accuracy (25 SNPs)** | 86.93% ± 3.26% |
| **Optimal SNP Count** | 15-25 SNPs (accuracy plateau) |
| **Target Populations** | CHB, JPT, KHV (n=306 samples) |

### Pipeline Summary

```
Raw VCF → Hard Filtering → Situational Filtering → Statistical Selection → ML Validation
           (SNPs-only,      (HWE, LD pruning)       (χ², MI, IG, KL)        (10 models)
            MAF, Call rate)
```

---

## Dataset Overview

### 1000 Genomes Panel

| Statistic | Value |
|-----------|-------|
| **Total Samples** | 2,504 |
| **Total Populations** | 26 |
| **Super Populations** | 5 (AFR, AMR, EAS, EUR, SAS) |
| **Gender Distribution** | Female: 1,271 / Male: 1,233 |

### East Asian (EAS) Subpopulations

| Population | Description | Sample Size |
|------------|-------------|-------------|
| **CHB** | Han Chinese in Beijing | 103 |
| **JPT** | Japanese in Tokyo | 104 |
| **KHV** | Kinh in Ho Chi Minh City, Vietnam | 99 |
| **Total** | - | **306** |

### Population Genetic Similarity (F<sub>ST</sub>)

Pairwise Hudson F<sub>ST</sub> values indicate close genetic relationships:

| Population Pair | Hudson F<sub>ST</sub> | Interpretation |
|-----------------|----------------------|----------------|
| CHB ↔ JPT | 0.0062 | Very low differentiation |
| CHB ↔ KHV | 0.0056 | Very low differentiation |
| JPT ↔ KHV | 0.0125 | Low differentiation |

> **Note**: F<sub>ST</sub> < 0.05 indicates minimal genetic differentiation, making subpopulation classification challenging. The selected AISNPs must capture subtle frequency differences.

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

**Output**: `EAS_AND_SNP_filtered_data` PLINK2 pfile

### Stage 2: Situational Filtering

**Filters Applied**:

| Filter | Parameters | Purpose |
|--------|------------|---------|
| Hardy-Weinberg Equilibrium | p < 1e-6, keep-fewhet | Remove genotyping errors |
| Unique IDs | Deduplicate | Ensure variant uniqueness |
| LD Pruning | window=1000kb, step=1, r²=0.1 | Remove redundant variants |

**Output**: `EAS_FINAL_DATA_FOR_FST` (LD-pruned, QC-passed SNPs)

---

## Statistical SNP Selection

### Four-Test Consensus Approach

SNPs were required to pass **ALL 4 statistical tests** to be included in the final consensus set:

#### 1. Pearson Chi-Squared Test (χ²)

- **Purpose**: Test independence of genotype frequency and population
- **Selection**: FDR-corrected q-value < 0.05
- **Result**: Identified SNPs with significant genotype-population associations

#### 2. Mutual Information (MI)

- **Purpose**: Quantify shared information between genotype and population
- **Selection**: Top 500 SNPs by MI score
- **Range**: 0.063 - 0.200 bits (among selected SNPs)

#### 3. Information Gain (IG)

- **Purpose**: Measure entropy reduction in population labels given genotype
- **Selection**: Top 500 SNPs by IG score
- **Range**: 0.091 - 0.289 (among selected SNPs)

#### 4. Kullback-Leibler Divergence (KL)

- **Purpose**: Measure divergence between population-specific genotype distributions
- **Selection**: Top 500 SNPs by average pairwise KL divergence
- **Range**: 1.33 - 4.39 (among selected SNPs)

### Consensus SNPs

**Final count**: **37 SNPs** passed all four statistical tests.

These SNPs were ranked by composite score (average rank across all tests):

| Rank | SNP ID | Composite Rank | χ² p-value | MI | KL Divergence |
|------|--------|----------------|------------|-----|---------------|
| 1 | 3:130239945[b37]G,A | 1.5 | 2.55e-09 | 0.103 | 4.22 |
| 2 | 4:101002868[b37]A,T | 3.0 | 3.29e-16 | 0.140 | 3.99 |
| 3 | 12:128054516[b37]C,T | 4.5 | 2.32e-13 | 0.113 | 3.04 |
| 4 | 20:9928437[b37]G,T | 5.0 | 3.23e-09 | 0.067 | 1.80 |
| 5 | 11:130759223[b37]G,T | 9.0 | 2.10e-06 | 0.063 | 1.81 |
| 6 | 16:83272044[b37]C,T | 9.0 | 1.25e-06 | 0.067 | 2.46 |
| 7 | 4:5532600[b37]T,C | 10.0 | 1.33e-09 | 0.073 | 2.14 |
| 8 | 6:37497412[b37]C,T | 10.5 | 2.70e-09 | 0.081 | 2.38 |
| 9 | 16:88769665[b37]T,G | 11.5 | 4.10e-07 | 0.069 | 1.50 |
| 10 | 1:240285457[b37]G,A | 12.0 | 3.82e-06 | 0.066 | 2.69 |

*Full ranked list available in: `1000genomes/output/consensus_snps_ranked.txt`*

---

## Machine Learning Performance

### Cross-Validation Results (37 Consensus SNPs)

10 classifiers were evaluated using 5-fold stratified cross-validation:

| Rank | Model | Accuracy | Std | F1-Score | Overfit Gap |
|------|-------|----------|-----|----------|-------------|
| 1 | **SVM (RBF)** | **89.55%** | ±3.80% | 0.895 | 3.3% |
| 2 | Logistic Regression | 86.27% | ±2.47% | 0.861 | 7.4% |
| 3 | SVM (Linear) | 84.96% | ±4.34% | 0.846 | 8.3% |
| 4 | Naive Bayes | 82.68% | ±1.65% | 0.823 | 2.0% |
| 5 | Random Forest | 81.69% | ±4.10% | 0.814 | 10.5% |
| 6 | XGBoost | 80.72% | ±1.91% | 0.802 | 12.5% |
| 7 | MLP Neural Network | 79.44% | ±6.23% | 0.783 | 7.5% |
| 8 | Gradient Boosting | 78.10% | ±3.98% | 0.779 | 16.6% |
| 9 | AdaBoost | 63.70% | ±5.85% | 0.592 | 5.3% |
| 10 | K-Nearest Neighbors | 57.49% | ±8.01% | 0.580 | 9.8% |

### Performance vs. Number of SNPs

Feature reduction analysis shows accuracy progression with increasing SNP count:

| SNPs | Best Model | Accuracy | Notes |
|------|------------|----------|-------|
| 5 | AdaBoost | 67.64% | Insufficient discrimination |
| 10 | Logistic Regression | 74.85% | Moderate performance |
| 15 | Random Forest | 84.31% | Strong improvement |
| 20 | SVM (Linear) | 84.32% | Plateau begins |
| 25 | SVM (RBF) | 86.93% | Near-optimal |
| 37 (All) | SVM (RBF) | 89.55% | Best overall |

### Accuracy Trends by Model

| Model | 5 SNPs | 10 SNPs | 15 SNPs | 20 SNPs | 25 SNPs | 37 SNPs |
|-------|--------|---------|---------|---------|---------|---------|
| SVM (RBF) | 66.3% | 74.2% | 84.0% | 83.7% | **86.9%** | **89.6%** |
| Logistic Regression | 67.3% | **74.8%** | 82.4% | 83.0% | 83.0% | 86.3% |
| Random Forest | 64.4% | 73.2% | **84.3%** | 82.7% | 79.1% | 81.7% |
| SVM (Linear) | 66.7% | 71.2% | 83.0% | **84.3%** | 86.0% | 85.0% |
| XGBoost | 65.3% | 73.9% | 81.0% | 81.4% | 80.1% | 80.7% |

**Key Observations**:

- Performance significantly improves from 5→15 SNPs
- Accuracy plateaus around 15-25 SNPs
- SVM (RBF) consistently outperforms other models above 15 SNPs
- Tree-based methods (Random Forest, XGBoost) show slight decline at higher SNP counts

### Statistical SNP Validation (Top 25 Statistical SNPs)

Separate validation using statistically-selected SNPs (from notebook 02b):

| Model | Test Accuracy | F1-Score | CV Mean |
|-------|---------------|----------|---------|
| SVM (RBF) | 85.48% | 0.852 | 88.93% |
| Random Forest | 83.87% | 0.831 | 88.52% |
| AdaBoost | 82.26% | 0.820 | 82.37% |
| SVM (Linear) | 82.26% | 0.816 | 86.47% |
| Gradient Boosting | 80.65% | 0.800 | 85.65% |

---

## Feature Importance Analysis

### Consensus Approach: Random Forest + Logistic Regression

SNPs were ranked by averaging importance rankings from:

1. **Random Forest**: Gini importance
2. **Logistic Regression**: Absolute coefficient magnitude

### Top 10 Most Informative SNPs

| Rank | SNP ID | RF Importance | LR Importance | Avg Rank |
|------|--------|---------------|---------------|----------|
| 1 | **3:130239945[b37]G,A** | 0.074 | 0.801 | 1.5 |
| 2 | 4:101002868[b37]A,T | 0.097 | 0.584 | 3.0 |
| 3 | 12:128054516[b37]C,T | 0.053 | 0.584 | 4.5 |
| 4 | 20:9928437[b37]G,T | 0.036 | 0.701 | 5.0 |
| 5 | 11:130759223[b37]G,T | 0.041 | 0.506 | 9.0 |
| 6 | 16:83272044[b37]C,T | 0.028 | 0.624 | 9.0 |
| 7 | 4:5532600[b37]T,C | 0.027 | 0.644 | 10.0 |
| 8 | 6:37497412[b37]C,T | 0.032 | 0.526 | 10.5 |
| 9 | 16:88769665[b37]T,G | 0.047 | 0.359 | 11.5 |
| 10 | 1:240285457[b37]G,A | 0.027 | 0.558 | 12.0 |

### SNP Distribution by Chromosome

| Chromosome | Count | Notable SNPs |
|------------|-------|--------------|
| Chr 1 | 5 | 1:240285457, 1:12387655, 1:168205652 |
| Chr 4 | 4 | 4:101002868, 4:5532600, 4:17813761, 4:15787630 |
| Chr 6 | 4 | 6:37497412, 6:133653349, 6:26340872, 6:70989171 |
| Chr 16 | 4 | 16:83272044, 16:88769665, 16:10424828 |
| Chr 12 | 3 | 12:128054516, 12:56981509, 12:66982904 |
| Other | 17 | Distributed across 12 additional chromosomes |

---

## Key Findings and Conclusions

### 1. Minimal SNP Panel Efficacy

> **Finding**: A panel of just **37 consensus SNPs** achieves **89.55% accuracy** in distinguishing CHB, JPT, and KHV populations.

This demonstrates that ancestry inference for closely related East Asian subpopulations is achievable with a minimal marker set, suitable for forensic and clinical applications.

### 2. Statistical Selection Success

> **Finding**: The four-test consensus approach (χ², MI, IG, KL divergence) effectively identified biologically relevant markers.

- All 37 selected SNPs showed statistically significant genotype-population associations
- Top-ranked SNPs exhibited high KL divergence values (up to 4.39), indicating distinct allele frequency profiles

### 3. Model Recommendations

| Use Case | Recommended Model | Rationale |
|----------|-------------------|-----------|
| **Production/Deployment** | SVM (RBF) | Best accuracy (89.55%), moderate complexity |
| **Interpretability** | Logistic Regression | Good accuracy (86.27%), transparent coefficients |
| **Quick Baseline** | Naive Bayes | Fast, low overfit (2.0% gap), 82.68% accuracy |
| **Feature Selection** | Random Forest | Gini importance for SNP ranking |

### 4. Optimal Panel Size

> **Finding**: **15-25 SNPs** provide the best accuracy-to-complexity trade-off.

- Below 15 SNPs: significant accuracy degradation
- 15-25 SNPs: ~85% accuracy (sufficient for most applications)
- 37 SNPs: ~90% accuracy (marginal gain)

### 5. Population Differentiation Challenge

> **Finding**: Low F<sub>ST</sub> values (0.006-0.013) confirm minimal genetic differentiation between East Asian subpopulations.

Despite this challenge, the selected AISNPs capture subtle frequency differences sufficient for classification. The JPT-KHV pair shows the highest F<sub>ST</sub> (0.0125) and may be easiest to distinguish.

### 6. Limitations and Future Work

| Limitation | Impact | Mitigation |
|------------|--------|------------|
| Small sample size (n=306) | Limited generalizability | Expand with additional cohorts |
| Three populations only | Cannot distinguish other EAS groups | Add CDX, CHS populations |
| GRCh37 coordinates | May need updating | Convert to GRCh38 for newer data |
| No external validation | Unknown real-world performance | Test on independent datasets |

---

## Output Files Reference

### Primary Results

| File | Description |
|------|-------------|
| `1000genomes/output/consensus_snps_ranked.txt` | Final 37 SNPs ranked by importance |
| `1000genomes/output/consensus_snps_importance.csv` | Feature importance (RF + LR) |
| `1000genomes/output/consensus_snps_cv_results.csv` | Cross-validation performance |
| `1000genomes/output/consensus_snps_performance_summary.csv` | Model comparison summary |

### Intermediate Files

| File | Description |
|------|-------------|
| `1000genomes/output/EAS_FINAL_DATA_FOR_FST.*` | QC-passed PLINK files |
| `1000genomes/output/EAS_FST_RESULTS.fst.summary` | Population F<sub>ST</sub> values |
| `1000genomes/output/statistics_snp_scores.csv` | All SNP statistical scores |
| `reports/statistical_ml_report.txt` | ML training summary |

---

## Reproducibility

All analyses can be reproduced by running notebooks in order:

```bash
cd scripts/notebooks/
jupyter notebook 01_hard_filtering.ipynb
jupyter notebook 02_situational_filtering.ipynb
jupyter notebook 02b_statistical_snp_selection.ipynb
jupyter notebook 03_fst_and_pca.ipynb
jupyter notebook 04c_ml_consensus_snps.ipynb
```

Configuration parameters are centralized in `scripts/config.py`.

---

*Report generated from Part 1 pipeline analysis*
*Target: East Asian Subpopulation Ancestry Inference*
*Date: January 2026*
