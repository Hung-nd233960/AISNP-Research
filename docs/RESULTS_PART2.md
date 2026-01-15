# Part 2 Results: Known AISNP Panel Comparison

Comprehensive comparison of known AISNP panels from research publications and commercial products for East Asian subpopulation classification (CHB, JPT, KHV).

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Dataset Overview](#dataset-overview)
3. [SNP Panel Sources](#snp-panel-sources)
4. [Machine Learning Performance](#machine-learning-performance)
5. [Cross-Source Comparison](#cross-source-comparison)
6. [Key Findings](#key-findings)
7. [Conclusions and Recommendations](#conclusions-and-recommendations)

---

## Executive Summary

### Objective

Compare classification performance of known AISNP panels (published research and commercial products) against our statistically-selected SNPs for distinguishing East Asian subpopulations.

### Key Results

| Metric | Value |
|--------|-------|
| **AISNP Sources Compared** | 6 |
| **Best Performing Source** | cal_et_al (52 SNPs) |
| **Best Accuracy** | 91.50% ± 4.33% (XGBoost) |
| **Best Statistical SNPs** | 89.55% ± 3.80% (SVM-RBF) |
| **Classifiers Tested** | 7 |
| **Cross-Validation** | 5-Fold Stratified |
| **Target Populations** | CHB, JPT, KHV (n=306) |

### Performance Ranking by Source

| Rank | Source | SNPs | Best Accuracy | Best Model |
|------|--------|------|---------------|------------|
| 1 | **cal_et_al** | 52 | **91.50%** | XGBoost |
| 2 | **statistical_all4** | 37 | **89.55%** | SVM (RBF) |
| 3 | hsiao_lin_hwa | 125 | 64.38% | SVM (RBF) |
| 4 | seldin_128 | 124 | 64.05% | SVM (RBF) |
| 5 | forenseq | 55 | 60.80% | SVM (RBF) |
| 6 | kidd_55 | 53 | 60.78% | Logistic Regression |

---

## Dataset Overview

### Sample Information

All AISNP panels were tested on the same sample set:

| Statistic | Value |
|-----------|-------|
| **Total Samples** | 306 |
| **Populations** | 3 (CHB, JPT, KHV) |
| **Source Data** | 1000 Genomes Project Phase 3 |
| **Genome Build** | GRCh37 |

### Population Distribution

| Population | Description | Sample Size |
|------------|-------------|-------------|
| **CHB** | Han Chinese in Beijing | 103 |
| **JPT** | Japanese in Tokyo | 104 |
| **KHV** | Kinh in Ho Chi Minh City, Vietnam | 99 |

---

## SNP Panel Sources

### Overview Table

| Source | SNPs | Samples | Populations | Description |
|--------|------|---------|-------------|-------------|
| cal_et_al | 52 | 306 | CHB, JPT, KHV | Cal_et_al (2025) AISNP panel |
| seldin_128 | 124 | 306 | CHB, JPT, KHV | Seldin et al. 128-marker panel |
| forenseq | 55 | 306 | CHB, JPT, KHV | ForenSeq commercial panel |
| kidd_55 | 53 | 306 | CHB, JPT, KHV | Kidd lab 55-marker panel |
| hsiao_lin_hwa | 125 | 306 | CHB, JPT, KHV | Hsiao/Lin/Hwa research panel |
| statistical_all4 | 37 | 306 | CHB, JPT, KHV | Our statistical selection (Part 1) |

### Panel Characteristics

#### 1. cal_et_al (52 SNPs)

- **Origin**: Calibration study for East Asian ancestry
- **Design**: Specifically selected for East Asian subpopulation differentiation
- **Performance**: **Best overall** (91.50% accuracy)

#### 2. statistical_all4 (37 SNPs)

- **Origin**: Our Part 1 statistical pipeline
- **Design**: Four-test consensus (χ², MI, IG, KL divergence)
- **Performance**: **Second best** (89.55% accuracy)

#### 3. seldin_128 (124 SNPs)

- **Origin**: Seldin et al. global ancestry panel
- **Design**: Continental-level ancestry differentiation
- **Performance**: Moderate (64.05% accuracy)

#### 4. hsiao_lin_hwa (125 SNPs)

- **Origin**: Research study panel
- **Design**: General ancestry inference
- **Performance**: Moderate (64.38% accuracy)

#### 5. forenseq (55 SNPs)

- **Origin**: Verogen ForenSeq commercial kit
- **Design**: Forensic ancestry inference
- **Performance**: Low for EAS subpopulations (60.80%)

#### 6. kidd_55 (53 SNPs)

- **Origin**: Kidd Lab 55 AISNP panel
- **Design**: Global biogeographic ancestry
- **Performance**: Low for EAS subpopulations (60.78%)

---

## Machine Learning Performance

### Classifiers Tested

| Classifier | Configuration |
|------------|---------------|
| Random Forest | n_estimators=100, max_depth=10 |
| XGBoost | n_estimators=100, max_depth=5, lr=0.1 |
| Logistic Regression | max_iter=1000, multinomial |
| SVM (RBF) | kernel=rbf, probability=True |
| K-Nearest Neighbors | n_neighbors=5 |
| Gradient Boosting | n_estimators=100, max_depth=5 |
| MLP Neural Network | hidden_layers=(100,50), early_stopping |

### Cross-Validation Results by Source

#### cal_et_al (52 SNPs) - **Best Performer**

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **XGBoost** | **91.50%** | ±4.33% | 0.914 | 8.5% |
| SVM (RBF) | 90.84% | ±4.36% | 0.909 | 8.5% |
| Random Forest | 90.85% | ±3.21% | 0.908 | 9.1% |
| Logistic Regression | 90.21% | ±4.24% | 0.901 | 9.8% |
| Gradient Boosting | 89.20% | ±3.07% | 0.891 | 10.8% |
| K-Nearest Neighbors | 84.31% | ±2.47% | 0.843 | 4.3% |
| MLP Neural Network | 83.64% | ±5.42% | 0.836 | 7.1% |

#### statistical_all4 (37 SNPs) - **Our Statistical Selection**

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **SVM (RBF)** | **89.55%** | ±3.80% | 0.895 | 3.3% |
| Logistic Regression | 86.27% | ±2.47% | 0.861 | 7.4% |
| Random Forest | 81.69% | ±4.10% | 0.814 | 10.5% |
| XGBoost | 80.72% | ±1.91% | 0.802 | 12.5% |
| MLP Neural Network | 79.44% | ±6.23% | 0.783 | 7.5% |
| Gradient Boosting | 78.10% | ±3.98% | 0.779 | 16.6% |
| K-Nearest Neighbors | 57.49% | ±8.01% | 0.580 | 9.8% |

#### seldin_128 (124 SNPs)

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **SVM (RBF)** | **64.05%** | ±8.12% | 0.640 | 36.0% |
| Logistic Regression | 55.88% | ±4.69% | 0.558 | 44.1% |
| XGBoost | 52.93% | ±4.13% | 0.527 | 47.1% |
| Random Forest | 50.97% | ±4.75% | 0.506 | 49.0% |
| Gradient Boosting | 49.03% | ±6.18% | 0.489 | 50.9% |
| MLP Neural Network | 48.69% | ±3.79% | 0.482 | 35.1% |
| K-Nearest Neighbors | 46.72% | ±3.73% | 0.456 | 17.3% |

#### hsiao_lin_hwa (125 SNPs)

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **SVM (RBF)** | **64.38%** | ±5.19% | 0.645 | 35.0% |
| Random Forest | 57.19% | ±4.47% | 0.568 | 42.8% |
| Gradient Boosting | 55.86% | ±3.62% | 0.554 | 44.1% |
| XGBoost | 54.56% | ±3.96% | 0.543 | 45.4% |
| Logistic Regression | 52.59% | ±5.36% | 0.526 | 47.4% |
| MLP Neural Network | 48.08% | ±9.23% | 0.477 | 42.6% |
| K-Nearest Neighbors | 38.57% | ±9.09% | 0.358 | 20.2% |

#### forenseq (55 SNPs)

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **SVM (RBF)** | **60.80%** | ±3.40% | 0.611 | 31.2% |
| Logistic Regression | 60.45% | ±6.68% | 0.605 | 20.0% |
| Random Forest | 60.13% | ±1.74% | 0.600 | 39.9% |
| Gradient Boosting | 57.53% | ±4.36% | 0.577 | 42.5% |
| XGBoost | 55.90% | ±3.83% | 0.560 | 44.1% |
| MLP Neural Network | 51.27% | ±10.04% | 0.492 | 14.3% |
| K-Nearest Neighbors | 49.69% | ±3.61% | 0.499 | 18.6% |

#### kidd_55 (53 SNPs)

| Model | Accuracy | Std | F1 | Overfit Gap |
|-------|----------|-----|-----|-------------|
| **Logistic Regression** | **60.78%** | ±7.62% | 0.607 | 19.2% |
| Random Forest | 59.48% | ±2.82% | 0.592 | 40.4% |
| SVM (RBF) | 59.49% | ±3.36% | 0.599 | 32.6% |
| Gradient Boosting | 57.53% | ±2.87% | 0.573 | 42.5% |
| XGBoost | 55.89% | ±2.77% | 0.560 | 44.1% |
| K-Nearest Neighbors | 48.39% | ±4.42% | 0.491 | 19.3% |
| MLP Neural Network | 48.34% | ±7.87% | 0.479 | 24.1% |

---

## Cross-Source Comparison

### Summary Statistics by Source

| Source | N_SNPs | Mean Acc | Max Acc | Std | Mean F1 | Max F1 | Overfit Gap |
|--------|--------|----------|---------|-----|---------|--------|-------------|
| **cal_et_al** | 52 | 88.65% | **91.50%** | 3.28% | 0.886 | 0.914 | 8.3% |
| **statistical_all4** | 37 | 79.04% | **89.55%** | 10.31% | 0.787 | 0.895 | 9.7% |
| hsiao_lin_hwa | 125 | 53.03% | 64.38% | 8.07% | 0.525 | 0.645 | 39.6% |
| forenseq | 55 | 56.54% | 60.80% | 4.52% | 0.563 | 0.611 | 30.1% |
| kidd_55 | 53 | 55.70% | 60.78% | 5.25% | 0.557 | 0.607 | 31.8% |
| seldin_128 | 124 | 52.61% | 64.05% | 5.88% | 0.523 | 0.640 | 39.9% |

### SNP Count vs Accuracy

| Source | SNPs | Best Accuracy | Efficiency (Acc/SNP) |
|--------|------|---------------|----------------------|
| statistical_all4 | 37 | 89.55% | **2.42%** |
| cal_et_al | 52 | 91.50% | 1.76% |
| forenseq | 55 | 60.80% | 1.11% |
| kidd_55 | 53 | 60.78% | 1.15% |
| seldin_128 | 124 | 64.05% | 0.52% |
| hsiao_lin_hwa | 125 | 64.38% | 0.51% |

> **Key Insight**: Our statistical selection (37 SNPs) achieves the highest efficiency, obtaining near-best accuracy with the fewest markers.

---

## Key Findings

### 1. Panel Design Matters More Than Size

> **Finding**: The cal_et_al (52 SNPs) and statistical_all4 (37 SNPs) panels dramatically outperform larger panels (124-125 SNPs).

Panels designed specifically for East Asian subpopulation differentiation vastly outperform general ancestry panels, even with fewer markers.

### 2. Our Statistical Selection is Highly Effective

> **Finding**: With only 37 SNPs, our statistical selection achieves 89.55% accuracy, comparable to the best published panel.

The four-test consensus approach (χ², MI, IG, KL divergence) successfully identified markers that differentiate CHB, JPT, and KHV populations.

### 3. General Ancestry Panels Underperform for Subpopulations

> **Finding**: Panels designed for continental-level ancestry (seldin_128, hsiao_lin_hwa, kidd_55, forenseq) perform poorly (~60-64%) for East Asian subpopulation classification.

These panels were designed to distinguish major continental groups (Africa, Europe, Asia, Americas) and lack the resolution for within-region differentiation.

### 4. Model Consistency

> **Finding**: SVM (RBF) is the most consistent top performer across multiple sources.

| Source | Best Model |
|--------|------------|
| cal_et_al | XGBoost |
| statistical_all4 | SVM (RBF) |
| seldin_128 | SVM (RBF) |
| hsiao_lin_hwa | SVM (RBF) |
| forenseq | SVM (RBF) |
| kidd_55 | Logistic Regression |

### 5. Overfitting Concerns with Larger Panels

> **Finding**: Larger panels (124-125 SNPs) show severe overfitting with gaps of 35-50% between train and test accuracy.

The high-performing panels (cal_et_al, statistical_all4) show reasonable overfit gaps (3-10%), indicating better generalization.

---

## Conclusions and Recommendations

### Primary Conclusions

1. **cal_et_al (52 SNPs)** is the best performing known AISNP panel for East Asian subpopulation classification at **91.50% accuracy**.

2. **Our statistical selection (37 SNPs)** achieves comparable performance at **89.55% accuracy** with fewer markers, demonstrating the effectiveness of our methodology.

3. **General ancestry panels are inadequate** for within-region subpopulation differentiation, achieving only 60-65% accuracy.

4. **Fewer, well-selected SNPs outperform** larger generic panels, emphasizing the importance of population-specific marker selection.

### Recommendations

| Use Case | Recommended Panel | Accuracy | SNPs |
|----------|-------------------|----------|------|
| **Best Accuracy** | cal_et_al | 91.50% | 52 |
| **Minimal Markers** | statistical_all4 | 89.55% | 37 |
| **Production Model** | cal_et_al + XGBoost | 91.50% | 52 |
| **Interpretable Model** | statistical_all4 + Logistic Regression | 86.27% | 37 |

### Future Directions

1. **Combine panels**: Merge cal_et_al and statistical_all4 to identify overlapping markers
2. **Expand populations**: Include CDX (Dai Chinese) and CHS (Southern Han Chinese)
3. **External validation**: Test on independent cohorts outside 1000 Genomes
4. **Marker annotation**: Investigate biological significance of top-performing SNPs

---

## Output Files

### Results Files

| File | Description |
|------|-------------|
| `output/part2/ml_comparison_results.csv` | Full cross-validation results |
| `output/part2/accuracy_by_source.csv` | Accuracy pivot table |
| `output/part2/f1_by_source.csv` | F1 score pivot table |
| `output/part2/ml_comparison_report.txt` | Summary report |

### Visualization Files

| File | Description |
|------|-------------|
| `graphs/part2/snp_counts_by_source.png` | SNP count bar chart |
| `graphs/part2/performance_heatmaps.png` | Accuracy/F1 heatmaps |
| `graphs/part2/best_models_comparison.png` | Best model comparison |
| `graphs/part2/confusion_matrices.png` | Confusion matrices |
| `graphs/part2/snps_vs_accuracy.png` | SNPs vs accuracy scatter |

---

*Report generated from Part 2 pipeline analysis*  
*Known AISNP Panel Comparison*  
*Date: January 2026*
