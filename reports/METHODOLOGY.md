# Methodology — AISNP Panel Selection for Within-East-Asian Population Discrimination

---

## 2.1 Data and Quality Filtering

We used whole-genome SNP data from the 1000 Genomes Project Phase 3, selecting 504 East Asian individuals across five subpopulations: CHB (Han Chinese Beijing), CHS (Han Chinese South), JPT (Japanese), KHV (Kinh Vietnamese), and CDX (Dai Chinese Xishuangbanna). For this study, subpopulations were merged into three target groups: **Han** (CHB + CHS, n=208), **JPT** (JPT, n=104), and **SEA** (KHV + CDX, n=192).

Raw genotype data was processed through a sequential quality filtering pipeline using plink2. The following filters were applied in order:

1. **SNP-only, biallelic**: indels, CNVs, and multi-allelic sites excluded
2. **Minor allele frequency (MAF)**: SNPs with MAF < 1/(2×504) ≈ 0.001 removed (near-monomorphic sites carry no discriminatory power)
3. **Call rate**: SNPs with genotyping completeness < 95% across samples removed
4. **Hardy-Weinberg Equilibrium (HWE)**: SNPs with p < 1×10⁻⁶ removed (mode: keep-fewhet)
5. **LD pruning**: pairwise r² < 0.1 within 1,000 kb windows, step size 1 variant

After filtering, **614,759 SNPs** were retained for downstream analysis.

---

## 2.2 Candidate Set Construction

Directly applying machine learning to 614,759 SNPs is computationally intractable and statistically noisy. We therefore constructed three candidate SNP sets using complementary nomination strategies, each capturing a different aspect of population informativeness.

### 2.2.1 Statistical Block

We evaluated four population-differentiation statistics across all 614,759 SNPs to identify those with significant allele frequency differences among Han, JPT, and SEA:

- **Chi-squared test (χ²)**: tests whether allele count distributions differ significantly across the three groups
- **Jensen-Shannon Divergence (JSD)**: a symmetric, bounded measure of divergence between population allele frequency distributions
- **Allele Frequency Difference (AFD)**: maximum pairwise allele frequency difference across all three population pairs

A preliminary evaluation (`test_evaluation` notebook) compared all candidate statistics for their ability to rank genuinely discriminative SNPs. Based on this evaluation, χ², JSD, and AFD were retained; other tests provided redundant signal. For each retained test, the top-500 highest-ranking SNPs were selected. The union of these three lists formed the **statistical block: 1,005 SNPs**.

### 2.2.2 FST Block

Wright's fixation index (F_ST) measures the proportion of total genetic variance attributable to differences between subpopulations. We computed pairwise Hudson F_ST (Bhatia et al. 2013) between all three population pairs (Han vs. JPT, Han vs. SEA, JPT vs. SEA) using plink2. The top-1,000 SNPs from each pairwise comparison were pooled into a union set, forming the **FST block: 2,508 SNPs**.

### 2.2.3 Combined Block (fst_stat)

The intersection of the statistical block and the FST block — SNPs independently nominated by both approaches — formed the **combined block: 1,003 SNPs**. This consensus set reduces method-specific noise: a SNP present in both blocks has been identified as discriminative by fundamentally different criteria (linear statistical tests vs. population-level differentiation), suggesting more robust signal.

The three candidate sets carried forward to machine learning are summarised below:

| Block | Construction | Size |
|---|---|---|
| stat | Union of top-500 per statistical test (χ², JSD, AFD) | 1,005 |
| FST | Union of top-1,000 per pairwise FST (3 pairs) | 2,508 |
| fst_stat | Intersection of stat and FST blocks | 1,003 |

---

## 2.3 ML-Based Panel Selection

SNP panel selection was framed as a feature selection problem: find the smallest subset of SNPs from the candidate set such that a classifier achieves maximum accuracy on the 3-class Han/JPT/SEA task. We implemented a three-stage pipeline designed to minimise evaluation leakage while systematically comparing all candidate configurations.

### 2.3.1 Stage 1 — Reductor and Candidate Set Selection

For each target panel size N ∈ {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100}, we compared all 9 combinations of 3 candidate sets × 3 reductors using 5-fold nested cross-validation. A **Random Forest** classifier (n_estimators=200) served as the fixed evaluation classifier at this stage for speed.

**Reductors** rank SNPs within the candidate set; the top-N by score are selected:

- **L1 Logistic Regression (L1-LR)**: SAGA solver with L1 penalty; produces sparse coefficient vectors that zero out non-discriminative SNPs. Encourages maximum sparsity.
- **ElasticNet**: SAGA solver with L1+L2 penalty (l1_ratio=0.5); balances sparsity (L1) with coefficient stability (L2). Particularly well-suited to correlated feature sets.
- **Random Forest importance (RF)**: Gini impurity-based mean decrease in impurity across all trees; captures non-linear interactions between SNPs.

To prevent leakage, the reductor was re-fit inside each training fold of the cross-validation — the SNP ranking was never computed using held-out samples. The output of Stage 1 is the **best (candidate set, reductor) pair per N**.

### 2.3.2 Stage 2 — Classifier Evaluation

Using the Stage-1 winning configuration for each N, we evaluated six classifiers via 5-fold nested cross-validation with the same leak-free ranker design:

| Classifier | Key hyperparameters |
|---|---|
| **Random Forest (RF)** | n_estimators=100, max_depth=10 |
| **XGBoost (XGB)** | n_estimators=200, max_depth=6, learning_rate=0.1, subsample=0.8 |
| **Logistic Regression (LR)** | LBFGS solver, max_iter=1000, L2 penalty |
| **SVM with RBF kernel (SVM-RBF)** | Default C and γ (sklearn defaults) |
| **SVM with Linear kernel (SVM-Lin)** | Default C |
| **Gradient Boosting Machine (GBM)** | n_estimators=100, max_depth=5 |

Performance was measured across four metrics per fold: **accuracy**, **weighted F1-score**, **Matthews Correlation Coefficient (MCC)**, and **ROC-AUC** (one-vs-rest, macro-averaged). Stage 2 accuracy is the **primary reported metric** throughout this study.

*Leakage disclosure*: Stage 1 configuration selection and Stage 2 classifier selection both use all 504 samples via cross-validation, introducing mild double-selection optimism. This is unavoidable at n=504 and is consistent with standard practice in small-cohort studies; we disclose it explicitly and address it in the Stage 3 hold-out.

### 2.3.3 Stage 3 — Panel Commitment

For each N, the full dataset was split 80/20 (stratified by population). The Stage-1 winning reductor was fit on the 80% training set to rank SNPs; the top-N were selected as the **committed panel**. The Stage-2 winning classifier was then trained on the 80% subset and evaluated on the 20% hold-out.

Stage 3 accuracy is a single-split confirmation score. It is expected to fall 4–8% below Stage 2 CV accuracy due to (a) single-split variance at small n and (b) the double-selection optimism accumulated across Stages 1 and 2. Stage 3 is not the primary metric; its purpose is to confirm that selected panels generalise beyond the CV folds and to produce the committed SNP lists.

Primary panels committed: **N = 35, N = 50, N = 70**.
Panel files saved as `panels/panel_N{N:03d}.csv` (rank, SNP ID, reductor score).

---

## 2.4 Benchmarking Against Published Panels

### 2.4.1 Classifier Evaluation of Published Panels

Three published AISNP panels were collected for benchmarking: Shi et al. 2019 (four nested panels: 36/59/98/142 SNPs), Cai et al. 2024 (34 EAS-specific SNPs), and Cao et al. 2022 (19 SNPs). Published panel SNPs were mapped to genomic coordinates via the MyVariant.info batch API and extracted from our MAF-filtered pre-LD-pruning genotype data; match rates are reported alongside nominal panel sizes throughout.

Each matched panel was evaluated using the same six-classifier suite and 5-fold stratified cross-validation protocol described in Section 2.3.2, enabling direct metric comparison. Two caveats apply: (a) the PAANDA-EA panels carry in-sample optimism as they were selected from the same 504 individuals, while published panels are external; (b) published panels evaluated on fewer than their nominal SNP count (due to match rate < 100%) may underperform relative to their original reported accuracy.

### 2.4.2 Panel Overlap Analysis

The PAANDA-EA committed panels (N=35, N=50, N=70) were converted from internal SNP IDs to rsIDs via MyVariant.info. Set intersections were computed between all panels (PAANDA-EA and the three published sources). Shared SNP loci were annotated using GeneCards and OMIM to characterise any known trait or disease associations.
