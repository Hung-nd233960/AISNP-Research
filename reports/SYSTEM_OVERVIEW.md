# AISNP Selection Pipeline — System Overview

## Problem Statement

Given whole-genome sequencing data from 504 individuals across three East Asian population groups, identify the **smallest possible set of Single Nucleotide Polymorphisms (SNPs)** that allows a machine learning classifier to accurately assign any individual to their ancestral population.

### Why it matters
- Forensic identification and missing-persons cases
- Pharmacogenomics (drug response varies by ancestry)
- Population genetics research
- Medical risk stratification

### The scale challenge
After standard quality filtering, **614,759 SNPs** remain. Machine learning classifiers cannot operate at this scale directly — training a Random Forest on 600k features per sample is computationally intractable and statistically noisy. A principled dimensionality reduction step is required before ML.

---

## Data

| Property | Value |
|---|---|
| Source | 1000 Genomes Project Phase 3 |
| Total samples used | 504 |
| Super-population | EAS (East Asian) |
| SNPs after pre-filtering | 614,759 |

### Population Design — 5 raw → 3 target groups

The five 1000 Genomes East Asian sub-populations are merged into three classification targets:

| Target Group | Raw Sub-populations | Description |
|---|---|---|
| **CN** | CHB + CHS | Han Chinese (Beijing + South China) |
| **JPT** | JPT | Japanese in Tokyo |
| **SEA** | KHV + CDX | Kinh (Vietnam) + Chinese Dai (Yunnan) |

| Group | N |
|---|---|
| CN | 208 |
| SEA | 192 |
| JPT | 104 |
| **Total** | **504** |

This is a deliberately hard classification problem: all three groups belong to the same super-population (EAS) and are genetically close. Simple continental-level ancestry markers are insufficient.

---

## System Architecture — Four Stages

```
1000 Genomes VCF (all chromosomes, 2,504 individuals)
        │
        ▼
┌───────────────────────────────────┐
│  STAGE 1 — PRE-FILTERING          │  614,759 SNPs remain
│  Quality control                  │
└───────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────┐
│  STAGE 2 — SNP SCORING            │
│                                   │
│  ┌──────────────┐  ┌───────────┐  │
│  │ FST-based    │  │Statistical│  │  Two competing approaches
│  │ (classical)  │  │(new)      │  │  → both reduce to ~500-1000 SNPs
│  └──────────────┘  └───────────┘  │
└───────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────┐
│  STAGE 3 — ML TRAINING            │  7 classifiers, 5-fold CV
│  Train on candidate panels        │
└───────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────┐
│  STAGE 4 — ML-BASED REDUCTION     │  ~50–150 SNPs remain
│  Panel size elbow curve           │
└───────────────────────────────────┘
```

---

## Stage 1 — Pre-filtering

Standard genomic quality control applied uniformly to all analyses. No methodological novelty here — this stage is fixed across all approaches.

### Hard Filters (always applied)

| Filter | Threshold | Rationale |
|---|---|---|
| SNP-only | Exclude indels, CNVs | Consistency across arrays and sequencing |
| Biallelic | Max 2 alleles | Required for standard population genetics statistics |
| MAF | ≥ 1/(504×2) ≈ 0.001 | Remove singletons (one allele across all samples) |
| Call rate | ≥ 95% | Remove poorly genotyped variants |

Input: ~78.4M raw variants → after hard filters: retained for next step.

### Situational Filters (applied before scoring)

| Filter | Setting | Rationale |
|---|---|---|
| HWE test | p < 1×10⁻⁶, keep-fewhet | Remove variants likely due to genotyping error |
| LD pruning | 1000 kb window, R² < 0.10 | Remove correlated variants; retain independent signals |

After LD pruning: **614,759 independent SNPs** pass to Stage 2.

### Tool
PLINK2 (v2.0), 16 threads, 16 GB RAM.

---

## Stage 2A — FST-Based Filtering (Baseline Approach)

### What is FST?

FST (Wright's Fixation Index) measures the proportion of genetic variance that falls **between** populations relative to total genetic variance:

```
FST = (H_T − H_S) / H_T
```

Where H_T = total heterozygosity, H_S = average within-population heterozygosity.

- FST = 0 → allele frequencies identical across populations
- FST = 1 → allele frequencies completely diverged (fixed in each population)

FST is computed **pairwise** across all three population pairs (CN↔JPT, CN↔SEA, JPT↔SEA) and aggregated. This is the classical approach in population genetics for identifying ancestry-informative markers.

### Implementation
- PLINK2 `--fst` with Hudson estimator
- Top **1,000 SNPs** by FST score passed to downstream ML
- This is the baseline against which the statistical approach is compared

### Key property
FST normalizes by heterozygosity — a rare allele with a 20% frequency difference gets a higher FST than a common allele with the same absolute difference. This makes FST sensitive to private alleles and drift-fixed variants.

---

## Stage 2B — Statistical Filtering (New Approach)

### Core motivation

Before ML can run, we need to go from 614,759 SNPs to a manageable panel. The question is: **can statistical tests do this independently of FST**, identifying a different (and complementary) class of ancestry-informative SNPs?

### The approach

Seven statistical tests are computed in parallel across all 614,759 SNPs using vectorized NumPy operations (no Python loop over SNPs — the entire genome is processed in seconds via a shared `(n_populations × 3_genotypes × n_SNPs)` tensor).

| Test | What it measures | Type |
|---|---|---|
| **χ² (Chi-squared)** | Independence of genotype from population label | Categorical |
| **Cramér's V** | Effect size derived from χ² | Effect size |
| **ANOVA F** | Difference in mean dosage (0/1/2) across populations | Parametric |
| **Kruskal-Wallis H** | Non-parametric analog of ANOVA | Non-parametric |
| **δAF (Max AF Delta)** | Max pairwise allele frequency difference | Direct |
| **JSD (Jensen-Shannon Divergence)** | Distributional divergence across populations | Information theory |
| **MI (Mutual Information)** | Information shared between genotype and population | Information theory |

### Key finding: most tests are redundant

Metric correlation analysis (Pearson R on 614,759 SNPs):

| Pair | Correlation | Conclusion |
|---|---|---|
| ANOVA F ↔ KW H | **1.00** | Completely redundant — keep one |
| JSD ↔ MI | **0.99** | Essentially identical — keep one |
| χ² ↔ ANOVA F | 0.98 | χ² subsumes ANOVA F signal |
| χ² ↔ JSD | 0.96 | Moderate redundancy |
| χ² ↔ δAF | **0.68** | Most independent pair |
| δAF ↔ JSD | 0.71 | Most independent after χ² |

### Chosen metrics for panel construction: χ², δAF, JSD

These three represent three distinct mathematical perspectives:
- **χ²** — categorical independence test (most statistically established)
- **δAF** — raw allele frequency divergence (most orthogonal to all others)
- **JSD** — information-theoretic divergence of genotype distributions

For each metric, the top-N SNPs by score are extracted. The **union** of all three lists is passed to ML. This was tested at N ∈ {100, 200, 500}; the union of top-500 from each (1,005 unique SNPs) was the best-performing input.

### δAF vs FST — are they the same?

No. Both measure allele frequency divergence but differ in normalization:
- **δAF** = max |AF_i − AF_j| — absolute, unweighted
- **FST** = Var(AF) / (p̄ × (1−p̄)) — normalized by heterozygosity

A rare allele (p̄ = 0.1) with δAF = 0.20 has FST ≈ 0.15. A common allele (p̄ = 0.5) with the same δAF = 0.20 has FST ≈ 0.04. **FST upweights rare variants; δAF treats all allele frequencies equally.** In practice, top-500 δAF overlaps with top-500 χ² at only 14% — confirming δAF captures a distinct signal.

---

## Stage 3 — ML Training

Seven classifiers are evaluated per approach using 5-fold stratified cross-validation:

| Classifier | Notes |
|---|---|
| Random Forest | 200 trees, no max depth |
| XGBoost | 200 trees, depth 6, lr 0.1 |
| Logistic Regression | Multinomial, lbfgs, max_iter 1000 |
| SVM (RBF kernel) | |
| SVM (Linear kernel) | |
| K-Nearest Neighbors | |
| Gradient Boosting | |
| MLP Neural Network | |
| Naive Bayes | |
| AdaBoost | |

Split: 80% train / 20% test, random state fixed at 42.

---

## Stage 4 — ML-Based Panel Reduction

### The core insight

Even after Stage 2 narrows the pool to ~1,000 SNPs, many are partially redundant. ML feature importance (Random Forest + Logistic Regression combined ranking) can identify the **minimal effective panel**.

### Method

Starting from the best-performing condition's SNP list (union top-500 = 1,005 SNPs), SNPs are ranked by average RF + LR importance rank. All 10 classifiers are then re-evaluated at progressively smaller panel sizes:

`[25, 30, 40, 50, 55, 60, 100, 150, 200, 300, 500, 750, 1005]`

### Elbow curve results (F3 V2 Statistical approach)

| Panel size | Best accuracy | Best classifier | Drop from baseline |
|---|---|---|---|
| 25 SNPs | 0.867 | Random Forest | −12.3% |
| 50 SNPs | **0.944** | SVM (RBF) | −4.6% |
| 100 SNPs | 0.978 | SVM (RBF) | −1.2% ← **elbow** |
| 150 SNPs | 0.990 | SVM (RBF) | ~0% |
| 200+ SNPs | ~0.990 | — | saturated |

**The elbow point is ~100 SNPs**: accuracy within 1.2% of the 1,005-SNP baseline. A 50-SNP panel retains 95.4% of the full-set performance.

---

## Performance Comparison — All Approaches at Top-50 SNPs

Fair comparison: all approaches evaluated at exactly 50 SNPs.

| Approach | Best Accuracy | Top Model |
|---|---|---|
| F1: FST-only | **0.9484** | LR / SVM-RBF (tied) |
| F3: V2 Statistical | **0.9444** | SVM (RBF) |
| F2/F3: FST + Statistical | 0.9206 | Random Forest |
| F2: V1 Statistical | 0.8849 | Logistic Regression |

### Key observations

1. **FST-only and Statistical V2 are essentially tied at top-50** (0.9484 vs 0.9444). Neither is definitively better at this scale.
2. **Statistical V2 surpasses FST at larger panels** — at 150 SNPs, statistical achieves 0.990 vs FST's plateau.
3. **Statistical V1 underperforms** (0.8849) — the V2 improvement (new metrics + union strategy) was substantial.
4. **FST + Statistical combination does not improve over either alone** at top-50. The union approach does not benefit from redundancy.

---

## The 34-SNP Core

34 SNPs appear in the top-50 of **every single approach** (FST-only, V1 FST+RF, V2 FST+RF). These are the highest-confidence ancestry-informative markers — independently confirmed by three methodologically distinct pipelines.

These 34 SNPs represent the practical minimum panel for CN/JPT/SEA classification. They are likely to be the most biologically robust markers and the best candidates for validation against published panels.

### SNP set membership summary (across all 5 sub-approaches)

| Category | Count |
|---|---|
| Unique to one approach only | 106 |
| Shared by exactly 2 approaches | 21 |
| Shared by all 3 factions | **34** |
| Total unique SNPs | 161 |

---

## Metric Correlation Summary

From the 7-test evaluation on 614,759 SNPs, the effective dimensionality is **3 independent signals**:

```
Cluster A (r ≥ 0.89):  χ²  ≈  Cramér's V  ≈  ANOVA F  ≈  KW H  ≈  JSD  ≈  MI
Cluster B (r ≤ 0.71):  δAF  (most orthogonal to everything)
```

The 7 tests reduce to 3 truly independent measures: **χ², δAF, JSD**.

---

## Limitations and Future Work

| Item | Status |
|---|---|
| Comparison against published AISNP panels | Not yet run — data collection in progress |
| RF vs LR for feature importance in reduction | Not yet systematically verified |
| Cross-validation of the 34-SNP core on held-out populations | Not done |
| Extension to broader EAS or global populations | Not explored |
| Statistical V2 on larger panel sizes with XGBoost sweep | Partial |

---

## Summary

| Stage | Input | Output | Key tool |
|---|---|---|---|
| Pre-filter | ~78M raw SNPs | 614,759 clean SNPs | PLINK2 |
| FST scoring | 614,759 SNPs | Top 1,000 by FST | PLINK2 `--fst` |
| Statistical scoring | 614,759 SNPs | Union top-500×3 = ~1,005 SNPs | NumPy vectorized |
| ML training | ~1,000–2,500 SNPs | 7-classifier CV results | scikit-learn, XGBoost |
| ML reduction | ~1,005 SNPs | **50–150 SNP final panel** | RF + LR importance |

**Main methodological claim**: Statistical filtering (χ², δAF, JSD union) is a valid, fast, and complementary alternative to FST for identifying ancestry-informative SNP panels. It converges to comparable accuracy (~94% at 50 SNPs, ~99% at 150 SNPs) while capturing a partly orthogonal subset of the genome.
