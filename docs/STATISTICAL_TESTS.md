# Statistical Tests Reference

Detailed explanation of statistical methods used for AISNP selection in notebook 02b.

## Table of Contents

1. [Overview](#overview)
2. [Pearson Chi-Squared Test](#1-pearson-chi-squared-test)
3. [Mutual Information](#2-mutual-information)
4. [Information Gain](#3-information-gain)
5. [KL Divergence](#4-kl-divergence)
6. [Multiple Testing Correction](#5-multiple-testing-correction)
7. [Consensus Selection](#6-consensus-selection)

---

## Overview

The statistical SNP selection process uses four complementary tests to identify variants that differentiate between populations:

| Test | Type | Measures | Selection |
|------|------|----------|-----------|
| χ² | Frequentist | Independence | FDR significant |
| Mutual Information | Information-theoretic | Association | Top 500 |
| Information Gain | Information-theoretic | Discrimination | Top 500 |
| KL Divergence | Information-theoretic | Distribution divergence | Top 500 |

**Consensus**: SNPs must pass **ALL 4 tests** to be selected.

---

## 1. Pearson Chi-Squared Test

### Purpose

Tests whether genotype frequencies are independent of population membership.

### Hypothesis

- **H₀**: Genotype distribution is the same across all populations
- **H₁**: Genotype distribution differs between populations

### Method

For each SNP, create a 3×3 contingency table:

|  | CHB | JPT | KHV | Total |
|--|-----|-----|-----|-------|
| Genotype 0 (0/0) | n₀₁ | n₀₂ | n₀₃ | n₀. |
| Genotype 1 (0/1) | n₁₁ | n₁₂ | n₁₃ | n₁. |
| Genotype 2 (1/1) | n₂₁ | n₂₂ | n₂₃ | n₂. |
| Total | n.₁ | n.₂ | n.₃ | n |

Calculate:

```
χ² = Σᵢⱼ (Oᵢⱼ - Eᵢⱼ)² / Eᵢⱼ

where:
  Oᵢⱼ = observed count in cell (i,j)
  Eᵢⱼ = expected count = (row total × column total) / grand total
```

### Interpretation

- Higher χ² → stronger evidence of population-genotype association
- p-value indicates statistical significance
- FDR correction applied for multiple testing

### Code

```python
from scipy.stats import chi2_contingency

contingency = pd.crosstab(population_series, genotype_series)
chi2, pvalue, dof, expected = chi2_contingency(contingency)
```

---

## 2. Mutual Information

### Purpose

Quantifies how much information about population is gained by knowing genotype.

### Formula

```
MI(G; P) = H(P) - H(P|G)
         = Σₚ Σₐ P(g,p) × log₂[P(g,p) / (P(g) × P(p))]
```

Where:

- G = genotype random variable (0, 1, 2)
- P = population random variable (CHB, JPT, KHV)
- H(P) = entropy of population
- H(P|G) = conditional entropy of population given genotype

### Interpretation

- MI = 0 → genotype and population are independent
- Higher MI → stronger association (more informative SNP)
- Measured in bits (base-2 logarithm)

### Properties

- Symmetric: MI(G; P) = MI(P; G)
- Non-negative: MI ≥ 0
- Bounded: MI ≤ min(H(G), H(P))

### Code

```python
from sklearn.metrics import mutual_info_score

mi = mutual_info_score(population_labels, genotype_values)
```

---

## 3. Information Gain

### Purpose

Measures the reduction in population entropy when genotype is known.

### Formula

```
IG(P, G) = H(P) - H(P|G)
         = H(P) - Σₐ P(G=g) × H(P|G=g)
```

Where:

- H(P) = -Σₚ P(p) × log₂ P(p) (population entropy)
- H(P|G=g) = entropy of population among samples with genotype g

### Interpretation

- Quantifies how much knowing the genotype "explains" population membership
- Higher IG → genotype is more predictive of population
- Related to decision tree feature selection

### Difference from MI

- IG is asymmetric (specifically measures P explained by G)
- In practice, IG ≈ MI for this application
- IG is directional: "How much does G tell us about P?"

### Code

```python
def compute_entropy(series):
    probs = series.value_counts(normalize=True)
    return -np.sum(probs * np.log2(probs + 1e-10))

def information_gain(genotype, population):
    h_pop = compute_entropy(population)
    h_pop_given_geno = 0
    for g in genotype.unique():
        mask = genotype == g
        p_g = mask.sum() / len(genotype)
        h_pop_given_geno += p_g * compute_entropy(population[mask])
    return h_pop - h_pop_given_geno
```

---

## 4. KL Divergence

### Purpose

Measures how different genotype distributions are between population pairs.

### Formula

For two populations P and Q with genotype distributions:

```
KL(P||Q) = Σₐ P(g) × log₂[P(g) / Q(g)]
```

For multi-population comparison, compute mean pairwise KL:

```
KL_mean = (1/C) × Σᵢ<ⱼ KL(Pᵢ||Pⱼ)

where C = number of population pairs = n(n-1)/2
```

### Interpretation

- KL = 0 → identical distributions
- Higher KL → more different distributions
- Asymmetric: KL(P||Q) ≠ KL(Q||P) generally

### Properties

- Non-negative: KL ≥ 0
- Not a true metric (asymmetric, no triangle inequality)
- Sensitive to zero probabilities (add small ε)

### Implementation Notes

```python
from scipy.stats import entropy

def kl_divergence(genotype, population):
    pops = population.unique()
    distributions = {}
    
    for p in pops:
        geno_pop = genotype[population == p]
        # Count genotypes 0, 1, 2
        dist = np.array([
            (geno_pop == 0).sum(),
            (geno_pop == 1).sum(),
            (geno_pop == 2).sum()
        ], dtype=float)
        # Normalize and add smoothing
        dist = dist / (dist.sum() + 1e-10) + 1e-10
        distributions[p] = dist
    
    # Mean pairwise KL
    kl_values = []
    for i, p1 in enumerate(pops):
        for j, p2 in enumerate(pops):
            if i < j:
                kl = entropy(distributions[p1], distributions[p2])
                kl_values.append(kl)
    
    return np.mean(kl_values)
```

---

## 5. Multiple Testing Correction

### Problem

When testing thousands of SNPs, many will appear significant by chance.

- 10,000 SNPs × α=0.05 → ~500 false positives expected

### Solutions

#### Bonferroni Correction

```
α_corrected = α / n_tests

Example: 0.05 / 10000 = 5×10⁻⁶
```

- Very conservative
- Controls family-wise error rate (FWER)
- May miss true positives

#### Benjamini-Hochberg FDR (Used in Pipeline)

```
1. Sort p-values: p₁ ≤ p₂ ≤ ... ≤ pₙ
2. Find largest k where pₖ ≤ (k/n) × α
3. Reject all hypotheses 1, 2, ..., k
```

- Less conservative than Bonferroni
- Controls false discovery rate (FDR)
- Better power to detect true positives

### Implementation

```python
from statsmodels.stats.multitest import fdrcorrection

alpha = 0.05
pvalues = results_df['chi2_pvalue'].values
reject, qvalues = fdrcorrection(pvalues, alpha=alpha)

# reject[i] = True means SNP i is significant
# qvalues[i] = adjusted p-value (q-value)
```

---

## 6. Consensus Selection

### Rationale

Each statistical test captures different aspects of population differentiation:

| Test | Captures |
|------|----------|
| χ² | Overall genotype-population dependence |
| MI | Information shared between genotype and population |
| IG | Predictive power of genotype for population |
| KL | Divergence of genotype distributions |

### Selection Criteria

```
Significant SNPs:
  χ²:  FDR-significant (q < 0.05)
  MI:  Top 500 by mutual information
  IG:  Top 500 by information gain
  KL:  Top 500 by mean KL divergence

Consensus = χ²_sig ∩ MI_top500 ∩ IG_top500 ∩ KL_top500
```

### Why All 4?

- Reduces false positives from any single test
- SNPs must be robust across different statistical frameworks
- Ensures both frequentist (χ²) and information-theoretic validation
- Captures SNPs that are informative by multiple criteria

### Typical Results

```
Example from EAS analysis:
  χ² FDR significant: 2,847 SNPs
  MI top 500: 500 SNPs
  IG top 500: 500 SNPs
  KL top 500: 500 SNPs
  
  All 4 tests: ~150-300 SNPs (highly confident AISNPs)
```

---

## Summary Table

| Metric | Formula | Range | Best Value | Interpretation |
|--------|---------|-------|------------|----------------|
| χ² | Σ(O-E)²/E | [0, ∞) | Higher | Stronger association |
| p-value | P(χ² ≥ obs \| H₀) | [0, 1] | Lower | More significant |
| MI | H(P) - H(P\|G) | [0, H(P)] | Higher | More informative |
| IG | H(P) - Σ P(g)H(P\|g) | [0, H(P)] | Higher | More discriminating |
| KL | Σ P log(P/Q) | [0, ∞) | Higher | More divergent |

---

*See also: [PIPELINE.md](PIPELINE.md), [ML_MODELS.md](ML_MODELS.md)*
