---
marp: true
theme: default
paginate: true
style: |
  section {
    font-family: 'Segoe UI', sans-serif;
    font-size: 22px;
  }
  h1 { color: #1a3a5c; font-size: 1.8em; }
  h2 { color: #1a3a5c; border-bottom: 2px solid #3a7bd5; padding-bottom: 6px; }
  table { font-size: 0.85em; width: 100%; }
  th { background: #1a3a5c; color: white; }
  code { background: #f0f4f8; padding: 2px 6px; border-radius: 3px; }
  .highlight { color: #d63031; font-weight: bold; }
  blockquote { border-left: 4px solid #3a7bd5; background: #f0f6ff; padding: 10px 16px; }
---

# Ancestry-Informative SNP Selection
## A Statistical Approach Beyond FST

**CN · JPT · SEA — East Asian Population Classification**

---

## Why Does Ancestry Matter?

- **Forensics** — ancestry inference from biological samples
- **Medicine** — drug response and disease risk vary by population
- **Population genetics** — tracing human migration and history

> The goal: find the **smallest set of DNA positions** (SNPs) that reliably identifies where a person's ancestors came from.

---

## The Data

**1000 Genomes Project** — 504 individuals, 3 East Asian groups

| Group | Composition | N |
|---|---|---|
| **CN** | Han Chinese (Beijing + South) | 208 |
| **JPT** | Japanese | 104 |
| **SEA** | Kinh (Vietnam) + Chinese Dai | 192 |

> These groups are **genetically close** — all East Asian.  
> This is a deliberately hard problem.

---

## The Scale Challenge

```
Whole genome sequencing
        ↓
~78,000,000 raw SNPs
        ↓  quality filters
  614,759 candidate SNPs
        ↓  ???
      ~50–150 SNPs  ← what we want
```

**Machine learning cannot handle 600k features directly.**  
We need a principled way to narrow the search first.

---

## Our Pipeline — 4 Stages

```
┌─────────────────────────────────────┐
│  STAGE 1   Pre-filtering            │  78M → 614,759 SNPs
│            Quality control (QC)     │
└──────────────────┬──────────────────┘
                   │
        ┌──────────┴──────────┐
        ▼                     ▼
┌──────────────┐     ┌─────────────────┐
│  STAGE 2A    │     │  STAGE 2B       │
│  FST-based   │     │  Statistical    │  ← NEW
│  (classical) │     │  (this work)    │
└──────┬───────┘     └───────┬─────────┘
       └──────────┬──────────┘
                  ▼
┌─────────────────────────────────────┐
│  STAGE 3   ML Training              │  7 classifiers, 5-fold CV
└──────────────────┬──────────────────┘
                   ▼
┌─────────────────────────────────────┐
│  STAGE 4   ML-based Reduction       │  → 50–150 SNP final panel
└─────────────────────────────────────┘
```

---

## Stage 1 — Pre-filtering (QC)

Standard genomics quality control. **No changes from classical practice.**

| Filter | Setting |
|---|---|
| SNP-only (no indels) | ✓ |
| Biallelic variants only | max 2 alleles |
| Min allele frequency | ≥ 0.1% |
| Genotyping call rate | ≥ 95% |
| Hardy-Weinberg Equilibrium | p > 1×10⁻⁶ |
| LD pruning | R² < 0.10, 1000 kb window |

**78,000,000 → 614,759 independent SNPs**

---

## Stage 2A — FST (Baseline)

**Fixation Index (FST)** — the classical population genetics approach

```
FST = (Total variance − Within-population variance) / Total variance
```

- Computed **pairwise** (CN↔JPT, CN↔SEA, JPT↔SEA) then aggregated
- Ranks SNPs by how "differentiated" they are across populations
- FST = 0 → same frequencies everywhere
- FST = 1 → completely diverged

> **Take the top 1,000 SNPs by FST → feed to ML**

Well-established, interpretable, widely used in population genetics.

---

## Stage 2B — Statistical Filtering (New)

**Question:** Can statistical tests find equally good panels — using a *different lens* than FST?

We evaluate **7 metrics** across all 614,759 SNPs simultaneously:

| # | Metric | Captures |
|---|---|---|
| 1 | **χ² (Chi-squared)** | Genotype-population independence |
| 2 | Cramér's V | Effect size (derived from χ²) |
| 3 | ANOVA F | Mean dosage difference |
| 4 | Kruskal-Wallis H | Non-parametric mean difference |
| 5 | **δAF (Max AF Delta)** | Largest raw frequency gap |
| 6 | **JSD (Jensen-Shannon Div.)** | Full distribution divergence |
| 7 | Mutual Information | Information overlap |

---

## Finding the Independent Signals

Most metrics are **highly correlated** — they say the same thing.

| Pair | r | Verdict |
|---|---|---|
| ANOVA F ↔ KW H | **1.00** | Identical — drop one |
| JSD ↔ MI | **0.99** | Identical — drop one |
| χ² ↔ δAF | **0.68** | Most independent pair |

**Three truly independent signals:**

> **χ²** (statistical test)  ·  **δAF** (frequency gap)  ·  **JSD** (distribution shape)

Union of top-500 from each → **~1,005 candidate SNPs** for ML

---

## Note: δAF ≠ FST

Both measure allele frequency divergence — but differently:

| | δAF | FST |
|---|---|---|
| Formula | max \|AF_i − AF_j\| | Var(AF) / (p̄ × (1−p̄)) |
| Normalization | None — absolute | By heterozygosity |
| Rare alleles | Same weight | **Upweighted** |

A rare allele and a common allele with the **same 20% frequency gap** get the same δAF but very different FST.

**Top-500 δAF overlaps top-500 χ² at only 14%** — it finds fundamentally different SNPs.

---

## Stage 4 — ML-Based Reduction

Starting from ~1,005 SNPs, rank by combined RF + LR feature importance.  
Evaluate accuracy at shrinking panel sizes:

| Panel size | Accuracy | Drop from full |
|---|---|---|
| 25 SNPs | 0.867 | −12.3% |
| 50 SNPs | **0.944** | −4.6% |
| 100 SNPs | 0.978 | −1.2% ← **elbow** |
| 150 SNPs | **0.990** | ~0% |
| 200+ SNPs | ~0.990 | saturated |

> **50 SNPs → 94% accuracy.**  
> **150 SNPs → near-perfect (~99%).**

---

## Results — Fair Comparison at Top-50 SNPs

All approaches evaluated at exactly **50 SNPs**:

| Approach | Best Accuracy | Best Model |
|---|---|---|
| **FST-only** | **0.948** | LR / SVM-RBF |
| **Statistical V2** | **0.944** | SVM (RBF) |
| FST + Statistical | 0.921 | Random Forest |
| Statistical V1 | 0.885 | Logistic Regression |

**Statistical V2 matches FST at top-50.**  
At larger panels (150+ SNPs), statistical surpasses FST.

---

## What Do Both Methods Agree On?

**34 SNPs** appear in the top-50 of **every approach** tested.

```
FST-only  ∩  V1 FST+RF  ∩  V2 FST+RF  =  34 SNPs
```

These are the **highest-confidence ancestry-informative markers** — independently confirmed by both classical and statistical methods from completely different mathematical frameworks.

> These 34 are the best candidates for a minimal validated panel.

---

## Key Takeaways

1. **Statistical filtering works.** χ², δAF, and JSD together rival FST at top-50 and surpass it at larger panels.

2. **Most metrics are redundant.** 7 tests collapse to 3 independent signals — good news for computational efficiency.

3. **50 SNPs is a practical floor** — ~94% accuracy. 150 SNPs gives near-perfect (~99%) performance.

4. **The two approaches are complementary**, not substitutes. Their SNP sets overlap at only ~34 high-confidence markers — they explore different genomic signal.

5. **Methodology over results** — the contribution is a fast, principled statistical framework for AISNP discovery that does not rely on FST.

---

## Future Work

| Next step | Description |
|---|---|
| Benchmark vs published panels | Compare against known AISNP sets from literature and commercial products |
| Validate 34-SNP core | Test on held-out or external populations |
| RF vs LR reduction | Systematically compare feature importance methods |
| Global populations | Extend beyond East Asian |
| Array compatibility | Check if top SNPs are present on common genotyping arrays |
