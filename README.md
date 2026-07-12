# PAANDA-EA

**P**anel for **A**ISNP-based **A**Ncestry **D**iscovery & **A**ssessment — **E**ast **A**sian

Repository: <https://github.com/Hung-nd233960/PAANDA-EA>

Automated selection and validation of Ancestry-Informative SNPs (AISNPs) for within-East-Asian population discrimination (Han / JPT / SEA) using 1000 Genomes Project Phase 3 data.

```bash
git clone https://github.com/Hung-nd233960/PAANDA-EA.git
```

## Overview

Given 504 EAS individuals across three subpopulations, the pipeline finds the smallest SNP panel that allows a machine-learning classifier to accurately assign ancestry at the Han Chinese / Japanese / Southeast Asian level — a deliberately hard classification task since all three groups belong to the same EAS super-population.

### Target Populations

| Group | Source sub-populations | N |
|-------|----------------------|---|
| **Han** | CHB (Han Chinese Beijing) + CHS (Han Chinese South) | 208 |
| **JPT** | JPT (Japanese in Tokyo) | 104 |
| **SEA** | KHV (Kinh Vietnamese) + CDX (Dai Chinese Xishuangbanna) | 192 |

## Pipeline Architecture

```
1000 Genomes Phase 3 (~22M SNPs, 504 EAS samples)
        │
        ▼
┌──────────────────────────┐
│  1. QUALITY FILTERING    │  614,759 SNPs retained
│  MAF · call rate ·       │
│  HWE · LD pruning        │
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  2. CANDIDATE SET        │  3 pools:
│     CONSTRUCTION         │  stat (1,005) · FST (2,508) · fst_stat (1,003)
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  3. THREE-STAGE ML SWEEP │  Stage 1: pool × reductor selection (9 configs)
│                          │  Stage 2: 6-classifier evaluation
│                          │  Stage 3: 80/20 panel commitment
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  4. BENCHMARKING         │  vs Cai 2024 · Shi 2019 · Cao 2022
└──────────────────────────┘
```

### Key Results

| Panel | N | Leaky CV¹ | Leak-free CV² |
|-------|---|-----------|---------------|
| PAANDA-EA (stat + ElasticNet) | 35 | ~92% | ~86% |
| PAANDA-EA (stat + ElasticNet) | 50 | ~93% | ~92% |
| PAANDA-EA (stat + ElasticNet) | 70 | ~95% | ~92% |
| Cai et al. 2024 (EAS-specific) | 34 | ~94.6% | ~94.6% |
| Cao et al. 2022 | 14 matched / 19 | ~82% | ~82% |
| Shi et al. 2019 | 116 matched / 142 | ~79% | ~79% |

¹ Candidate pool selected on all 504 samples (optimistic — see `reports/METHODOLOGY.md` §2.3.1).
² Nested CV with the whole selection chain rebuilt inside each fold (notebook 08b / `make compare-nested`). Published panels are external, so their numbers are unchanged and directly comparable to the leak-free column. Full leak-free analysis, winner-consistency stats, and recommendations: **`reports/NESTED_CV_MASTER_REPORT.md`**. Note the FST pool overtakes the stat pool for large panels (N≥65): FST+EN+SVM-RBF reaches ~94% at N=70, ~level with Cai-34.

## Project Structure

```
AISNP_Research/
├── notebooks/sea_jpt_cn/
│   ├── pre-filtering/
│   │   ├── 01_hard_filtering.ipynb       # MAF, call rate, biallelic filter
│   │   ├── 02_situational_filtering.ipynb # HWE, LD pruning → 614,759 SNPs
│   │   └── 03_vcf_to_matrix.ipynb        # Build genotype matrix cache
│   │
│   ├── statistical/
│   │   ├── 04a_snp_selection.ipynb       # χ², δAF, JSD ranking
│   │   ├── 05a_stat_training.ipynb       # stat candidate set (1,005 SNPs)
│   │   ├── 05b_reduction.ipynb           # stat/FST correlation analysis
│   │   ├── 05c_fst_stat_training.ipynb   # fst_stat candidate set (1,003 SNPs)
│   │   └── analysis/test_evaluation.ipynb
│   │
│   ├── fst/
│   │   ├── 04b_fst_and_pca.ipynb         # Hudson FST, PCA
│   │   └── 05b_fst_only_training.ipynb   # FST candidate set (2,508 SNPs)
│   │
│   └── self_evaluation/
│       ├── 08_unified_panel_sweep.ipynb  # ★ 3-stage ML sweep (main)
│       ├── 09_published_panel_comparison.ipynb # Benchmark vs published panels
│       ├── 10_panel_overlap.ipynb        # rsID conversion + overlap
│       └── 11_results.ipynb              # All figures and tables
│
├── data/
│   └── published_panels/                 # Cai, Cao, Shi rsID lists + coord maps
│
├── scripts/
│   ├── config.py                         # Centralized paths and parameters
│   └── notebook_init.py                  # Shared notebook setup
│
└── reports/
    ├── SYSTEM_OVERVIEW.md                # Full pipeline description
    └── METHODOLOGY.md                    # Methods section (paper-ready)
```

## Quick Start

### Prerequisites

```bash
# System
Python 3.10+
PLINK2  https://www.cog-genomics.org/plink/2.0/

# Python packages
pip install pandas numpy scipy scikit-learn xgboost matplotlib seaborn \
            tqdm statsmodels requests
```

### Conda Environment

```bash
conda create -n aisnp --file spec-file.txt
conda activate aisnp
```

### Data Preparation

Before running any notebooks, raw VCF files from the 1000 Genomes Project must be
normalized and merged into a single PLINK2 pfile. The script handles this automatically:

```bash
# Usage: scripts/vcf_preprocessing.sh <vcf_dir> <output_prefix> [threads] [memory_mb]
bash scripts/vcf_preprocessing.sh /mnt/data/1000genomes allchr 16 32000
```

What it does:
1. **Normalize** — splits multi-allelic sites, left-aligns indels (bcftools norm), autosomes only
2. **Convert** — each chromosome VCF → PLINK2 pgen/pvar/psam
3. **Merge** — all chromosomes into a single `allchr.{pgen,pvar,psam}` file

The output path must match `PATHS.PLINK_MERGED` in `paths.yaml` / `paths.local.yaml`.
This step only needs to run once; all notebooks assume the merged pfile already exists.

### Run Order

```bash
# 0. Data preparation (once, before everything else)
bash scripts/vcf_preprocessing.sh <vcf_dir> <output_prefix> [threads] [memory_mb]

# 1. Quality filtering
01_hard_filtering → 02_situational_filtering → 03_vcf_to_matrix

# 2. Candidate set construction (can run in parallel)
04a_snp_selection → 05a_stat_training → 05c_fst_stat_training
04b_fst_and_pca   → 05b_fst_only_training

# 3. Main sweep (requires step 2)
08_unified_panel_sweep        # ~30 min, uses cache on re-run

# 4. Benchmarking
09_published_panel_comparison # requires 08
10_panel_overlap              # requires 08

# 5. Results
11_results                    # requires 08 + 09
```

## Methods Summary

### Quality Filtering

| Filter | Threshold |
|--------|-----------|
| SNP-only, biallelic | — |
| MAF | ≥ 1/(2×504) ≈ 0.001 |
| Call rate | ≥ 95% |
| HWE | p ≥ 1×10⁻⁶ (keep-fewhet) |
| LD pruning | r² < 0.10, 1,000 kb window |

### Candidate Sets

| Pool | Construction | Size |
|------|-------------|------|
| **stat** | Union of top-500 per test (χ², δAF, JSD) | 1,005 |
| **FST** | Union of top-1,000 per pairwise Hudson FST | 2,508 |
| **fst_stat** | Intersection of stat and FST | 1,003 |

### Three-Stage ML Sweep (notebook 08)

**Stage 1** — For each N ∈ {5, 10, … 100}: compare 3 pools × 3 reductors (L1-LR, ElasticNet, RF) via 5-fold nested CV. Winner: **stat + ElasticNet** at nearly all N ≥ 15.

**Stage 2** — Evaluate 6 classifiers on Stage 1 winner config:

| Classifier | Key hyperparameters |
|------------|-------------------|
| Random Forest (RF) | n_estimators=100, max_depth=10 |
| XGBoost (XGB) | n_estimators=200, max_depth=6, lr=0.1, subsample=0.8 |
| Logistic Regression (LR) | lbfgs, max_iter=1000 |
| SVM-RBF | RBF kernel, scaled |
| SVM-Linear | Linear kernel, scaled |
| GBM | n_estimators=100, max_depth=5 |

**Stage 3** — 80/20 stratified hold-out to commit final panels. Committed panels: **N = 35, N = 50, N = 70**.

### Benchmarking (notebook 09)

Published panels evaluated on the same 504-sample cohort with identical 6-classifier + 5-fold CV protocol. Match rates: Cai 34/34 (100%), Cao 14/19 (74%), Shi 116/142 (82%).

## References

- 1000 Genomes Project: <https://www.internationalgenome.org/>
- PLINK2: <https://www.cog-genomics.org/plink/2.0/>
- Cai et al. 2024 — EAS-specific 34-SNP panel
- Cao et al. 2022 — 19-SNP panel
- Shi et al. 2019 — 36/59/98/142-SNP nested panels

## License

For research and educational purposes only.
