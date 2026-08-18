# PAANDA-EA

**P**anel for **A**ISNP-based **A**Ncestry **D**iscovery & **A**ssessment — **E**ast **A**sian

Repository: <https://github.com/Hung-nd233960/PAANDA-EA>

Automated, leak-free selection and validation of Ancestry-Informative SNPs (AISNPs) for within-East-Asian population discrimination (Han / JPT / SEA) using 1000 Genomes Project Phase 3 data.

```bash
git clone https://github.com/Hung-nd233960/PAANDA-EA.git
```

**Full write-up:** [`reports/REPORT.md`](reports/REPORT.md) (complete methodology, derivations, statistical tests, marker documentation) · [`reports/REPORT_SHORT.md`](reports/REPORT_SHORT.md) (condensed version, includes a point-by-point reviewer-response section).

## Overview

Given 504 EAS individuals across three subpopulations, the pipeline finds the smallest SNP panel that allows a machine-learning classifier to accurately assign ancestry at the Han Chinese / Japanese / Southeast Asian level — a deliberately hard classification task since all three groups belong to the same EAS super-population.

Because candidate-pool construction, feature ranking, and classification all use the population label, every one of those steps is rebuilt independently inside each cross-validation fold — nothing that touches a label is fit on a sample before that sample's own held-out score is computed. Two accuracy estimates are reported throughout: **Blind Selection**, in which an inner CV loop picks the (pool, reducer, classifier) configuration per fold without ever seeing that fold's test data (zero selection optimism), and **Best Fixed Configuration**, the single best-scoring, named, deployable panel across the full grid (carries a small, disclosed "winner's-curse" optimism). See `reports/REPORT.md` §IV for the full nested-CV design.

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
│  1. QUALITY FILTERING    │  614,759 SNPs retained (label-agnostic; safe to
│  MAF · call rate ·       │  compute once, outside cross-validation)
│  HWE · LD pruning        │
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  2. CANDIDATE POOLS      │  3 pools, REBUILT PER FOLD from training samples only:
│     (label-using)        │  stat (1,005) · FST (2,508) · fst_stat (965)
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  3. LEAK-FREE NESTED CV  │  Outer: 5-fold. Inner: 5-fold, picks (pool, reducer,
│                          │  classifier) blind to the outer-test fold.
│                          │  → Blind Selection (bias-free) +
│                          │    Best Fixed Configuration (named, deployable)
└──────────────────────────┘
        │
        ▼
┌──────────────────────────┐
│  4. BENCHMARKING         │  vs Cai 2024 · Shi 2019 · Cao 2022, re-genotyped
│                          │  on this cohort, same protocol; paired significance
│                          │  testing against Cai-34 at matched size and N=70
└──────────────────────────┘
```

`fst_stat` is **not** a set intersection — it applies the `stat` pool's own ranking (χ²/JSD/AFD) to `FST`'s 2,508 candidates only, i.e. a cascade of one criterion's ranking restricted to another criterion's scope. An earlier version of this pipeline built `fst_stat` as a full-cohort intersection *before* the CV split, which leaked label information into candidate eligibility; that bug is fixed in the current architecture (`scripts/nested_selection.py:fst_stat_pool_indices`) and is why `fst_stat`'s size changed from 1,003 (old, leaky) to 965 (current, leak-free).

### Key Results

| Panel | N | Blind Selection | Best Fixed Configuration |
|-------|---|---|---|
| PAANDA-EA (fst_stat + EN + SVM-RBF) | 35 | 86.91% [85.34, 88.48] | 87.90% [86.08, 89.72] |
| PAANDA-EA (stat + EN + SVM-RBF) | 50 | 88.28% [84.37, 92.20] | 91.67% [90.00, 93.33] |
| PAANDA-EA (FST + EN + SVM-RBF) | 70 | 92.86% [89.68, 96.04] | 94.25% [92.25, 96.26] |
| Cai et al. 2024 (EAS-specific) | 34 | — | 94.64% |
| Cao et al. 2022 | 14 matched / 19 | — | 81.94% |
| Shi et al. 2019 | 116 matched / 142 | — | 78.57% |

Both columns are unbiased estimates from the same 5 outer folds (95% CI in brackets, t-distribution, 4 df) — Blind Selection carries no selection optimism at all; Best Fixed Configuration is the single named panel, with its residual optimism explicitly measured rather than assumed away. ElasticNet is the dominant reducer at every evaluated panel size (18/18); candidate-pool choice is the least trustworthy part of a hindsight read of the full grid (§IV-D3, §VI-A of the full report).

**"Parity" with Cai-34 is conditional on panel size, not unconditional.** A paired *t*-test on the 5 matched outer folds finds no significant difference from Cai-34 at N=70 (mean diff −0.39 pts, 95% CI [−3.72, +2.94], *p*=0.76), but a significant, sizeable deficit at the more size-matched N=35 (mean diff −6.74 pts, 95% CI [−9.10, −4.39], *p*=0.0014). Read together: this pipeline reaches parity with expert curation only once given roughly double Cai-34's marker count; at matched size, expert curation still wins clearly. Full statistical detail: `reports/REPORT.md` §V-D, `scripts/ci_and_paired_test.py`.

## Project Structure

```
AISNP_Research/
├── notebooks/sea_jpt_cn/
│   ├── pre-filtering/
│   │   ├── 01_hard_filtering.ipynb        # MAF, call rate, biallelic filter
│   │   ├── 02_situational_filtering.ipynb # HWE, LD pruning → 614,759 SNPs
│   │   └── 03_vcf_to_matrix.ipynb         # Build genotype matrix cache
│   │
│   ├── statistical/
│   │   ├── 04a_snp_selection.ipynb        # χ², δAF, JSD ranking
│   │   ├── 05a_stat_training.ipynb        # stat candidate set (1,005 SNPs)
│   │   ├── 05b_reduction.ipynb            # stat/FST correlation analysis
│   │   ├── 05c_fst_stat_training.ipynb    # fst_stat candidate set (full-cohort; superseded by the per-fold version)
│   │   └── analysis/test_evaluation.ipynb
│   │
│   ├── fst/
│   │   ├── 04b_fst_and_pca.ipynb          # Hudson FST, PCA
│   │   └── 05b_fst_only_training.ipynb    # FST candidate set (2,508 SNPs)
│   │
│   └── self_evaluation/
│       ├── 08_unified_panel_sweep.ipynb        # Frozen leaky baseline (candidate pools built on all 504 samples) — kept only as the "how much did leakage inflate accuracy" reference point, not the headline pipeline
│       ├── 09_published_panel_comparison.ipynb # Benchmark vs published panels
│       ├── 10_panel_overlap.ipynb              # rsID conversion + overlap
│       └── 11_results.ipynb                    # Figures/tables for the leaky-baseline notebook chain
│
├── data/
│   └── published_panels/                  # Cai, Cao, Shi rsID lists + coord maps
│
├── scripts/
│   ├── config.py, notebook_init.py        # Centralized paths, shared setup
│   ├── nested_selection.py, nested_fst.py # Per-fold candidate-pool construction (leak-free)
│   ├── nested_cv_sweep.py                 # ★ Leak-free nested CV sweep — the current headline pipeline
│   │                                        (Tier B: isolates the leakage; Tier A: Blind Selection)
│   ├── nested_comparison.py               # Leak-free vs. frozen-leaky comparison
│   ├── ci_and_paired_test.py              # 95% CIs + paired t-test/Wilcoxon vs. Cai-34
│   ├── committed_panel_marker_table.py    # Per-marker rsID/AF/FST/fold-selection-frequency table
│   ├── subpopulation_breakdown.py         # Recall broken out to the 5 original 1000G subpopulations
│   ├── fetch_rsids.py                     # rsID lookup via Ensembl GRCh37 VEP API
│   └── report_figures.py, marker_stability_figure.py, system_diagram_*.py, build_slides.py, build_poster.py
│                                          # Figure/deliverable generation for reports/
│
└── reports/
    ├── REPORT.md                          # ★ Canonical write-up: full methodology, results, statistics, marker documentation
    ├── REPORT_SHORT.md                    # ★ Condensed write-up + point-by-point reviewer-response section
    ├── figures/                           # All report figures (PNG/SVG) + committed_panel_markers.csv
    ├── METHODOLOGY.md, SYSTEM_OVERVIEW.md, NESTED_CV_MASTER_REPORT.md
    │                                       # Earlier-phase working documents; superseded by REPORT.md
    └── sea_jpt_cn/                         # Per-notebook rendered reports
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

# 1. Quality filtering (label-agnostic, run once)
make prefilter          # 01_hard_filtering → 02_situational_filtering → 03_vcf_to_matrix

# 2. Leak-free nested CV — the current headline pipeline
make nested              # scripts/nested_cv_sweep.py --tier both (Tier A: Blind Selection; Tier B: leakage isolation)
make compare-nested       # vs. the frozen leaky baseline below

# 3. Frozen leaky baseline (optional, for the leakage-magnitude comparison only)
make sweep benchmark overlap results   # notebooks 08 → 09 → 10 → 11

# 4. Retroactive statistics + supplementary tables (no retraining; reuses archived per-fold CSVs)
conda run -n aisnp python scripts/ci_and_paired_test.py
conda run -n aisnp python scripts/committed_panel_marker_table.py
conda run -n aisnp python scripts/fetch_rsids.py
conda run -n aisnp python scripts/subpopulation_breakdown.py
```

`make help` lists every target; `make clean-nested` / `make clean-sweep` / `make clean` clear cached intermediates for a given stage.

## Methods Summary

### Quality Filtering

| Filter | Threshold |
|--------|-----------|
| SNP-only, biallelic | — |
| MAF | ≥ 1/(2×504) ≈ 0.001 |
| Call rate | ≥ 95% |
| HWE | p ≥ 1×10⁻⁶ (keep-fewhet) |
| LD pruning | r² < 0.10, 1,000 kb window |

### Candidate Pools (rebuilt inside every CV fold)

| Pool | Construction | Size |
|------|-------------|------|
| **stat** | Union of top-500 SNPs per test (χ², JSD, AFD) | 1,005 |
| **FST** | Union of top-1,000 SNPs per pairwise Hudson $F_{ST}$ comparison | 2,508 |
| **fst_stat** | stat's own ranking, cascaded onto FST's 2,508 candidates only | 965 |

### Leak-Free Nested CV (`scripts/nested_cv_sweep.py`)

- **Outer loop:** 5-fold `StratifiedKFold` (seed 42) over all 504 samples. Pools, reducer, and classifier are rebuilt from each outer-training fold only.
- **Inner loop (Blind Selection):** a further 5-fold split of each outer-training fold, used in two scored stages — Stage 1 picks (pool, reducer) via a fixed Random Forest probe; Stage 2, given that pair, picks the classifier from all 6 candidates. Both stages are scored by mean accuracy across the 5 inner folds.
- **Final estimate:** the winning configuration is refit once on the full outer-training fold and scored once on the untouched outer-test fold. The headline Blind Selection accuracy at each N is the unweighted mean of the 5 outer-fold scores.
- **Best Fixed Configuration:** at a given N, the single (pool, reducer, classifier) combination scoring highest across the full grid (18 N × 3 pools × 3 reducers × 6 classifiers × 5 folds, 4,860 fits) — named as the deployable panel, with its winner's-curse gap over Blind Selection explicitly reported.

| Classifier | Key hyperparameters |
|------------|-------------------|
| Random Forest (RF) | n_estimators=100, max_depth=10 |
| XGBoost (XGB) | n_estimators=200, max_depth=6, lr=0.1, subsample=0.8 |
| Logistic Regression (LR) | lbfgs, max_iter=1000 |
| SVM-RBF | RBF kernel, scaled |
| SVM-Linear | Linear kernel, scaled |
| GBM | n_estimators=100, max_depth=5 |

Committed panels: **N = 35, 50, 70**. Full derivation, precision/winner's-curse analysis, and component-level trust hierarchy (reducer > classifier > pool): `reports/REPORT.md` §IV–§VI.

### Benchmarking

Published panels re-genotyped and evaluated on the same 504-sample cohort with the identical 6-classifier, 5-fold CV protocol. Match rates: Cai 34/34 (100%), Cao 14/19 (74%), Shi 116/142 (82%). Paired significance testing against Cai-34 (5 matched outer folds): not distinguishable at N=70 (*p*=0.76); significant deficit at the size-matched N=35 (*p*=0.0014) — see `reports/REPORT.md` §V-D and `reports/REPORT_SHORT.md` §8 for the full statistical treatment and caveats (only 5 paired folds; both tests are underpowered by conventional standards).

## References

- 1000 Genomes Project: <https://www.internationalgenome.org/>
- PLINK2: <https://www.cog-genomics.org/plink/2.0/>
- Cai et al. 2024 — EAS-specific 34-SNP panel
- Cao et al. 2022 — 19-SNP panel
- Shi et al. 2019 — 36/59/98/142-SNP nested panels

Full citation list with DOIs: `reports/references.bib` (if present) or the References section of `reports/REPORT.md`.

## License

For research and educational purposes only.
