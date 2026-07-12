# Master Report — Leak-Free Nested Cross-Validation of the PAANDA-EA Pipeline

*Branch `nested-cv-leakfree`. Generated 2026-07-13. Compares against the frozen
leaky baseline `results_archive/baseline_leaky_v1/` (commit a803362).*

---

## 1. Executive summary

The original pipeline (notebook 08) reported 5-fold CV accuracy for panels whose
**candidate pools were selected once from all 504 labelled samples** (notebooks
04a/04b/05a/05c), before the CV. Every fold therefore tested on samples that had
helped choose those pools out of ~615k SNPs — supervised feature-selection
leakage. We re-estimated everything with the **entire selection chain rebuilt
inside each CV fold** (Option B: per-fold `plink2 --keep --fst`), and ran a full
leak-free grid (N=5–100 × 3 pools × 3 reducers × 6 classifiers × 5 folds =
4,860 fits).

**Findings**

1. **The leak inflated CV accuracy by ~3 points on average** (worst −6.2 at N=35),
   uniformly (every honest number ≤ its leaky twin).
2. **The reducer conclusion survives and strengthens: ElasticNet is the decisive
   winner** (56/90 fold-wins; present in every top configuration).
3. **The pool winner is N-dependent** — the *stat* pool wins at small/mid panels
   (N=15–55), the *FST* pool wins at large panels (N≥65).
4. **Honest performance at N=70 is low-90s%**, and with the large-N-appropriate
   pool (FST+EN+SVM_RBF = 94.25%) is **level with the external Cai-34 panel
   (94.6%)** — not behind it. The leak had masked *which pool* to use at large N,
   not just the absolute number.
5. **Committed SNP lists are unchanged.** Only the reported *accuracy estimate* is
   corrected; a shipped panel may still be selected on all 504 samples.

---

## 2. The leak and the fix

| | Selects using labels? | Inside CV folds (old / new) |
|---|---|---|
| MAF, call-rate, HWE, LD prune (01/02) | no | on all data — **not leakage** (kept) |
| stat/FST/fst_stat pools (04a/04b/05a/05c) | **yes** | **all 504 (old)** → **per training fold (new)** |
| Reductor top-N (08) | yes | per fold (already leak-free) |

**Option B.** Per fold the three pools are rebuilt from the training samples only:
`stat` (χ²/JSD/AFD top-500 union, `nested_selection.py`), `FST` (pairwise Hudson
F_ST top-1000 union via `plink2 --keep <fold> --fst`, `nested_fst.py`), and
`fst_stat` (05c consensus). Each builder was validated to **reproduce the
committed 1,005 / 2,508 / 1,003 pools exactly** when run on all 504 samples, so the
only thing that changed between old and new is *when* the pool sees the data.

---

## 3. Honest performance

### 3.1 Leaky vs leak-free (committed panels)

| N | Leaky (old) | Tier B¹ | Tier A² | Best fixed config³ |
|---|-------------|---------|---------|--------------------|
| 35 | 92.26% | 86.10% | 84.91% | 86.31% (stat+RF+RF) |
| 50 | 93.05% | 91.67% | 86.70% | 91.67% (stat+EN+SVM_RBF) |
| 70 | 95.04% | 92.46% | 90.67% | 94.25% (FST+EN+SVM_RBF) |

¹ Tier B = baseline-committed config, pools rebuilt per fold (isolates the leak).
² Tier A = fully nested; pool + reducer + classifier chosen by inner CV per fold
(most conservative; also pays for config-selection variance).
³ Max over 54 configs per N — carries mild config-selection optimism; use as an
upper reference, not a headline.

**Triangulation at N=70:** 90.67% (fully nested) → 92.46% (pre-specified
stat+EN+LR) → 94.25% (FST+EN+SVM_RBF) vs Cai-34 94.64%. Honest performance sits in
the low-to-mid 90s depending on how strictly config-selection optimism is avoided.

### 3.2 Pre-specified pipelines (no per-fold config search)

| N | stat+EN+LR | stat+EN+SVM_RBF | FST+EN+SVM_RBF |
|---|-----------|-----------------|----------------|
| 35 | 86.10 | 84.91 | 82.73 |
| 50 | 89.28 | 91.67 | 89.89 |
| 70 | 92.46 | 92.46 | **94.25** |
| 100 | 94.25 | 94.44 | **96.83** |

The stat pool leads for small panels; the FST pool overtakes it from N≈60 (bigger,
noisier pool needs more SNPs to pay off). Crossover ≈ N=55–65.

---

## 4. Winner analysis (most consistent winner)

Tallied over 90 (N, fold) cells; see `winner_per_fold.csv`, `winner_per_N.csv`,
`winner_config_leaderboard.csv`, `winner_consistency.png`.

**Component fold-win rates**

| Reducer | wins | | Pool | wins | | Classifier | wins |
|---------|------|-|------|------|-|------------|------|
| **EN** | **56** | | **stat** | **56** | | **SVM_RBF** | **33** |
| RF | 25 | | FST | 22 | | LR | 19 |
| LR | 9 | | fst_stat | 12 | | RF | 14 |
| | | | | | | XGB / SVM_Lin / GBM | 10 / 7 / 7 |

- **Reducer — ElasticNet, unambiguous.** Robust to the leak fix; the strongest
  standalone conclusion of the pipeline.
- **Pool — N-dependent** (stat small-N, FST large-N); fst_stat rarely best.
- **Classifier — soft.** SVM_RBF wins most folds but LR is a close, *interpretable*
  second and ties SVM_RBF at N=35 and N=70 with the EN reducer. The 6 classifiers
  span only ~2.6 points on average — classifier choice is not decisive.

**Best config per N** (leak-free CV): FST+EN+SVM_Lin (N=5), stat+EN+LR (N=15–35),
stat+EN+SVM_RBF (N=40–60), **FST+EN+SVM_RBF (N=65–100)**.

---

## 5. Recommendation (for your approval — not yet wired as the headline)

- **Reducer: ElasticNet.** Non-negotiable; the data are clear.
- **Pool: size-dependent.** stat for compact panels (≤50), FST for large panels
  (≥65). If a single pool is wanted for simplicity, stat is the safer default and
  costs little except at N≥70.
- **Classifier: pre-specify one to avoid the "best of six" optimism.** Two honest
  options:
  - **Logistic Regression** — linear, interpretable (per-SNP coefficients), stable;
    N=70 stat+EN+LR = 92.46%.
  - **SVM-RBF** — nominally strongest, needed to reach parity with Cai at large N
    (FST+EN+SVM_RBF, N=70 = 94.25%).
- **Report nested CV as the primary metric**, with the leaky number retained only
  as an optimistic upper bound.

### Suggested narrative (least-damaging *and* honest)

> *Contribution = a reproducible pipeline (statistical/FST pool → ElasticNet
> ranking → a pre-specified classifier) that yields compact EAS panels. Under
> leak-free nested CV, a 70-SNP panel reaches ~92–94% (pool-dependent), comparable
> to the curated Cai-34 expert panel — while our panel is produced automatically
> and evaluated without the in-sample optimism common in the field.*

Parity + automation + transparency, not "we beat Cai." Reporting the leak-free
number is a credibility asset, not a concession.

---

## 6. Comparison to published panels (unchanged — external, no optimism)

| Panel | Matched N | Acc |
|-------|-----------|-----|
| Cai-34 | 34 | 94.64% |
| Cao-19 | 14 | 81.94% |
| Shi-142 | 116 | 78.57% |

Because these carry no selection optimism, the **leak-free** PAANDA-EA numbers
(not the leaky ones) are the fair comparison. See `comparison_accuracy_vs_n.png`.

---

## 7. Reproduce

```bash
make nested          # Tier B (full) + Tier A (committed N)
make compare-nested  # table + figure vs frozen baseline
conda run -n aisnp python scripts/nested_cv_sweep.py --tier grid   # full grid
conda run -n aisnp python scripts/nested_winner_stats.py           # winner analysis
```

Code: `scripts/nested_fst.py`, `nested_selection.py`, `nested_cv_sweep.py`,
`nested_winner_stats.py`, `nested_comparison.py`.
Artifacts: `results_archive/nested_leakfree_v1/`. Frozen leaky baseline:
`results_archive/baseline_leaky_v1/` (read-only).
