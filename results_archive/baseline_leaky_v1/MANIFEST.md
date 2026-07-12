# Baseline snapshot — `baseline_leaky_v1`

**Frozen:** 2026-07-13
**Git commit at snapshot:** `a803362` ("Change README")
**Source on disk:** `/mnt/data/aisnp_data/1000genomes/outputs/self_evaluation/`
**Status:** READ-ONLY. Do not edit or regenerate in place. This is the reference
we compare the leak-free (nested-CV) version against.

## Why this exists

The current pipeline selects the statistical / FST candidate pools (stat = 1,005
SNPs, FST = 2,508, fst_stat = 1,003) **once, on all 504 labelled samples**, in
notebooks `04a/04b/05a/05c` — *before* the notebook-08 cross-validation. Every CV
fold in stages 1–3 therefore tests on samples that helped choose those pools out
of ~615k SNPs. This is supervised feature-selection leakage (Ambroise & McLachlan
2002). The notebook-08 *reductor* step is leak-free (re-fit per fold), but the
upstream statistical filter is not.

These numbers are the **optimistic / leaky** baseline. The forthcoming nested-CV
version moves the entire supervised chain (χ²/JSD/AFD/FST → pool → reductor →
top-N → classifier) inside the outer training fold. The gap between the two = the
leakage magnitude.

## Key baseline numbers (to beat / compare against)

### Stage 2 — primary reported CV accuracy (5-fold), leaky
| N | panel | reductor | classifier | acc_mean | mcc_mean |
|---|-------|----------|-----------|----------|----------|
| 35 | stat | EN | LR      | 0.9226 | 0.8798 |
| 50 | stat | EN | SVM_RBF | 0.9305 | 0.8925 |
| 70 | stat | EN | SVM_RBF | 0.9504 | 0.9234 |

Full curve N=5..100 in `self_evaluation/08_unified_panel_sweep/stage2_best.csv`.

### Stage 3 — single 80/20 hold-out, leaky pools
| N | classifier | acc | mcc |
|---|-----------|------|------|
| 35 | LR      | 0.8812 | 0.8176 |
| 50 | SVM_RBF | 0.9010 | 0.8474 |
| 70 | SVM_RBF | 0.9010 | 0.8491 |

### Benchmark vs published (from `11_results/table2_published_panels.csv`)
| Panel | Matched N | Best clf | Accuracy |
|-------|-----------|----------|----------|
| PAANDA-EA N=35 | 35 | LR | 0.9226 |
| PAANDA-EA N=50 | 50 | SVM_RBF | 0.9305 |
| PAANDA-EA N=70 | 70 | SVM_RBF | 0.9504 |
| cai_eas34 (external) | 34 | SVM_RBF | 0.9464 |
| cao_19 (external) | 14 | LR | 0.8194 |
| shi_142 (external) | 116 | SVM_RBF | 0.7857 |

Note: the published panels (Cai/Cao/Shi) are **external** — no selection optimism.
The PAANDA-EA rows above carry the leak; the nested-CV rerun is the fair
comparison to the external panels.

## Contents

- `self_evaluation/08_unified_panel_sweep/` — stage1/2/3 CSVs, per-N panel lists, OOF preds, figures
- `self_evaluation/09_published_panel_comparison/` — published panel CV results
- `self_evaluation/10_panel_overlap/` — rsID overlap heatmaps
- `self_evaluation/11_results/` — final figures + `table1_our_pipeline.csv`, `table2_published_panels.csv`
- `CHECKSUMS.md5` — md5 of every archived file; verify with `md5sum -c CHECKSUMS.md5`

## Verify integrity

```bash
cd results_archive/baseline_leaky_v1 && md5sum -c CHECKSUMS.md5
```
