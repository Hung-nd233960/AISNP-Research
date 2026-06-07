# TODO — Unified Panel Sweep

## STATUS: complete — results in `outputs/self_evaluation/08_unified_panel_sweep/`

System compares **3 candidate panels** via a single **3-stage unified notebook**.

| path      | candidate set |
|-----------|---------------|
| FST       | 2508 FST SNPs (`PATHS.ML_DATA`) |
| stat      | 1005 union top-500 SNPs (`05a` feature importance file) |
| FST+stat  | 1003 consensus SNPs (`05c` rf_feature_importances.csv) |

## Key findings

- **stat + ElasticNet** wins Stage 1 at almost every N (≥15).
- **Best classifier** varies: LR at most N, SVM_RBF at higher N (≥50).
- Performance range: ~69% acc at N=5 → ~95% acc at N=55–75 (Stage 2 CV).
- Stage 3 single-split confirmation is 5–8% below Stage 2 CV — expected due to
  double selection optimism (Stage 1 + Stage 2 both selected by CV on all 504 samples).

## Architecture: 3 stages in `08_unified_panel_sweep.ipynb`

**Stage 1 — Reductor × Panel selection (nested CV)**
- For each N in fine grid [5…100]:
  - 3 panels × 3 reductors (RF, LR, ElasticNet) = 9 configs
  - 5-fold nested CV; ranker re-fit inside each training fold → no leakage
  - Fixed RF eval classifier (fast)
  - → best (panel, reductor) per N

**Stage 2 — Classifier evaluation (nested CV, fixed config)**
- For each N: use Stage-1 winner config, run 6 classifiers (RF, XGB, LR, SVM_RBF, SVM_Lin, GBM)
- Ranker re-fit inside each fold (same leak-free design)
- Reports acc / F1 / MCC / ROC-AUC per classifier
- → best classifier per N (Stage 2 numbers = primary reported performance)

**Stage 3 — Final panel commitment**
- For each N: 80/20 stratified split
- Fit reductor on 80% → top-N SNPs → **final panel list**
- Train best Stage-2 classifier on 80%, confirm on 20%
- Saves `panels/panel_N{N:03d}.csv` per N
- Stage 3 acc = single-split confirmation, NOT the primary metric

## Leakage disclosure

| Step | Uses labels on | Type | Severity |
|------|----------------|------|----------|
| FST / chi² / JSD candidate set | All 504 samples | Filter-method, not model-fit | Mild, field-standard |
| Stage 1 config selection | All 504 samples (via CV) | Model-selection bias | Mild, unavoidable at n=504 |
| Ranker inside each CV fold | Training fold only | None | Clean |
| Classifier training | Training fold only | None | Clean |

## Run order (manual)

Prerequisites (must exist before running 08):
- `PATHS.ML_DATA` (from `04b_fst_and_pca`)
- `statistical/05a_stat_training/05a_stat_training_feature_importance.csv` (from `05a`)
- `statistical/05c_fst_stat_training/rf_feature_importances.csv` (from `05c`)

Run: `05a → 05c → 08_unified_panel_sweep`

## Output files (from 08)

- `stage1_results.csv` / `stage1_best.csv` — best (panel, reductor) per N
- `stage2_results.csv` / `stage2_agg.csv` / `stage2_best.csv` — classifier performance
- `stage3_results.csv` — final confirmation scores
- `panels/panel_N{N:03d}.csv` — final SNP list per N (rank, snp_id, reductor_score)
- `stage1_overview.png` / `stage2_performance.png` — plots

## Verification after run

- `stage1_best.csv` has 18 rows (one per N), `panel` and `reductor` columns.
- `stage2_best.csv` has 18 rows with acc/F1/MCC/ROC-AUC.
- `panels/` contains 18 CSV files (panel_N005.csv … panel_N100.csv).
- Summary table printed at end shows all 3 stages side-by-side.

## Separate track (unchanged)

`08_known_aisnps_ml` still reads `PATHS.STAT_ML_DATA` (published-panel comparison).
`STAT_SCORES / STAT_ML_DATA / STAT_ALL4` config props kept for it.
See `COMPARISON.md` for collection status of published panels.

## Old 05x notebooks (kept for exploratory / elbow analysis)

`05b_fst_only_training`, `05b_reduction`, `05c_fst_stat_training` remain intact
for their elbow plots, PCA, and single-path sweeps. They are NOT the primary
comparison path — `08_unified_panel_sweep` supersedes them for cross-path
comparison and final panel selection.

`06_self_evaluation.ipynb` — deleted (superseded by 08 Stage 2 summary).
`07_reduction_methods.ipynb` — never created; Stage 1 of 08 covers the same question.
