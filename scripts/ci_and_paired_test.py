"""Retroactive confidence intervals and paired significance testing --
purely post-hoc on already-archived/cached results, no retraining.

Part A: 95% CI on Blind Selection and Best Fixed Configuration accuracy at
the three committed panel sizes (N=35/50/70), from already-archived
per-fold values (nested_tierA_results.csv, grid_results.csv).

Part B: paired comparison at N=70 between this pipeline's committed panel
(FST+EN+SVM-RBF) and the external Cai-34 panel, on the *same* 5 outer
folds. Cai-34's cached out-of-fold predictions were built with the same
StratifiedKFold(seed=42) as everything else in this codebase, so fold
membership is reconstructed by re-running just the splitter (no model
fitting) and slicing the cached prediction arrays -- verified with a
sanity check against the known aggregate 94.64%.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nested_cv_sweep import N_CV, RANDOM_STATE, Data

REPO = Path(__file__).resolve().parent.parent
ARCHIVE = REPO / "results_archive/nested_leakfree_v1/08b_nested_cv_sweep"
OOFPRED_PKL = (
    REPO / "results_archive/baseline_leaky_v1/self_evaluation/09_published_panel_comparison"
    / "panel_oofpreds.pkl"
)

T_CRIT_DF4 = stats.t.ppf(0.975, df=4)  # 95% CI, n=5, df=n-1=4


def ci95(values: np.ndarray) -> tuple[float, float, float]:
    mean = values.mean()
    sd = values.std(ddof=1)
    se = sd / np.sqrt(len(values))
    return mean, mean - T_CRIT_DF4 * se, mean + T_CRIT_DF4 * se


print("=" * 70)
print("PART A: 95% CI at each committed panel size")
print("=" * 70)

tier_a = pd.read_csv(ARCHIVE / "nested_tierA_results.csv")
grid = pd.read_csv(ARCHIVE / "grid_results.csv")

BFC_WINNERS = {
    35: ("fst_stat", "EN", "SVM_RBF"),
    50: ("stat", "EN", "SVM_RBF"),
    70: ("FST", "EN", "SVM_RBF"),
}

for N in [35, 50, 70]:
    row = tier_a[tier_a.n_snps == N].iloc[0]
    # tier_a stores mean+std (not raw per-fold values) -- CI from those directly
    se = row.acc_std / np.sqrt(N_CV)
    lo, hi = row.acc_mean - T_CRIT_DF4 * se, row.acc_mean + T_CRIT_DF4 * se
    print(f"BS  N={N:3d}: mean={row.acc_mean*100:.2f}%  sd={row.acc_std*100:.2f}  "
          f"95% CI [{lo*100:.2f}%, {hi*100:.2f}%]")

    pool, red, clf = BFC_WINNERS[N]
    sub = grid[(grid.n_snps == N) & (grid.panel == pool) & (grid.reductor == red) & (grid.classifier == clf)]
    assert len(sub) == 5, f"expected 5 folds, got {len(sub)}"
    vals = sub.sort_values("fold")["acc"].values
    mean, lo, hi = ci95(vals)
    print(f"BFC N={N:3d}: mean={mean*100:.2f}%  sd={vals.std(ddof=1)*100:.2f}  "
          f"95% CI [{lo*100:.2f}%, {hi*100:.2f}%]  folds={np.round(vals*100,2)}")
    print()

print("=" * 70)
print("PART B: Paired comparison at N=70, FST+EN+SVM-RBF vs. Cai-34")
print("=" * 70)

data = Data()
class_order = list(data.le.classes_)

# Reconstruct the exact fold split used throughout this codebase.
skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
fold_of = np.empty(len(data.y), dtype=int)
for fold_idx, (_, te) in enumerate(skf.split(data.G, data.y)):
    fold_of[te] = fold_idx

# Our N=70 committed panel's per-fold accuracy (already archived).
sub70 = grid[(grid.n_snps == 70) & (grid.panel == "FST") & (grid.reductor == "EN") & (grid.classifier == "SVM_RBF")]
ours = sub70.sort_values("fold")["acc"].values
print("Ours  (FST+EN+SVM-RBF) per-fold accuracy:", np.round(ours * 100, 2))

# Cai-34's cached OOF predictions, sliced into the same 5 folds.
with open(OOFPRED_PKL, "rb") as fh:
    d = pickle.load(fh)
pkl_classes = [str(c) for c in d["classes"]]
assert pkl_classes == class_order, (pkl_classes, class_order)
y_true_cai, y_pred_cai = d["oofpreds"]["cai_eas34"]["SVM_RBF"]

# Sanity check: does the cached array's sample order match data.y's order?
# If so, slicing by fold_of (computed from data.y) reproduces the known 94.64%.
overall_acc = (y_true_cai == y_pred_cai).mean()
assert abs(overall_acc - 0.9464) < 0.001, f"cache mismatch: {overall_acc}"
assert np.array_equal(y_true_cai, data.y), "Cai-34 cache sample order does NOT match Data.y order -- cannot slice by fold_of directly"

cai_fold_acc = np.array([
    (y_true_cai[fold_of == f] == y_pred_cai[fold_of == f]).mean()
    for f in range(N_CV)
])
print("Cai-34 (external, SVM-RBF) per-fold accuracy:", np.round(cai_fold_acc * 100, 2))
reconstructed_overall = (y_true_cai == y_pred_cai).mean()
print(f"Reconstructed overall Cai-34 accuracy: {reconstructed_overall*100:.2f}% (known: 94.64%) -- {'MATCH' if abs(reconstructed_overall-0.9464)<0.001 else 'MISMATCH'}")
print()

diff = ours - cai_fold_acc
print("Per-fold difference (ours - Cai-34):", np.round(diff * 100, 2))
mean_diff, lo_diff, hi_diff = ci95(diff)
print(f"Mean paired difference: {mean_diff*100:+.2f} pts, 95% CI [{lo_diff*100:+.2f}, {hi_diff*100:+.2f}]")

t_stat, t_p = stats.ttest_rel(ours, cai_fold_acc)
print(f"Paired t-test: t={t_stat:.3f}, p={t_p:.3f}")

w_stat, w_p = stats.wilcoxon(ours, cai_fold_acc)
print(f"Wilcoxon signed-rank: W={w_stat:.3f}, p={w_p:.3f}")
