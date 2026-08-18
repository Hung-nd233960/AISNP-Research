"""Retroactive breakdown of out-of-fold predictions by original 1000
Genomes subpopulation code (CHB/CHS/JPT/KHV/CDX), not just the merged
Han/JPT/SEA label used for training. No retraining -- reuses the exact
per-sample out-of-fold predictions already computed for the committed
panels and Cai-34, joined against the original panel file's population
codes purely for post-hoc reporting.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from confusion_matrices import COMMITTED, OOFPRED_PKL, Data, out_of_fold_predictions

REPO = Path(__file__).resolve().parent.parent
PANEL_FILE = Path("/mnt/data/aisnp_data/1000genomes/raw_downloads/integrated_call_samples_v3.20130502.ALL.panel")

SUBPOP_TO_MERGED = {"CHB": "Han", "CHS": "Han", "JPT": "JPT", "KHV": "SEA", "CDX": "SEA"}

data = Data()
class_order = list(data.le.classes_)  # merged 3-class encoding order
merged_label = {i: c for i, c in enumerate(class_order)}
DISPLAY = {"CN": "Han", "JPT": "JPT", "SEA": "SEA"}

panel_df = pd.read_csv(PANEL_FILE, sep="\t")
subpop_lut = dict(zip(panel_df["sample"], panel_df["pop"]))
subpops = np.array([subpop_lut[s] for s in data.samples])
assert set(np.unique(subpops)) == set(SUBPOP_TO_MERGED), f"unexpected codes: {set(np.unique(subpops))}"

# sanity check: does the subpop-derived merged label match the label actually used?
derived_merged = np.array([SUBPOP_TO_MERGED[s] for s in subpops])
actual_merged = np.array([DISPLAY[c] for c in data.y_str])
assert np.array_equal(derived_merged, actual_merged), "subpopulation-to-merged mapping does not match the labels actually used"
print("Sanity check passed: subpopulation codes correctly reconstruct the merged Han/JPT/SEA labels.\n")

subpop_n = pd.Series(subpops).value_counts()
print("Subpopulation sample counts:", dict(subpop_n), "\n")


def per_subpop_recall(y_true_idx: np.ndarray, y_pred_idx: np.ndarray) -> pd.DataFrame:
    y_true_lab = np.array([merged_label[i] for i in y_true_idx])
    y_pred_lab = np.array([merged_label[i] for i in y_pred_idx])
    rows = []
    for sp in ["CHB", "CHS", "JPT", "KHV", "CDX"]:
        mask = subpops == sp
        n = mask.sum()
        correct = (y_pred_lab[mask] == y_true_lab[mask]).sum()
        rows.append({"subpop": sp, "merged_into": SUBPOP_TO_MERGED[sp], "n": n,
                     "correct": correct, "recall_%": round(100 * correct / n, 2)})
    return pd.DataFrame(rows)


print("=" * 70)
print("Committed panels")
print("=" * 70)
for label, n, pool, reductor, clf in COMMITTED:
    y_true, y_pred = out_of_fold_predictions(data, pool, reductor, clf, n)
    df = per_subpop_recall(y_true, y_pred)
    print(f"\n{label}")
    print(df.to_string(index=False))

print("\n" + "=" * 70)
print("Cai-34 (external)")
print("=" * 70)
with open(OOFPRED_PKL, "rb") as fh:
    d = pickle.load(fh)
pkl_classes = [str(c) for c in d["classes"]]
assert pkl_classes == class_order
y_true_cai, y_pred_cai = d["oofpreds"]["cai_eas34"]["SVM_RBF"]
df_cai = per_subpop_recall(y_true_cai, y_pred_cai)
print(df_cai.to_string(index=False))
