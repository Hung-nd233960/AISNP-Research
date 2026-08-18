"""Out-of-fold confusion matrices for the three committed panels vs. Cai-34.

Compares the Best Fixed Configuration at each committed panel size (N=35:
fst_stat+EN+SVM-RBF, N=50: stat+EN+SVM-RBF, N=70: FST+EN+SVM-RBF — matching
Table X) against the external, expert-curated Cai-34 panel (N=34, SVM-RBF —
matching Table XII), so error patterns can be read side by side across panel
sizes and against the external benchmark.

Leak-free for the three committed panels: each pool is rebuilt from the outer
training fold only, exactly as in scripts/nested_cv_sweep.py. Concatenating
predictions across the 5 outer test folds gives one out-of-fold prediction per
sample, covering all 504 samples exactly once.

Cai-34 needs no new computation: its out-of-fold SVM-RBF predictions are
already cached in
results_archive/baseline_leaky_v1/self_evaluation/09_published_panel_comparison/
panel_oofpreds.pkl (published panels carry no selection optimism, so this
cache is valid regardless of the leak fix — verified here to reproduce
Table XII's 94.64% exactly before use).

Writes: reports/figures/fig7_confusion_matrices.png
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nested_cv_sweep import N_CV, RANDOM_STATE, Data, fit_reductor, make_classifiers

REPO = Path(__file__).resolve().parent.parent
OOFPRED_PKL = (
    REPO
    / "results_archive/baseline_leaky_v1/self_evaluation/09_published_panel_comparison"
    / "panel_oofpreds.pkl"
)
OUT = REPO / "reports/figures"
OUT.mkdir(parents=True, exist_ok=True)

# (label, N, pool, reductor, classifier) — the four to plot, in display order.
COMMITTED = [
    ("N=35 (fst_stat+EN+SVM-RBF)", 35, "fst_stat", "EN", "SVM_RBF"),
    ("N=50 (stat+EN+SVM-RBF)", 50, "stat", "EN", "SVM_RBF"),
    ("N=70 (FST+EN+SVM-RBF)", 70, "FST", "EN", "SVM_RBF"),
]
LABEL_DISPLAY = {"CN": "Han", "JPT": "JPT", "SEA": "SEA"}
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"


def out_of_fold_predictions(
    data: Data, pool_name: str, reductor: str, classifier: str, n_target: int
) -> tuple[np.ndarray, np.ndarray]:
    skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
    ctmpl, needs_scale = make_classifiers()[classifier]

    y_true = np.empty(len(data.y), dtype=int)
    y_pred = np.empty(len(data.y), dtype=int)

    for tr, te in skf.split(data.G, data.y):
        pools = data.build_pools(tr, which=[pool_name])
        col_idx = pools[pool_name]
        Xtr_pool = data.G[tr][:, col_idx].astype(np.float32)
        Xte_pool = data.G[te][:, col_idx].astype(np.float32)
        y_tr = data.y[tr]

        imp, scaler = fit_reductor(reductor, Xtr_pool, y_tr)
        top = np.argsort(imp)[::-1][:n_target]

        if needs_scale:
            Xa = scaler.transform(Xtr_pool)[:, top]
            Xb = scaler.transform(Xte_pool)[:, top]
        else:
            Xa = Xtr_pool[:, top]
            Xb = Xte_pool[:, top]

        clf = clone(ctmpl)
        clf.fit(Xa, y_tr)
        y_pred[te] = clf.predict(Xb)
        y_true[te] = data.y[te]

    return y_true, y_pred


def confusion(y_true: np.ndarray, y_pred: np.ndarray, n_classes: int) -> np.ndarray:
    cm = np.zeros((n_classes, n_classes), dtype=int)
    for t, p in zip(y_true, y_pred):
        cm[t, p] += 1
    return cm


def load_cai34_confusion(class_order: list[str]) -> np.ndarray:
    with open(OOFPRED_PKL, "rb") as fh:
        d = pickle.load(fh)
    pkl_classes = [str(c) for c in d["classes"]]
    assert pkl_classes == class_order, (pkl_classes, class_order)
    y_true, y_pred = d["oofpreds"]["cai_eas34"]["SVM_RBF"]
    acc = (y_true == y_pred).mean()
    assert abs(acc - 0.9464) < 0.001, f"Cai-34 cache mismatch: acc={acc}"
    return confusion(y_true, y_pred, len(class_order))


def plot(cms: dict[str, np.ndarray], class_order: list[str]) -> None:
    labels = [LABEL_DISPLAY[c] for c in class_order]
    vmax = max(cm.max() for cm in cms.values())
    titles = list(cms.keys())

    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.4), dpi=200)
    for ax, title in zip(axes, titles):
        cm = cms[title]
        ax.imshow(cm, cmap="Blues", vmin=0, vmax=vmax)
        for i in range(len(labels)):
            for j in range(len(labels)):
                val = cm[i, j]
                frac = val / vmax if vmax else 0
                ax.text(
                    j, i, str(val), ha="center", va="center",
                    color="white" if frac > 0.55 else INK,
                    fontsize=11, fontweight="bold" if i == j else "normal",
                )
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xlabel("Predicted")
        if title == titles[0]:
            ax.set_ylabel("True")
        acc = np.trace(cm) / cm.sum()
        ax.set_title(f"{title}\nacc {acc * 100:.1f}%", fontsize=10, color=INK)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)

    fig.suptitle(
        "Out-of-fold confusion matrices: committed panels vs. external Cai-34",
        fontsize=11, color=INK_SECONDARY,
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig7_confusion_matrices.png", bbox_inches="tight")
    fig.savefig(OUT / "fig7_confusion_matrices.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = Data()
    class_order = list(data.le.classes_)

    cms: dict[str, np.ndarray] = {}
    for label, n, pool, reductor, clf in COMMITTED:
        y_true, y_pred = out_of_fold_predictions(data, pool, reductor, clf, n)
        cm = confusion(y_true, y_pred, len(class_order))
        cms[label] = cm
        acc = np.trace(cm) / cm.sum()
        print(f"{label}: acc={acc * 100:.2f}%\n{cm}")

    cms["Cai-34 (external, SVM-RBF)"] = load_cai34_confusion(class_order)
    acc = np.trace(cms["Cai-34 (external, SVM-RBF)"]) / cms["Cai-34 (external, SVM-RBF)"].sum()
    print(f"Cai-34: acc={acc * 100:.2f}%\n{cms['Cai-34 (external, SVM-RBF)']}")

    plot(cms, class_order)
    print(f"Wrote {OUT / 'fig7_confusion_matrices.png'}")


if __name__ == "__main__":
    main()
