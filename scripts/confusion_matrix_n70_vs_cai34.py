"""Special two-panel confusion matrix: N=70 (this pipeline's committed
panel) on the left, Cai-34 (external expert-curated panel) on the right.

For manual insertion into the poster (not referenced by REPORT.md). Reuses
the exact data-computation functions from confusion_matrices.py -- no new
model fitting beyond what that script already does, same leak-free
out-of-fold predictions.

Writes: reports/figures/poster_confusion_n70_vs_cai34.png (+ .svg)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from confusion_matrices import (
    Data,
    INK,
    INK_SECONDARY,
    LABEL_DISPLAY,
    confusion,
    load_cai34_confusion,
    out_of_fold_predictions,
)

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "reports/figures"
OUT.mkdir(parents=True, exist_ok=True)


def plot(cms: dict[str, np.ndarray], class_order: list[str]) -> None:
    labels = [LABEL_DISPLAY[c] for c in class_order]
    vmax = max(cm.max() for cm in cms.values())
    titles = list(cms.keys())

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.6), dpi=320)
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
                    fontsize=13, fontweight="bold" if i == j else "normal",
                )
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xlabel("Predicted")
        if title == titles[0]:
            ax.set_ylabel("True")
        acc = np.trace(cm) / cm.sum()
        ax.set_title(f"{title}\nacc {acc * 100:.1f}%", fontsize=12, color=INK)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)

    fig.suptitle(
        "N = 70 vs. Cai-34: out-of-fold confusion matrices",
        fontsize=12.5, color=INK_SECONDARY,
    )
    fig.tight_layout()
    fig.savefig(OUT / "poster_confusion_n70_vs_cai34.png", bbox_inches="tight")
    fig.savefig(OUT / "poster_confusion_n70_vs_cai34.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = Data()
    class_order = list(data.le.classes_)

    y_true, y_pred = out_of_fold_predictions(data, "FST", "EN", "SVM_RBF", 70)
    cm_n70 = confusion(y_true, y_pred, len(class_order))
    acc = np.trace(cm_n70) / cm_n70.sum()
    print(f"N=70 (FST+EN+SVM-RBF): acc={acc * 100:.2f}%\n{cm_n70}")

    cm_cai = load_cai34_confusion(class_order)
    acc_cai = np.trace(cm_cai) / cm_cai.sum()
    print(f"Cai-34: acc={acc_cai * 100:.2f}%\n{cm_cai}")

    cms = {"N = 70 (FST+EN+SVM-RBF)": cm_n70, "Cai-34 (external, SVM-RBF)": cm_cai}
    plot(cms, class_order)
    print(f"Wrote {OUT / 'poster_confusion_n70_vs_cai34.png'}")


if __name__ == "__main__":
    main()
