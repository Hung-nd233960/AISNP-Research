"""Quick 3-pie-chart figure: distribution of Best Fixed Configuration's
per-panel-size winner, aggregated across all 18 panel sizes (one winner per
N, picked from the 5-fold CV mean) -- for side-by-side comparison against
Fig. 3 (fig3_winner_consistency.png), which tallies the grid's fold-cell
winner across all 90 panel-size x fold cells instead.

Same color mapping and flat-design conventions as report_figures.py's
fig3_winner_consistency(), so the two figures read as one system.

Writes: reports/figures/bfc_18n_winner_pies.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
ARCHIVE = REPO / "results_archive/nested_leakfree_v1/08b_nested_cv_sweep"
OUT = REPO / "reports/figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#2a78d6"
GREEN = "#008300"
MAGENTA = "#e87ba4"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
SURFACE = "#fcfcfb"

REDUCTOR_LABELS = {"EN": "ElasticNet", "RF": "RF importance", "LR": "L1-LR"}
CLASSIFIER_LABELS = {
    "SVM_RBF": "SVM-RBF", "SVM_Lin": "SVM-Lin", "LR": "LR",
    "RF": "RF", "XGB": "XGB", "GBM": "GBM",
}
classifier_color = {
    "SVM_RBF": BLUE, "LR": GREEN, "RF": MAGENTA,
    "XGB": "#eda100", "SVM_Lin": "#1baf7a", "GBM": "#eb6834",
}
reductor_color = {"EN": BLUE, "RF": GREEN, "LR": MAGENTA}
pool_color = {"stat": BLUE, "FST": GREEN, "fst_stat": MAGENTA}

grid = pd.read_csv(ARCHIVE / "grid_results.csv")
fold_mean = grid.groupby(["n_snps", "panel", "reductor", "classifier"], as_index=False)["acc"].mean()
Ns = sorted(grid.n_snps.unique())
winners = pd.DataFrame([
    fold_mean[fold_mean.n_snps == N].sort_values("acc", ascending=False).iloc[0]
    for N in Ns
])
total = len(winners)

panels = [
    ("Pool", winners["panel"].value_counts(), {}, pool_color),
    ("Reducer", winners["reductor"].value_counts(), REDUCTOR_LABELS, reductor_color),
    ("Classifier", winners["classifier"].value_counts(), CLASSIFIER_LABELS, classifier_color),
]

fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.6), dpi=220)
fig.patch.set_facecolor(SURFACE)
for ax, (title, counts, relabel, colors) in zip(axes, panels):
    counts = counts.sort_values(ascending=False)
    labels = [relabel.get(k, k) for k in counts.index]
    wedge_colors = [colors[k] for k in counts.index]
    wedges, _, autotexts = ax.pie(
        counts.values,
        colors=wedge_colors,
        autopct=lambda p: f"{p:.0f}%",
        pctdistance=0.75,
        startangle=90,
        counterclock=False,
        wedgeprops={"linewidth": 2, "edgecolor": SURFACE},
        textprops={"color": "white", "fontsize": 9, "fontweight": "bold"},
    )
    ax.legend(
        wedges, [f"{l} ({c}/{total})" for l, c in zip(labels, counts.values)],
        loc="upper center", bbox_to_anchor=(0.5, -0.02), frameon=False,
        fontsize=8.5, ncol=1,
    )
    ax.set_title(title, fontsize=11, color=INK)

fig.suptitle(
    "Best Fixed Configuration's per-N winner (18 panel sizes, 5-fold mean)",
    fontsize=11.5, color=INK_SECONDARY,
)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(OUT / "bfc_18n_winner_pies.png", bbox_inches="tight", facecolor=SURFACE)
print(f"Wrote {OUT / 'bfc_18n_winner_pies.png'}")
print("Pool:", dict(winners["panel"].value_counts()))
print("Reducer:", dict(winners["reductor"].value_counts()))
print("Classifier:", dict(winners["classifier"].value_counts()))
