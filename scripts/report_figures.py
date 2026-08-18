"""Generate the figures embedded in reports/REPORT.md from archived nested-CV results.

No new model fitting — purely reads already-computed results and plots them.

Reads:
  results_archive/nested_leakfree_v1/08b_nested_cv_sweep/{grid_results,winner_per_N,
    winner_per_fold,nested_tierA_results}.csv
  results_archive/baseline_leaky_v1/self_evaluation/09_published_panel_comparison/
    published_panel_results.csv

Writes (PNG, 200 dpi, Word/LaTeX-embeddable; numbered to match their figure
number in reports/REPORT.md, not the order generated below):
  reports/figures/fig3_winner_consistency.png
  reports/figures/fig4_pool_crossover.png
  reports/figures/fig5_pool_divergence.png
  reports/figures/fig6_accuracy_vs_n.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
ARCHIVE = REPO / "results_archive/nested_leakfree_v1/08b_nested_cv_sweep"
PUBLISHED = (
    REPO
    / "results_archive/baseline_leaky_v1/self_evaluation/09_published_panel_comparison"
    / "published_panel_results.csv"
)
OUT = REPO / "reports/figures"
OUT.mkdir(parents=True, exist_ok=True)

# Validated categorical palette (dataviz skill reference), fixed slot order.
BLUE = "#2a78d6"
GREEN = "#008300"
MAGENTA = "#e87ba4"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID_COLOR = "#e1e0d9"
SURFACE = "#fcfcfb"

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "axes.edgecolor": INK_MUTED,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 10,
    }
)

PUBLISHED_LABELS = {"cai_eas34": "Cai-34", "cao_19": "Cao-19", "shi_142": "Shi-142"}
REDUCTOR_LABELS = {"EN": "ElasticNet", "RF": "RF importance", "LR": "L1-LR"}
CLASSIFIER_LABELS = {
    "SVM_RBF": "SVM-RBF", "SVM_Lin": "SVM-Lin", "LR": "LR",
    "RF": "RF", "XGB": "XGB", "GBM": "GBM",
}


def _clean_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", color=GRID_COLOR, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def fig6_accuracy_vs_n() -> None:
    best = pd.read_csv(ARCHIVE / "winner_per_N.csv").sort_values("n_snps")
    tier_a = pd.read_csv(ARCHIVE / "nested_tierA_results.csv")
    pub = pd.read_csv(PUBLISHED)
    pub = pub[pub["panel"].isin(PUBLISHED_LABELS)]
    pub_best = pub.loc[pub.groupby("panel")["acc"].idxmax()]

    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=200)
    ax.plot(
        best["n_snps"], best["best_acc"], color=BLUE, linewidth=2, marker="o",
        markersize=6, markerfacecolor=BLUE, markeredgecolor=SURFACE,
        markeredgewidth=1.2, label="Best Fixed Configuration", zorder=3,
    )
    tier_a = tier_a.sort_values("n_snps")
    if len(tier_a) > 3:
        # Full curve now available (all 18 N) — plot as a line like the other series.
        ax.plot(
            tier_a["n_snps"], tier_a["acc_mean"] * 100, color=GREEN, linewidth=2,
            marker="D", markersize=6, markerfacecolor=GREEN, markeredgecolor=SURFACE,
            markeredgewidth=1.2, label="Blind Selection", zorder=4,
        )
    else:
        # Sparse fallback: only the 3 committed panel sizes have been run.
        ax.scatter(
            tier_a["n_snps"], tier_a["acc_mean"] * 100, color=GREEN, s=90, marker="D",
            edgecolor=SURFACE, linewidth=1.2, label="Blind Selection", zorder=4,
        )
    for _, row in pub_best.iterrows():
        n, acc = row["n_snps"], row["acc"] * 100
        ax.scatter([n], [acc], color=INK_MUTED, marker="*", s=150,
                   edgecolor=SURFACE, linewidth=0.8, zorder=5)
        ax.annotate(
            PUBLISHED_LABELS.get(row["panel"], row["panel"]), (n, acc),
            textcoords="offset points", xytext=(6, -3), fontsize=8.5,
            color=INK_SECONDARY,
        )

    _clean_axes(ax)
    ax.set_xlabel("Panel size N")
    ax.set_ylabel("Cross-validated accuracy (%)")
    ax.set_ylim(60, 100)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig6_accuracy_vs_n.png")
    fig.savefig(OUT / "fig6_accuracy_vs_n.svg")
    plt.close(fig)


def fig4_pool_crossover() -> None:
    grid = pd.read_csv(ARCHIVE / "grid_results.csv")
    fold_mean = grid.groupby(
        ["n_snps", "panel", "reductor", "classifier"], as_index=False
    )["acc"].mean()
    best_per_pool = fold_mean.loc[fold_mean.groupby(["n_snps", "panel"])["acc"].idxmax()]

    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=200)
    colors = {"stat": BLUE, "FST": GREEN, "fst_stat": MAGENTA}
    for pool in ["stat", "FST", "fst_stat"]:
        d = best_per_pool[best_per_pool["panel"] == pool].sort_values("n_snps")
        ax.plot(
            d["n_snps"], d["acc"] * 100, color=colors[pool], linewidth=2, marker="o",
            markersize=6, markerfacecolor=colors[pool], markeredgecolor=SURFACE,
            markeredgewidth=1.0, label=pool, zorder=3,
        )

    _clean_axes(ax)
    ax.set_xlabel("Panel size N")
    ax.set_ylabel("Best cross-validated accuracy (%), by pool")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_pool_crossover.png")
    fig.savefig(OUT / "fig4_pool_crossover.svg")
    plt.close(fig)

    # Report the actual crossover N so the figure and Methods/Results prose agree.
    pivot = best_per_pool.pivot(index="n_snps", columns="panel", values="acc").sort_index()
    leader = pivot.idxmax(axis=1)
    print("Per-N pool leader:", dict(zip(pivot.index, leader)))


def fig3_winner_consistency() -> None:
    wf = pd.read_csv(ARCHIVE / "winner_per_fold.csv")
    total = len(wf)

    # Fixed categorical slot order (validated palette), assigned per pie
    # independently — reducer and classifier both reuse the codes "LR"/"RF"
    # for unrelated things, so a single shared color map would clash.
    classifier_color = {
        "SVM_RBF": BLUE, "LR": GREEN, "RF": MAGENTA,
        "XGB": "#eda100", "SVM_Lin": "#1baf7a", "GBM": "#eb6834",
    }
    reductor_color = {"EN": BLUE, "RF": GREEN, "LR": MAGENTA}
    pool_color = {"stat": BLUE, "FST": GREEN, "fst_stat": MAGENTA}

    panels = [
        ("Pool", wf["panel"].value_counts(), {}, pool_color),
        ("Reducer", wf["reductor"].value_counts(), REDUCTOR_LABELS, reductor_color),
        ("Classifier", wf["classifier"].value_counts(), CLASSIFIER_LABELS, classifier_color),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.6), dpi=200)
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

    fig.suptitle("Fold-win share by component (% of 90 panel-size × fold cells)",
                 fontsize=10.5, color=INK_SECONDARY, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_winner_consistency.png", bbox_inches="tight")
    fig.savefig(OUT / "fig3_winner_consistency.svg", bbox_inches="tight")
    plt.close(fig)

    for title, counts, _, _ in panels:
        print(f"{title} win counts:", counts.to_dict())


def fig5_pool_divergence() -> None:
    """Grid (hindsight) pool winner vs. Blind Selection's own blind pool pick,
    per N and in aggregate — the evidence that the grid's pool-crossover story
    is not corroborated by the bias-free estimate."""
    grid = pd.read_csv(ARCHIVE / "grid_results.csv")
    folds = pd.read_csv(ARCHIVE / "nested_tierA_folds.csv")
    pool_color = {"stat": BLUE, "FST": GREEN, "fst_stat": MAGENTA}

    fold_mean = grid.groupby(["n_snps", "panel", "reductor", "classifier"], as_index=False)["acc"].mean()
    Ns = sorted(grid.n_snps.unique())

    grid_winner_pool = {
        N: fold_mean[fold_mean.n_snps == N].sort_values("acc", ascending=False).iloc[0]["panel"]
        for N in Ns
    }
    bs_modal_pool = {N: folds[folds.n_snps == N]["panel"].value_counts().index[0] for N in Ns}

    fig = plt.figure(figsize=(13, 8), dpi=200)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.3], hspace=0.55)

    ax1 = fig.add_subplot(gs[0])
    for i, N in enumerate(Ns):
        gw, bm = grid_winner_pool[N], bs_modal_pool[N]
        ax1.add_patch(plt.Rectangle((i, 1), 0.9, 0.8, facecolor=pool_color[gw], edgecolor="none"))
        ax1.add_patch(plt.Rectangle((i, 0), 0.9, 0.8, facecolor=pool_color[bm], edgecolor="none"))
        match = gw == bm
        ax1.text(i + 0.45, 2.05, "✓" if match else "✗", ha="center", va="bottom",
                  fontsize=11, color=(INK_SECONDARY if match else "#c0392b"), fontweight="bold")
        ax1.text(i + 0.45, -0.35, str(N), ha="center", va="top", fontsize=9, color=INK_SECONDARY)
    ax1.set_xlim(-0.3, len(Ns))
    ax1.set_ylim(-0.9, 2.4)
    ax1.set_yticks([0.4, 1.4])
    ax1.set_yticklabels(["Blind Selection\n(blind, modal)", "Grid\n(hindsight, winner)"], fontsize=9.5)
    ax1.set_xticks([])
    ax1.set_title("Pool choice per panel size N: grid hindsight winner vs. Blind Selection's blind modal pick",
                   fontsize=11.5, color=INK, loc="left", pad=14)
    for spine in ax1.spines.values():
        spine.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c) for c in [BLUE, GREEN, MAGENTA]]
    ax1.legend(handles, ["stat", "FST", "fst_stat"], loc="upper center",
               bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=9.5)

    ax2 = fig.add_subplot(gs[1])
    grid_tally = pd.Series(grid_winner_pool).value_counts()
    bs_tally = folds["panel"].value_counts()
    pools = ["stat", "FST", "fst_stat"]
    grid_pct = [grid_tally.get(p, 0) / len(Ns) * 100 for p in pools]
    bs_pct = [bs_tally.get(p, 0) / len(folds) * 100 for p in pools]

    x = range(len(pools))
    w = 0.32
    bars1 = ax2.bar([xi - w / 2 - 0.02 for xi in x], grid_pct, width=w,
                     color=[pool_color[p] for p in pools])
    bars2 = ax2.bar([xi + w / 2 + 0.02 for xi in x], bs_pct, width=w,
                     color=[pool_color[p] for p in pools], alpha=0.45,
                     hatch="////", edgecolor=SURFACE, linewidth=1)
    for b in bars1:
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, f"{b.get_height():.0f}%",
                  ha="center", fontsize=10, color=INK, fontweight="bold")
    for b in bars2:
        ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, f"{b.get_height():.0f}%",
                  ha="center", fontsize=10, color=INK_SECONDARY)
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(pools, fontsize=11)
    ax2.set_ylabel("Share of choices (%)")
    ax2.set_ylim(0, 75)
    ax2.grid(True, axis="y", color=GRID_COLOR, linewidth=0.7, zorder=0)
    ax2.set_axisbelow(True)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    ax2.set_title("Overall pool preference: hindsight (grid) vs. blind (Blind Selection)",
                   fontsize=11.5, color=INK, loc="left")
    from matplotlib.patches import Patch
    leg_handles = [
        Patch(facecolor=INK_MUTED, label="Grid — % of 18 N sizes won"),
        Patch(facecolor=INK_MUTED, alpha=0.45, hatch="////", edgecolor=SURFACE,
              label="Blind Selection — % of 90 blind picks"),
    ]
    ax2.legend(handles=leg_handles, loc="upper right", frameon=False, fontsize=9.5)

    fig.tight_layout()
    fig.savefig(OUT / "fig5_pool_divergence.png", bbox_inches="tight")
    fig.savefig(OUT / "fig5_pool_divergence.svg", bbox_inches="tight")
    plt.close(fig)

    n_match = sum(1 for N in Ns if grid_winner_pool[N] == bs_modal_pool[N])
    print(f"Pool per-N agreement: {n_match}/{len(Ns)}")
    print("Grid pool tally:", dict(grid_tally))
    print("Blind Selection pool tally:", dict(bs_tally))


if __name__ == "__main__":
    fig3_winner_consistency()
    fig4_pool_crossover()
    fig5_pool_divergence()
    fig6_accuracy_vs_n()
    print(f"Wrote figures to {OUT}")
