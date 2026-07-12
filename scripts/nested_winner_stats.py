"""Winner analysis over the full leak-free grid (grid_results.csv).

Answers: across N=5..100 and 5 CV folds, which pool / reductor / classifier
wins most consistently? Reports fold-level winners, per-N tied-best sets
(within 1 SE), and config/component leaderboards. Writes CSVs, a figure, and a
markdown fragment for the master report.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from notebook_init import setup


def main():
    _, PATHS, *_ = setup()
    out = Path(str(PATHS.outputs_dir("self_evaluation/08b_nested_cv_sweep")))
    g = pd.read_csv(out / "grid_results.csv")
    g["config"] = g["panel"] + "+" + g["reductor"] + "+" + g["classifier"]

    # ---- (1) per-fold winners: argmax config per (N, fold) ----
    idx = g.groupby(["n_snps", "fold"])["acc"].idxmax()
    winners = g.loc[idx, ["n_snps", "fold", "panel", "reductor", "classifier", "config", "acc"]]
    winners.to_csv(out / "winner_per_fold.csv", index=False)
    n_cells = len(winners)

    def rate(col):
        vc = winners[col].value_counts()
        return (vc / n_cells * 100).round(1)

    # ---- (2) per-N mean + tied-best (within 1 SE of the best config) ----
    agg = (g.groupby(["n_snps", "config", "panel", "reductor", "classifier"])
             .agg(acc_mean=("acc", "mean"), acc_se=("acc", lambda x: x.std() / np.sqrt(len(x))))
             .reset_index())
    per_n_best, tied_counter = [], {}
    for N, sub in agg.groupby("n_snps"):
        b = sub.loc[sub["acc_mean"].idxmax()]
        thr = b["acc_mean"] - b["acc_se"]
        tied = sub[sub["acc_mean"] >= thr]
        for c in tied["config"]:
            tied_counter[c] = tied_counter.get(c, 0) + 1
        per_n_best.append({
            "n_snps": N, "best_config": b["config"],
            "best_acc": round(b["acc_mean"] * 100, 2),
            "n_tied_within_1se": len(tied),
            "tied_configs": "; ".join(sorted(tied["config"])),
        })
    per_n_best_df = pd.DataFrame(per_n_best)
    per_n_best_df.to_csv(out / "winner_per_N.csv", index=False)

    # ---- (3) leaderboards ----
    n_N = g["n_snps"].nunique()
    top1 = winners["config"].value_counts()
    config_board = (agg.groupby(["config", "panel", "reductor", "classifier"])["acc_mean"]
                       .mean().reset_index().rename(columns={"acc_mean": "mean_acc_over_N"}))
    config_board["mean_acc_over_N"] = (config_board["mean_acc_over_N"] * 100).round(2)
    config_board["top1_folds"] = config_board["config"].map(top1).fillna(0).astype(int)
    config_board["top1_pct"] = (config_board["top1_folds"] / n_cells * 100).round(1)
    config_board["tied_best_N"] = config_board["config"].map(tied_counter).fillna(0).astype(int)
    config_board["tied_best_pct"] = (config_board["tied_best_N"] / n_N * 100).round(1)
    config_board = config_board.sort_values(
        ["top1_folds", "mean_acc_over_N"], ascending=False).reset_index(drop=True)
    config_board.to_csv(out / "winner_config_leaderboard.csv", index=False)

    # component-level mean acc (marginal)
    comp = {}
    for col in ["panel", "reductor", "classifier"]:
        comp[col] = pd.DataFrame({
            "mean_acc": (g.groupby(col)["acc"].mean() * 100).round(2),
            "fold_win_rate_pct": rate(col),
        }).sort_values("mean_acc", ascending=False)

    # ---- (4) figure ----
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, col, title in zip(axes, ["panel", "reductor", "classifier"],
                              ["Pool", "Reductor", "Classifier"]):
        d = comp[col]
        ax.barh(d.index[::-1], d["fold_win_rate_pct"][::-1], color="#377eb8", alpha=0.85)
        for i, (name, row) in enumerate(d[::-1].iterrows()):
            ax.text(row["fold_win_rate_pct"] + 0.5, i,
                    f"{row['fold_win_rate_pct']:.0f}%  ({row['mean_acc']:.1f})", va="center", fontsize=8)
        ax.set_title(f"{title}: fold-win rate\n(label = win% (mean acc))")
        ax.set_xlabel(f"% of {n_cells} (N,fold) cells won")
        ax.grid(alpha=0.3, axis="x")
    fig.suptitle("Most consistent winner across N=5..100 x 5 folds (leak-free grid)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "winner_consistency.png", dpi=150, bbox_inches="tight")

    # ---- (5) markdown fragment ----
    lines = ["## Winner consistency (leak-free grid, N=5..100 x 5 folds)", "",
             f"Winners tallied over {n_cells} (N, fold) cells; mean acc = mean over all folds & N.", ""]
    for col, label in [("panel", "Pool"), ("reductor", "Reductor"), ("classifier", "Classifier")]:
        lines.append(f"**{label}** — fold-win rate (mean acc):")
        for name, row in comp[col].iterrows():
            lines.append(f"- {name}: {row['fold_win_rate_pct']:.0f}%  ({row['mean_acc']:.1f}%)")
        lines.append("")
    lines.append("**Top configurations** (by fold-win count):")
    lines.append("")
    lines.append("| config | mean acc over N | top-1 folds | within-1SE of best (of N) |")
    lines.append("|--------|-----------------|-------------|----------------------------|")
    for _, r in config_board.head(8).iterrows():
        lines.append(f"| {r['config']} | {r['mean_acc_over_N']}% | "
                     f"{r['top1_folds']}/{n_cells} ({r['top1_pct']}%) | {r['tied_best_N']}/{n_N} ({r['tied_best_pct']}%) |")
    (out / "winner_stats_fragment.md").write_text("\n".join(lines))

    print("Saved winner analysis to:", out)
    for col in ["panel", "reductor", "classifier"]:
        print(f"\n== {col} ==\n{comp[col].to_string()}")
    print("\n== top configs ==")
    print(config_board.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
