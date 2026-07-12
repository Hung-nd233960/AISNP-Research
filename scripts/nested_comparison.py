"""Compare leak-free nested CV against the frozen leaky baseline.

Reads:
  results_archive/baseline_leaky_v1/...          (leaky Stage-2 numbers, published panels)
  outputs/self_evaluation/08b_nested_cv_sweep/   (nested Tier B / Tier A results)

Writes (to the 08b output dir):
  comparison_leaky_vs_nested.csv     N | leaky | nested | Delta
  comparison_accuracy_vs_n.png       leaky vs nested curve + published panels
  comparison_summary.md              human-readable summary
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

REPO = Path(__file__).resolve().parent.parent
BASE = REPO / "results_archive/baseline_leaky_v1/self_evaluation"


def main():
    _, PATHS, *_ = setup()
    out = Path(str(PATHS.outputs_dir("self_evaluation/08b_nested_cv_sweep")))

    leaky = pd.read_csv(BASE / "08_unified_panel_sweep/stage2_best.csv")[
        ["n_snps", "panel", "reductor", "classifier", "acc_mean", "mcc_mean"]
    ].rename(columns={"acc_mean": "leaky_acc", "mcc_mean": "leaky_mcc",
                      "panel": "leaky_panel", "reductor": "leaky_reductor",
                      "classifier": "leaky_clf"})

    tb = pd.read_csv(out / "nested_tierB_results.csv")[
        ["n_snps", "acc_mean", "acc_std", "mcc_mean"]
    ].rename(columns={"acc_mean": "nestedB_acc", "acc_std": "nestedB_std",
                      "mcc_mean": "nestedB_mcc"})

    cmp = leaky.merge(tb, on="n_snps", how="outer")
    cmp["delta_acc"] = cmp["nestedB_acc"] - cmp["leaky_acc"]

    ta_path = out / "nested_tierA_results.csv"
    if ta_path.exists():
        ta = pd.read_csv(ta_path)[["n_snps", "acc_mean", "mcc_mean"]].rename(
            columns={"acc_mean": "nestedA_acc", "mcc_mean": "nestedA_mcc"})
        cmp = cmp.merge(ta, on="n_snps", how="left")

    cmp = cmp.sort_values("n_snps").reset_index(drop=True)
    cmp.to_csv(out / "comparison_leaky_vs_nested.csv", index=False)

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(cmp["n_snps"], cmp["leaky_acc"] * 100, "o-", color="#e41a1c",
            lw=2, ms=6, label="Leaky (pool selected on all 504) — old")
    ax.errorbar(cmp["n_snps"], cmp["nestedB_acc"] * 100,
                yerr=cmp["nestedB_std"] * 100, fmt="s-", color="#377eb8",
                lw=2, ms=6, capsize=3, label="Nested Tier B (in-fold pools) — leak-free")
    if "nestedA_acc" in cmp:
        sub = cmp.dropna(subset=["nestedA_acc"])
        ax.plot(sub["n_snps"], sub["nestedA_acc"] * 100, "D", color="#4daf4a",
                ms=10, label="Nested Tier A (fully nested)")

    # published external panels (already leak-free) as reference markers
    try:
        pub = pd.read_csv(BASE / "11_results/table2_published_panels.csv")
        pub = pub[pub["Panel"].str.startswith(("cai", "cao", "shi"))]
        for _, r in pub.iterrows():
            ax.scatter(r["Matched N"], r["Accuracy"] * 100, marker="*",
                       s=160, color="gray", zorder=5)
            ax.annotate(r["Panel"], (r["Matched N"], r["Accuracy"] * 100),
                        fontsize=7, color="dimgray",
                        xytext=(3, -9), textcoords="offset points")
    except Exception as e:
        print(f"(published panels not plotted: {e})")

    ax.set_xlabel("Panel size N")
    ax.set_ylabel("CV accuracy (%)")
    ax.set_title("PAANDA-EA: leaky vs leak-free nested CV\n"
                 "(gray ★ = external published panels, no selection optimism)")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(out / "comparison_accuracy_vs_n.png", dpi=150, bbox_inches="tight")

    # ---- markdown summary ----
    key = cmp[cmp["n_snps"].isin([35, 50, 70])]
    lines = ["# Leaky vs leak-free nested CV — comparison", "",
             "Accuracy is 5-fold CV (seed 42). Δ = nested − leaky (negative = leakage removed).",
             "", "| N | config | leaky acc | nested Tier B | Δ | nested Tier A |",
             "|---|--------|-----------|---------------|---|---------------|"]
    for _, r in cmp.iterrows():
        ta_txt = (f"{r['nestedA_acc']*100:.2f}%" if "nestedA_acc" in cmp
                  and pd.notna(r.get("nestedA_acc")) else "—")
        lines.append(
            f"| {int(r['n_snps'])} | {r['leaky_panel']}+{r['leaky_reductor']}+{r['leaky_clf']} "
            f"| {r['leaky_acc']*100:.2f}% | {r['nestedB_acc']*100:.2f}% ± {r['nestedB_std']*100:.2f} "
            f"| {r['delta_acc']*100:+.2f} | {ta_txt} |")
    mean_drop = cmp["delta_acc"].mean() * 100
    lines += ["",
              f"**Mean Δ across N:** {mean_drop:+.2f} points.",
              f"**Committed panels (N=35/50/70) mean Δ:** {key['delta_acc'].mean()*100:+.2f} points.",
              "",
              "External panels carry no selection optimism, so the nested numbers are the "
              "fair comparison to Cai/Cao/Shi. See `comparison_accuracy_vs_n.png`."]
    (out / "comparison_summary.md").write_text("\n".join(lines))

    print("Saved comparison to:", out)
    print(cmp[["n_snps", "leaky_acc", "nestedB_acc", "delta_acc"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
