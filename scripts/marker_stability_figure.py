"""Marker-selection-stability figure (Reviewer 2, minor comment 4): how often
each committed panel's markers reappear across this pipeline's own 5 leak-free
outer folds, versus the committed (full-504-sample) panel used for deployment.

Reads:  reports/figures/committed_panel_markers.csv (fold_selection_freq column,
        built by scripts/committed_panel_marker_table.py)
Writes: reports/figures/fig8_marker_stability.png / .svg
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CSV = REPO / "reports/figures/committed_panel_markers.csv"
OUT = REPO / "reports/figures"

INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID_COLOR = "#e1e0d9"
SURFACE = "#fcfcfb"

# Sequential single-hue ramp (blue), light -> dark, for the ordinal N=35/50/70 series.
N_COLORS = {35: "#a8c8ec", 50: "#5892d6", 70: "#154a86"}

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

df = pd.read_csv(CSV)
buckets = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
Ns = [35, 50, 70]

fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=200)
ax.grid(True, axis="y", color=GRID_COLOR, linewidth=0.7, zorder=0)

width = 0.24
x = np.arange(len(buckets))
for i, N in enumerate(Ns):
    sub = df[df.N == N]
    pct = [(sub.fold_selection_freq == b).mean() * 100 for b in buckets]
    offset = (i - 1) * width
    bars = ax.bar(
        x + offset, pct, width=width * 0.9, color=N_COLORS[N],
        edgecolor=SURFACE, linewidth=0.8, zorder=3,
        label=f"N = {N} (n={len(sub)} markers)",
    )
    for rect, v in zip(bars, pct):
        if v > 0:
            ax.annotate(
                f"{v:.0f}%", (rect.get_x() + rect.get_width() / 2, rect.get_height()),
                ha="center", va="bottom", fontsize=7.5, color=INK_SECONDARY,
                textcoords="offset points", xytext=(0, 2),
            )

ax.set_xticks(x)
ax.set_xticklabels([f"{int(b*5)}/5" for b in buckets])
ax.set_xlabel("Times the marker also appears in the 5 leak-free per-fold panels")
ax.set_ylabel("Share of committed panel's markers (%)")
ax.set_title(
    "Marker-selection stability: committed panel vs. its own 5 outer folds",
    fontsize=11, color=INK, pad=12,
)
ax.legend(frameon=False, loc="upper left", fontsize=9)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
fig.tight_layout()

fig.savefig(OUT / "fig8_marker_stability.png")
fig.savefig(OUT / "fig8_marker_stability.svg")
print(f"Wrote {OUT / 'fig8_marker_stability.png'}")

for N in Ns:
    sub = df[df.N == N]
    print(
        f"N={N}: mean freq={sub.fold_selection_freq.mean():.2f}, "
        f"5/5={  (sub.fold_selection_freq == 1.0).mean()*100:.1f}%, "
        f"0/5={(sub.fold_selection_freq == 0.0).mean()*100:.1f}%"
    )
