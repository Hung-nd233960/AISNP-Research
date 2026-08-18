"""Diagram describing Cai et al. (2024)'s own pipeline -- for comparison
against this project's system diagrams (fig1_system_overview,
fig2_panelfinder_detail). Same flat-design conventions and palette, limited
strictly to what the paper itself describes: a single embedded
Random-Forest selection pass, not a leak-free nested-CV search.

Writes: reports/figures/fig_cai_pipeline.png (+ .svg)
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
SURFACE = "#fdfcf9"
BLUE = "#6fa7e8"
GREEN = "#7fbf7f"
MAGENTA = "#e6a7c3"
GOLD = "#dfb978"
BLOCK_TINT = "#f7f5f1"

OUT = Path(__file__).resolve().parent.parent / "reports" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

DPI = 320

fig, ax = plt.subplots(figsize=(14.2, 6.2), dpi=DPI)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ax.set_xlim(0, 14.5)
ax.set_ylim(0, 6.4)
ax.axis("off")
ax.set_aspect("equal")


def block(x, y, w, h, title, accent, fontsize=14):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.28",
        linewidth=2.2, edgecolor=accent, facecolor=BLOCK_TINT, zorder=2,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h - 0.4, title, ha="center", va="center",
             fontsize=fontsize, color=INK, weight="bold", zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def chip(x, y, w, h, text, accent=INK_MUTED, fontsize=10.5, fill="white"):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.015,rounding_size=0.13",
        linewidth=1.4, edgecolor=accent, facecolor=fill, zorder=3,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=fontsize, color=INK, linespacing=1.35, zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def node(x, y, w, h, text, fontsize=11, fill="white", accent=INK):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.26",
        linewidth=2.4, edgecolor=accent, facecolor=fill, zorder=3,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=fontsize, color=INK, weight="bold", linespacing=1.4, zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def connect(p1, p2, y=None, label=None, color=INK_SECONDARY, label_dy=0.2, label_fs=9):
    yy = y if y is not None else (p1["cy"] + p2["cy"]) / 2
    ax.add_patch(FancyArrowPatch(
        (p1["right"], yy), (p2["left"], yy), arrowstyle="-|>", mutation_scale=17,
        linewidth=1.8, color=color, zorder=2, shrinkA=0, shrinkB=0,
    ))
    if label:
        ax.text((p1["right"] + p2["left"]) / 2, yy + label_dy, label,
                 ha="center", va="bottom", fontsize=label_fs, color=INK_SECONDARY, zorder=4)


# ---------------------------------------------------------------------------
# Stage 1: Input node
# ---------------------------------------------------------------------------
p_input = node(0.3, 2.2, 2.15, 1.6, "1000 Genomes\nPhase 3\n\n26 populations,\ngenome-wide",
               fontsize=10.5, fill="white", accent=INK)

# ---------------------------------------------------------------------------
# Stage 2: Feature Selection block
# ---------------------------------------------------------------------------
fs_x = 2.95
p_fs = block(fs_x, 1.0, 2.9, 4.2, "Feature Selection", GOLD)
chip(fs_x + 0.3, 3.55, 2.3, 0.85, "Random Forest\nclassifier", accent=GOLD, fontsize=10.5)
chip(fs_x + 0.3, 2.5, 2.3, 0.85, "One-vs-rest\nscreen", accent=GOLD, fontsize=10.5)
chip(fs_x + 0.3, 1.3, 2.3, 0.9, "Embedded selection,\ngenome-wide", accent=GOLD, fontsize=10)

connect(p_input, p_fs, y=(p_input["cy"] + p_fs["cy"]) / 2)

# ---------------------------------------------------------------------------
# Stage 3: 58-SNP Table block -- panel sizes + Cai-34 highlighted
# ---------------------------------------------------------------------------
out_x = 6.35
p_out = block(out_x, 1.0, 3.2, 4.2, "58-SNP Table", BLUE)
chip(out_x + 0.3, 3.55, 2.6, 0.85, "N of six sizes,\n50 to 2,000 SNPs", accent=BLUE, fontsize=10.5)
chip(out_x + 0.3, 1.3, 2.6, 2.0, "34 within-EAS\n(\"Cai-34\")", accent=MAGENTA, fontsize=12)

connect(p_fs, p_out, y=(p_fs["cy"] + p_out["cy"]) / 2)

# ---------------------------------------------------------------------------
# Stage 4: Evaluation block
# ---------------------------------------------------------------------------
ev_x = 10.15
p_ev = block(ev_x, 1.3, 3.2, 3.6, "Evaluation", MAGENTA, fontsize=13)
chip(ev_x + 0.3, 3.15, 2.6, 0.75, "XGBoost\nclassifier", accent=MAGENTA, fontsize=10.5)
chip(ev_x + 0.3, 2.2, 2.6, 0.75, "0.94 continental\n(5-population)", accent=MAGENTA, fontsize=10)
chip(ev_x + 0.3, 1.5, 2.6, 0.55, "0.92 intra-EAS", accent=MAGENTA, fontsize=10)

connect(p_out, p_ev, y=(p_out["cy"] + p_ev["cy"]) / 2)

fig.suptitle("Cai et al. (2024) — Pipeline", fontsize=16, color=INK, weight="bold", y=0.97)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT / "fig_cai_pipeline.png", bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT / "fig_cai_pipeline.svg", bbox_inches="tight", facecolor=SURFACE)
print(f"Wrote {OUT / 'fig_cai_pipeline.png'} at {DPI} dpi")
