"""Detail diagram for the "Panel Finder" block: candidate pools -> reducers
-> (top-N) -> committed panels -> evaluation.

Flat design, pastel palette, same conventions as system_diagram_overview.py.
Connectors are computed from actual box edges.

Writes: reports/figures/fig2_panelfinder_detail.png
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

fig, ax = plt.subplots(figsize=(16, 7.6), dpi=DPI)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ax.set_xlim(0, 16.2)
ax.set_ylim(0, 7.8)
ax.axis("off")
ax.set_aspect("equal")


def block(x, y, w, h, title, accent, sub=None, fontsize=14.5):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.28",
        linewidth=2.2, edgecolor=accent, facecolor=BLOCK_TINT, zorder=2,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h - 0.42, title, ha="center", va="center",
             fontsize=fontsize, color=INK, weight="bold", zorder=4)
    if sub:
        ax.text(x + w / 2, y + h - 0.78, sub, ha="center", va="center",
                 fontsize=9.5, color=INK_SECONDARY, style="italic", zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def chip(x, y, w, h, text, accent=INK_MUTED, fontsize=11, fill="white"):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.015,rounding_size=0.15",
        linewidth=1.5, edgecolor=accent, facecolor=fill, zorder=3,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=fontsize, color=INK, linespacing=1.35, zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def node(x, y, w, h, text, fontsize=11.5, fill="white", accent=INK):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.26",
        linewidth=2.4, edgecolor=accent, facecolor=fill, zorder=3,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
             fontsize=fontsize, color=INK, weight="bold", linespacing=1.4, zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def connect(p1, p2, y=None, label=None, color=INK_SECONDARY, label_dy=0.22, label_fs=9.5):
    yy = y if y is not None else (p1["cy"] + p2["cy"]) / 2
    ax.add_patch(FancyArrowPatch(
        (p1["right"], yy), (p2["left"], yy), arrowstyle="-|>", mutation_scale=18,
        linewidth=1.9, color=color, zorder=2, shrinkA=0, shrinkB=0,
    ))
    if label:
        ax.text((p1["right"] + p2["left"]) / 2, yy + label_dy, label,
                 ha="center", va="bottom", fontsize=label_fs, color=INK_SECONDARY, zorder=4)


# ---------------------------------------------------------------------------
# Input node
# ---------------------------------------------------------------------------
p_in = node(0.25, 3.35, 2.05, 1.5, "614,759\nSNPs\n\n(Hard Filter\noutput)", fontsize=10.5,
            fill="white", accent=INK)

# ---------------------------------------------------------------------------
# Column 1: Candidate Pools
# ---------------------------------------------------------------------------
c1_x = 2.85
p_c1 = block(c1_x, 0.9, 3.1, 6.0, "Candidate Pools", BLUE, sub="population-differentiation criteria")
inner_x = c1_x + 0.32
inner_w = 3.1 - 0.64
chip(inner_x, 4.75, inner_w, 1.15, "stat\n\nχ² · JSD · AFD\n1,005 SNPs", accent=BLUE, fontsize=10.5)
chip(inner_x, 3.35, inner_w, 1.15, "FST\n\nHudson $F_{ST}$\n2,508 SNPs", accent=GREEN, fontsize=10.5)
chip(inner_x, 1.35, inner_w, 1.75, "fst_stat\n\nstat cascaded\nonto FST\n965 SNPs", accent=MAGENTA, fontsize=10.5)

connect(p_in, p_c1, y=(p_in["cy"] + p_c1["cy"]) / 2)

# ---------------------------------------------------------------------------
# Column 2: Reducers
# ---------------------------------------------------------------------------
c2_x = 6.75
p_c2 = block(c2_x, 0.9, 2.8, 6.0, "Reducers", GOLD, sub="rank SNPs within pool")
inner_x2 = c2_x + 0.28
inner_w2 = 2.8 - 0.56
chip(inner_x2, 4.75, inner_w2, 1.15, "L1-LR", accent=GOLD, fontsize=12)
chip(inner_x2, 3.35, inner_w2, 1.15, "RF\nimportance", accent=GOLD, fontsize=12)
chip(inner_x2, 1.95, inner_w2, 1.15, "EN\nElasticNet", accent=GOLD, fontsize=12)

connect(p_c1, p_c2, y=(p_c1["cy"] + p_c2["cy"]) / 2)

# ---------------------------------------------------------------------------
# Output node: committed panels, fed by a "top-N" labeled arrow
# ---------------------------------------------------------------------------
out_x = 10.55
p_out = node(out_x, 3.4, 2.05, 1.4, "Panels\n\nN = 35, 50, 70", fontsize=11.5,
             fill="#e7effa", accent=BLUE)

connect(p_c2, p_out, y=(p_c2["cy"] + p_out["cy"]) / 2, label="by N = 35, 50, 70", label_fs=9.5)

# ---------------------------------------------------------------------------
# Column 3: Evaluation
# ---------------------------------------------------------------------------
c3_x = 13.5
p_c3 = block(c3_x, 0.7, 2.45, 6.4, "Evaluation", MAGENTA, sub="score every panel", fontsize=13.5)
labels = ["RF", "XGB", "LR", "SVM-RBF", "SVM-Lin", "GBM"]
inner_x3 = c3_x + 0.22
inner_w3 = 2.45 - 0.44
ev_h, ev_gap = 0.68, 0.16
ev_top = p_c3["top"] - 0.96
for i, lab in enumerate(labels):
    cy_top = ev_top - i * (ev_h + ev_gap)
    chip(inner_x3, cy_top - ev_h, inner_w3, ev_h, lab, accent=MAGENTA, fontsize=10.5)

connect(p_out, p_c3, y=(p_out["cy"] + p_c3["cy"]) / 2)

fig.suptitle('"Panel Finder" -- Detail', fontsize=17, color=INK, weight="bold", y=0.975)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT / "fig2_panelfinder_detail.png", bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT / "fig2_panelfinder_detail.svg", bbox_inches="tight", facecolor=SURFACE)
print(f"Wrote {OUT / 'fig2_panelfinder_detail.png'} at {DPI} dpi")
