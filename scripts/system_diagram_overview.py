"""Overall system diagram: 1000 Genomes -> Hard Filter -> Panel Finder ->
committed panels -> comparison to published panels.

Flat design, pastel palette, matching the report's established figure
conventions (same box style as bs_vs_bfc_flowchart.py and report_figures.py).
All connectors are computed from actual box edges (never eyeballed), so nothing
goes visually disjointed if a box size changes. High-DPI PNG for print/slide use.

Writes: reports/figures/fig1_system_overview.png
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

fig, ax = plt.subplots(figsize=(16.4, 6.8), dpi=DPI)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ax.set_xlim(0, 16.8)
ax.set_ylim(0, 7.0)
ax.axis("off")
ax.set_aspect("equal")


def block(x, y, w, h, title, accent, sub=None, fontsize=14):
    """A big rounded-rectangle column/stage container with a bold title.
    Returns (cx, top_y, bottom_y, left_x, right_x) for downstream connectors."""
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.28",
        linewidth=2.2, edgecolor=accent, facecolor=BLOCK_TINT, zorder=2,
    )
    ax.add_patch(b)
    ax.text(x + w / 2, y + h - 0.4, title, ha="center", va="center",
             fontsize=fontsize, color=INK, weight="bold", zorder=4)
    if sub:
        ax.text(x + w / 2, y + h - 0.74, sub, ha="center", va="center",
                 fontsize=9.5, color=INK_SECONDARY, style="italic", zorder=4)
    return dict(cx=x + w / 2, cy=y + h / 2, top=y + h, bot=y, left=x, right=x + w)


def chip(x, y, w, h, text, accent=INK_MUTED, fontsize=10.5, fill="white"):
    """A small rounded-rectangle component box inside a block."""
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
    """A data-artifact node (input/output): solid fill, sharp accent."""
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
    """Horizontal connector between the right edge of p1 and the left edge of
    p2, at height y (defaults to the average of both centers)."""
    yy = y if y is not None else (p1["cy"] + p2["cy"]) / 2
    ax.add_patch(FancyArrowPatch(
        (p1["right"], yy), (p2["left"], yy), arrowstyle="-|>", mutation_scale=17,
        linewidth=1.8, color=color, zorder=2, shrinkA=0, shrinkB=0,
    ))
    if label:
        ax.text((p1["right"] + p2["left"]) / 2, yy + label_dy, label,
                 ha="center", va="bottom", fontsize=label_fs, color=INK_SECONDARY, zorder=4)


def vconnect(p1, p2, x=None, color=INK_MUTED):
    """Vertical connector from the bottom edge of p1 down to the top edge of p2."""
    xx = x if x is not None else p1["cx"]
    ax.add_patch(FancyArrowPatch(
        (xx, p1["bot"]), (xx, p2["top"]), arrowstyle="-|>", mutation_scale=12,
        linewidth=1.4, color=color, zorder=2, shrinkA=0, shrinkB=0,
    ))


# ---------------------------------------------------------------------------
# Stage 1: Input node
# ---------------------------------------------------------------------------
p_input = node(0.3, 2.7, 2.05, 1.6, "1000 Genomes\nPhase 3\n\n504 EAS samples",
               fontsize=11, fill="white", accent=INK)

# ---------------------------------------------------------------------------
# Stage 2: Hard Filter block
# ---------------------------------------------------------------------------
hf_x = 2.85
p_hf = block(hf_x, 0.9, 2.75, 4.9, "Hard Filter", GOLD, sub="label-agnostic QC")
hf_items = ["SNP-only, biallelic", "MAF ≥ 0.001", "Call rate ≥ 95%",
            "HWE p ≥ 1×10⁻⁶", "LD pruning"]
hf_chip_h, hf_gap = 0.62, 0.18
hf_top = p_hf["top"] - 0.98
for i, item in enumerate(hf_items):
    cy = hf_top - i * (hf_chip_h + hf_gap)
    chip(hf_x + 0.28, cy - hf_chip_h, 2.75 - 0.56, hf_chip_h, item, accent=GOLD, fontsize=10)

connect(p_input, p_hf, y=(p_input["cy"] + p_hf["cy"]) / 2)

# ---------------------------------------------------------------------------
# Stage 3: Panel Finder block (collapsed -- see detail diagram)
# ---------------------------------------------------------------------------
pf_x = 6.35
p_pf = block(pf_x, 0.5, 4.4, 5.9, "Panel Finder", BLUE, sub="leak-free, rebuilt per CV fold")

pf_inner_x = pf_x + 0.3
pf_inner_w = 4.4 - 0.6
c1 = chip(pf_inner_x, 4.3, pf_inner_w, 0.85,
          "Candidate Pools\nstat · FST · fst_stat", accent=BLUE, fontsize=10.5)
c2 = chip(pf_inner_x, 3.05, pf_inner_w, 0.85,
          "Reducers\nL1-LR · RF · EN", accent=BLUE, fontsize=10.5)
c3 = chip(pf_inner_x, 1.75, pf_inner_w, 0.85,
          "Evaluation\nRF · XGB · LR · SVM-RBF · SVM-Lin · GBM",
          accent=BLUE, fontsize=9.5)
c4 = chip(pf_inner_x, 0.75, pf_inner_w, 0.6,
          "Blind Selection + Best Fixed Configuration", accent=INK_SECONDARY, fontsize=9)

vconnect(c1, c2)
ax.text(pf_x + 4.4 / 2 + 1.35, (c1["bot"] + c2["top"]) / 2, "top-N\n(35,50,70)",
        ha="left", va="center", fontsize=7.3, color=INK_SECONDARY, style="italic")
vconnect(c2, c3)
vconnect(c3, c4)

connect(p_hf, p_pf, y=(p_hf["cy"] + p_pf["cy"]) / 2, label="614,759 SNPs")

# ---------------------------------------------------------------------------
# Stage 4: Output node
# ---------------------------------------------------------------------------
out_x = 11.55
p_out = node(out_x, 2.75, 2.15, 1.6, "Committed\npanels\n\nN = 35, 50, 70", fontsize=11,
             fill="#e7effa", accent=BLUE)

connect(p_pf, p_out, y=(p_pf["cy"] + p_out["cy"]) / 2)

# ---------------------------------------------------------------------------
# Stage 5: Comparison block
# ---------------------------------------------------------------------------
cmp_x = 14.5
p_cmp = block(cmp_x, 1.15, 2.0, 4.4, "Comparison", MAGENTA, sub="to published panels", fontsize=12)
cmp_items = ["Cai-34", "Cao-19", "Shi-142"]
cmp_chip_h, cmp_gap = 0.55, 0.2
cmp_top = p_cmp["top"] - 0.98
for i, item in enumerate(cmp_items):
    cy = cmp_top - i * (cmp_chip_h + cmp_gap)
    chip(cmp_x + 0.15, cy - cmp_chip_h, 2.0 - 0.3, cmp_chip_h, item, accent=MAGENTA, fontsize=10)
last_top = cmp_top - len(cmp_items) * (cmp_chip_h + cmp_gap) + cmp_gap
chip(cmp_x + 0.15, last_top - 0.66, 2.0 - 0.3, 0.66,
     "same 5-fold\nCV protocol", accent=INK_MUTED, fontsize=8.4)

connect(p_out, p_cmp, y=(p_out["cy"] + p_cmp["cy"]) / 2)

fig.suptitle("Overall System Diagram", fontsize=16, color=INK, weight="bold", y=0.97)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUT / "fig1_system_overview.png", bbox_inches="tight", facecolor=SURFACE)
fig.savefig(OUT / "fig1_system_overview.svg", bbox_inches="tight", facecolor=SURFACE)
print(f"Wrote {OUT / 'fig1_system_overview.png'} at {DPI} dpi")
