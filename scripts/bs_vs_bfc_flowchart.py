"""Infographic-style explainer: Blind Selection vs Best Fixed Configuration.

Flat design throughout (no gradients/shadows) with small hand-drawn icons
(fold badges, an eye for "can it see the test answers yet", a mini
scoreboard for the 54-config comparison) to make the timing distinction —
the entire point of the diagram — visually concrete, not just textual.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
SURFACE = "#fcfcfb"
BLUE = "#2a78d6"
GREEN = "#008300"
BLUE_TINT = "#eef4fc"
GREEN_TINT = "#eaf5ea"
GOLD = "#c98500"

fig, ax = plt.subplots(figsize=(12.5, 13.5), dpi=200)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)
ax.set_xlim(0, 10.6)
ax.set_ylim(0, 16.2)
ax.axis("off")
ax.set_aspect("equal")


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

def box(x, y, w, h, text, accent=None, fontsize=10, weight="normal", fill="white"):
    b = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.1",
        linewidth=1.4, edgecolor=accent if accent else INK_MUTED,
        facecolor=fill, zorder=3,
    )
    ax.add_patch(b)
    ax.text(
        x + w / 2, y + h / 2, text, ha="center", va="center",
        fontsize=fontsize, color=INK, weight=weight, linespacing=1.35, zorder=4,
    )


def arrow(x, y1, y2, color=INK_MUTED):
    ax.add_patch(FancyArrowPatch(
        (x, y1), (x, y2), arrowstyle="-|>", mutation_scale=13,
        linewidth=1.3, color=color, zorder=2,
    ))


def fold_badges(x_center, y, n=5, filled_idx=None, r=0.16, gap=0.42):
    """A row of n small numbered circles representing the 5 CV folds."""
    total_w = (n - 1) * gap
    x0 = x_center - total_w / 2
    for i in range(n):
        cx = x0 + i * gap
        is_on = filled_idx is not None and i in filled_idx
        c = Circle((cx, y), r, linewidth=1.3,
                   edgecolor=INK_SECONDARY,
                   facecolor=INK_SECONDARY if is_on else "white", zorder=5)
        ax.add_patch(c)
        ax.text(cx, y, str(i + 1), ha="center", va="center", fontsize=7.5,
                color="white" if is_on else INK_SECONDARY, weight="bold", zorder=6)


def eye_icon(cx, cy, r=0.24, open_=True, color=INK):
    e = Ellipse((cx, cy), r * 2.1, r * 1.25, linewidth=1.5,
                edgecolor=color, facecolor="white", zorder=5)
    ax.add_patch(e)
    if open_:
        ax.add_patch(Circle((cx, cy), r * 0.34, facecolor=color, edgecolor="none", zorder=6))
    else:
        ax.plot([cx - r * 1.05, cx + r * 1.05], [cy - r * 0.62, cy + r * 0.62],
                color=color, linewidth=1.8, zorder=7, solid_capstyle="round")


def scoreboard_icon(x, y, w, h, n=12, winner_idx=9, color=INK_MUTED, winner_color=GOLD):
    """A tiny bar-chart glyph: many short bars, one taller + starred = the pick."""
    heights = [0.35, 0.55, 0.42, 0.68, 0.3, 0.5, 0.6, 0.4, 0.72, 1.0, 0.48, 0.58]
    heights = heights[:n]
    bw = w / n * 0.68
    gap = w / n
    for i, hf in enumerate(heights):
        bx = x + i * gap
        bh = h * hf
        is_win = i == winner_idx
        ax.add_patch(Rectangle((bx, y), bw, bh,
                                facecolor=winner_color if is_win else color,
                                edgecolor="none", zorder=5))
    if winner_idx is not None:
        wx = x + winner_idx * gap + bw / 2
        ax.plot(wx, y + h * heights[winner_idx] + 0.14, marker="*",
                markersize=11, color=winner_color, zorder=6,
                markeredgecolor="none")


def prop_bar(x, y, w, h, frac, color, tint, label):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                 linewidth=1.2, edgecolor=INK_MUTED, facecolor=tint, zorder=2))
    ax.add_patch(FancyBboxPatch((x, y), w * frac, h, boxstyle="round,pad=0,rounding_size=0.05",
                                 linewidth=0, facecolor=color, zorder=3))
    ax.text(x + w + 0.15, y + h / 2, label, ha="left", va="center",
            fontsize=11, color=INK, weight="bold", zorder=4)


# ---------------------------------------------------------------------------
# Column headers
# ---------------------------------------------------------------------------

box(0.4, 15.2, 4.6, 0.8, "BLIND SELECTION", accent=GREEN, fill=GREEN_TINT,
    fontsize=13, weight="bold")
box(5.6, 15.2, 4.6, 0.8, "BEST FIXED CONFIGURATION", accent=BLUE, fill=BLUE_TINT,
    fontsize=12.5, weight="bold")

# ---------------------------------------------------------------------------
# Left column: Blind Selection
# ---------------------------------------------------------------------------
lx, lw = 0.4, 4.6
lcx = lx + lw / 2

ax.text(lcx, 14.75, "same 5 folds, used one at a time", ha="center", fontsize=8.7,
        color=INK_SECONDARY, style="italic")
fold_badges(lcx, 14.35, filled_idx=[0])
arrow(lcx, 14.1, 13.95)

box(lx, 13.0, lw, 0.85, "For fold #1: set its test slice aside", fontsize=10)
arrow(lcx, 13.0, 12.9)

eye_icon(lx + 0.35, 12.15, r=0.22, open_=False, color=GREEN)
box(lx + 0.75, 11.55, lw - 0.75, 1.2,
    "Using ONLY the training data,\nguess the best of 54 configs\n(this fold's test slice: not seen)",
    fontsize=9.3)
arrow(lcx, 11.55, 11.45)

box(lx, 10.2, lw, 1.0, "Run the guessed config on the\nheld-out test slice → record its score", fontsize=9.5)
arrow(lcx, 10.2, 10.1)

fold_badges(lcx, 9.65, filled_idx=[0, 1, 2, 3, 4])
box(lx, 8.75, lw, 0.75, "Repeat for folds 2–5, average the 5 scores", fontsize=9.5)
arrow(lcx, 8.75, 8.65)

box(lx, 7.55, lw, 0.95, "Blind Selection = 90.67%\n(at N = 70)", accent=GREEN, fill=GREEN_TINT,
    fontsize=12, weight="bold")

box(lx, 6.35, lw, 1.05,
    "Every guess was made blind —\nno fold's config was chosen after\nseeing that fold's test answers", fontsize=9.2)

# ---------------------------------------------------------------------------
# Right column: Best Fixed Configuration
# ---------------------------------------------------------------------------
rx, rw = 5.6, 4.6
rcx = rx + rw / 2

ax.text(rcx, 14.75, "same 5 folds, all used already", ha="center", fontsize=8.7,
        color=INK_SECONDARY, style="italic")
fold_badges(rcx, 14.35, filled_idx=[0, 1, 2, 3, 4])
arrow(rcx, 14.1, 13.95)

box(rx, 13.0, rw, 0.85, "Run all 54 configs on all 5 folds", fontsize=10)
arrow(rcx, 13.0, 12.9)

eye_icon(rx + 0.35, 12.2, r=0.22, open_=True, color=BLUE)
box(rx + 0.75, 11.65, rw - 0.75, 1.05,
    "Every config's average score is\nnow known — test answers already\nseen, for all 54 of them",
    fontsize=9.3)
arrow(rcx, 11.65, 11.55)

box(rx, 10.2, rw, 1.35, "", fill="white")
scoreboard_icon(rx + 0.55, 10.85, rw - 1.1, 0.45, winner_idx=9)
ax.text(rcx, 10.45, "Look at the scoreboard — pick\nwhichever config has the best average",
        ha="center", va="center", fontsize=9.3, color=INK, linespacing=1.35, zorder=4)
arrow(rcx, 10.2, 10.1)

box(rx, 8.75, rw, 0.75, "That single config is now “the winner”", fontsize=9.5)
arrow(rcx, 8.75, 8.65)

box(rx, 7.55, rw, 0.95, "Best Fixed Configuration = 94.25%\n(FST+EN+SVM-RBF, N = 70)",
    accent=BLUE, fill=BLUE_TINT, fontsize=12, weight="bold")

box(rx, 6.35, rw, 1.05,
    "The winner was chosen AFTER\nseeing every config's test score —\nflattered by that peek", fontsize=9.2)

# ---------------------------------------------------------------------------
# Shared takeaway + proportional comparison bar
# ---------------------------------------------------------------------------
box(0.4, 4.55, 9.8, 1.35,
    "Same 5 folds. Same 54 candidate configs. The only difference is WHEN the winning\n"
    "config gets picked relative to seeing the test scores — before (blindfolded, left) or\n"
    "after (eyes open, right). That timing difference is the entire gap below.",
    fontsize=10.4, fill="#f5f4f0")

ax.text(0.4, 4.05, "90.67%  vs  94.25%", fontsize=11.5, color=INK, weight="bold")
prop_bar(0.4, 3.15, 7.6, 0.55, 90.67 / 100, GREEN, GREEN_TINT, "Blind Selection — 90.67%")
prop_bar(0.4, 2.35, 7.6, 0.55, 94.25 / 100, BLUE, BLUE_TINT, "Best Fixed Configuration — 94.25%")

ax.text(5.3, 1.35,
        "Both numbers come from the exact same underlying grid — the gap is entirely\n"
        "about when the peek happens, not about different data or different models.",
        ha="center", fontsize=9, color=INK_SECONDARY, style="italic")

fig.tight_layout()
OUT = Path(__file__).resolve().parent.parent / "reports/figures"
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "bs_vs_bfc_flowchart.png", bbox_inches="tight")
print(f"Wrote {OUT / 'bs_vs_bfc_flowchart.png'}")
