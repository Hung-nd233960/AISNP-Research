"""Supplementary marker table for the three committed panels (N=35/50/70):
genomic coordinates, reference build, alleles, per-population allele
frequency, max pairwise Hudson FST, and selection frequency across the 5
outer CV folds -- answering reviewer comment #6.

Panel SNP identities are the "shipped" panel: pool + EN reductor fit on
all 504 samples (matching how a deployed panel would actually be chosen,
per this report's own Limitations §VI-E -- only the accuracy *estimate* is
nested, not the final committed SNP list). Selection-frequency-across-folds
is a separate, honest cross-check: how often does this same SNP also show
up in the leak-free, per-fold-rebuilt top-N panel?

No rsIDs: the genotype matrix's own column identifiers are chr:pos[b37]
ref,alt -- no rsID annotation is embedded, and none is fetched here
(would need an external dbSNP/Ensembl lookup, not done in this pass).

Writes: reports/figures/committed_panel_markers.csv (one file, N column
distinguishes the three panels)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nested_cv_sweep import N_CV, RANDOM_STATE, Data, fit_reductor

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "reports/figures/committed_panel_markers.csv"

COMMITTED = [(35, "fst_stat"), (50, "stat"), (70, "FST")]

COL_RE = re.compile(r"^(\d+|X|Y):(\d+)\[(b\d+)\](\w),(\w)$")


def parse_col(col: str) -> tuple[str, int, str, str, str]:
    m = COL_RE.match(col)
    assert m, f"unparseable column id: {col}"
    chrom, pos, build, ref, alt = m.groups()
    return chrom, int(pos), build, ref, alt


def allele_freq(G_col: np.ndarray) -> float:
    """G_col in {0,1,2} dosage of the alt allele; frequency = mean/2."""
    return float(np.nanmean(G_col) / 2.0)


def hudson_fst(p1: float, n1: int, p2: float, n2: int) -> float:
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return float(num / den) if den > 0 else np.nan


data = Data()
class_order = list(data.le.classes_)
DISPLAY = {"CN": "Han", "JPT": "JPT", "SEA": "SEA"}
label_of = np.array([DISPLAY[s] for s in data.y_str])

all_idx = np.arange(len(data.y))

# --- Selection frequency across the 5 outer folds (leak-free, per-fold panels)
skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
fold_panels = {N: [] for N, _ in COMMITTED}
for tr, _ in skf.split(data.G, data.y):
    pools = data.build_pools(tr, which=list({p for _, p in COMMITTED}))
    for N, pool in COMMITTED:
        col_idx = pools[pool]
        X = data.G[tr][:, col_idx].astype(np.float32)
        imp, _ = fit_reductor("EN", X, data.y[tr])
        top_local = np.argsort(imp)[::-1][:N]
        fold_panels[N].append(set(col_idx[top_local].tolist()))

# --- Final "shipped" panel: pool + EN fit on all 504 samples
pools_full = data.build_pools(all_idx, which=list({p for _, p in COMMITTED}))

rows = []
for N, pool in COMMITTED:
    col_idx = pools_full[pool]
    X = data.G[:, col_idx].astype(np.float32)
    imp, _ = fit_reductor("EN", X, data.y)
    top_local = np.argsort(imp)[::-1][:N]
    top_global = col_idx[top_local]

    for rank, gi in enumerate(top_global, start=1):
        col = data.cols[gi]
        chrom, pos, build, ref, alt = parse_col(col)
        Gcol = data.G[:, gi]

        afs = {}
        for grp in ["Han", "JPT", "SEA"]:
            mask = label_of == grp
            afs[grp] = allele_freq(Gcol[mask])

        ns = {grp: int((label_of == grp).sum()) for grp in ["Han", "JPT", "SEA"]}
        pairs = [("Han", "JPT"), ("Han", "SEA"), ("JPT", "SEA")]
        fsts = [hudson_fst(afs[a], ns[a], afs[b], ns[b]) for a, b in pairs]
        max_fst = float(np.nanmax(fsts))

        sel_freq = sum(1 for s in fold_panels[N] if gi in s) / N_CV

        rows.append({
            "N": N, "rank": rank, "chr": chrom, "pos": pos, "build": build,
            "ref": ref, "alt": alt,
            "AF_Han": round(afs["Han"], 4), "AF_JPT": round(afs["JPT"], 4), "AF_SEA": round(afs["SEA"], 4),
            "max_pairwise_FST": round(max_fst, 4),
            "fold_selection_freq": sel_freq,
        })

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"Wrote {OUT} ({len(df)} rows)")
for N, _ in COMMITTED:
    sub = df[df.N == N]
    print(f"\nN={N}: chr distribution:", dict(sub["chr"].value_counts()))
    print(f"N={N}: mean fold_selection_freq = {sub.fold_selection_freq.mean():.2f}, "
          f"fraction always selected (5/5) = {(sub.fold_selection_freq == 1.0).mean():.2f}")
    print(sub.head(5).to_string(index=False))
