"""Train-only candidate-pool construction for leak-free nested CV.

Rebuilds the `stat` and `fst_stat` candidate pools from a given set of samples
(a CV training fold) using the *exact* statistics the original pipeline used, so
that on the full cohort `stat_pool_indices` reproduces the committed 1,005-SNP
stat pool (04a/05a) bit-for-bit.

`fst_stat_pool_indices` intentionally does NOT reproduce the original notebook
05c's consensus pool: 05c filtered the FST block using a different statistical
trio (chi-squared, mutual information, pairwise KL divergence) than the stat
pool uses (chi-squared, JSD, allele-frequency delta), with no documented reason
for the mismatch beyond 05c's own note that running on the smaller FST-filtered
subset is "much faster than running on the full 600k set" — a runtime remark,
not a rationale for *which* statistics to use. That mismatched trio is treated
here as an artifact of an older architecture, not a deliberate design choice.
This module instead cascades the stat pool's own construction (`stat_pool_indices`)
onto the FST block's candidates, so `fst_stat` differs from `stat` only in which
614,759-SNP subset it starts from, not in which statistics it applies.

The FST pool itself comes from `nested_fst.fst_pool_for_samples` (plink2 --keep).

Everything here operates on the numpy genotype matrix
(G: n_samples x n_snps int8, missing = -1) and returns *column indices into G*.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# stat pool — replicates notebook 04a exactly (chi2_stat / jsd / af_max_delta)
# ---------------------------------------------------------------------------

def _geno_count_tensor(G: np.ndarray, pop_labels: np.ndarray):
    """Return (pop_geno_counts [n_pops,3,n_snps], n_sizes [n_pops], pop_names)."""
    pop_names = sorted(np.unique(pop_labels))
    n_pops = len(pop_names)
    n_snps = G.shape[1]
    pop_geno_counts = np.zeros((n_pops, 3, n_snps), dtype=np.float32)
    n_sizes = np.zeros(n_pops, dtype=np.float32)
    for i, p in enumerate(pop_names):
        mask = pop_labels == p
        n_sizes[i] = mask.sum()
        Gp = G[mask]
        for v in range(3):
            pop_geno_counts[i, v] = (Gp == v).sum(axis=0)
    return pop_geno_counts, n_sizes, pop_names


def stat_scores(G: np.ndarray, pop_labels: np.ndarray) -> pd.DataFrame:
    """Per-SNP chi2_stat, jsd, af_max_delta — identical to 04a_snp_selection."""
    pop_geno_counts, n_sizes, _ = _geno_count_tensor(G, pop_labels)
    n_total = G.shape[0]

    # chi-squared (Pearson) over (n_pops x 3) contingency table
    row_sums = pop_geno_counts.sum(axis=1)               # (n_pops, n_snps)
    col_sums = pop_geno_counts.sum(axis=0)               # (3, n_snps)
    expected = (row_sums[:, None, :] * col_sums[None, :, :]) / n_total
    safe_expected = np.where(expected > 0.5, expected, np.inf)
    chi2_stat = ((pop_geno_counts - expected) ** 2 / safe_expected).sum(axis=(0, 1))

    # Jensen-Shannon divergence across the 3 genotype distributions
    pop_geno_freq = (pop_geno_counts + 1e-10) / (n_sizes[:, None, None] + 3e-10)
    M = pop_geno_freq.mean(axis=0)                       # (3, n_snps)
    kl_from_M = (pop_geno_freq * np.log(pop_geno_freq / (M[None, :, :] + 1e-10))).sum(axis=1)
    jsd = np.maximum(kl_from_M.mean(axis=0), 0)

    # max pairwise allele-frequency delta
    af = (pop_geno_counts[:, 1, :] + 2 * pop_geno_counts[:, 2, :]) / (2 * n_sizes[:, None] + 1e-10)
    n_pops = af.shape[0]
    af_max_delta = np.zeros(G.shape[1], dtype=np.float32)
    for i in range(n_pops):
        for j in range(i + 1, n_pops):
            af_max_delta = np.maximum(af_max_delta, np.abs(af[i] - af[j]))

    return pd.DataFrame({
        "col": np.arange(G.shape[1]),
        "chi2_stat": chi2_stat,
        "jsd": jsd,
        "af_max_delta": af_max_delta,
    })


def stat_pool_indices(G: np.ndarray, pop_labels: np.ndarray, top_n: int = 500) -> np.ndarray:
    """Union of top-`top_n` by chi2_stat, jsd, af_max_delta → sorted column indices."""
    s = stat_scores(G, pop_labels)
    union = set(s.nlargest(top_n, "chi2_stat")["col"])
    union |= set(s.nlargest(top_n, "jsd")["col"])
    union |= set(s.nlargest(top_n, "af_max_delta")["col"])
    return np.array(sorted(union), dtype=np.int64)


# ---------------------------------------------------------------------------
# fst_stat pool — stat-pool statistics (chi2/JSD/AFD), cascaded onto the FST
# block's own candidates. See module docstring: this replaces an earlier
# MI/KL-based consensus filter (notebook 05c) that used a different, and
# undocumented-as-intentional, statistical trio than the stat pool itself.
# ---------------------------------------------------------------------------

def fst_stat_pool_indices(
    G: np.ndarray,
    pop_labels: np.ndarray,
    fst_pool_idx: Sequence[int],
    top_n: int = 500,
) -> np.ndarray:
    """Cascade stat_pool_indices onto the FST block's candidates only.

    Identical mechanism to stat_pool_indices (union of top-`top_n` by
    chi2_stat, jsd, af_max_delta) — the only difference from the stat pool is
    that it screens the FST block's SNPs (fst_pool_idx) instead of all
    614,759, so `fst_stat` and `stat` use the same statistics and differ only
    in candidate scope. Returns sorted column indices into G.
    """
    fst_pool_idx = np.asarray(sorted(fst_pool_idx), dtype=np.int64)
    G_fst = G[:, fst_pool_idx]
    local_idx = stat_pool_indices(G_fst, pop_labels, top_n=top_n)
    return np.sort(fst_pool_idx[local_idx])
