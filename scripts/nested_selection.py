"""Train-only candidate-pool construction for leak-free nested CV.

Rebuilds the `stat` and `fst_stat` candidate pools from a given set of samples
(a CV training fold) using the *exact* statistics the original pipeline used, so
that on the full cohort these functions reproduce the committed 1,005-SNP stat
pool (04a/05a) and 1,003-SNP consensus pool (05c) bit-for-bit.

The FST pool itself comes from `nested_fst.fst_pool_for_samples` (plink2 --keep).

Everything here operates on the numpy genotype matrix
(G: n_samples x n_snps int8, missing = -1) and returns *column indices into G*.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

from stats import compute_snp_stats_vectorized


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
# fst_stat pool — replicates notebook 05c consensus filter over the FST pool
# ---------------------------------------------------------------------------

def fst_stat_pool_indices(
    G: np.ndarray,
    pop_labels: np.ndarray,
    fst_pool_idx: Sequence[int],
    snp_ids: Sequence[str],
) -> np.ndarray:
    """05c consensus: from the FST pool, keep chi²-FDR-sig ∧ top-50% MI ∧ top-50% KL.

    Falls back to >=2 tests, then chi²-only, matching 05c's cascade.
    Returns sorted column indices into G.
    """
    fst_pool_idx = list(fst_pool_idx)
    sub = G[:, fst_pool_idx]
    col_names = [snp_ids[i] for i in fst_pool_idx]
    df = pd.DataFrame(sub, columns=col_names)
    df.insert(0, "pop", pop_labels)
    df.insert(0, "sample", np.arange(len(df)))

    stats_df = compute_snp_stats_vectorized(df, pd.Series(pop_labels), verbose=False)
    # stats_df.snp_id is in col_names order → map back to G columns
    id_to_gcol = {snp_ids[i]: i for i in fst_pool_idx}

    chi2_pvals = stats_df["chi2_pvalue"].fillna(1).values
    chi2_reject, _ = fdrcorrection(chi2_pvals, alpha=0.05)
    top_n = max(50, len(stats_df) // 2)

    sig_chi2 = set(stats_df.loc[chi2_reject, "snp_id"])
    sig_mi = set(stats_df.nlargest(top_n, "mutual_information")["snp_id"])
    sig_kl = set(stats_df.nlargest(top_n, "kl_divergence")["snp_id"])

    from collections import Counter
    counts = Counter(list(sig_chi2) + list(sig_mi) + list(sig_kl))
    all3 = [s for s, c in counts.items() if c == 3]
    ge2 = [s for s, c in counts.items() if c >= 2]

    if len(all3) >= 25:
        consensus = all3
    elif len(ge2) >= 25:
        consensus = ge2
    else:
        consensus = list(sig_chi2)

    idx = sorted(id_to_gcol[s] for s in consensus if s in id_to_gcol)
    return np.array(idx, dtype=np.int64)
