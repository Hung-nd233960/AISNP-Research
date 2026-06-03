"""Statistical QC and FST analysis utilities (reads PLINK2 output files)."""

import os
from itertools import combinations
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from sklearn.preprocessing import LabelEncoder

from config import POPULATIONS


def compute_snp_stats_vectorized(
    genotype_df: pd.DataFrame,
    populations: pd.Series,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compute chi², mutual information, and KL divergence for every SNP at once.

    Builds a (n_pops × n_genotypes × n_snps) contingency array in one pass using
    matrix operations, then derives all three statistics from it — no Python loop
    over SNPs.

    Args:
        genotype_df: DataFrame with columns ['sample','pop', snp1, snp2, ...]
                     Genotype values: 0/1/2, missing encoded as -1.
        populations:  Series of population labels aligned with genotype_df rows.

    Returns:
        DataFrame with columns: snp_id, chi2, chi2_pvalue,
                                 mutual_information, kl_divergence
    """
    snp_cols = [c for c in genotype_df.columns if c not in ("sample", "pop")]
    geno_matrix = genotype_df[snp_cols].to_numpy(dtype=np.int8)  # (n_samples, n_snps)
    n_samples, n_snps = geno_matrix.shape

    le = LabelEncoder()
    pop_enc = le.fit_transform(populations.values)
    n_pops = len(le.classes_)
    n_genos = 3  # 0, 1, 2

    if verbose:
        print(f"Vectorised stats: {n_snps:,} SNPs × {n_samples} samples × "
              f"{n_pops} populations")

    # ------------------------------------------------------------------
    # Build counts[p, g, snp] = # samples in population p with genotype g
    # Subset rows first, then sum — avoids bool/int matmul accumulation bugs.
    # numpy .sum(axis=0) on bool always accumulates in int64: safe, no overflow.
    # ------------------------------------------------------------------
    counts = np.zeros((n_pops, n_genos, n_snps), dtype=np.float32)
    for p in range(n_pops):
        pop_genos = geno_matrix[pop_enc == p]        # (n_i, n_snps) int8
        for g in range(n_genos):
            counts[p, g] = (pop_genos == g).sum(axis=0)

    total = counts.sum(axis=(0, 1))                  # (n_snps,) — valid calls only

    # ------------------------------------------------------------------
    # Chi-squared (Pearson)
    # ------------------------------------------------------------------
    row_sums = counts.sum(axis=1)                    # (n_pops,  n_snps)
    col_sums = counts.sum(axis=0)                    # (n_genos, n_snps)
    expected = (
        row_sums[:, np.newaxis, :]                   # (n_pops,  1,       n_snps)
        * col_sums[np.newaxis, :, :]                 # (1,       n_genos, n_snps)
        / np.maximum(total, 1)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2_stats = np.nansum(
            np.where(expected > 0, (counts - expected) ** 2 / expected, 0.0),
            axis=(0, 1),
        )
    dof = (n_pops - 1) * (n_genos - 1)
    chi2_pvals = chi2_dist.sf(chi2_stats, df=dof)

    # ------------------------------------------------------------------
    # Mutual information
    # MI = Σ p(x,y) log( p(x,y) / (p(x) p(y)) )
    # ------------------------------------------------------------------
    p_joint = counts / np.maximum(total, 1)
    p_pop  = p_joint.sum(axis=1, keepdims=True)      # (n_pops,  1,       n_snps)
    p_geno = p_joint.sum(axis=0, keepdims=True)      # (1,       n_genos, n_snps)
    outer  = p_pop * p_geno
    with np.errstate(divide="ignore", invalid="ignore"):
        mi = np.nansum(
            np.where(
                (p_joint > 0) & (outer > 0),
                p_joint * np.log(p_joint / outer),
                0.0,
            ),
            axis=(0, 1),
        )

    # ------------------------------------------------------------------
    # Symmetric KL divergence (mean over all population pairs)
    # KL(P||Q) = Σ P(g) log( P(g) / Q(g) )
    # ------------------------------------------------------------------
    p_cond = counts / np.maximum(row_sums[:, np.newaxis, :], 1)  # (n_pops, n_genos, n_snps)
    kl_sum = np.zeros(n_snps, dtype=np.float32)
    n_pairs = 0
    for i, j in combinations(range(n_pops), 2):
        p = p_cond[i]   # (n_genos, n_snps)
        q = p_cond[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            kl_fwd = np.nansum(
                np.where(p > 0, p * np.log(p / np.maximum(q, 1e-10)), 0.0), axis=0
            )
            kl_rev = np.nansum(
                np.where(q > 0, q * np.log(q / np.maximum(p, 1e-10)), 0.0), axis=0
            )
        kl_sum += (kl_fwd + kl_rev) / 2
        n_pairs += 1
    kl_div = kl_sum / n_pairs

    if verbose:
        print("Done.")

    return pd.DataFrame({
        "snp_id":              snp_cols,
        "chi2":                chi2_stats.astype(np.float32),
        "chi2_pvalue":         chi2_pvals,
        "mutual_information":  mi.astype(np.float32),
        "kl_divergence":       kl_div,
    })


def analyze_afreq(
    afreq_path: str,
    num_samples: int = POPULATIONS.NUM_SAMPLES,
    maf_threshold: float = 0.01,
    cr_threshold: float = 0.95,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Load a PLINK2 .afreq file and compute MAF + call-rate QC summary."""
    if not os.path.exists(afreq_path):
        raise FileNotFoundError(f".afreq file not found: {afreq_path}")

    if verbose:
        print(f"Loading: {afreq_path}")

    an_total = num_samples * 2
    df = pd.read_csv(
        afreq_path,
        sep=r"\s+",
        usecols=["ALT_FREQS", "OBS_CT"],
        dtype={"ALT_FREQS": "float32", "OBS_CT": "uint32"},
    )
    df["MAF"] = df["ALT_FREQS"].apply(lambda f: min(float(f), 1.0 - float(f)))
    df["CALL_RATE"] = df["OBS_CT"] / an_total

    low_maf = int((df["MAF"] < maf_threshold).sum())
    low_cr = int((df["CALL_RATE"] < cr_threshold).sum())
    summary = {
        "num_variants": len(df),
        "low_maf_variants": low_maf,
        "low_cr_variants": low_cr,
        "low_maf_pct": low_maf / len(df) * 100,
        "low_cr_pct": low_cr / len(df) * 100,
        "maf_describe": df["MAF"].describe(),
        "cr_describe": df["CALL_RATE"].describe(),
    }

    if verbose:
        print(f"  {len(df)} variants | low MAF: {low_maf} | low CR: {low_cr}")

    return df, summary


def analyze_hardy(
    hardy_path: str,
    p_threshold: float = 1e-6,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Load a PLINK2 .hardy file and compute HWE deviation summary."""
    if not os.path.exists(hardy_path):
        raise FileNotFoundError(f".hardy file not found: {hardy_path}")

    if verbose:
        print(f"Loading: {hardy_path}")

    cols = ["#CHROM", "ID", "A1", "AX", "HOM_A1_CT", "HET_A1_CT",
            "TWO_AX_CT", "O(HET_A1)", "E(HET_A1)", "P"]
    dtypes = {
        "#CHROM": "category", "ID": "category", "A1": "category", "AX": "category",
        "HOM_A1_CT": "uint32", "HET_A1_CT": "uint32", "TWO_AX_CT": "uint32",
        "O(HET_A1)": "float32", "E(HET_A1)": "float32", "P": "float64",
    }
    df = pd.read_csv(hardy_path, sep=r"\s+", usecols=cols, dtype=dtypes)

    hwe_fail = int((df["P"] < p_threshold).sum())
    summary = {
        "num_variants": len(df),
        "hwe_fail_variants": hwe_fail,
        "hwe_fail_pct": hwe_fail / len(df) * 100,
        "pvalue_describe": df["P"].describe(),
        "observed_het_describe": df["O(HET_A1)"].describe(),
        "expected_het_describe": df["E(HET_A1)"].describe(),
    }

    if verbose:
        print(f"  {len(df)} variants | HWE fail (P<{p_threshold}): {hwe_fail} ({summary['hwe_fail_pct']:.4f}%)")

    return df, summary


def analyze_fst(
    fst_path: str,
    top_n: int = 1000,
    plot: bool = False,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load a PLINK2 .fst.var file and return (full_df, top_n_df)."""
    if verbose:
        print(f"Loading FST file: {fst_path}")

    df = pd.read_csv(fst_path, sep=r"\s+")

    if verbose:
        print(f"  {len(df)} variants | NaN: {df['HUDSON_FST'].isna().sum()} "
              f"| pos: {(df['HUDSON_FST'] > 0).sum()}")

    fst_valid = df[df["HUDSON_FST"].notna()]
    top_df = fst_valid.nlargest(top_n, "HUDSON_FST")

    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 4))
        fst_valid["HUDSON_FST"].hist(bins=50)
        plt.title("Distribution of per-variant Hudson FST")
        plt.xlabel("FST")
        plt.ylabel("Count")
        plt.show()

    return df, top_df


def merge_fst_top_variants(*dfs: pd.DataFrame, id_col: str = "ID") -> pd.DataFrame:
    """Concatenate FST top-variant DataFrames and deduplicate by variant ID."""
    merged = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=[id_col])
    print(f"Merged {len(dfs)} FST tables → {len(merged)} unique variants")
    return merged
