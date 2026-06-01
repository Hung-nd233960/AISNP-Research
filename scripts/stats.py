"""Statistical QC and FST analysis utilities (reads PLINK2 output files)."""

import os
from typing import Any, Dict, Tuple

import pandas as pd

from config import POPULATIONS


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
        plt.xlabel("FST"); plt.ylabel("Count")
        plt.show()

    return df, top_df


def merge_fst_top_variants(*dfs: pd.DataFrame, id_col: str = "ID") -> pd.DataFrame:
    """Concatenate FST top-variant DataFrames and deduplicate by variant ID."""
    merged = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=[id_col])
    print(f"Merged {len(dfs)} FST tables → {len(merged)} unique variants")
    return merged
