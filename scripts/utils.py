"""
Utility functions for data processing, quality control, and analysis.
Shared across filtering, FST analysis, and ML pipelines.
"""

import os
import csv
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional, Union
import pandas as pd
import numpy as np

# Import configuration
try:
    from config import PATHS, HARD_FILTERS, SITUATIONAL_FILTERS, PLINK, POPULATIONS
except ImportError:
    # Allow running from different directories
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from config import PATHS, HARD_FILTERS, SITUATIONAL_FILTERS, PLINK, POPULATIONS


# =============================================================================
# File and Directory Utilities
# =============================================================================


def ensure_dir(path: Union[str, Path]) -> Path:
    """Create directory if it doesn't exist."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_project_root() -> Path:
    """Get the project root directory."""
    return PATHS.ROOT


# =============================================================================
# PLINK2 Command Utilities
# =============================================================================


def run_plink2_command(
    args: List[str],
    check: bool = True,
    capture_output: bool = False,
    verbose: bool = True,
) -> subprocess.CompletedProcess:
    """
    Execute a PLINK2 command with error handling.

    Args:
        args: List of arguments to pass to PLINK2
        check: Raise exception on non-zero exit code
        capture_output: Capture stdout/stderr
        verbose: Print command before execution

    Returns:
        CompletedProcess instance
    """
    cmd = [PLINK.EXECUTABLE, "--threads", str(PLINK.THREADS)] + args

    if verbose:
        print(f"Running: {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=True if capture_output else None,
    )

    return result


def run_bcftools_command(
    args: List[str],
    check: bool = True,
    capture_output: bool = False,
    verbose: bool = True,
) -> subprocess.CompletedProcess:
    """Execute a bcftools command with error handling."""
    cmd = ["bcftools"] + args

    if verbose:
        print(f"Running: {' '.join(cmd)}")

    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=True if capture_output else None,
    )


# =============================================================================
# Panel and Sample Utilities
# =============================================================================


def read_panel(path: Union[str, Path]) -> pd.DataFrame:
    """Read PLINK panel file (tab-delimited)."""
    return pd.read_csv(path, sep="\t")


def extract_population_samples(
    panel_df: pd.DataFrame,
    populations: List[str],
    output_csv: Union[str, Path],
    output_list: Union[str, Path],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Extract samples from specific populations and save to files.

    Args:
        panel_df: DataFrame with sample, pop, super_pop columns
        populations: List of population codes to extract
        output_csv: Path for full metadata CSV
        output_list: Path for sample ID list only

    Returns:
        Tuple of (filtered DataFrame, list of sample IDs)
    """
    filtered_df = panel_df[panel_df["pop"].isin(populations)].copy()

    # Save full metadata (no header for PLINK compatibility)
    filtered_df.to_csv(output_csv, index=False, header=False)

    # Save sample list only
    sample_list = filtered_df["sample"].tolist()
    sample_df = pd.DataFrame({"sample": sample_list})
    sample_df.to_csv(output_list, index=False, header=False)

    print(f"Extracted {len(sample_list)} samples from populations: {populations}")

    return filtered_df, sample_list


# =============================================================================
# Frequency Analysis Utilities
# =============================================================================


def analyze_afreq(
    afreq_path: str,
    num_samples: int = POPULATIONS.NUM_SAMPLES,
    maf_threshold: float = 0.01,
    cr_threshold: float = 0.95,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load and analyze a PLINK2 .afreq file.
    Computes MAF, Call Rate, and QC failure counts.

    Args:
        afreq_path: Path to .afreq file
        num_samples: Number of diploid samples
        maf_threshold: MAF threshold for QC summary
        cr_threshold: Call rate threshold for QC summary
        verbose: Print summary

    Returns:
        df: DataFrame with MAF and CALL_RATE added
        summary: dict with QC counts + descriptive stats
    """
    if not os.path.exists(afreq_path):
        raise FileNotFoundError(f".afreq file not found: {afreq_path}")

    if verbose:
        print(f"Loading data: {afreq_path}")

    # Alleles = 2 per sample (diploid)
    an_total = num_samples * 2

    required_cols = ["ALT_FREQS", "OBS_CT"]
    dtype_map = {
        "ALT_FREQS": "float32",
        "OBS_CT": "uint32",
    }

    df = pd.read_csv(
        afreq_path,
        sep=r"\s+",
        usecols=required_cols,
        dtype=dtype_map,
    )

    # Compute MAF and Call Rate
    df["MAF"] = df["ALT_FREQS"].apply(lambda f: min(float(f), 1.0 - float(f)))
    df["CALL_RATE"] = df["OBS_CT"] / an_total

    if verbose:
        print(f"Loaded {len(df)} variants.")
        print(df.head())

    # QC summary
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
        print("\n--- QC Summary ---")
        print(
            f"Variants failing MAF < {maf_threshold}: {low_maf} ({summary['low_maf_pct']:.2f}%)"
        )
        print(
            f"Variants failing Call Rate < {cr_threshold}: {low_cr} ({summary['low_cr_pct']:.2f}%)"
        )

    return df, summary


# =============================================================================
# Hardy-Weinberg Analysis Utilities
# =============================================================================


def analyze_hardy(
    hardy_path: str,
    p_threshold: float = 1e-6,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Load and analyze a PLINK2 .hardy file.

    Args:
        hardy_path: Path to .hardy file
        p_threshold: P-value threshold for HWE failure
        verbose: Print summary

    Returns:
        df: Hardy dataframe with optimized dtypes
        summary: Dict of summary statistics
    """
    if not os.path.exists(hardy_path):
        raise FileNotFoundError(f".hardy file not found: {hardy_path}")

    if verbose:
        print(f"Loading HWE results: {hardy_path}")

    required_cols = [
        "#CHROM",
        "ID",
        "A1",
        "AX",
        "HOM_A1_CT",
        "HET_A1_CT",
        "TWO_AX_CT",
        "O(HET_A1)",
        "E(HET_A1)",
        "P",
    ]

    dtype_map = {
        "#CHROM": "category",
        "ID": "category",
        "A1": "category",
        "AX": "category",
        "HOM_A1_CT": "uint32",
        "HET_A1_CT": "uint32",
        "TWO_AX_CT": "uint32",
        "O(HET_A1)": "float32",
        "E(HET_A1)": "float32",
        "P": "float64",
    }

    df = pd.read_csv(hardy_path, sep=r"\s+", usecols=required_cols, dtype=dtype_map)

    if verbose:
        print(f"Loaded {len(df)} variants.")
        print(df.head())

    hwe_fail = (df["P"] < p_threshold).sum()

    summary = {
        "num_variants": len(df),
        "hwe_fail_variants": int(hwe_fail),
        "hwe_fail_pct": float(hwe_fail / len(df) * 100),
        "pvalue_describe": df["P"].describe(),
        "observed_het_describe": df["O(HET_A1)"].describe(),
        "expected_het_describe": df["E(HET_A1)"].describe(),
    }

    if verbose:
        print("\n--- HWE QC Summary ---")
        print(
            f"Variants failing HWE (P < {p_threshold}): "
            f"{hwe_fail} ({summary['hwe_fail_pct']:.4f}%)"
        )

    return df, summary


# =============================================================================
# FST Analysis Utilities
# =============================================================================


def analyze_fst(
    fst_path: str,
    top_n: int = 1000,
    plot: bool = False,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load and analyze PLINK2 FST file.

    Args:
        fst_path: Path to .fst.var file
        top_n: Number of top variants to extract
        plot: Whether to generate plots
        verbose: Print summary

    Returns:
        df: Full FST dataframe
        top_df: Top-N variants by FST
    """
    if verbose:
        print(f"Loading FST file: {fst_path}")

    df = pd.read_csv(fst_path, sep=r"\s+")

    if verbose:
        print(df.head())
        print(df.describe())
        print(f"Total variants: {len(df)}")
        print(f"NaN FST: {df['HUDSON_FST'].isna().sum()}")
        print(f"Negative FST: {(df['HUDSON_FST'] < 0).sum()}")
        print(f"Zero FST: {(df['HUDSON_FST'] == 0).sum()}")
        print(f"Positive FST: {(df['HUDSON_FST'] > 0).sum()}")

    fst_filtered = df[df["HUDSON_FST"].notna()]
    top_df = fst_filtered.nlargest(top_n, "HUDSON_FST")

    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 4))
        fst_filtered["HUDSON_FST"].hist(bins=50)
        plt.title("Distribution of per-variant Hudson FST")
        plt.xlabel("FST")
        plt.ylabel("Count")
        plt.show()

    if verbose:
        print(f"\nTop {top_n} variants:")
        print(top_df.head())

    return df, top_df


def merge_fst_top_variants(*dfs: pd.DataFrame, id_col: str = "ID") -> pd.DataFrame:
    """
    Merge multiple FST top-variant DataFrames keeping unique variants.

    Args:
        dfs: Variable number of DataFrames to merge
        id_col: Column name for variant IDs

    Returns:
        Merged DataFrame with unique variants
    """
    merged = pd.concat(dfs, ignore_index=True)
    merged_unique = merged.drop_duplicates(subset=[id_col])
    print(
        f"Merged {len(dfs)} dataframes: {merged.shape[0]} total -> {merged_unique.shape[0]} unique variants"
    )
    return merged_unique


# =============================================================================
# Variant ID and BED Conversion Utilities
# =============================================================================


def variant_to_bed(variant: str) -> Tuple[str, int, int, str]:
    """
    Convert variant string (e.g., '12:58124534[b37]C,G') to BED format.

    Args:
        variant: Variant ID string in format "chrom:pos[build]ref,alt"

    Returns:
        Tuple of (chrom, start, end, name)
    """
    try:
        chrom_pos, allele_info = variant.split("]")
        chrom_pos = chrom_pos.replace("[b37", "")
        chrom, pos_str = chrom_pos.split(":")
        pos = int(pos_str)
        start = pos - 1  # BED is 0-based
        end = pos  # BED end is 1-based
        name = allele_info
        return (chrom, start, end, name)
    except Exception as e:
        raise ValueError(f"Failed to parse variant '{variant}': {e}")


def variants_to_bed_file(
    snp_file: Union[str, Path],
    bed_file: Union[str, Path],
    verbose: bool = True,
) -> int:
    """
    Convert a file of variant IDs to BED format.

    Args:
        snp_file: Input file with one variant ID per line
        bed_file: Output BED file path
        verbose: Print progress

    Returns:
        Number of variants converted
    """
    count = 0
    with open(snp_file, "r") as sf, open(bed_file, "w", newline="") as bf:
        bed_writer = csv.writer(bf, delimiter="\t")
        for line in sf:
            variant = line.strip()
            if variant:
                chrom, start, end, name = variant_to_bed(variant)
                bed_writer.writerow([chrom, start, end, name])
                count += 1

    if verbose:
        print(f"Converted {count} variants to BED format: {bed_file}")

    return count


# =============================================================================
# VCF to Numeric Matrix Utilities
# =============================================================================


def genotype_to_numeric(gt: str) -> Optional[int]:
    """
    Convert genotype string to numeric value.

    Args:
        gt: Genotype string (e.g., "0|0", "0/1", "1|1")

    Returns:
        0 for hom_ref, 1 for het, 2 for hom_alt, None for missing
    """
    gt_base = str(gt).split(":")[0]
    if gt_base in {"0|0", "0/0"}:
        return 0
    elif gt_base in {"0|1", "1|0", "0/1", "1/0"}:
        return 1
    elif gt_base in {"1|1", "1/1"}:
        return 2
    else:
        return None  # Handle missing or unexpected genotypes


def vcf_to_numeric_matrix(
    vcf_file: Union[str, Path],
    output_csv: Union[str, Path],
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Convert VCF file to numeric genotype matrix.

    Args:
        vcf_file: Path to VCF file (without ## header lines)
        output_csv: Output CSV path
        verbose: Print progress

    Returns:
        DataFrame with samples as rows and variants as columns
    """
    if verbose:
        print(f"Reading VCF: {vcf_file}")

    # Read VCF header
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                break

    # Load data
    df = pd.read_csv(vcf_file, comment="#", sep="\t", names=header, dtype=str)

    # Create variant ID
    df["pop_id"] = df["#CHROM"] + ":" + df["POS"]

    sample_cols = header[9:]
    df_samples = df[["pop_id"] + sample_cols].copy()

    # Convert genotypes
    for col in sample_cols:
        df_samples[col] = df_samples[col].apply(genotype_to_numeric)

    # Transpose so samples are rows
    df_final = df_samples.set_index("pop_id").T.reset_index()
    df_final = df_final.rename(columns={"index": "sample"})

    df_final.to_csv(output_csv, index=False)

    if verbose:
        print(f"Saved numeric matrix to {output_csv}")
        print(f"Shape: {df_final.shape}")
        print(df_final.head())

    return df_final


def add_population_labels(
    numeric_csv: Union[str, Path],
    population_csv: Union[str, Path],
    output_csv: Union[str, Path],
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Merge numeric matrix with population labels.

    Args:
        numeric_csv: Path to numeric genotype matrix
        population_csv: Path to population metadata
        output_csv: Output path for merged data
        verbose: Print progress

    Returns:
        Merged DataFrame
    """
    population_df = pd.read_csv(population_csv, header=None)
    population_df.columns = ["sample", "pop", "super_pop"]

    numeric_df = pd.read_csv(numeric_csv)

    df_final = pd.merge(numeric_df, population_df, on="sample", how="left")
    df_final = df_final.drop(columns=["super_pop"])

    df_final.to_csv(output_csv, index=False)

    if verbose:
        print(f"Merged with population labels: {output_csv}")
        print(f"Shape: {df_final.shape}")
        print(df_final.head())

    return df_final


# =============================================================================
# Reporting Utilities
# =============================================================================


def save_report(
    report_text: str,
    output_path: Union[str, Path],
    verbose: bool = True,
) -> None:
    """Save a text report to file."""
    ensure_dir(Path(output_path).parent)
    with open(output_path, "w") as f:
        f.write(report_text)
    if verbose:
        print(f"Report saved to: {output_path}")


def count_variants(pfile_prefix: Union[str, Path]) -> int:
    """Count number of variants in a PLINK2 pfile."""
    pvar_path = f"{pfile_prefix}.pvar"
    if os.path.exists(pvar_path):
        with open(pvar_path) as f:
            # Skip header lines starting with #
            count = sum(1 for line in f if not line.startswith("#"))
        return count
    return -1


def count_samples(pfile_prefix: Union[str, Path]) -> int:
    """Count number of samples in a PLINK2 pfile."""
    psam_path = f"{pfile_prefix}.psam"
    if os.path.exists(psam_path):
        with open(psam_path) as f:
            # Skip header line
            count = sum(1 for line in f if not line.startswith("#")) - 1
        return count
    return -1


if __name__ == "__main__":
    # Test utilities
    print("Testing utility functions...")
    print(f"Project root: {get_project_root()}")
