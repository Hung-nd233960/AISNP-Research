"""bcftools wrappers and VCF-to-numeric-matrix conversion."""

import subprocess
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd


def run_bcftools_command(
    args: List[str],
    check: bool = True,
    capture_output: bool = False,
    verbose: bool = True,
) -> subprocess.CompletedProcess:
    """Execute a bcftools command."""
    cmd = ["bcftools"] + args
    if verbose:
        print(f"Running: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=True if capture_output else None,
    )


def genotype_to_numeric(gt: str) -> Optional[int]:
    """
    Convert a genotype string to 0 / 1 / 2 / None.

    0 = hom-ref, 1 = het, 2 = hom-alt, None = missing/unknown.
    """
    gt_base = str(gt).split(":")[0]
    if gt_base in {"0|0", "0/0"}:
        return 0
    elif gt_base in {"0|1", "1|0", "0/1", "1/0"}:
        return 1
    elif gt_base in {"1|1", "1/1"}:
        return 2
    return None


def vcf_to_numeric_matrix(
    vcf_file: Union[str, Path],
    output_csv: Union[str, Path],
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Convert a VCF file to a numeric genotype matrix (samples × variants).
    Saves to output_csv and returns the DataFrame.
    """
    if verbose:
        print(f"Reading VCF: {vcf_file}")

    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                break

    df = pd.read_csv(vcf_file, comment="#", sep="\t", names=header, dtype=str)
    df["pop_id"] = df["#CHROM"] + ":" + df["POS"]

    sample_cols = header[9:]
    df_samples = df[["pop_id"] + sample_cols].copy()
    for col in sample_cols:
        df_samples[col] = df_samples[col].apply(genotype_to_numeric)

    df_final = df_samples.set_index("pop_id").T.reset_index()
    df_final = df_final.rename(columns={"index": "sample"})
    df_final.to_csv(output_csv, index=False)

    if verbose:
        print(f"Saved numeric matrix to {output_csv}  shape={df_final.shape}")

    return df_final


def add_population_labels(
    numeric_csv: Union[str, Path],
    population_csv: Union[str, Path],
    output_csv: Union[str, Path],
    verbose: bool = True,
) -> pd.DataFrame:
    """Merge a numeric genotype matrix with population labels."""
    pop_df = pd.read_csv(population_csv, header=None, names=["sample", "pop", "super_pop"])
    num_df = pd.read_csv(numeric_csv)

    df_final = pd.merge(num_df, pop_df[["sample", "pop"]], on="sample", how="left")
    df_final.to_csv(output_csv, index=False)

    if verbose:
        print(f"Merged with population labels → {output_csv}  shape={df_final.shape}")

    return df_final
