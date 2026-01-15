#!/usr/bin/env python3
"""
Add population panel information to genotyped samples CSV.

Merges sample-level population and super_pop labels from 1000 Genomes panel file.
"""

import argparse
from pathlib import Path
import pandas as pd


def add_panels(
    genotype_csv: Path,
    panel_file: Path,
    output_file: Path | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Add population panels to genotyped samples.

    Args:
        genotype_csv: Path to CSV with genotyped samples (must have 'sample' column).
        panel_file: Path to 1000 Genomes panel file.
        output_file: Output path (optional). If None, returns DataFrame only.
        verbose: Print progress.

    Returns:
        Merged DataFrame with pop and super_pop columns.
    """
    # Load panel file
    panel_df = pd.read_csv(panel_file, sep=r"\s+")
    if verbose:
        print(f"Loaded panel: {len(panel_df)} samples")

    # Load genotype file
    genotype_df = pd.read_csv(genotype_csv)
    if verbose:
        print(
            f"Loaded genotypes: {len(genotype_df)} samples, {len(genotype_df.columns)-1} variants"
        )

    # Merge
    merged_df = pd.merge(
        genotype_df,
        panel_df[["sample", "pop", "super_pop"]],
        on="sample",
        how="left",
    )

    # Check for missing panels
    missing = merged_df["pop"].isna().sum()
    if missing > 0 and verbose:
        print(f"Warning: {missing} samples have no panel information")

    if verbose:
        print(f"Merged: {len(merged_df)} samples")
        print(f"Populations: {merged_df['pop'].nunique()}")
        print(f"Super populations: {merged_df['super_pop'].nunique()}")

    # Save if output specified
    if output_file:
        merged_df.to_csv(output_file, index=False)
        if verbose:
            print(f"Saved to: {output_file}")

    return merged_df


def main():
    parser = argparse.ArgumentParser(
        description="Add population panel information to genotyped samples."
    )
    parser.add_argument(
        "genotype_csv",
        type=Path,
        help="Input CSV with genotyped samples (must have 'sample' column)",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("1000genomes/integrated_call_samples_v3.20130502.ALL.panel"),
        help="1000 Genomes panel file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file path (default: <input>_with_panels.csv)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress output",
    )

    args = parser.parse_args()

    # Default output name
    if args.output is None:
        stem = args.genotype_csv.stem
        args.output = args.genotype_csv.parent / f"{stem}_with_panels.csv"

    add_panels(
        genotype_csv=args.genotype_csv,
        panel_file=args.panel,
        output_file=args.output,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
