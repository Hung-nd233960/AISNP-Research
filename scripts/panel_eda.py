#!/usr/bin/env python3
"""
Panel EDA - Exploratory Data Analysis for 1000 Genomes panel data.

Provides summary statistics and visualizations for population structure.
"""

import argparse
from pathlib import Path

import pandas as pd


def basic_eda(df: pd.DataFrame, name: str = "DataFrame") -> dict:
    """
    Perform basic exploratory data analysis on a DataFrame.

    Args:
        df: Input DataFrame.
        name: Name for display.

    Returns:
        Dictionary with EDA results.
    """
    print("=" * 80)
    print(f"📊 BASIC EDA REPORT FOR: {name}")
    print("=" * 80)

    results = {}

    # 1. Shape
    print(f"\n🔎 **Shape**: {df.shape[0]} rows × {df.shape[1]} columns")
    results["shape"] = df.shape

    # 2. Column types
    print("\n🔎 **Column Data Types**")
    print(df.dtypes.to_string())
    results["dtypes"] = df.dtypes.to_dict()

    # 3. First rows
    print("\n🔎 **First 5 Rows**")
    print(df.head().to_string())

    # 4. Missing values
    missing = df.isna().sum()
    if missing.sum() > 0:
        print("\n🔎 **Missing Values**")
        print(missing[missing > 0].to_string())
    else:
        print("\n🔎 **Missing Values**: None")
    results["missing"] = missing.to_dict()

    # 5. Unique counts
    print("\n🔎 **Unique Values per Column**")
    print(df.nunique().to_string())
    results["nunique"] = df.nunique().to_dict()

    # 6. Categorical distributions
    cat_cols = df.select_dtypes(include=["object", "category"]).columns
    if len(cat_cols) > 0:
        print("\n🔎 **Categorical Column Distributions**")
        for col in cat_cols:
            print(f"\n➡ {col}:")
            print(df[col].value_counts().to_string())
            results[f"{col}_counts"] = df[col].value_counts().to_dict()

    # 7. Numeric summary
    num_cols = df.select_dtypes(include=["number"]).columns
    if len(num_cols) > 0:
        print("\n🔎 **Numeric Summary Statistics**")
        print(df[num_cols].describe().to_string())
        results["numeric_summary"] = df[num_cols].describe().to_dict()

    return results


def panel_eda(
    panel_file: Path,
    super_pop: str | None = None,
    populations: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Perform EDA on 1000 Genomes panel file.

    Args:
        panel_file: Path to panel file.
        super_pop: Filter by super population.
        populations: Filter by specific populations.
        verbose: Print results.

    Returns:
        Filtered panel DataFrame.
    """
    # Load panel
    panel_df = pd.read_csv(panel_file, sep=r"\s+")

    if verbose:
        basic_eda(panel_df, name=f"1000 Genomes Panel ({panel_file.name})")

    # Filter if requested
    if super_pop:
        panel_df = panel_df[panel_df["super_pop"] == super_pop]
        if verbose:
            print(f"\n🔎 **Filtered to {super_pop}**: {len(panel_df)} samples")
            print(panel_df["pop"].value_counts().to_string())

    if populations:
        panel_df = panel_df[panel_df["pop"].isin(populations)]
        if verbose:
            print(f"\n🔎 **Filtered to {populations}**: {len(panel_df)} samples")

    return panel_df


def main():
    parser = argparse.ArgumentParser(
        description="Exploratory Data Analysis for 1000 Genomes panel data."
    )
    parser.add_argument(
        "panel_file",
        type=Path,
        help="Path to 1000 Genomes panel file",
    )
    parser.add_argument(
        "--super-pop",
        type=str,
        default=None,
        help="Filter by super population (e.g., EAS, EUR, AFR)",
    )
    parser.add_argument(
        "--populations",
        type=str,
        nargs="+",
        default=None,
        help="Filter by specific populations (e.g., CHB JPT KHV)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Save filtered panel to CSV",
    )

    args = parser.parse_args()

    df = panel_eda(
        panel_file=args.panel_file,
        super_pop=args.super_pop,
        populations=args.populations,
    )

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved filtered panel to: {args.output}")


if __name__ == "__main__":
    main()
