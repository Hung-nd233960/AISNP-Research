#!/usr/bin/env python3
"""
PCA visualization for population genetics data.

Creates publication-ready PCA scatter plots colored by population.
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")


def load_pca_data(
    eigenvec_file: Path,
    panel_file: Path,
    population_col: str = "pop",
) -> pd.DataFrame:
    """
    Load PCA results and merge with population labels.

    Args:
        eigenvec_file: PLINK2 .eigenvec file.
        panel_file: 1000 Genomes panel file.
        population_col: Column for population labels ('pop' or 'super_pop').

    Returns:
        DataFrame with PCA coordinates and population labels.
    """
    # Load eigenvector file
    pca_df = pd.read_csv(eigenvec_file, sep=r"\s+", header=None)

    # PLINK2 format: first column is sample ID (or #IID)
    # Check if first row is header
    if pca_df.iloc[0, 0] == "#IID" or pca_df.iloc[0, 0] == "IID":
        pca_df = pd.read_csv(eigenvec_file, sep=r"\s+")
        pca_df = pca_df.rename(columns={"#IID": "IID"})
    else:
        n_pcs = pca_df.shape[1] - 1
        pca_df.columns = ["IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]

    # Load panel
    panel_df = pd.read_csv(panel_file, sep=r"\s+")

    # Merge
    df = pca_df.merge(
        panel_df[["sample", population_col]],
        left_on="IID",
        right_on="sample",
        how="left",
    )

    return df


def plot_pca(
    df: pd.DataFrame,
    population_col: str = "pop",
    pc_x: str = "PC1",
    pc_y: str = "PC2",
    output_file: Path | None = None,
    title: str = "PCA Plot",
    figsize: tuple = (10, 8),
    alpha: float = 0.7,
    point_size: int = 20,
) -> None:
    """
    Create PCA scatter plot colored by population.

    Args:
        df: DataFrame with PCA coordinates and population labels.
        population_col: Column with population labels.
        pc_x: X-axis PC column.
        pc_y: Y-axis PC column.
        output_file: Output file path (PNG/PDF).
        title: Plot title.
        figsize: Figure size.
        alpha: Point transparency.
        point_size: Point size.
    """
    plt.figure(figsize=figsize)

    populations = df[population_col].dropna().unique()
    colors = plt.cm.tab10(range(len(populations)))

    for pop, color in zip(sorted(populations), colors):
        subset = df[df[population_col] == pop]
        plt.scatter(
            subset[pc_x],
            subset[pc_y],
            c=[color],
            label=pop,
            s=point_size,
            alpha=alpha,
            edgecolors="none",
        )

    plt.xlabel(pc_x)
    plt.ylabel(pc_y)
    plt.title(title)
    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        borderaxespad=0,
    )
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved: {output_file}")

    plt.close()


def plot_pca_grid(
    df: pd.DataFrame,
    population_col: str = "pop",
    n_pcs: int = 4,
    output_file: Path | None = None,
    title: str = "PCA Grid",
) -> None:
    """
    Create grid of PC pairs (PC1 vs PC2, PC1 vs PC3, etc.).

    Args:
        df: DataFrame with PCA coordinates.
        population_col: Column with population labels.
        n_pcs: Number of PCs to include.
        output_file: Output file path.
        title: Plot title.
    """
    import itertools

    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    pc_pairs = list(itertools.combinations(pc_cols, 2))

    n_plots = len(pc_pairs)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = axes.flatten() if n_plots > 1 else [axes]

    populations = sorted(df[population_col].dropna().unique())
    colors = plt.cm.tab10(range(len(populations)))

    for ax, (pc_x, pc_y) in zip(axes, pc_pairs):
        for pop, color in zip(populations, colors):
            subset = df[df[population_col] == pop]
            ax.scatter(
                subset[pc_x],
                subset[pc_y],
                c=[color],
                label=pop,
                s=15,
                alpha=0.6,
                edgecolors="none",
            )
        ax.set_xlabel(pc_x)
        ax.set_ylabel(pc_y)

    # Hide empty subplots
    for ax in axes[n_plots:]:
        ax.set_visible(False)

    # Single legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.1, 0.5))

    fig.suptitle(title)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved: {output_file}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Create PCA plots from PLINK2 eigenvector files."
    )
    parser.add_argument(
        "eigenvec_file",
        type=Path,
        help="PLINK2 .eigenvec file",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("docs/integrated_call_samples_v3.20130502.ALL.panel"),
        help="1000 Genomes panel file",
    )
    parser.add_argument(
        "--pop-col",
        type=str,
        default="pop",
        choices=["pop", "super_pop"],
        help="Population column to color by",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file (PNG or PDF)",
    )
    parser.add_argument(
        "--grid",
        action="store_true",
        help="Create grid of PC pairs instead of single plot",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="PCA Plot",
        help="Plot title",
    )

    args = parser.parse_args()

    # Default output
    if args.output is None:
        args.output = args.eigenvec_file.with_suffix(".png")

    # Load data
    df = load_pca_data(
        eigenvec_file=args.eigenvec_file,
        panel_file=args.panel,
        population_col=args.pop_col,
    )

    print(f"Loaded {len(df)} samples with {df[args.pop_col].nunique()} populations")

    # Create plot
    if args.grid:
        plot_pca_grid(
            df=df,
            population_col=args.pop_col,
            output_file=args.output,
            title=args.title,
        )
    else:
        plot_pca(
            df=df,
            population_col=args.pop_col,
            output_file=args.output,
            title=args.title,
        )


if __name__ == "__main__":
    main()
