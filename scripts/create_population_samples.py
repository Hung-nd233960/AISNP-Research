#!/usr/bin/env python3
"""
Create population sample CSV and sample list files for a given configuration.

- Filters the 1000 Genomes panel to RAW_SUBPOPS.
- Applies POP_MERGE_MAP (if any) to create merged population labels.
- Writes:
  * <prefix>_subpopulation_samples.csv (no header)
  * <prefix>_subpopulation_samples_list.csv (sample IDs only)

Example:
  python scripts/create_population_samples.py --config sea_jpt_cn
"""

import argparse
from pathlib import Path
import pandas as pd

import config as cfg
from utils import read_panel


def create_population_samples(panel_file: Path) -> tuple[Path, Path]:
    """
    Create population sample files for the active config.

    Args:
        panel_file: Path to integrated panel file.

    Returns:
        Tuple of (samples_csv_path, samples_list_path)
    """
    panel_df = read_panel(panel_file)

    # Filter to raw subpopulations in the panel
    raw_pops = list(cfg.POPULATIONS.RAW_SUBPOPS)
    filtered_df = panel_df[panel_df["pop"].isin(raw_pops)].copy()

    # Apply merge map if present
    if cfg.POPULATIONS.POP_MERGE_MAP:
        filtered_df["pop"] = filtered_df["pop"].map(cfg.POPULATIONS.POP_MERGE_MAP)

    # Prepare output paths
    samples_csv = cfg.PATHS.get_absolute(cfg.PATHS.EAS_SAMPLES_CSV)
    samples_list = cfg.PATHS.get_absolute(cfg.PATHS.EAS_SAMPLES_LIST)
    samples_csv.parent.mkdir(parents=True, exist_ok=True)
    samples_list.parent.mkdir(parents=True, exist_ok=True)

    # Save full metadata (no header for PLINK compatibility)
    filtered_df[["sample", "pop", "super_pop"]].to_csv(
        samples_csv, index=False, header=False
    )

    # Save sample list only
    filtered_df[["sample"]].to_csv(samples_list, index=False, header=False)

    print(f"Saved samples CSV: {samples_csv}")
    print(f"Saved samples list: {samples_list}")
    print(f"Total samples: {len(filtered_df)}")
    print(f"Population counts:\n{filtered_df['pop'].value_counts().to_string()}")

    return samples_csv, samples_list


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create population sample files for a given config."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="sea_jpt_cn",
        choices=["khv_jpt_chb", "sea_jpt_cn"],
        help="Population config name",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=None,
        help="Override panel file path (optional)",
    )

    args = parser.parse_args()

    # Activate config (updates PATHS and POPULATIONS)
    cfg.set_population_config(args.config)

    panel_file = (
        args.panel if args.panel else cfg.PATHS.get_absolute(cfg.PATHS.PANEL_FILE)
    )

    create_population_samples(panel_file)


if __name__ == "__main__":
    main()
