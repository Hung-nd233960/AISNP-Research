"""File I/O helpers: directories, panels, BED conversion, reports."""

import csv
from pathlib import Path
from typing import List, Tuple, Union

import pandas as pd

from config import PATHS


def ensure_dir(path: Union[str, Path]) -> Path:
    """Create directory (and parents) if it doesn't exist, then return it."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_project_root() -> Path:
    return PATHS.ROOT


def save_report(report_text: str, output_path: Union[str, Path], verbose: bool = True) -> None:
    """Write a text report to disk, creating parent directories as needed."""
    ensure_dir(Path(output_path).parent)
    with open(output_path, "w") as f:
        f.write(report_text)
    if verbose:
        print(f"Report saved: {output_path}")


def read_panel(path: Union[str, Path]) -> pd.DataFrame:
    """Read a tab-delimited 1000 Genomes panel file."""
    return pd.read_csv(path, sep="\t")


def extract_population_samples(
    panel_df: pd.DataFrame,
    populations: List[str],
    output_csv: Union[str, Path],
    output_list: Union[str, Path],
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Filter panel to the given populations and write two files:
      output_csv  — full metadata (no header, PLINK-compatible)
      output_list — sample IDs only
    """
    filtered = panel_df[panel_df["pop"].isin(populations)].copy()
    filtered.to_csv(output_csv, index=False, header=False)
    samples = filtered["sample"].tolist()
    pd.DataFrame({"sample": samples}).to_csv(output_list, index=False, header=False)
    print(f"Extracted {len(samples)} samples from {populations}")
    return filtered, samples


def variant_to_bed(variant: str) -> Tuple[str, int, int, str]:
    """
    Parse a variant ID like '12:58124534[b37]C,G' into a BED tuple
    (chrom, 0-based-start, end, name).
    """
    try:
        chrom_pos, allele_info = variant.split("]")
        chrom_pos = chrom_pos.replace("[b37", "")
        chrom, pos_str = chrom_pos.split(":")
        pos = int(pos_str)
        return chrom, pos - 1, pos, allele_info
    except Exception as e:
        raise ValueError(f"Cannot parse variant '{variant}': {e}")


def variants_to_bed_file(
    snp_file: Union[str, Path],
    bed_file: Union[str, Path],
    verbose: bool = True,
) -> int:
    """Convert a one-variant-per-line SNP file to BED format. Returns count."""
    count = 0
    with open(snp_file) as sf, open(bed_file, "w", newline="") as bf:
        writer = csv.writer(bf, delimiter="\t")
        for line in sf:
            v = line.strip()
            if v:
                writer.writerow(variant_to_bed(v))
                count += 1
    if verbose:
        print(f"Wrote {count} variants → {bed_file}")
    return count
