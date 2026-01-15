#!/usr/bin/env python3
"""
Convert CSV SNP data to BED format for bcftools region extraction.

Input CSV must have columns: Chromosome, AmpliconStartPosition, AmpliconEndPosition, Locus
"""

import argparse
import csv
from pathlib import Path


def csv_to_bed(
    input_csv: Path,
    output_bed: Path | None = None,
    chrom_prefix: str = "",
    verbose: bool = True,
) -> Path:
    """
    Convert CSV SNP coordinates to BED format.

    Args:
        input_csv: Input CSV file with SNP coordinates.
        output_bed: Output BED file (default: same name with .bed extension).
        chrom_prefix: Prefix to add to chromosome (e.g., "chr").
        verbose: Print progress.

    Returns:
        Path to output BED file.
    """
    if output_bed is None:
        output_bed = input_csv.with_suffix(".bed")

    count = 0
    with open(input_csv, "r") as csvfile, open(output_bed, "w", newline="") as bedfile:
        reader = csv.DictReader(csvfile)
        writer = csv.writer(bedfile, delimiter="\t")

        for row in reader:
            chrom = f"{chrom_prefix}{row['Chromosome']}"
            # BED format: 0-based start, 1-based end
            start = int(row["AmpliconStartPosition"]) - 1
            end = int(row["AmpliconEndPosition"])
            name = row.get("Locus", row.get("rsid", f"snp_{count}"))

            writer.writerow([chrom, start, end, name])
            count += 1

    if verbose:
        print(f"Converted {count} SNPs from {input_csv} to {output_bed}")

    return output_bed


def rsid_csv_to_bed(
    input_csv: Path,
    output_bed: Path | None = None,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    rsid_col: str = "rsid",
    verbose: bool = True,
) -> Path:
    """
    Convert CSV with rsIDs and positions to BED format.

    Args:
        input_csv: Input CSV file.
        output_bed: Output BED file.
        chrom_col: Column name for chromosome.
        pos_col: Column name for position.
        rsid_col: Column name for rsID.
        verbose: Print progress.

    Returns:
        Path to output BED file.
    """
    if output_bed is None:
        output_bed = input_csv.with_suffix(".bed")

    count = 0
    with open(input_csv, "r") as csvfile, open(output_bed, "w", newline="") as bedfile:
        reader = csv.DictReader(csvfile)
        writer = csv.writer(bedfile, delimiter="\t")

        for row in reader:
            chrom = row[chrom_col]
            pos = int(row[pos_col])
            rsid = row.get(rsid_col, f"snp_{count}")

            # BED: 0-based start, 1-based end (for single position)
            writer.writerow([chrom, pos - 1, pos, rsid])
            count += 1

    if verbose:
        print(f"Converted {count} SNPs to {output_bed}")

    return output_bed


def main():
    parser = argparse.ArgumentParser(description="Convert CSV SNP data to BED format.")
    parser.add_argument(
        "input_csv",
        type=Path,
        help="Input CSV file with SNP coordinates",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output BED file (default: <input>.bed)",
    )
    parser.add_argument(
        "--chrom-prefix",
        type=str,
        default="",
        help="Prefix to add to chromosome (e.g., 'chr')",
    )
    parser.add_argument(
        "--format",
        choices=["amplicon", "rsid"],
        default="amplicon",
        help="Input format: 'amplicon' (with start/end) or 'rsid' (with single position)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress output",
    )

    args = parser.parse_args()

    if args.format == "amplicon":
        csv_to_bed(
            input_csv=args.input_csv,
            output_bed=args.output,
            chrom_prefix=args.chrom_prefix,
            verbose=not args.quiet,
        )
    else:
        rsid_csv_to_bed(
            input_csv=args.input_csv,
            output_bed=args.output,
            verbose=not args.quiet,
        )


if __name__ == "__main__":
    main()
