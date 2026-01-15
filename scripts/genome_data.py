#!/usr/bin/env python3
"""
Extract genotypes from 1000 Genomes VCF files for specific SNP loci.

Uses cyvcf2 for efficient VCF parsing and outputs a sample x SNP genotype matrix.
"""

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd
from cyvcf2 import VCF


def load_loci(bed_file: Path) -> pd.DataFrame:
    """
    Load BED loci into a DataFrame.

    Args:
        bed_file: Path to BED file with columns: chrom, start, end, rsid

    Returns:
        DataFrame with loci information.
    """
    loci = []
    with open(bed_file, "r") as bed:
        for line in bed:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom, start, end, rsid = parts[:4]
                loci.append((chrom, int(start), int(end), rsid))

    return pd.DataFrame(loci, columns=["chrom", "start", "end", "rsid"])


def load_samples_from_vcf(vcf_path: str) -> List[str]:
    """
    Extract sample names from VCF file.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of sample names.
    """
    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    vcf.close()
    return samples


def process_variant_window(
    vcf: VCF,
    chrom: str,
    start: int,
    end: int,
    rsid: str,
    df: pd.DataFrame,
    samples: List[str],
    verbose: bool = True,
) -> bool:
    """
    Process one variant window and extract genotypes.

    Genotype encoding: 0 = hom ref, 1 = het, 2 = hom alt

    Args:
        vcf: cyvcf2 VCF object.
        chrom: Chromosome identifier.
        start: Start position (0-based).
        end: End position.
        rsid: SNP identifier.
        df: DataFrame to update with genotypes.
        samples: List of sample names.
        verbose: Print progress.

    Returns:
        True if variant was found and processed.
    """
    if verbose:
        print(f"  Looking for {rsid} at {chrom}:{start}-{end}...")

    found = False
    region = f"{chrom}:{start}-{end}"

    for variant in vcf(region):
        if variant is None:
            continue

        # Check position matches expected
        if variant.POS != start + 1:
            if verbose:
                print(f"    Skipped: variant at {variant.POS} != expected {start+1}")
            continue

        # Skip multi-allelic or complex variants
        if len(variant.ALT) > 1 or len(variant.ALT[0]) > 1 or len(variant.REF) > 1:
            if verbose:
                print(
                    f"    Skipped: complex variant REF={variant.REF}, ALT={variant.ALT}"
                )
            continue

        found = True
        if verbose:
            print(
                f"    Found: {rsid} at {variant.CHROM}:{variant.POS} (REF={variant.REF}, ALT={variant.ALT[0]})"
            )

        # Extract genotypes
        genotypes = variant.genotypes

        for sample, gt in zip(samples, genotypes):
            allele1, allele2 = gt[0], gt[1]

            # Handle missing genotypes
            if allele1 < 0 or allele2 < 0:
                df.loc[df["sample"] == sample, rsid] = -1
            else:
                # Sum of alt alleles: 0=hom ref, 1=het, 2=hom alt
                df.loc[df["sample"] == sample, rsid] = allele1 + allele2

        break  # Only process first matching variant

    if not found and verbose:
        print(f"    Not found: {rsid}")

    return found


def process_chromosome(
    chrom: str,
    loci_df: pd.DataFrame,
    vcf_template: str,
    df: pd.DataFrame,
    samples: List[str],
    verbose: bool = True,
) -> int:
    """
    Process all loci in a single chromosome.

    Args:
        chrom: Chromosome identifier.
        loci_df: DataFrame of loci for this chromosome.
        vcf_template: Template string for VCF file paths (use {chrom} placeholder).
        df: DataFrame to update with genotypes.
        samples: List of sample names.
        verbose: Print progress.

    Returns:
        Number of variants found.
    """
    vcf_path = vcf_template.format(chrom=chrom)

    if not Path(vcf_path).exists():
        if verbose:
            print(f"Warning: VCF not found: {vcf_path}")
        return 0

    if verbose:
        print(f"Processing chromosome {chrom} ({len(loci_df)} loci)...")

    vcf = VCF(vcf_path)
    found_count = 0

    for _, locus in loci_df.iterrows():
        ok = process_variant_window(
            vcf=vcf,
            chrom=chrom,
            start=int(locus["start"]),
            end=int(locus["end"]),
            rsid=str(locus["rsid"]),
            df=df,
            samples=samples,
            verbose=verbose,
        )
        if ok:
            found_count += 1

    vcf.close()
    return found_count


def extract_genotypes(
    bed_file: Path,
    vcf_template: str,
    output_file: Path,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Extract genotypes from VCF files for loci in BED file.

    Args:
        bed_file: Path to BED file with target loci.
        vcf_template: Template for VCF paths (use {chrom} placeholder).
        output_file: Output CSV path.
        verbose: Print progress.

    Returns:
        DataFrame with sample x SNP genotype matrix.
    """
    # Load loci
    loci = load_loci(bed_file)
    if verbose:
        print(f"Loaded {len(loci)} loci from {bed_file}")

    # Get sample names from first chromosome VCF
    first_chrom = str(loci.iloc[0]["chrom"])
    first_vcf = vcf_template.format(chrom=first_chrom)
    samples = load_samples_from_vcf(first_vcf)
    if verbose:
        print(f"Found {len(samples)} samples in VCF")

    # Initialize DataFrame
    df = pd.DataFrame({"sample": samples})

    # Group loci by chromosome
    loci_by_chrom: Dict[str, pd.DataFrame] = {
        chrom: loci[loci["chrom"] == chrom] for chrom in loci["chrom"].unique()
    }

    # Process each chromosome
    total_found = 0
    for chrom, chrom_loci in loci_by_chrom.items():
        total_found += process_chromosome(
            chrom=str(chrom),
            loci_df=chrom_loci,
            vcf_template=vcf_template,
            df=df,
            samples=samples,
            verbose=verbose,
        )

    if verbose:
        print(f"\nFound {total_found}/{len(loci)} variants")

    # Save output
    df.to_csv(output_file, index=False)
    if verbose:
        print(f"Saved to: {output_file}")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Extract genotypes from 1000 Genomes VCF files."
    )
    parser.add_argument(
        "bed_file",
        type=Path,
        help="BED file with target loci (chrom, start, end, rsid)",
    )
    parser.add_argument(
        "--vcf-template",
        type=str,
        default="1000genomes/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        help="Template for VCF paths (use {chrom} placeholder)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output CSV file",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress output",
    )

    args = parser.parse_args()

    extract_genotypes(
        bed_file=args.bed_file,
        vcf_template=args.vcf_template,
        output_file=args.output,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
