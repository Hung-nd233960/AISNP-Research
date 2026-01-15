"""
HARD FILTERS: Strict quality control filters.

These filters are applied consistently across all datasets and should
not be modified without strong justification.

Filter Order:
1. SNP-only + Biallelic filter
2. Sample subsetting (optional)
3. MAF (Minor Allele Frequency) filter
4. Call rate filter

Rationale:
- SNP-only: Focuses analysis on single nucleotide polymorphisms
- Biallelic: Removes multi-allelic sites that complicate analysis
- MAF: Removes very rare variants (singletons/doubletons)
- Call rate: Ensures high genotyping quality
"""

import os
from pathlib import Path
from typing import Optional, List, Union

# Import configuration and utilities
try:
    from config import PATHS, HARD_FILTERS, PLINK
    from utils import run_plink2_command, count_variants, count_samples
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from config import PATHS, HARD_FILTERS, PLINK
    from utils import run_plink2_command, count_variants, count_samples


# =============================================================================
# HARD FILTER 1: SNP-only and Biallelic
# =============================================================================


def filter_snp_biallelic(
    input_vcf: Optional[str] = None,
    input_pfile: Optional[str] = None,
    output_pfile: str = None,
    keep_samples: Optional[str] = None,
    max_alleles: int = HARD_FILTERS.MAX_ALLELES,
    verbose: bool = True,
) -> str:
    """
    HARD FILTER 1: SNP-only and biallelic filter.

    Keeps only:
    - Single nucleotide polymorphisms (removes indels, CNVs, etc.)
    - Variants with at most 2 alleles (biallelic)

    Args:
        input_vcf: Path to input VCF file (mutually exclusive with input_pfile)
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        keep_samples: Path to sample list file (optional)
        max_alleles: Maximum number of alleles (default: 2)
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if input_vcf is None and input_pfile is None:
        raise ValueError("Must provide either input_vcf or input_pfile")

    if output_pfile is None:
        output_pfile = str(PATHS.PLINK_SNP_FILTERED)

    args = []

    # Input specification
    if input_vcf:
        args.extend(["--vcf", input_vcf])
    else:
        args.extend(["--pfile", input_pfile])

    # SNP-only filter
    args.append("--snps-only")

    # Biallelic filter
    args.extend(["--max-alleles", str(max_alleles)])

    # Sample subsetting (optional)
    if keep_samples:
        args.extend(["--keep", keep_samples])

    # Output
    args.extend(["--make-pgen", "--out", output_pfile])

    if verbose:
        print("=" * 60)
        print("HARD FILTER 1: SNP-only and Biallelic")
        print("=" * 60)
        print(f"  SNP-only: True")
        print(f"  Max alleles: {max_alleles}")
        if keep_samples:
            print(f"  Sample filter: {keep_samples}")

    run_plink2_command(args)

    if verbose:
        n_variants = count_variants(output_pfile)
        n_samples = count_samples(output_pfile)
        print(f"\nOutput: {output_pfile}")
        print(f"  Variants: {n_variants}")
        print(f"  Samples: {n_samples}")

    return output_pfile


# =============================================================================
# HARD FILTER 2: Minor Allele Frequency (MAF)
# =============================================================================


def filter_maf(
    input_pfile: str,
    output_pfile: str = None,
    min_af: float = HARD_FILTERS.MIN_AF,
    verbose: bool = True,
) -> str:
    """
    HARD FILTER 2: Minor Allele Frequency filter.

    Removes variants with very low allele frequency.
    Default threshold: 1/612 = 0.0016 (at least 1 allele in 306 samples)

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        min_af: Minimum allele frequency threshold
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = str(PATHS.PLINK_MAF_FILTERED)

    args = [
        "--pfile",
        input_pfile,
        "--min-af",
        str(min_af),
        "--make-pgen",
        "--out",
        output_pfile,
    ]

    if verbose:
        print("=" * 60)
        print("HARD FILTER 2: Minor Allele Frequency (MAF)")
        print("=" * 60)
        print(f"  Min AF threshold: {min_af}")

        # Count input variants
        n_input = count_variants(input_pfile)
        print(f"  Input variants: {n_input}")

    run_plink2_command(args)

    if verbose:
        n_output = count_variants(output_pfile)
        n_removed = n_input - n_output if n_input > 0 else 0
        print(f"\nOutput: {output_pfile}")
        print(f"  Output variants: {n_output}")
        print(
            f"  Removed: {n_removed} ({n_removed/n_input*100:.2f}%)"
            if n_input > 0
            else ""
        )

    return output_pfile


# =============================================================================
# HARD FILTER 3: Call Rate (Genotyping Completeness)
# =============================================================================


def filter_call_rate(
    input_pfile: str,
    output_pfile: str = None,
    min_call_rate: float = HARD_FILTERS.MIN_CALL_RATE,
    verbose: bool = True,
) -> str:
    """
    HARD FILTER 3: Call rate filter.

    Removes variants with too many missing genotypes.
    Default: Keep variants with >= 95% genotyping rate.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        min_call_rate: Minimum call rate (0.95 = 95%)
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = input_pfile + "_CR_filtered"

    # PLINK2 uses --geno for max missing rate (inverse of call rate)
    max_missing = 1.0 - min_call_rate

    args = [
        "--pfile",
        input_pfile,
        "--geno",
        str(max_missing),
        "--make-pgen",
        "--out",
        output_pfile,
    ]

    if verbose:
        print("=" * 60)
        print("HARD FILTER 3: Call Rate")
        print("=" * 60)
        print(f"  Min call rate: {min_call_rate} ({min_call_rate*100:.1f}%)")
        print(f"  Max missing rate: {max_missing}")

        n_input = count_variants(input_pfile)
        print(f"  Input variants: {n_input}")

    run_plink2_command(args)

    if verbose:
        n_output = count_variants(output_pfile)
        n_removed = n_input - n_output if n_input > 0 else 0
        print(f"\nOutput: {output_pfile}")
        print(f"  Output variants: {n_output}")
        print(f"  Removed: {n_removed}" if n_input > 0 else "")

    return output_pfile


# =============================================================================
# Frequency Statistics Calculation
# =============================================================================


def calculate_frequencies(
    input_pfile: str,
    output_prefix: str = None,
    verbose: bool = True,
) -> str:
    """
    Calculate allele frequency statistics.
    Output goes to .afreq file for QC inspection.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        verbose: Print progress

    Returns:
        Output prefix (adds .afreq extension)
    """
    if output_prefix is None:
        output_prefix = input_pfile + "_info"

    args = [
        "--pfile",
        input_pfile,
        "--freq",
        "--out",
        output_prefix,
    ]

    if verbose:
        print(f"Calculating frequency statistics: {output_prefix}.afreq")

    run_plink2_command(args)

    return output_prefix


# =============================================================================
# Pipeline: Apply All Hard Filters
# =============================================================================


def apply_all_hard_filters(
    input_vcf: str = None,
    input_pfile: str = None,
    keep_samples: str = None,
    output_prefix: str = None,
    min_af: float = HARD_FILTERS.MIN_AF,
    min_call_rate: float = HARD_FILTERS.MIN_CALL_RATE,
    calculate_stats: bool = True,
    verbose: bool = True,
) -> str:
    """
    Apply all hard filters in sequence.

    Pipeline:
    1. SNP-only + Biallelic filter
    2. MAF filter
    3. Call rate filter (if data has missing genotypes)

    Args:
        input_vcf: Path to input VCF file
        input_pfile: Path to input PLINK2 pfile prefix
        keep_samples: Path to sample list file
        output_prefix: Base prefix for all outputs
        min_af: Minimum allele frequency
        min_call_rate: Minimum call rate
        calculate_stats: Calculate frequency stats after each step
        verbose: Print progress

    Returns:
        Final output pfile prefix
    """
    if output_prefix is None:
        output_prefix = str(PATHS.OUTPUT_DIR / "hard_filtered")

    print("=" * 60)
    print("APPLYING ALL HARD FILTERS")
    print("=" * 60)

    # Step 1: SNP-only + Biallelic
    step1_output = f"{output_prefix}_step1_snp_biallelic"
    filter_snp_biallelic(
        input_vcf=input_vcf,
        input_pfile=input_pfile,
        output_pfile=step1_output,
        keep_samples=keep_samples,
        verbose=verbose,
    )

    if calculate_stats:
        calculate_frequencies(step1_output, verbose=verbose)

    # Step 2: MAF filter
    step2_output = f"{output_prefix}_step2_maf"
    filter_maf(
        input_pfile=step1_output,
        output_pfile=step2_output,
        min_af=min_af,
        verbose=verbose,
    )

    if calculate_stats:
        calculate_frequencies(step2_output, verbose=verbose)

    # Step 3: Call rate filter (optional - skip if call rate is 100%)
    # For 1000 Genomes data, call rate is typically 100%
    # Uncomment below to enable:
    # step3_output = f"{output_prefix}_step3_callrate"
    # filter_call_rate(
    #     input_pfile=step2_output,
    #     output_pfile=step3_output,
    #     min_call_rate=min_call_rate,
    #     verbose=verbose,
    # )

    final_output = step2_output

    print("\n" + "=" * 60)
    print("HARD FILTERING COMPLETE")
    print("=" * 60)
    print(f"Final output: {final_output}")
    print(f"  Variants: {count_variants(final_output)}")
    print(f"  Samples: {count_samples(final_output)}")

    return final_output


# =============================================================================
# CLI Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Apply hard filters to genetic data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter VCF file with sample subsetting
  python hard_filters.py --vcf data.vcf.gz --samples samples.txt --output filtered
  
  # Filter existing pfile with custom MAF
  python hard_filters.py --pfile input --min-af 0.01 --output filtered
        """,
    )

    parser.add_argument("--vcf", help="Input VCF file")
    parser.add_argument("--pfile", help="Input PLINK2 pfile prefix")
    parser.add_argument("--samples", help="Sample list file")
    parser.add_argument("--output", "-o", required=True, help="Output prefix")
    parser.add_argument(
        "--min-af",
        type=float,
        default=HARD_FILTERS.MIN_AF,
        help=f"Minimum allele frequency (default: {HARD_FILTERS.MIN_AF})",
    )
    parser.add_argument(
        "--min-call-rate",
        type=float,
        default=HARD_FILTERS.MIN_CALL_RATE,
        help=f"Minimum call rate (default: {HARD_FILTERS.MIN_CALL_RATE})",
    )
    parser.add_argument(
        "--no-stats", action="store_true", help="Skip frequency statistics calculation"
    )

    args = parser.parse_args()

    apply_all_hard_filters(
        input_vcf=args.vcf,
        input_pfile=args.pfile,
        keep_samples=args.samples,
        output_prefix=args.output,
        min_af=args.min_af,
        min_call_rate=args.min_call_rate,
        calculate_stats=not args.no_stats,
    )
