"""
FST-based Variant Selection Pipeline.

FST (Fixation Index) measures genetic differentiation between populations.
High FST variants are "population-informative" and useful for:
- Ancestry inference
- Population classification
- Panel development

Workflow:
1. Calculate pairwise FST between populations
2. Analyze FST distribution
3. Select top-N variants per population pair
4. Merge and deduplicate variant lists
5. Extract selected variants for downstream analysis
"""

import os
import csv
from pathlib import Path
from typing import List, Tuple, Dict, Union, Optional
import pandas as pd
import numpy as np

# Import configuration and utilities
try:
    from config import PATHS, SITUATIONAL_FILTERS, POPULATIONS
    from utils import (
        run_plink2_command,
        run_bcftools_command,
        analyze_fst,
        merge_fst_top_variants,
        variant_to_bed,
        variants_to_bed_file,
        ensure_dir,
    )
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from config import PATHS, SITUATIONAL_FILTERS, POPULATIONS
    from utils import (
        run_plink2_command,
        run_bcftools_command,
        analyze_fst,
        merge_fst_top_variants,
        variant_to_bed,
        variants_to_bed_file,
        ensure_dir,
    )


# =============================================================================
# PSAM File Management
# =============================================================================


def add_population_to_psam(
    psam_path: Union[str, Path],
    samples_csv: Union[str, Path],
    output_path: Union[str, Path] = None,
    verbose: bool = True,
) -> str:
    """
    Add population labels to a PLINK2 .psam file.

    Args:
        psam_path: Path to input .psam file
        samples_csv: Path to CSV with sample,pop,super_pop columns
        output_path: Path for output .psam file (default: overwrite input)
        verbose: Print progress

    Returns:
        Output path
    """
    if output_path is None:
        output_path = psam_path

    # Read sample metadata
    samples_df = pd.read_csv(
        samples_csv,
        header=None,
        names=["sample", "pop", "super_pop"],
    )

    # Read psam file
    psam_df = pd.read_csv(psam_path, sep="\t")

    # Remove SEX column if present (we'll add pop instead)
    if "SEX" in psam_df.columns:
        psam_df = psam_df.drop(columns=["SEX"])

    # Remove super_pop from samples (keep only pop)
    samples_df = samples_df[["sample", "pop"]]

    # Merge on sample ID
    merged_df = pd.merge(
        psam_df,
        samples_df,
        left_on="#IID",
        right_on="sample",
        how="left",
    )
    merged_df = merged_df.drop(columns=["sample"])

    if verbose:
        print(f"Added population labels to psam file")
        print(f"  Samples: {len(merged_df)}")
        print(f"  Populations: {merged_df['pop'].value_counts().to_dict()}")

    # Save
    merged_df.to_csv(output_path, sep="\t", index=False)

    return str(output_path)


# =============================================================================
# FST Calculation
# =============================================================================


def calculate_pairwise_fst(
    input_pfile: str,
    output_prefix: str = None,
    pop_column: str = "pop",
    verbose: bool = True,
) -> str:
    """
    Calculate pairwise FST between all populations.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        pop_column: Column name in .psam containing population labels
        verbose: Print progress

    Returns:
        Output prefix (creates multiple .fst.var files)
    """
    if output_prefix is None:
        output_prefix = str(PATHS.FST_RESULTS)

    args = [
        "--pfile",
        input_pfile,
        "--fst",
        pop_column,
        "report-variants",
        "--out",
        output_prefix,
    ]

    if verbose:
        print("=" * 60)
        print("FST CALCULATION")
        print("=" * 60)
        print(f"  Input: {input_pfile}")
        print(f"  Population column: {pop_column}")

    run_plink2_command(args)

    if verbose:
        print(f"\nFST files created: {output_prefix}.*")

    return output_prefix


# =============================================================================
# FST Analysis and Variant Selection
# =============================================================================


def find_fst_files(
    fst_prefix: str,
    verbose: bool = True,
) -> List[str]:
    """
    Find all FST variant files from a prefix.

    Args:
        fst_prefix: FST output prefix
        verbose: Print progress

    Returns:
        List of FST file paths
    """
    base_dir = Path(fst_prefix).parent
    base_name = Path(fst_prefix).name

    fst_files = list(base_dir.glob(f"{base_name}.*.fst.var"))

    if verbose:
        print(f"Found {len(fst_files)} FST files:")
        for f in fst_files:
            print(f"  - {f.name}")

    return [str(f) for f in fst_files]


def analyze_all_fst_files(
    fst_prefix: str,
    top_n: int = SITUATIONAL_FILTERS.FST_TOP_N,
    plot: bool = False,
    verbose: bool = True,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Analyze all FST files from pairwise comparisons.

    Args:
        fst_prefix: FST output prefix
        top_n: Number of top variants per comparison
        plot: Generate plots
        verbose: Print progress

    Returns:
        all_fst: Dict of full FST dataframes keyed by comparison
        top_fst: Dict of top-N dataframes keyed by comparison
    """
    fst_files = find_fst_files(fst_prefix, verbose=verbose)

    all_fst = {}
    top_fst = {}

    for fst_file in fst_files:
        # Extract comparison name from filename
        # e.g., "EAS_FST_RESULTS.CHB.JPT.fst.var" -> "CHB.JPT"
        name_parts = Path(fst_file).stem.replace(".fst", "").split(".")
        comparison = ".".join(name_parts[-2:]) if len(name_parts) >= 2 else fst_file

        if verbose:
            print(f"\nAnalyzing: {comparison}")

        df, top_df = analyze_fst(fst_file, top_n=top_n, plot=plot, verbose=verbose)
        all_fst[comparison] = df
        top_fst[comparison] = top_df

    return all_fst, top_fst


def select_top_fst_variants(
    fst_prefix: str,
    output_file: str = None,
    top_n: int = SITUATIONAL_FILTERS.FST_TOP_N,
    verbose: bool = True,
) -> Tuple[pd.DataFrame, str]:
    """
    Select top FST variants from all pairwise comparisons.

    Args:
        fst_prefix: FST output prefix
        output_file: Path to output variant list file
        top_n: Number of top variants per comparison
        verbose: Print progress

    Returns:
        merged_df: DataFrame of unique top variants
        output_file: Path to saved variant list
    """
    if output_file is None:
        output_file = str(PATHS.TOP_SNPS_FILE)

    # Analyze all FST files
    all_fst, top_fst = analyze_all_fst_files(fst_prefix, top_n=top_n, verbose=verbose)

    # Merge top variants from all comparisons
    top_dfs = list(top_fst.values())
    merged_df = merge_fst_top_variants(*top_dfs, id_col="ID")

    if verbose:
        print(f"\n{'='*60}")
        print("FST VARIANT SELECTION SUMMARY")
        print(f"{'='*60}")
        print(f"  Top N per comparison: {top_n}")
        print(f"  Number of comparisons: {len(top_fst)}")
        print(f"  Total unique variants: {len(merged_df)}")

    # Save variant IDs
    ensure_dir(Path(output_file).parent)
    with open(output_file, "w") as f:
        for snp_id in merged_df["ID"]:
            f.write(f"{snp_id}\n")

    if verbose:
        print(f"\nVariant list saved: {output_file}")

    return merged_df, output_file


# =============================================================================
# Variant Extraction
# =============================================================================


def extract_fst_variants(
    input_pfile: str,
    variant_list: str,
    output_pfile: str = None,
    verbose: bool = True,
) -> str:
    """
    Extract selected FST variants from pfile.

    Args:
        input_pfile: Input pfile with unique variant IDs
        variant_list: Path to variant list file
        output_pfile: Output pfile prefix
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = str(PATHS.FST_FILTERED)

    args = [
        "--pfile",
        input_pfile,
        "--extract",
        variant_list,
        "--make-pgen",
        "--out",
        output_pfile,
    ]

    if verbose:
        print(f"Extracting FST-selected variants...")

    run_plink2_command(args)

    return output_pfile


def extract_variants_vcf(
    input_vcf: str,
    variant_list: str,
    sample_list: str,
    output_vcf: str,
    verbose: bool = True,
) -> str:
    """
    Extract variants from VCF using bcftools.

    Args:
        input_vcf: Input VCF file
        variant_list: Path to variant list in SNP ID format
        sample_list: Path to sample list
        output_vcf: Output VCF path
        verbose: Print progress

    Returns:
        Output VCF path
    """  # First convert variant list to BED format
    bed_file = variant_list.replace(".txt", ".bed")
    variants_to_bed_file(variant_list, bed_file, verbose=verbose)

    args = [
        "view",
        "-S",
        sample_list,
        "-R",
        bed_file,
        "-Oz",
        "-o",
        output_vcf,
        input_vcf,
    ]

    if verbose:
        print(f"Extracting variants from VCF...")

    run_bcftools_command(args)

    return output_vcf


# =============================================================================
# Full FST Selection Pipeline
# =============================================================================


def run_fst_selection_pipeline(
    input_pfile: str = None,
    samples_csv: str = None,
    top_n: int = SITUATIONAL_FILTERS.FST_TOP_N,
    output_dir: str = None,
    verbose: bool = True,
) -> Dict[str, str]:
    """
    Run complete FST-based variant selection pipeline.

    Pipeline:
    1. Add population labels to psam
    2. Calculate pairwise FST
    3. Analyze FST and select top variants
    4. Extract selected variants

    Args:
        input_pfile: Input pfile (should be LD-pruned)
        samples_csv: Path to sample metadata CSV
        top_n: Number of top FST variants per comparison
        output_dir: Output directory
        verbose: Print progress

    Returns:
        Dictionary with output file paths
    """
    if input_pfile is None:
        input_pfile = str(PATHS.PLINK_LD_PRUNED)
    if samples_csv is None:
        samples_csv = str(PATHS.EAS_SAMPLES_CSV)
    if output_dir is None:
        output_dir = str(PATHS.OUTPUT_DIR)

    print("=" * 60)
    print("FST SELECTION PIPELINE")
    print("=" * 60)

    outputs = {}

    # Step 1: Add population labels
    psam_path = f"{input_pfile}.psam"
    psam_with_pop = f"{input_pfile}_with_pop.psam"
    add_population_to_psam(psam_path, samples_csv, psam_with_pop, verbose=verbose)

    # Replace original psam with the one with pop labels
    os.replace(psam_with_pop, psam_path)
    outputs["psam"] = psam_path

    # Step 2: Calculate FST
    fst_prefix = f"{output_dir}/EAS_FST_RESULTS"
    calculate_pairwise_fst(input_pfile, fst_prefix, verbose=verbose)
    outputs["fst_prefix"] = fst_prefix

    # Step 3: Select top variants
    top_snps_file = f"{output_dir}/top_snps.txt"
    merged_df, top_snps_file = select_top_fst_variants(
        fst_prefix,
        top_snps_file,
        top_n=top_n,
        verbose=verbose,
    )
    outputs["top_snps"] = top_snps_file
    outputs["top_snps_df"] = merged_df

    # Step 4: Convert to BED format
    bed_file = f"{output_dir}/top_snps.bed"
    variants_to_bed_file(top_snps_file, bed_file, verbose=verbose)
    outputs["bed_file"] = bed_file

    # Step 5: Extract variants
    # Use pfile with unique IDs for extraction
    unique_ids_pfile = str(PATHS.PLINK_UNIQUE_IDS)
    fst_filtered = f"{output_dir}/FST_FILTERED"
    extract_fst_variants(unique_ids_pfile, top_snps_file, fst_filtered, verbose=verbose)
    outputs["fst_filtered"] = fst_filtered

    print("\n" + "=" * 60)
    print("FST SELECTION COMPLETE")
    print("=" * 60)
    print(f"  Top SNPs file: {top_snps_file}")
    print(f"  BED file: {bed_file}")
    print(f"  Filtered pfile: {fst_filtered}")

    return outputs


# =============================================================================
# CLI Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="FST-based variant selection for population classification"
    )

    parser.add_argument("--pfile", help="Input PLINK2 pfile prefix")
    parser.add_argument("--samples", help="Sample metadata CSV")
    parser.add_argument(
        "--top-n",
        type=int,
        default=SITUATIONAL_FILTERS.FST_TOP_N,
        help=f"Top N variants per comparison (default: {SITUATIONAL_FILTERS.FST_TOP_N})",
    )
    parser.add_argument("--output-dir", "-o", help="Output directory")

    args = parser.parse_args()

    run_fst_selection_pipeline(
        input_pfile=args.pfile,
        samples_csv=args.samples,
        top_n=args.top_n,
        output_dir=args.output_dir,
    )
