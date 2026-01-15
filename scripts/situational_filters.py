"""
SITUATIONAL FILTERS: Context-dependent quality control filters.

These filters may vary based on:
- Study design (population genetics vs. GWAS vs. ancestry inference)
- Research questions
- Population structure
- Analysis goals

Filter Categories:
1. Hardy-Weinberg Equilibrium (HWE) - removes genotyping errors
2. Linkage Disequilibrium (LD) pruning - reduces redundancy
3. Variant ID standardization - for cross-study compatibility
4. FST-based selection - for population-informative markers

When to use each filter:
- HWE: Almost always, but threshold varies
- LD pruning: For PCA, structure analysis, ancestry; NOT for association
- FST selection: For ancestry panel development, population classification
"""

import os
from pathlib import Path
from typing import Optional, List, Union

# Import configuration and utilities
try:
    from config import PATHS, SITUATIONAL_FILTERS, PLINK
    from utils import run_plink2_command, count_variants, count_samples
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from config import PATHS, SITUATIONAL_FILTERS, PLINK
    from utils import run_plink2_command, count_variants, count_samples


# =============================================================================
# SITUATIONAL FILTER 1: Hardy-Weinberg Equilibrium (HWE)
# =============================================================================


def calculate_hwe_stats(
    input_pfile: str,
    output_prefix: str = None,
    verbose: bool = True,
) -> str:
    """
    Calculate Hardy-Weinberg statistics without filtering.
    Use this to inspect HWE distribution before applying filter.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        verbose: Print progress

    Returns:
        Output prefix (adds .hardy extension)
    """
    if output_prefix is None:
        output_prefix = input_pfile + "_hardy"

    # Use --hardy to compute HWE statistics without filtering (plink2 rejects 1.0 thresholds)
    args = [
        "--pfile",
        input_pfile,
        "--hardy",
        "--out",
        output_prefix,
    ]

    if verbose:
        print(f"Calculating HWE statistics: {output_prefix}.hardy")

    run_plink2_command(args)

    return output_prefix


def filter_hwe(
    input_pfile: str,
    output_pfile: str = None,
    p_threshold: float = SITUATIONAL_FILTERS.HWE_P_THRESHOLD,
    filter_mode: str = SITUATIONAL_FILTERS.HWE_FILTER_MODE,
    verbose: bool = True,
) -> str:
    """
    SITUATIONAL FILTER 1: Hardy-Weinberg Equilibrium filter.

    Removes variants that deviate significantly from HWE expectations.

    Context considerations:
    - Population genetics: Use relaxed threshold (1e-6)
    - GWAS controls: Use stricter threshold (1e-10)
    - Case-control: Apply to controls only
    - Admixed populations: May need to skip or be very relaxed

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        p_threshold: P-value threshold (variants with P < threshold are removed)
        filter_mode: "keep-fewhet" to keep excess homozygotes,
                     None for standard filtering
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = str(PATHS.PLINK_HWE_FILTERED)

    args = [
        "--pfile",
        input_pfile,
        "--hwe",
        str(p_threshold),
    ]

    if filter_mode:
        args.append(filter_mode)

    args.extend(["--make-pgen", "--out", output_pfile])

    if verbose:
        print("=" * 60)
        print("SITUATIONAL FILTER 1: Hardy-Weinberg Equilibrium")
        print("=" * 60)
        print(f"  P-value threshold: {p_threshold}")
        print(f"  Filter mode: {filter_mode or 'standard'}")
        print(f"  Note: This is SITUATIONAL - adjust threshold based on study design")

        n_input = count_variants(input_pfile)
        print(f"  Input variants: {n_input}")

    run_plink2_command(args)

    if verbose:
        n_output = count_variants(output_pfile)
        n_removed = n_input - n_output if n_input > 0 else 0
        print(f"\nOutput: {output_pfile}")
        print(f"  Output variants: {n_output}")
        print(f"  Removed for HWE deviation: {n_removed}" if n_input > 0 else "")

    return output_pfile


# =============================================================================
# SITUATIONAL FILTER 2: Variant ID Standardization
# =============================================================================


def set_variant_ids(
    input_pfile: str,
    output_pfile: str = None,
    id_format: str = "@:#[b37]$r,$a",
    verbose: bool = True,
) -> str:
    """
    Set standardized variant IDs.

    This is important for:
    - Cross-study comparisons
    - Merging datasets
    - Extracting specific variants

    Default format: "chrom:pos[build]ref,alt"
    Example: "12:58124534[b37]C,G"

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        id_format: PLINK2 variant ID format string
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = str(PATHS.PLINK_UNIQUE_IDS)

    args = [
        "--pfile",
        input_pfile,
        "--set-all-var-ids",
        id_format,
        "--make-pgen",
        "--out",
        output_pfile,
    ]

    if verbose:
        print("=" * 60)
        print("SITUATIONAL: Variant ID Standardization")
        print("=" * 60)
        print(f"  ID format: {id_format}")

    run_plink2_command(args)

    if verbose:
        print(f"\nOutput: {output_pfile}")

    return output_pfile


# =============================================================================
# SITUATIONAL FILTER 3: Linkage Disequilibrium (LD) Pruning
# =============================================================================


def calculate_ld_prune_list(
    input_pfile: str,
    output_prefix: str = None,
    window_kb: int = SITUATIONAL_FILTERS.LD_WINDOW_KB,
    step: int = SITUATIONAL_FILTERS.LD_STEP,
    r2_threshold: float = SITUATIONAL_FILTERS.LD_R2_THRESHOLD,
    verbose: bool = True,
) -> str:
    """
    Calculate LD pruning list (without applying it).

    Creates .prune.in (keep) and .prune.out (remove) files.

    Recommended settings by use case:
    - PCA/Structure: window=1000kb, step=1, r2=0.1 (aggressive)
    - Ancestry: window=500kb, step=5, r2=0.2 (moderate)
    - Association: Usually skip LD pruning

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        window_kb: Window size in kilobases
        step: Step size in SNPs
        r2_threshold: R-squared threshold for LD
        verbose: Print progress

    Returns:
        Output prefix (adds .prune.in and .prune.out)
    """
    if output_prefix is None:
        output_prefix = f"{input_pfile}_LD_pruned_{window_kb}kb_{step}_{r2_threshold}"

    args = [
        "--pfile",
        input_pfile,
        "--indep-pairwise",
        f"{window_kb}kb",
        str(step),
        str(r2_threshold),
        "--out",
        output_prefix,
    ]

    if verbose:
        print("=" * 60)
        print("SITUATIONAL FILTER 3: LD Pruning (Calculate List)")
        print("=" * 60)
        print(f"  Window: {window_kb} kb")
        print(f"  Step: {step} SNPs")
        print(f"  R² threshold: {r2_threshold}")
        print(f"  Note: This is SITUATIONAL - skip for association studies")

    run_plink2_command(args)

    if verbose:
        # Count variants in prune.in file
        prune_in = f"{output_prefix}.prune.in"
        if os.path.exists(prune_in):
            with open(prune_in) as f:
                n_keep = sum(1 for _ in f)
            print(f"\nVariants to keep: {n_keep}")
            print(f"Prune list: {prune_in}")

    return output_prefix


def apply_ld_prune_list(
    input_pfile: str,
    prune_file: str,
    output_pfile: str = None,
    verbose: bool = True,
) -> str:
    """
    Apply a pre-calculated LD pruning list.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        prune_file: Path to .prune.in file
        output_pfile: Path to output PLINK2 pfile prefix
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    if output_pfile is None:
        output_pfile = str(PATHS.PLINK_LD_PRUNED)

    args = [
        "--pfile",
        input_pfile,
        "--extract",
        prune_file,
        "--make-pgen",
        "--out",
        output_pfile,
    ]

    if verbose:
        print(f"Applying LD pruning from: {prune_file}")
        n_input = count_variants(input_pfile)
        print(f"  Input variants: {n_input}")

    run_plink2_command(args)

    if verbose:
        n_output = count_variants(output_pfile)
        print(f"\nOutput: {output_pfile}")
        print(f"  Output variants: {n_output}")

    return output_pfile


def filter_ld(
    input_pfile: str,
    output_pfile: str = None,
    window_kb: int = SITUATIONAL_FILTERS.LD_WINDOW_KB,
    step: int = SITUATIONAL_FILTERS.LD_STEP,
    r2_threshold: float = SITUATIONAL_FILTERS.LD_R2_THRESHOLD,
    verbose: bool = True,
) -> str:
    """
    Calculate and apply LD pruning in one step.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_pfile: Path to output PLINK2 pfile prefix
        window_kb: Window size in kilobases
        step: Step size in SNPs
        r2_threshold: R-squared threshold
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
    # Step 1: Calculate prune list
    prune_prefix = calculate_ld_prune_list(
        input_pfile=input_pfile,
        window_kb=window_kb,
        step=step,
        r2_threshold=r2_threshold,
        verbose=verbose,
    )

    # Step 2: Apply prune list
    prune_file = f"{prune_prefix}.prune.in"
    output = apply_ld_prune_list(
        input_pfile=input_pfile,
        prune_file=prune_file,
        output_pfile=output_pfile,
        verbose=verbose,
    )

    return output


# =============================================================================
# SITUATIONAL FILTER 4: FST Calculation
# =============================================================================


def calculate_fst(
    input_pfile: str,
    output_prefix: str = None,
    pop_column: str = "pop",
    verbose: bool = True,
) -> str:
    """
    Calculate FST statistics between populations.

    FST measures genetic differentiation between populations.
    High FST variants are "population-informative" - useful for
    ancestry inference and population classification.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        pop_column: Column name in .psam file containing population labels
        verbose: Print progress

    Returns:
        Output prefix (creates .fst and .fst.var files)
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
        print("SITUATIONAL: FST Calculation")
        print("=" * 60)
        print(f"  Population column: {pop_column}")
        print(f"  Note: Requires population labels in .psam file")

    run_plink2_command(args)

    if verbose:
        print(f"\nFST results: {output_prefix}.*")

    return output_prefix


def extract_variants(
    input_pfile: str,
    variant_list: str,
    output_pfile: str,
    verbose: bool = True,
) -> str:
    """
    Extract specific variants from a pfile.

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        variant_list: Path to file with variant IDs (one per line)
        output_pfile: Path to output PLINK2 pfile prefix
        verbose: Print progress

    Returns:
        Output pfile prefix
    """
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
        print(f"Extracting variants from: {variant_list}")

    run_plink2_command(args)

    if verbose:
        n_output = count_variants(output_pfile)
        print(f"\nOutput: {output_pfile}")
        print(f"  Variants: {n_output}")

    return output_pfile


# =============================================================================
# SITUATIONAL: PCA Calculation
# =============================================================================


def calculate_pca(
    input_pfile: str,
    output_prefix: str = None,
    n_components: int = 10,
    verbose: bool = True,
) -> str:
    """
    Calculate principal components from genetic data.

    Best practices:
    - Apply LD pruning before PCA
    - Remove relatives if present
    - Consider removing outliers iteratively

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Path to output file prefix
        n_components: Number of PCs to compute
        verbose: Print progress

    Returns:
        Output prefix (creates .eigenvec and .eigenval files)
    """
    if output_prefix is None:
        output_prefix = str(PATHS.PCA_FILE)

    args = [
        "--pfile",
        input_pfile,
        "--pca",
        str(n_components),
        "--out",
        output_prefix,
    ]

    if verbose:
        print("=" * 60)
        print("SITUATIONAL: PCA Calculation")
        print("=" * 60)
        print(f"  Components: {n_components}")

    run_plink2_command(args)

    if verbose:
        print(f"\nPCA results: {output_prefix}.eigenvec")
        print(f"Eigenvalues: {output_prefix}.eigenval")

    return output_prefix


# =============================================================================
# Pipeline: Apply All Situational Filters
# =============================================================================


def apply_all_situational_filters(
    input_pfile: str,
    output_prefix: str = None,
    apply_hwe: bool = True,
    apply_ld: bool = True,
    apply_fst: bool = True,
    apply_pca: bool = True,
    hwe_threshold: float = SITUATIONAL_FILTERS.HWE_P_THRESHOLD,
    ld_window: int = SITUATIONAL_FILTERS.LD_WINDOW_KB,
    ld_step: int = SITUATIONAL_FILTERS.LD_STEP,
    ld_r2: float = SITUATIONAL_FILTERS.LD_R2_THRESHOLD,
    pop_column: str = "pop",
    verbose: bool = True,
) -> dict:
    """
    Apply situational filters in sequence.

    Pipeline (each step is optional):
    1. HWE filter
    2. Set variant IDs
    3. LD pruning
    4. FST calculation
    5. PCA

    Args:
        input_pfile: Path to input PLINK2 pfile prefix
        output_prefix: Base prefix for all outputs
        apply_hwe: Apply HWE filter
        apply_ld: Apply LD pruning
        apply_fst: Calculate FST
        apply_pca: Calculate PCA
        hwe_threshold: HWE p-value threshold
        ld_window: LD pruning window size
        ld_step: LD pruning step size
        ld_r2: LD pruning R² threshold
        pop_column: Population column name
        verbose: Print progress

    Returns:
        Dictionary with output paths for each step
    """
    if output_prefix is None:
        output_prefix = str(PATHS.OUTPUT_DIR / "situational_filtered")

    print("=" * 60)
    print("APPLYING SITUATIONAL FILTERS")
    print("=" * 60)
    print(f"  HWE filter: {apply_hwe}")
    print(f"  LD pruning: {apply_ld}")
    print(f"  FST calculation: {apply_fst}")
    print(f"  PCA: {apply_pca}")

    outputs = {}
    current_pfile = input_pfile

    # Step 1: HWE filter
    if apply_hwe:
        hwe_output = f"{output_prefix}_hwe"
        calculate_hwe_stats(current_pfile, verbose=verbose)
        current_pfile = filter_hwe(
            current_pfile,
            output_pfile=hwe_output,
            p_threshold=hwe_threshold,
            verbose=verbose,
        )
        outputs["hwe"] = current_pfile

    # Step 2: Set variant IDs
    id_output = f"{output_prefix}_unique_ids"
    current_pfile = set_variant_ids(
        current_pfile,
        output_pfile=id_output,
        verbose=verbose,
    )
    outputs["unique_ids"] = current_pfile

    # Step 3: LD pruning
    if apply_ld:
        ld_output = f"{output_prefix}_ld_pruned"
        current_pfile = filter_ld(
            current_pfile,
            output_pfile=ld_output,
            window_kb=ld_window,
            step=ld_step,
            r2_threshold=ld_r2,
            verbose=verbose,
        )
        outputs["ld_pruned"] = current_pfile

    # Step 4: FST calculation
    if apply_fst:
        fst_output = f"{output_prefix}_fst"
        calculate_fst(
            current_pfile,
            output_prefix=fst_output,
            pop_column=pop_column,
            verbose=verbose,
        )
        outputs["fst"] = fst_output

    # Step 5: PCA
    if apply_pca:
        pca_output = f"{output_prefix}_pca"
        calculate_pca(
            current_pfile,
            output_prefix=pca_output,
            verbose=verbose,
        )
        outputs["pca"] = pca_output

    outputs["final"] = current_pfile

    print("\n" + "=" * 60)
    print("SITUATIONAL FILTERING COMPLETE")
    print("=" * 60)
    print(f"Final output: {current_pfile}")
    print(f"  Variants: {count_variants(current_pfile)}")

    return outputs


# =============================================================================
# CLI Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Apply situational filters to genetic data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--pfile", required=True, help="Input PLINK2 pfile prefix")
    parser.add_argument("--output", "-o", required=True, help="Output prefix")
    parser.add_argument("--no-hwe", action="store_true", help="Skip HWE filter")
    parser.add_argument("--no-ld", action="store_true", help="Skip LD pruning")
    parser.add_argument("--no-fst", action="store_true", help="Skip FST calculation")
    parser.add_argument("--no-pca", action="store_true", help="Skip PCA")
    parser.add_argument(
        "--hwe-threshold",
        type=float,
        default=SITUATIONAL_FILTERS.HWE_P_THRESHOLD,
        help="HWE p-value threshold",
    )
    parser.add_argument(
        "--ld-window",
        type=int,
        default=SITUATIONAL_FILTERS.LD_WINDOW_KB,
        help="LD window size (kb)",
    )
    parser.add_argument(
        "--ld-r2",
        type=float,
        default=SITUATIONAL_FILTERS.LD_R2_THRESHOLD,
        help="LD R² threshold",
    )
    parser.add_argument(
        "--pop-column", default="pop", help="Population column name in .psam"
    )

    args = parser.parse_args()

    apply_all_situational_filters(
        input_pfile=args.pfile,
        output_prefix=args.output,
        apply_hwe=not args.no_hwe,
        apply_ld=not args.no_ld,
        apply_fst=not args.no_fst,
        apply_pca=not args.no_pca,
        hwe_threshold=args.hwe_threshold,
        ld_window=args.ld_window,
        ld_r2=args.ld_r2,
        pop_column=args.pop_column,
    )
