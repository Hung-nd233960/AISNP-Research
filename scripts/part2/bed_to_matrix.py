"""
bed_to_matrix.py - Convert BED files to ML-ready genotype matrices

This module provides functions to:
1. Load SNP lists from BED files
2. Extract genotypes from VCF or PLINK pfiles
3. Create ML-ready matrices with sample labels
"""

import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from tqdm import tqdm
import tempfile
import os


def bed_to_snp_list(bed_path: str) -> Tuple[List[str], pd.DataFrame]:
    """
    Load SNP IDs and coordinates from BED file.

    Args:
        bed_path: Path to BED file (chr, start, end, rsid)

    Returns:
        Tuple of (list of rsIDs, DataFrame with full info)
    """
    # BED files have no header
    df = pd.read_csv(
        bed_path, sep="\t", header=None, names=["chr", "start", "end", "rsid"]
    )

    rsids = df["rsid"].tolist()
    print(f"Loaded {len(rsids)} SNPs from {bed_path}")

    return rsids, df


def extract_genotypes_from_vcf(
    vcf_path: str, snp_ids: List[str], output_prefix: str
) -> str:
    """
    Extract specific SNPs from VCF using bcftools.

    Args:
        vcf_path: Path to input VCF file
        snp_ids: List of SNP IDs to extract
        output_prefix: Output file prefix

    Returns:
        Path to extracted VCF
    """
    # Create temp file with SNP IDs
    snp_file = f"{output_prefix}_snp_ids.txt"
    with open(snp_file, "w") as f:
        for snp in snp_ids:
            f.write(f"{snp}\n")

    output_vcf = f"{output_prefix}_extracted.vcf"

    # Use bcftools to extract
    cmd = f"bcftools view -i 'ID=@{snp_file}' {vcf_path} -o {output_vcf}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Warning: bcftools error: {result.stderr}")
        # Try with plink2 as fallback
        return None

    print(f"Extracted VCF: {output_vcf}")
    return output_vcf


def extract_genotypes_from_pfile(
    pfile_prefix: str, snp_ids: List[str], output_prefix: str
) -> str:
    """
    Extract specific SNPs from PLINK2 pfile.

    Args:
        pfile_prefix: Path prefix for pfile (without .pgen/.pvar/.psam)
        snp_ids: List of SNP IDs to extract
        output_prefix: Output file prefix

    Returns:
        Path prefix to extracted pfile
    """
    # Create temp file with SNP IDs
    snp_file = f"{output_prefix}_snp_ids.txt"
    with open(snp_file, "w") as f:
        for snp in snp_ids:
            f.write(f"{snp}\n")

    # Use plink2 to extract
    cmd = (
        f"plink2 --pfile {pfile_prefix} "
        f"--extract {snp_file} "
        f"--make-pgen --out {output_prefix}"
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"PLINK2 extraction warning: {result.stderr}")

    # Count extracted variants
    pvar_file = f"{output_prefix}.pvar"
    if os.path.exists(pvar_file):
        with open(pvar_file) as f:
            n_variants = sum(1 for line in f if not line.startswith("#"))
        print(f"Extracted {n_variants} variants to {output_prefix}")

    return output_prefix


def pfile_to_genotype_matrix(
    pfile_prefix: str, samples_csv: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert PLINK2 pfile to genotype matrix DataFrame.

    Args:
        pfile_prefix: Path prefix for pfile
        samples_csv: Optional CSV with sample, pop columns

    Returns:
        DataFrame with samples as rows, SNPs as columns, plus pop label
    """
    # Export to VCF first (simpler parsing)
    vcf_path = f"{pfile_prefix}_temp.vcf"
    cmd = f"plink2 --pfile {pfile_prefix} --export vcf --out {pfile_prefix}_temp"
    subprocess.run(cmd, shell=True, capture_output=True)

    # Parse VCF
    genotypes = {}
    sample_names = []

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                sample_names = fields[9:]
                continue

            fields = line.strip().split("\t")
            snp_id = fields[2] if fields[2] != "." else f"{fields[0]}:{fields[1]}"

            genos = []
            for gt_field in fields[9:]:
                gt = gt_field.split(":")[0]
                if gt in ["0/0", "0|0"]:
                    genos.append(0)
                elif gt in ["0/1", "1/0", "0|1", "1|0"]:
                    genos.append(1)
                elif gt in ["1/1", "1|1"]:
                    genos.append(2)
                else:
                    genos.append(np.nan)

            genotypes[snp_id] = genos

    # Create DataFrame
    df = pd.DataFrame(genotypes, index=sample_names)
    df = df.reset_index().rename(columns={"index": "sample"})

    # Add population labels if provided
    if samples_csv:
        samples_df = pd.read_csv(samples_csv, header=None)
        if samples_df.shape[1] >= 2:
            samples_df.columns = ["sample", "pop"] + [
                f"col{i}" for i in range(2, samples_df.shape[1])
            ]
            df = df.merge(samples_df[["sample", "pop"]], on="sample", how="left")
            # Reorder columns
            cols = ["sample", "pop"] + [
                c for c in df.columns if c not in ["sample", "pop"]
            ]
            df = df[cols]

    # Cleanup temp VCF
    if os.path.exists(vcf_path):
        os.remove(vcf_path)
    log_file = f"{pfile_prefix}_temp.log"
    if os.path.exists(log_file):
        os.remove(log_file)

    return df


def create_ml_matrix(
    pfile_prefix: str,
    snp_ids: List[str],
    samples_csv: str,
    output_path: str,
    source_name: str = "known_aisnps",
) -> pd.DataFrame:
    """
    Create ML-ready matrix from pfile for specific SNPs.

    Args:
        pfile_prefix: Path to input pfile
        snp_ids: List of SNP IDs to include
        samples_csv: CSV with sample, pop columns
        output_path: Output CSV path
        source_name: Name for intermediate files

    Returns:
        DataFrame ready for ML (sample, pop, snp1, snp2, ...)
    """
    output_dir = Path(output_path).parent
    temp_prefix = str(output_dir / f"{source_name}_extracted")

    # Extract SNPs
    print(f"Extracting {len(snp_ids)} SNPs from pfile...")
    extract_genotypes_from_pfile(pfile_prefix, snp_ids, temp_prefix)

    # Convert to matrix
    print("Converting to genotype matrix...")
    df = pfile_to_genotype_matrix(temp_prefix, samples_csv)

    # Save
    df.to_csv(output_path, index=False)
    print(f"ML matrix saved: {output_path}")
    print(f"  Shape: {df.shape}")
    print(f"  Samples: {len(df)}")
    print(f"  SNPs: {len(df.columns) - 2}")  # Exclude sample, pop

    # Cleanup temp files
    for ext in [".pgen", ".pvar", ".psam", ".log", "_snp_ids.txt"]:
        temp_file = f"{temp_prefix}{ext}"
        if os.path.exists(temp_file):
            os.remove(temp_file)

    return df


def merge_ml_matrices(
    matrices: List[pd.DataFrame], source_names: List[str]
) -> pd.DataFrame:
    """
    Merge multiple ML matrices (from different SNP sets) for comparison.

    Args:
        matrices: List of DataFrames with sample, pop, snps...
        source_names: Names for each source

    Returns:
        Merged DataFrame with source column
    """
    merged = []
    for df, name in zip(matrices, source_names):
        df = df.copy()
        df["source"] = name
        merged.append(df)

    return pd.concat(merged, ignore_index=True)


if __name__ == "__main__":
    print("bed_to_matrix.py - BED to ML matrix conversion utilities")
    print("Import this module to use the functions.")
