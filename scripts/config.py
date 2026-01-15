"""
Configuration module for bioinformatics analysis pipeline.
Centralized settings for paths, thresholds, and parameters.

Filter Categories:
==================
1. HARD FILTERS: Strict quality control applied to all datasets
   - SNP-only (no indels, no CNVs)
   - Biallelic (max 2 alleles)
   - MAF (minimum allele frequency)
   - Call rate (genotyping completeness)

2. SITUATIONAL FILTERS: Context-dependent, vary by study design
   - Hardy-Weinberg Equilibrium (HWE)
   - Linkage Disequilibrium (LD) pruning
   - FST-based variant selection
"""

from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List


@dataclass(frozen=True)
class PathConfig:
    """Path configuration for input/output directories."""

    # Root directories
    ROOT: Path = Path("/home/Plutonium/Documents/BioinfoMidterm")
    DATA_DIR: Path = Path("1000genomes")
    OUTPUT_DIR: Path = Path("1000genomes/output")
    VCF_DIR: Path = Path("1000genomes/vcf")
    REPORTS_DIR: Path = Path("reports")

    # Input data
    VCF_FILE: Path = Path("1000genomes/main_vcf/ALL_merged.vcf.gz")
    PANEL_FILE: Path = Path(
        "1000genomes/main_vcf/integrated_call_samples_v3.20130502.ALL.panel"
    )

    # Sample lists
    EAS_SAMPLES_CSV: Path = Path("1000genomes/EAS_subpopulation_samples.csv")
    EAS_SAMPLES_LIST: Path = Path("1000genomes/EAS_subpopulation_samples_list.csv")

    # Intermediate outputs - Hard filtered
    PLINK_SNP_FILTERED: Path = Path("1000genomes/output/EAS_AND_SNP_filtered_data")
    PLINK_MAF_FILTERED: Path = Path(
        "1000genomes/output/EAS_AND_SNP_filtered_data_MAF_filtered"
    )

    # Intermediate outputs - Situational filtered
    PLINK_HWE_FILTERED: Path = Path("1000genomes/output/EAS_SNP_MAF_HWE_filtered")
    PLINK_UNIQUE_IDS: Path = Path(
        "1000genomes/output/EAS_SNP_MAF_HWE_filtered_unique_ids"
    )
    PLINK_LD_PRUNED: Path = Path("1000genomes/output/EAS_FINAL_DATA_FOR_FST")

    # FST and SNP selection
    FST_RESULTS: Path = Path("1000genomes/output/EAS_FST_RESULTS")
    TOP_SNPS_FILE: Path = Path("1000genomes/output/top_snps.txt")
    TOP_SNPS_BED: Path = Path("1000genomes/output/top_snps.bed")
    FST_FILTERED: Path = Path("1000genomes/output/FST_FILTERED")

    # PCA outputs
    PCA_FILE: Path = Path("1000genomes/output/FST_PCA")

    # ML data
    ML_DATA: Path = Path("1000genomes/vcf/vcf_numeric_transposed_with_population.csv")
    ML_MODELS_DIR: Path = Path("output/ml_models")

    def get_absolute(self, relative_path: Path) -> Path:
        """Convert relative path to absolute."""
        return self.ROOT / relative_path


@dataclass(frozen=True)
class HardFilterThresholds:
    """
    HARD FILTERS: Strict quality control parameters.
    These should be applied consistently across all datasets.

    Rationale:
    - SNP_ONLY: Focus on single nucleotide polymorphisms for ancestry inference
    - MAX_ALLELES: Biallelic variants are more reliable and interpretable
    - MIN_AF: Remove singletons/doubletons (1/612 alleles for 306 samples)
    - MIN_CALL_RATE: Ensure high genotyping quality
    """

    # SNP-only filter (HARD) - removes indels, CNVs, structural variants
    SNP_ONLY: bool = True

    # Biallelic filter (HARD) - maximum number of alleles
    MAX_ALLELES: int = 2

    # MAF filter (HARD) - minimum minor allele frequency
    # 1/612 = 0.00163 (at least 1 allele in 306 diploid samples)
    MIN_AF: float = 0.0016

    # Call rate filter (HARD) - minimum genotyping completeness
    MIN_CALL_RATE: float = 0.95

    # Minimum observed count (for frequency analysis)
    MIN_OBS_CT: int = 1


@dataclass(frozen=True)
class SituationalFilterThresholds:
    """
    SITUATIONAL FILTERS: Context-dependent parameters.
    These may vary based on study design, research questions, and populations.

    Rationale:
    - HWE: Removes genotyping errors, but strict filtering may remove
           variants under selection (context-dependent)
    - LD: Reduces redundancy for PCA/ancestry analysis, but keeps
           all info for association studies
    - FST: Selects population-informative variants for ancestry inference
    """

    # Hardy-Weinberg Equilibrium (SITUATIONAL)
    # More relaxed threshold (1e-6) for population genetics
    # Stricter threshold (1e-10) for GWAS
    HWE_P_THRESHOLD: float = 1e-6
    HWE_FILTER_MODE: str = "keep-fewhet"  # Options: "keep-fewhet", "keep-all", None

    # LD pruning parameters (SITUATIONAL)
    # Aggressive for PCA: window=1000kb, step=1, r2=0.1
    # Moderate for ancestry: window=500kb, step=5, r2=0.2
    # Light for association: window=50kb, step=5, r2=0.5
    LD_WINDOW_KB: int = 1000
    LD_STEP: int = 1
    LD_R2_THRESHOLD: float = 0.1

    # FST selection (SITUATIONAL)
    # Top N variants per population pair
    FST_TOP_N: int = 1000
    FST_MIN_VALUE: float = 0.0  # Minimum FST to consider


@dataclass(frozen=True)
class Plink2Config:
    """PLINK2 executable and common options."""

    EXECUTABLE: str = "plink2"
    THREADS: int = 8
    MEMORY_MB: int = 16000  # 16GB


@dataclass(frozen=True)
class PopulationConfig:
    """Population definitions for the study."""

    # East Asian subpopulations of interest
    EAS_SUBPOPS: tuple = ("CHB", "JPT", "KHV")

    # Super population
    SUPER_POP: str = "EAS"

    # Number of samples (for frequency calculations)
    NUM_SAMPLES: int = 306


@dataclass(frozen=True)
class MLConfig:
    """Machine learning configuration."""

    # Data splits
    TEST_SIZE: float = 0.2
    RANDOM_STATE: int = 42

    # Cross-validation
    CV_FOLDS: int = 5

    # Random Forest hyperparameters
    RF_N_ESTIMATORS: int = 200
    RF_MAX_DEPTH: Optional[int] = None
    RF_MIN_SAMPLES_SPLIT: int = 2
    RF_MIN_SAMPLES_LEAF: int = 1

    # XGBoost hyperparameters
    XGB_N_ESTIMATORS: int = 200
    XGB_MAX_DEPTH: int = 6
    XGB_LEARNING_RATE: float = 0.1
    XGB_SUBSAMPLE: float = 0.8

    # Logistic Regression
    LR_MAX_ITER: int = 1000
    LR_SOLVER: str = "lbfgs"

    # Feature selection
    TOP_N_FEATURES: int = 25


# =============================================================================
# Global configuration instances
# =============================================================================
PATHS = PathConfig()
HARD_FILTERS = HardFilterThresholds()
SITUATIONAL_FILTERS = SituationalFilterThresholds()
PLINK = Plink2Config()
POPULATIONS = PopulationConfig()
ML = MLConfig()


def print_config_summary():
    """Print a summary of current configuration."""
    print("=" * 60)
    print("CONFIGURATION SUMMARY")
    print("=" * 60)

    print("\n--- HARD FILTERS (Always Applied) ---")
    print(f"  SNP-only:      {HARD_FILTERS.SNP_ONLY}")
    print(f"  Max alleles:   {HARD_FILTERS.MAX_ALLELES}")
    print(f"  Min AF:        {HARD_FILTERS.MIN_AF}")
    print(f"  Min call rate: {HARD_FILTERS.MIN_CALL_RATE}")

    print("\n--- SITUATIONAL FILTERS (Context-Dependent) ---")
    print(f"  HWE threshold: {SITUATIONAL_FILTERS.HWE_P_THRESHOLD}")
    print(f"  HWE mode:      {SITUATIONAL_FILTERS.HWE_FILTER_MODE}")
    print(f"  LD window:     {SITUATIONAL_FILTERS.LD_WINDOW_KB}kb")
    print(f"  LD step:       {SITUATIONAL_FILTERS.LD_STEP}")
    print(f"  LD R² cutoff:  {SITUATIONAL_FILTERS.LD_R2_THRESHOLD}")
    print(f"  FST top N:     {SITUATIONAL_FILTERS.FST_TOP_N}")

    print("\n--- POPULATIONS ---")
    print(f"  Subpopulations: {POPULATIONS.EAS_SUBPOPS}")
    print(f"  Num samples:    {POPULATIONS.NUM_SAMPLES}")

    print("=" * 60)


if __name__ == "__main__":
    print_config_summary()
