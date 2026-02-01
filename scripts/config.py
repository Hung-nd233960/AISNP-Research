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


@dataclass
class PathConfig:
    """Path configuration for input/output directories.

    Paths are dynamically generated based on the population configuration.
    - khv_jpt_chb: Uses 1000genomes/output/ directory
    - sea_jpt_cn: Uses 1000genomes/output_sea_jpt_cn/ directory
    """

    # Configuration name (determines output directory)
    config_name: str = "khv_jpt_chb"

    # Root directories (constant)
    ROOT: Path = Path("/home/Plutonium/Documents/BioinfoMidterm")

    # Input data (shared between configs)
    VCF_FILE: Path = Path("1000genomes/main_vcf/ALL_merged.vcf.gz")
    PANEL_FILE: Path = Path(
        "1000genomes/main_vcf/integrated_call_samples_v3.20130502.ALL.panel"
    )

    @property
    def _prefix(self) -> str:
        """Get the prefix for output files based on config."""
        if self.config_name == "khv_jpt_chb":
            return "EAS"
        else:
            return "SEA_JPT_CN"

    @property
    def _output_subdir(self) -> str:
        """Get the output subdirectory based on config."""
        if self.config_name == "khv_jpt_chb":
            return "output"
        else:
            return "output_sea_jpt_cn"

    @property
    def DATA_DIR(self) -> Path:
        return Path("1000genomes")

    @property
    def OUTPUT_DIR(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}")

    @property
    def VCF_DIR(self) -> Path:
        return Path(f"1000genomes/vcf_{self.config_name}")

    @property
    def REPORTS_DIR(self) -> Path:
        return Path(f"reports/{self.config_name}")

    # Sample lists
    @property
    def EAS_SAMPLES_CSV(self) -> Path:
        return Path(f"1000genomes/{self._prefix}_subpopulation_samples.csv")

    @property
    def EAS_SAMPLES_LIST(self) -> Path:
        return Path(f"1000genomes/{self._prefix}_subpopulation_samples_list.csv")

    # Intermediate outputs - Hard filtered
    @property
    def PLINK_SNP_FILTERED(self) -> Path:
        return Path(
            f"1000genomes/{self._output_subdir}/{self._prefix}_AND_SNP_filtered_data"
        )

    @property
    def PLINK_MAF_FILTERED(self) -> Path:
        return Path(
            f"1000genomes/{self._output_subdir}/{self._prefix}_AND_SNP_filtered_data_MAF_filtered"
        )

    # Intermediate outputs - Situational filtered
    @property
    def PLINK_HWE_FILTERED(self) -> Path:
        return Path(
            f"1000genomes/{self._output_subdir}/{self._prefix}_SNP_MAF_HWE_filtered"
        )

    @property
    def PLINK_UNIQUE_IDS(self) -> Path:
        return Path(
            f"1000genomes/{self._output_subdir}/{self._prefix}_SNP_MAF_HWE_filtered_unique_ids"
        )

    @property
    def PLINK_LD_PRUNED(self) -> Path:
        return Path(
            f"1000genomes/{self._output_subdir}/{self._prefix}_FINAL_DATA_FOR_FST"
        )

    # FST and SNP selection
    @property
    def FST_RESULTS(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}/{self._prefix}_FST_RESULTS")

    @property
    def TOP_SNPS_FILE(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}/top_snps.txt")

    @property
    def TOP_SNPS_BED(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}/top_snps.bed")

    @property
    def FST_FILTERED(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}/FST_FILTERED")

    # PCA outputs
    @property
    def PCA_FILE(self) -> Path:
        return Path(f"1000genomes/{self._output_subdir}/FST_PCA")

    # ML data
    @property
    def ML_DATA(self) -> Path:
        return Path(
            f"1000genomes/vcf_{self.config_name}/vcf_numeric_transposed_with_population.csv"
        )

    @property
    def ML_MODELS_DIR(self) -> Path:
        return Path(f"output/ml_models/{self.config_name}")

    def get_absolute(self, relative_path: Path) -> Path:
        """Convert relative path to absolute."""
        return self.ROOT / relative_path

    def ensure_output_dirs(self) -> None:
        """Create output directories if they don't exist."""
        dirs_to_create = [
            self.OUTPUT_DIR,
            self.VCF_DIR,
            self.REPORTS_DIR,
            self.ML_MODELS_DIR,
        ]
        for dir_path in dirs_to_create:
            abs_path = self.get_absolute(dir_path)
            abs_path.mkdir(parents=True, exist_ok=True)


def get_path_config(config_name: str = "khv_jpt_chb") -> PathConfig:
    """
    Get path configuration for a specific population config.

    Args:
        config_name: One of "khv_jpt_chb" or "sea_jpt_cn"

    Returns:
        PathConfig instance with appropriate paths
    """
    valid_configs = ["khv_jpt_chb", "sea_jpt_cn"]
    if config_name not in valid_configs:
        raise ValueError(f"Unknown config: {config_name}. Available: {valid_configs}")
    return PathConfig(config_name=config_name)


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

    # Configuration name
    CONFIG_NAME: str = "khv_jpt_chb"

    # Super population
    SUPER_POP: str = "EAS"

    # East Asian subpopulations in the raw data
    RAW_SUBPOPS: tuple = ("CHB", "JPT", "KHV")

    # Target populations after any merging (same as RAW for khv_jpt_chb)
    TARGET_POPS: tuple = ("CHB", "JPT", "KHV")

    # Population merge mapping (empty dict means no merging)
    # Format: {raw_pop: merged_group}
    POP_MERGE_MAP: dict = field(default_factory=dict)

    # Number of samples (for frequency calculations)
    NUM_SAMPLES: int = 306

    @property
    def EAS_SUBPOPS(self) -> tuple:
        """Backward compatibility alias for TARGET_POPS."""
        return self.TARGET_POPS


# =============================================================================
# Predefined Population Configurations
# =============================================================================


def get_khv_jpt_chb_config() -> PopulationConfig:
    """
    KHV-JPT-CHB configuration: 3 populations, no merging.
    Uses: KHV, JPT, CHB
    """
    return PopulationConfig(
        CONFIG_NAME="khv_jpt_chb",
        SUPER_POP="EAS",
        RAW_SUBPOPS=("CHB", "JPT", "KHV"),
        TARGET_POPS=("CHB", "JPT", "KHV"),
        POP_MERGE_MAP={},
        NUM_SAMPLES=306,
    )


def get_sea_jpt_cn_config() -> PopulationConfig:
    """
    SEA-JPT-CN configuration: 5 raw populations merged into 3 groups.
    Raw populations: CHB, JPT, CHS, CDX, KHV
    Merged groups:
      - CN: CHB + CHS (Han Chinese)
      - SEA: KHV + CDX (Southeast Asian)
      - JPT: JPT (Japanese)
    """
    return PopulationConfig(
        CONFIG_NAME="sea_jpt_cn",
        SUPER_POP="EAS",
        RAW_SUBPOPS=("CHB", "JPT", "CHS", "CDX", "KHV"),
        TARGET_POPS=("CN", "JPT", "SEA"),
        POP_MERGE_MAP={
            "CHB": "CN",
            "CHS": "CN",
            "JPT": "JPT",
            "KHV": "SEA",
            "CDX": "SEA",
        },
        NUM_SAMPLES=504,  # All 5 EAS subpopulations
    )


def get_population_config(config_name: str = "khv_jpt_chb") -> PopulationConfig:
    """
    Get population configuration by name.

    Args:
        config_name: One of "khv_jpt_chb" or "sea_jpt_cn"

    Returns:
        PopulationConfig instance
    """
    configs = {
        "khv_jpt_chb": get_khv_jpt_chb_config,
        "sea_jpt_cn": get_sea_jpt_cn_config,
    }

    if config_name not in configs:
        raise ValueError(
            f"Unknown population config: {config_name}. "
            f"Available: {list(configs.keys())}"
        )

    return configs[config_name]()


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

# Default config name
_CURRENT_CONFIG_NAME = "khv_jpt_chb"

PATHS = get_path_config(_CURRENT_CONFIG_NAME)
HARD_FILTERS = HardFilterThresholds()
SITUATIONAL_FILTERS = SituationalFilterThresholds()
PLINK = Plink2Config()

# Default population config - can be changed by calling set_population_config()
POPULATIONS = get_khv_jpt_chb_config()
ML = MLConfig()


def set_population_config(config_name: str) -> None:
    """
    Set the global population and path configuration.

    Args:
        config_name: One of "khv_jpt_chb" or "sea_jpt_cn"

    This updates both POPULATIONS and PATHS global variables.
    """
    global POPULATIONS, PATHS, _CURRENT_CONFIG_NAME
    _CURRENT_CONFIG_NAME = config_name
    POPULATIONS = get_population_config(config_name)
    PATHS = get_path_config(config_name)

    # Ensure output directories exist
    PATHS.ensure_output_dirs()


def get_current_config_name() -> str:
    """Get the current configuration name."""
    return _CURRENT_CONFIG_NAME


def print_config_summary():
    """Print a summary of current configuration."""
    print("=" * 60)
    print("CONFIGURATION SUMMARY")
    print("=" * 60)

    print("\n--- PATHS ---")
    print(f"  Config name:    {PATHS.config_name}")
    print(f"  Output dir:     {PATHS.OUTPUT_DIR}")
    print(f"  VCF dir:        {PATHS.VCF_DIR}")
    print(f"  Samples CSV:    {PATHS.EAS_SAMPLES_CSV}")

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
    print(f"  Config name:    {POPULATIONS.CONFIG_NAME}")
    print(f"  Raw subpops:    {POPULATIONS.RAW_SUBPOPS}")
    print(f"  Target pops:    {POPULATIONS.TARGET_POPS}")
    if POPULATIONS.POP_MERGE_MAP:
        print(f"  Merge mapping:  {POPULATIONS.POP_MERGE_MAP}")
    print(f"  Num samples:    {POPULATIONS.NUM_SAMPLES}")

    print("=" * 60)


if __name__ == "__main__":
    print_config_summary()
