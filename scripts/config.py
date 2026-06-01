"""
Configuration module for the SEA-JPT-CN ancestry inference pipeline.

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
from typing import Optional

# Project root is two levels up from this file (scripts/config.py → project root)
PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _load_path_roots() -> dict:
    """
    Load path roots from paths.local.yaml (gitignored) or paths.yaml.
    Falls back to hardcoded defaults if neither file exists.
    """
    defaults = {
        "genomes_data": "data/1000genomes",
        "output": "output",
        "reports": "reports",
    }
    for candidate in ("paths.local.yaml", "paths.yaml"):
        p = PROJECT_ROOT / candidate
        if p.exists():
            try:
                import yaml  # type: ignore
                with open(p) as f:
                    data = yaml.safe_load(f) or {}
                return {**defaults, **{k: v for k, v in data.items() if k in defaults}}
            except Exception:
                pass
    return defaults


_PATH_ROOTS = _load_path_roots()


@dataclass
class PathConfig:
    """Path configuration for the SEA-JPT-CN pipeline.

    Roots (genomes_data, output, reports) are read from paths.local.yaml or
    paths.yaml at the project root; edit those files to change data locations.
    """

    ROOT: Path = field(default_factory=lambda: PROJECT_ROOT)

    _genomes_data: str = field(default_factory=lambda: _PATH_ROOTS["genomes_data"])
    _output_root: str = field(default_factory=lambda: _PATH_ROOTS["output"])
    _reports_root: str = field(default_factory=lambda: _PATH_ROOTS["reports"])

    def __post_init__(self):
        self.PLINK_MERGED = Path(f"{self._genomes_data}/plink/allchr")
        self.VCF_FILE = Path(f"{self._genomes_data}/main_vcf/ALL_merged.vcf.gz")
        self.PANEL_FILE = Path(
            f"{self._genomes_data}/integrated_call_samples_v3.20130502.ALL.panel"
        )

    @property
    def DATA_DIR(self) -> Path:
        return Path(self._genomes_data)

    @property
    def OUTPUT_DIR(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn")

    @property
    def VCF_DIR(self) -> Path:
        return Path(f"{self._genomes_data}/vcf_sea_jpt_cn")

    @property
    def REPORTS_DIR(self) -> Path:
        return Path(f"{self._reports_root}/sea_jpt_cn")

    @property
    def SAMPLES_CSV(self) -> Path:
        return Path(f"{self._genomes_data}/SEA_JPT_CN_subpopulation_samples.csv")

    @property
    def SAMPLES_LIST(self) -> Path:
        return Path(f"{self._genomes_data}/SEA_JPT_CN_subpopulation_samples_list.csv")

    # Intermediate outputs — hard filtered
    @property
    def PLINK_SNP_FILTERED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_AND_SNP_filtered_data")

    @property
    def PLINK_MAF_FILTERED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_AND_SNP_filtered_data_MAF_filtered")

    # Intermediate outputs — situational filtered
    @property
    def PLINK_HWE_FILTERED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_SNP_MAF_HWE_filtered")

    @property
    def PLINK_UNIQUE_IDS(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_SNP_MAF_HWE_filtered_unique_ids")

    @property
    def PLINK_LD_PRUNED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_FINAL_DATA_FOR_FST")

    # FST and SNP selection
    @property
    def FST_RESULTS(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/SEA_JPT_CN_FST_RESULTS")

    @property
    def TOP_SNPS_FILE(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/top_snps.txt")

    @property
    def TOP_SNPS_BED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/top_snps.bed")

    @property
    def FST_FILTERED(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/FST_FILTERED")

    @property
    def PCA_FILE(self) -> Path:
        return Path(f"{self._genomes_data}/output_sea_jpt_cn/FST_PCA")

    @property
    def ML_DATA(self) -> Path:
        return Path(f"{self._genomes_data}/vcf_sea_jpt_cn/vcf_numeric_transposed_with_population.csv")

    @property
    def ML_MODELS_DIR(self) -> Path:
        return Path(f"{self._output_root}/ml_models/sea_jpt_cn")

    def get_absolute(self, relative_path: Path) -> Path:
        return self.ROOT / relative_path

    def ensure_output_dirs(self) -> None:
        for dir_path in [self.OUTPUT_DIR, self.VCF_DIR, self.REPORTS_DIR, self.ML_MODELS_DIR]:
            self.get_absolute(dir_path).mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class HardFilterThresholds:
    """HARD FILTERS: applied consistently across all datasets."""

    SNP_ONLY: bool = True
    MAX_ALLELES: int = 2
    # 1/(504*2) — one allele across all 504 samples
    MIN_AF: float = 1 / (504 * 2)
    MIN_CALL_RATE: float = 0.95
    MIN_OBS_CT: int = 1


@dataclass(frozen=True)
class SituationalFilterThresholds:
    """SITUATIONAL FILTERS: context-dependent, vary by study design."""

    # More relaxed (1e-6) for population genetics; stricter (1e-10) for GWAS
    HWE_P_THRESHOLD: float = 1e-6
    HWE_FILTER_MODE: str = "keep-fewhet"

    # Aggressive LD pruning for PCA/ancestry
    LD_WINDOW_KB: int = 1000
    LD_STEP: int = 1
    LD_R2_THRESHOLD: float = 0.1

    FST_TOP_N: int = 1000
    FST_MIN_VALUE: float = 0.0


@dataclass(frozen=True)
class Plink2Config:
    EXECUTABLE: str = "plink2"
    THREADS: int = 16
    MEMORY_MB: int = 16000


@dataclass(frozen=True)
class PopulationConfig:
    """SEA-JPT-CN: 5 raw 1000 Genomes populations merged into 3 groups.

    Raw: CHB, JPT, CHS, CDX, KHV
    Merged:
      CN  = CHB + CHS (Han Chinese)
      JPT = JPT       (Japanese)
      SEA = KHV + CDX (Southeast Asian)
    """

    SUPER_POP: str = "EAS"
    RAW_SUBPOPS: tuple = ("CHB", "JPT", "CHS", "CDX", "KHV")
    TARGET_POPS: tuple = ("CN", "JPT", "SEA")
    POP_MERGE_MAP: dict = field(default_factory=lambda: {
        "CHB": "CN", "CHS": "CN",
        "JPT": "JPT",
        "KHV": "SEA", "CDX": "SEA",
    })
    NUM_SAMPLES: int = 504


@dataclass(frozen=True)
class MLConfig:
    TEST_SIZE: float = 0.2
    RANDOM_STATE: int = 42
    CV_FOLDS: int = 5

    RF_N_ESTIMATORS: int = 200
    RF_MAX_DEPTH: Optional[int] = None
    RF_MIN_SAMPLES_SPLIT: int = 2
    RF_MIN_SAMPLES_LEAF: int = 1

    XGB_N_ESTIMATORS: int = 200
    XGB_MAX_DEPTH: int = 6
    XGB_LEARNING_RATE: float = 0.1
    XGB_SUBSAMPLE: float = 0.8

    LR_MAX_ITER: int = 1000
    LR_SOLVER: str = "lbfgs"

    TOP_N_FEATURES: int = 25


# ---------------------------------------------------------------------------
# Module-level constants — import these directly, no setter needed
# ---------------------------------------------------------------------------

PATHS = PathConfig()
HARD_FILTERS = HardFilterThresholds()
SITUATIONAL_FILTERS = SituationalFilterThresholds()
PLINK = Plink2Config()
POPULATIONS = PopulationConfig()
ML = MLConfig()


def print_config_summary() -> None:
    print("=" * 60)
    print("CONFIGURATION SUMMARY  (sea_jpt_cn)")
    print("=" * 60)
    print(f"\n  genomes_data:  {PATHS._genomes_data}")
    print(f"  output root:   {PATHS._output_root}")
    print(f"  OUTPUT_DIR:    {PATHS.OUTPUT_DIR}")
    print(f"  ML_MODELS_DIR: {PATHS.ML_MODELS_DIR}")
    print(f"\n  Hard filters:  MAF≥{HARD_FILTERS.MIN_AF:.4f}  call_rate≥{HARD_FILTERS.MIN_CALL_RATE}")
    print(f"  LD pruning:    {SITUATIONAL_FILTERS.LD_WINDOW_KB}kb / r²<{SITUATIONAL_FILTERS.LD_R2_THRESHOLD}")
    print(f"  FST top-N:     {SITUATIONAL_FILTERS.FST_TOP_N}")
    print(f"\n  Populations:   {POPULATIONS.TARGET_POPS}  (n={POPULATIONS.NUM_SAMPLES})")
    print("=" * 60)


if __name__ == "__main__":
    print_config_summary()
