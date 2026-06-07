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
                with open(p, encoding="utf-8") as f:
                    data = yaml.safe_load(f) or {}
                return {**defaults, **{k: v for k, v in data.items() if k in defaults}}
            except Exception:
                pass
    return defaults


_PATH_ROOTS = _load_path_roots()


@dataclass
class PathConfig:
    """Path configuration for the SEA-JPT-CN pipeline.

    Directory layout inside genomes_data:
      raw_downloads/   downloaded 1000 Genomes files — never written by pipeline
      cache/           intermediate PLINK/VCF files per notebook (01…)
      outputs/         primary outputs per notebook (01…)

    Roots are read from paths.local.yaml or paths.yaml at the project root.
    """

    ROOT: Path = field(default_factory=lambda: PROJECT_ROOT)

    _genomes_data: str = field(default_factory=lambda: _PATH_ROOTS["genomes_data"])
    _output_root: str = field(default_factory=lambda: _PATH_ROOTS["output"])
    _reports_root: str = field(default_factory=lambda: _PATH_ROOTS["reports"])

    def __post_init__(self):
        # Raw inputs — never written by the pipeline
        self.PLINK_MERGED = Path(f"{self._genomes_data}/plink/allchr")
        self.VCF_FILE     = Path(f"{self._genomes_data}/main_vcf/ALL_merged.vcf.gz")
        self.PANEL_FILE   = Path(f"{self._genomes_data}/raw_downloads/integrated_call_samples_v3.20130502.ALL.panel")

    # ------------------------------------------------------------------
    # Public accessors for the configured roots
    # ------------------------------------------------------------------

    @property
    def genomes_data(self) -> str:
        return self._genomes_data

    @property
    def output_root(self) -> str:
        return self._output_root

    # ------------------------------------------------------------------
    # Directory helpers
    # ------------------------------------------------------------------

    def cache_dir(self, notebook: str) -> Path:
        """Intermediate files for a notebook, e.g. PATHS.cache_dir('01_hard_filtering')."""
        return Path(f"{self._genomes_data}/cache/{notebook}")

    def outputs_dir(self, notebook: str) -> Path:
        """Primary outputs for a notebook, e.g. PATHS.outputs_dir('02_situational_filtering')."""
        return Path(f"{self._genomes_data}/outputs/{notebook}")

    # ------------------------------------------------------------------
    # Raw data
    # ------------------------------------------------------------------

    @property
    def DATA_DIR(self) -> Path:
        return Path(f"{self._genomes_data}/raw_downloads")

    @property
    def SAMPLES_CSV(self) -> Path:
        return Path(f"{self._genomes_data}/SEA_JPT_CN_subpopulation_samples.csv")

    @property
    def SAMPLES_LIST(self) -> Path:
        return Path(f"{self._genomes_data}/SEA_JPT_CN_subpopulation_samples_list.csv")

    # Backward-compat aliases
    @property
    def EAS_SAMPLES_CSV(self) -> Path:
        return self.SAMPLES_CSV

    @property
    def EAS_SAMPLES_LIST(self) -> Path:
        return self.SAMPLES_LIST

    @property
    def REPORTS_DIR(self) -> Path:
        return Path(f"{self._reports_root}/sea_jpt_cn")

    # ------------------------------------------------------------------
    # 01_hard_filtering
    # ------------------------------------------------------------------

    @property
    def PLINK_SNP_FILTERED(self) -> Path:
        """cache — SNP+biallelic filtered pfile (intermediate within 01)."""
        return self.cache_dir("01_hard_filtering") / "SEA_JPT_CN_SNP_filtered"

    @property
    def PLINK_MAF_FILTERED(self) -> Path:
        """output — MAF-filtered pfile (final product of 01, input to 02)."""
        return self.outputs_dir("01_hard_filtering") / "SEA_JPT_CN_MAF_filtered"

    # ------------------------------------------------------------------
    # 02_situational_filtering
    # ------------------------------------------------------------------

    @property
    def PLINK_HWE_FILTERED(self) -> Path:
        """cache — HWE-filtered pfile."""
        return self.cache_dir("02_situational_filtering") / "SEA_JPT_CN_HWE_filtered"

    @property
    def PLINK_UNIQUE_IDS(self) -> Path:
        """cache — variant-ID-standardised pfile."""
        return self.cache_dir("02_situational_filtering") / "SEA_JPT_CN_unique_ids"

    @property
    def PLINK_LD_PRUNED(self) -> Path:
        """output — LD-pruned pfile (final product of 02, input to 03+04)."""
        return self.outputs_dir("02_situational_filtering") / "SEA_JPT_CN_LD_pruned"

    # ------------------------------------------------------------------
    # 03_vcf_to_matrix
    # ------------------------------------------------------------------

    @property
    def GENOTYPE_MATRIX(self) -> Path:
        """output — int8 genotype matrix, numpy .npz (arrays: G, samples)."""
        return self.outputs_dir("03_vcf_to_matrix") / "genotype_matrix.npz"

    @property
    def GENOTYPE_MATRIX_COLS(self) -> Path:
        """output — SNP IDs matching columns of GENOTYPE_MATRIX, one per line."""
        return self.outputs_dir("03_vcf_to_matrix") / "genotype_matrix_cols.txt"

    # ------------------------------------------------------------------
    # 04a_statistical_snp_selection  (cross-notebook files named explicitly)
    # ------------------------------------------------------------------

    @property
    def STAT_SCORES(self) -> Path:
        """output — per-SNP statistical test scores (read by 05b, 05c, 08)."""
        return self.outputs_dir("04a_statistical_snp_selection") / "statistical_snp_scores.csv"

    @property
    def STAT_ML_DATA(self) -> Path:
        """output — stat-selected SNP genotype matrix (read by 05c, 08)."""
        return self.outputs_dir("04a_statistical_snp_selection") / "statistical_ml_data.csv"

    @property
    def STAT_ALL4_SNPS(self) -> Path:
        """output — SNPs significant in all 4 tests (read by 05b, 08b)."""
        return self.outputs_dir("04a_statistical_snp_selection") / "statistical_all4_snps.csv"

    # ------------------------------------------------------------------
    # fst_only/04b_fst_and_pca
    # ------------------------------------------------------------------

    @property
    def FST_RESULTS(self) -> Path:
        """cache — raw .fst.var files from PLINK2."""
        return self.cache_dir("fst_only/04b_fst_and_pca") / "SEA_JPT_CN_FST_RESULTS"

    @property
    def TOP_SNPS_FILE(self) -> Path:
        return self.outputs_dir("fst_only/04b_fst_and_pca") / "top_snps.txt"

    @property
    def TOP_SNPS_BED(self) -> Path:
        return self.outputs_dir("fst_only/04b_fst_and_pca") / "top_snps.bed"

    @property
    def FST_FILTERED(self) -> Path:
        """output — FST-selected pfile."""
        return self.outputs_dir("fst_only/04b_fst_and_pca") / "FST_FILTERED"

    @property
    def PCA_FILE(self) -> Path:
        return self.outputs_dir("fst_only/04b_fst_and_pca") / "FST_PCA"

    @property
    def ML_DATA(self) -> Path:
        """output — FST-filtered genotype matrix CSV (read by 05a, 05b, 06, 08b)."""
        return self.outputs_dir("fst_only/04b_fst_and_pca") / "ml_data_with_pop.csv"

    # ------------------------------------------------------------------
    # 05a (statistical_only) / 05b (FST) / 05c (FST+stat)  ML training
    # ------------------------------------------------------------------

    @property
    def ML_MODELS_FST(self) -> Path:
        return self.outputs_dir("fst/05b_fst_only_training")

    @property
    def ML_MODELS_DIR(self) -> Path:
        """Backward-compat alias → FST-only models dir."""
        return self.ML_MODELS_FST

    # ------------------------------------------------------------------
    # evaluation/
    # ------------------------------------------------------------------

    @property
    def EVAL_DIR(self) -> Path:
        return self.outputs_dir("evaluation")

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def get_absolute(self, relative_path: Path) -> Path:
        return self.ROOT / relative_path

    def ensure_output_dirs(self) -> None:
        dirs = [
            self.cache_dir("01_hard_filtering"),
            self.outputs_dir("01_hard_filtering"),
            self.cache_dir("02_situational_filtering"),
            self.outputs_dir("02_situational_filtering"),
            self.outputs_dir("03_vcf_to_matrix"),
            self.outputs_dir("04a_statistical_snp_selection"),
            self.cache_dir("fst_only/04b_fst_and_pca"),
            self.outputs_dir("fst_only/04b_fst_and_pca"),
            self.outputs_dir("fst/05b_fst_only_training"),
            self.outputs_dir("statistical/04a_snp_selection"),
            self.outputs_dir("statistical/05a_stat_training"),
            self.outputs_dir("statistical/05b_reduction"),
            self.outputs_dir("statistical/05c_fst_stat_training"),
            self.outputs_dir("evaluation"),
            self.REPORTS_DIR,
        ]
        for d in dirs:
            self.get_absolute(d).mkdir(parents=True, exist_ok=True)


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
    print(f"\n  genomes_data : {PATHS.genomes_data}")
    print(f"\n  cache/")
    print(f"    01         : {PATHS.cache_dir('01_hard_filtering')}")
    print(f"    02         : {PATHS.cache_dir('02_situational_filtering')}")
    print(f"    04         : {PATHS.cache_dir('fst_only/04b_fst_and_pca')}")
    print(f"\n  outputs/")
    print(f"    01         : {PATHS.outputs_dir('01_hard_filtering')}")
    print(f"    02         : {PATHS.outputs_dir('02_situational_filtering')}")
    print(f"    03         : {PATHS.outputs_dir('03_vcf_to_matrix')}")
    print(f"    03b        : {PATHS.outputs_dir('04a_statistical_snp_selection')}")
    print(f"    04         : {PATHS.outputs_dir('fst_only/04b_fst_and_pca')}")
    print(f"    05b FST    : {PATHS.ML_MODELS_FST}")
    print(f"\n  Hard filters : MAF≥{HARD_FILTERS.MIN_AF:.4f}  call_rate≥{HARD_FILTERS.MIN_CALL_RATE}")
    print(f"  LD pruning   : {SITUATIONAL_FILTERS.LD_WINDOW_KB}kb / r²<{SITUATIONAL_FILTERS.LD_R2_THRESHOLD}")
    print(f"  FST top-N    : {SITUATIONAL_FILTERS.FST_TOP_N}")
    print(f"\n  Populations  : {POPULATIONS.TARGET_POPS}  (n={POPULATIONS.NUM_SAMPLES})")
    print("=" * 60)


if __name__ == "__main__":
    print_config_summary()
