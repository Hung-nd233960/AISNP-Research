# Configuration Reference

Complete reference for all configurable parameters in `scripts/config.py`.

## Table of Contents

1. [Path Configuration](#path-configuration)
2. [Hard Filter Thresholds](#hard-filter-thresholds)
3. [Situational Filter Thresholds](#situational-filter-thresholds)
4. [PLINK2 Configuration](#plink2-configuration)
5. [Population Configuration](#population-configuration)
6. [ML Configuration](#ml-configuration)
7. [Modifying Configuration](#modifying-configuration)

---

## Path Configuration

**Class**: `PathConfig`

Defines all input/output file paths.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ROOT` | `/home/Plutonium/Documents/BioinfoMidterm` | Project root directory |
| `DATA_DIR` | `1000genomes` | Main data directory |
| `OUTPUT_DIR` | `1000genomes/output` | PLINK output files |
| `VCF_DIR` | `1000genomes/vcf` | VCF file directory |
| `REPORTS_DIR` | `reports` | Text reports |
| `VCF_FILE` | `1000genomes/main_vcf/ALL_merged.vcf.gz` | Input VCF |
| `PANEL_FILE` | `1000genomes/main_vcf/integrated_call_samples_v3.20130502.ALL.panel` | Sample panel |
| `EAS_SAMPLES_CSV` | `1000genomes/EAS_subpopulation_samples.csv` | EAS sample list |
| `PLINK_LD_PRUNED` | `1000genomes/output/EAS_FINAL_DATA_FOR_FST` | Final LD-pruned pfile |
| `ML_DATA` | `1000genomes/vcf/vcf_numeric_transposed_with_population.csv` | ML data file |

---

## Hard Filter Thresholds

**Class**: `HardFilterThresholds`

Strict quality control parameters applied to all datasets.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `SNP_ONLY` | `True` | bool | Remove indels, CNVs, structural variants |
| `MAX_ALLELES` | `2` | int | Maximum alleles per variant (biallelic) |
| `MIN_AF` | `0.0016` | float | Minimum minor allele frequency |
| `MIN_CALL_RATE` | `0.95` | float | Minimum genotyping completeness |
| `MIN_OBS_CT` | `1` | int | Minimum observed allele count |

### Rationale

**SNP_ONLY**:

- Indels have higher error rates
- SNPs are more reliable for ancestry inference
- Simplifies downstream analysis

**MAX_ALLELES = 2**:

- Biallelic variants are more reliable
- Triallelic+ variants often indicate errors
- Simplifies genotype encoding (0, 1, 2)

**MIN_AF = 0.0016**:

- Calculated as 1/(2×306) = 1/612 ≈ 0.0016
- Removes singletons (alleles appearing only once)
- Reduces noise from rare variants

**MIN_CALL_RATE = 0.95**:

- Ensures high-quality genotyping
- Reduces bias from missing data
- Industry standard threshold

---

## Situational Filter Thresholds

**Class**: `SituationalFilterThresholds`

Context-dependent parameters that may vary by study design.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `HWE_P_THRESHOLD` | `1e-6` | float | Hardy-Weinberg p-value threshold |
| `HWE_FILTER_MODE` | `"keep-fewhet"` | str | HWE filter mode |
| `LD_WINDOW_KB` | `1000` | int | LD pruning window size (kb) |
| `LD_STEP` | `1` | int | LD pruning step size (variants) |
| `LD_R2_THRESHOLD` | `0.1` | float | LD R² correlation cutoff |
| `FST_TOP_N` | `1000` | int | Top FST variants to select |
| `FST_MIN_VALUE` | `0.0` | float | Minimum FST value |

### HWE Filter Modes

| Mode | Description |
|------|-------------|
| `"keep-fewhet"` | Keep variants with fewer heterozygotes than expected (default) |
| `"keep-all"` | Keep all variants passing threshold |
| `None` | Disable HWE filtering |

### LD Pruning Presets

| Use Case | Window | Step | R² | Rationale |
|----------|--------|------|-----|-----------|
| **Ancestry/PCA** | 1000kb | 1 | 0.1 | Aggressive pruning, independent SNPs |
| Moderate | 500kb | 5 | 0.2 | Balanced approach |
| Association | 50kb | 5 | 0.5 | Light pruning, preserve info |

### Choosing Parameters

**For ancestry inference** (current settings):

- Aggressive LD pruning (R² = 0.1) removes correlated SNPs
- Large window (1000kb) captures long-range LD
- Reduces redundancy for PCA and ML

**For GWAS/association studies**:

- Light LD pruning (R² = 0.5) preserves more variants
- Smaller window (50kb) for local LD patterns
- Keeps tag SNPs for causal variant detection

---

## PLINK2 Configuration

**Class**: `Plink2Config`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EXECUTABLE` | `"plink2"` | PLINK2 executable name |
| `THREADS` | `8` | Number of CPU threads |
| `MEMORY_MB` | `16000` | Memory allocation (MB) |

Adjust based on your system resources:

```python
# For smaller systems
THREADS = 4
MEMORY_MB = 8000

# For larger servers
THREADS = 16
MEMORY_MB = 32000
```

---

## Population Configuration

**Class**: `PopulationConfig`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EAS_SUBPOPS` | `("CHB", "JPT", "KHV")` | Target populations |
| `SUPER_POP` | `"EAS"` | Super-population code |
| `NUM_SAMPLES` | `306` | Total sample count |

### Population Codes (1000 Genomes)

| Code | Population | Region |
|------|------------|--------|
| CHB | Han Chinese in Beijing | China |
| JPT | Japanese in Tokyo | Japan |
| KHV | Kinh in Ho Chi Minh City | Vietnam |
| CHS | Southern Han Chinese | China |
| CDX | Chinese Dai in Xishuangbanna | China |

---

## ML Configuration

**Class**: `MLConfig`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TEST_SIZE` | `0.2` | Test set fraction |
| `RANDOM_STATE` | `42` | Random seed |
| `CV_FOLDS` | `5` | Cross-validation folds |
| `RF_N_ESTIMATORS` | `200` | Random Forest trees |
| `RF_MAX_DEPTH` | `None` | RF max tree depth |
| `XGB_N_ESTIMATORS` | `200` | XGBoost trees |
| `XGB_MAX_DEPTH` | `6` | XGBoost max depth |
| `XGB_LEARNING_RATE` | `0.1` | XGBoost learning rate |
| `LR_MAX_ITER` | `1000` | Logistic Regression iterations |
| `TOP_N_FEATURES` | `25` | Top features for analysis |

---

## Modifying Configuration

### Option 1: Edit config.py directly

```python
# scripts/config.py

@dataclass(frozen=True)
class SituationalFilterThresholds:
    LD_R2_THRESHOLD: float = 0.2  # Changed from 0.1
```

### Option 2: Override in notebook

```python
from config import SITUATIONAL_FILTERS
from dataclasses import replace

# Create modified config
custom_filters = replace(SITUATIONAL_FILTERS, LD_R2_THRESHOLD=0.2)
```

### Option 3: Environment variables (advanced)

```python
import os

@dataclass
class PathConfig:
    ROOT: Path = Path(os.getenv('BIOINFO_ROOT', '/home/Plutonium/Documents/BioinfoMidterm'))
```

---

## Printing Current Configuration

```python
from config import print_config_summary

print_config_summary()
```

Output:

```
============================================================
CONFIGURATION SUMMARY
============================================================

--- HARD FILTERS (Always Applied) ---
  SNP-only:      True
  Max alleles:   2
  Min AF:        0.0016
  Min call rate: 0.95

--- SITUATIONAL FILTERS (Context-Dependent) ---
  HWE threshold: 1e-06
  HWE mode:      keep-fewhet
  LD window:     1000kb
  LD step:       1
  LD R² cutoff:  0.1
  FST top N:     1000

--- POPULATIONS ---
  Subpopulations: ('CHB', 'JPT', 'KHV')
  Num samples:    306
============================================================
```

---

*See also: [PIPELINE.md](PIPELINE.md), [STATISTICAL_TESTS.md](STATISTICAL_TESTS.md)*
