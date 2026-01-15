# Configuration Guide

Complete reference for configuring scripts and notebooks. All configuration is centralized in `scripts/config.py`.

## Overview

The `config.py` module uses frozen dataclasses to define:

- File paths
- Filter thresholds
- Population definitions
- ML parameters
- PLINK2 commands

Changes to `config.py` affect all scripts and notebooks.

---

## File Paths

Edit `PathConfig` dataclass to change input/output locations.

### Default Configuration

```python
# scripts/config.py

@dataclass(frozen=True)
class PathConfig:
    # Input directories
    VCF_DIR = "1000genomes"
    PLINK_DIR = "1000genomes/output"
    
    # VCF file
    VCF_FILE = "1000genomes/snps_only.vcf"
    
    # Panel information
    PANEL_FILE = "docs/integrated_call_samples_v3.20130502.ALL.panel"
    
    # Output directories
    PLINK_SNP_FILTERED = "1000genomes/output/EAS_AND_SNP_filtered_data"
    PLINK_MAF_FILTERED = "1000genomes/output/EAS_AND_SNP_filtered_data_MAF_filtered"
    PLINK_HWE_FILTERED = "1000genomes/output/EAS_AND_SNP_filtered_data_MAF_filtered_hardy"
    
    # Sample lists
    EAS_SAMPLES_CSV = "1000genomes/EAS_subpopulation_samples.csv"
    EAS_SAMPLES_LIST = "1000genomes/EAS_subpopulation_samples_list.csv"
    
    # ML data
    ML_DATA = "vcf_numeric_transposed.csv"
    ML_MODELS_DIR = Path("models")
```

### Changing Paths

To use different data:

```python
# Option 1: Edit config.py directly
PATHS = PathConfig(
    VCF_FILE="my_custom_vcf.vcf.gz",
    PLINK_DIR="my_output_dir",
    # ...
)

# Option 2: Override in notebook
import os
from config import PATHS

# Create new paths object with custom values
custom_paths = {
    "VCF_FILE": "new_vcf.vcf.gz",
    "PLINK_DIR": "new_output",
}

# Use with scripts
```

---

## Hard Filter Thresholds

Adjust SNP quality control parameters.

### Default Configuration

```python
@dataclass(frozen=True)
class HardFilterThresholds:
    """Strict QC filters applied to all datasets."""
    
    # Minimum allele frequency (removes rare variants)
    # 1/612 alleles = ~0.0016 (singletons in 306 samples)
    MIN_AF = 0.0016
    
    # Maximum genotype missingness rate
    # Allow up to 2% missing genotypes
    MAX_GENO = 0.02
    
    # PLINK2 SNP filter
    SNPS_ONLY = True
    
    # PLINK2 biallelic filter
    MAX_ALLELES = 2
```

### When to Change

| Parameter | Default | When to Change |
|-----------|---------|---|
| `MIN_AF` | 0.0016 | Rare variant studies (lower) vs common variants (higher) |
| `MAX_GENO` | 0.02 | Strict QC (0.01) vs lenient (0.05) |
| `SNPS_ONLY` | True | Set False to include indels |
| `MAX_ALLELES` | 2 | Increase for multiallelic SNPs |

### Example: Strict Quality Control

```python
# config.py
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.05,       # Only common variants (>5%)
    MAX_GENO=0.01,     # Strict: <1% missing
    SNPS_ONLY=True,
    MAX_ALLELES=2,
)
```

### Example: Rare Variant Study

```python
# config.py
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.001,      # Rare variants
    MAX_GENO=0.05,     # Lenient: 5% missing OK
    SNPS_ONLY=False,   # Include indels
    MAX_ALLELES=2,
)
```

---

## Situational Filter Thresholds

Context-dependent quality control parameters.

### Default Configuration

```python
@dataclass(frozen=True)
class SituationalFilterThresholds:
    """Optional filters applied based on study design."""
    
    # Hardy-Weinberg Equilibrium
    HWE_P_THRESHOLD = 1e-6
    HWE_MODE = "keep-fewhet"  # Conservative for population studies
    
    # LD Pruning (linkage disequilibrium)
    LD_WINDOW_KB = 1000        # 1000 kb window
    LD_STEP = 100              # Step every 100 variants
    LD_R2_THRESHOLD = 0.1      # Remove if r² > 0.1
    
    # Variant ID standardization
    SET_VAR_ID_FORMAT = "@:#:\$r:\$a"
```

### When to Apply Each Filter

| Filter | Use For | Threshold |
|--------|---------|-----------|
| **HWE** | Population structure, GWAS | 1e-6 (stringent) to 0.05 (liberal) |
| **LD Pruning** | PCA, association tests | r² = 0.1 to 0.5 |
| **Variant ID** | Merging datasets | Standard: @:#:$r:$a |

### Example: PCA Analysis

```python
# config.py
SITUATIONAL_FILTERS = SituationalFilterThresholds(
    HWE_P_THRESHOLD=1e-6,      # Remove HWE violations
    LD_WINDOW_KB=1000,
    LD_STEP=100,
    LD_R2_THRESHOLD=0.1,       # Aggressive LD pruning for PCA
)
```

### Example: Ancestry Inference

```python
# config.py
SITUATIONAL_FILTERS = SituationalFilterThresholds(
    HWE_P_THRESHOLD=1e-6,
    LD_WINDOW_KB=0,            # No LD pruning
    LD_STEP=0,
    LD_R2_THRESHOLD=0,
    # Focus on informative SNPs instead
)
```

---

## Population Configuration

Define population groups for analysis.

### Default Configuration

```python
@dataclass(frozen=True)
class PopulationConfig:
    """Population definitions."""
    
    # EAS subpopulations (East Asian)
    EAS_SUBPOPS = ["CHB", "JPT", "KHV"]  # Chinese, Japanese, Kinh
    NUM_SAMPLES = 306                     # 306 samples total
    
    # Sample counts per population
    EAS_POP_COUNTS = {
        "CHB": 103,  # Han Chinese in Beijing
        "JPT": 104,  # Japanese in Tokyo
        "KHV": 99,   # Kinh in Ho Chi Minh City
    }
```

### Changing Populations

```python
# Study only CHB and JPT
POPULATIONS = PopulationConfig(
    EAS_SUBPOPS=["CHB", "JPT"],
    NUM_SAMPLES=207,
    EAS_POP_COUNTS={"CHB": 103, "JPT": 104},
)
```

---

## Machine Learning Configuration

Model training and evaluation parameters.

### Default Configuration

```python
@dataclass(frozen=True)
class MLConfig:
    """Machine learning parameters."""
    
    # Data split
    TEST_SIZE = 0.2
    RANDOM_STATE = 42
    
    # Random Forest
    RF_N_ESTIMATORS = 100
    RF_MAX_DEPTH = None
    RF_MIN_SAMPLES_SPLIT = 2
    
    # XGBoost
    XGBOOST_N_ESTIMATORS = 300
    XGBOOST_LEARNING_RATE = 0.05
    XGBOOST_MAX_DEPTH = 6
    XGBOOST_SUBSAMPLE = 0.8
    
    # Logistic Regression
    LR_C = 1.0
    LR_SOLVER = "lbfgs"
    LR_MAX_ITER = 1000
    
    # Cross-validation
    CV_FOLDS = 5
```

### Tuning Models

#### Random Forest

```python
# For interpretability (shallow trees)
RF_N_ESTIMATORS = 100
RF_MAX_DEPTH = 10

# For accuracy (deeper trees)
RF_N_ESTIMATORS = 500
RF_MAX_DEPTH = None
```

#### XGBoost

```python
# Conservative (less overfit risk)
XGBOOST_N_ESTIMATORS = 100
XGBOOST_LEARNING_RATE = 0.1
XGBOOST_MAX_DEPTH = 4

# Aggressive (higher accuracy potential)
XGBOOST_N_ESTIMATORS = 500
XGBOOST_LEARNING_RATE = 0.05
XGBOOST_MAX_DEPTH = 8
```

---

## PLINK2 Configuration

Control PLINK2 command execution.

### Default Configuration

```python
@dataclass(frozen=True)
class Plink2Config:
    """PLINK2 execution settings."""
    
    # Path to plink2 executable
    PLINK2_CMD = "plink2"
    
    # Number of threads
    N_THREADS = -1  # Auto-detect (all cores)
    
    # Memory limit (GB)
    MEMORY_GB = None  # No limit
    
    # Temporary files
    KEEP_TEMP = False
```

### Customization

```python
# Single thread (for consistency)
PLINK2 = Plink2Config(
    PLINK2_CMD="plink2",
    N_THREADS=1,
)

# Limited memory
PLINK2 = Plink2Config(
    PLINK2_CMD="plink2",
    N_THREADS=4,
    MEMORY_GB=8,
)
```

---

## How to Apply Configuration Changes

### Method 1: Edit config.py Directly

```python
# scripts/config.py

# Change hard filter thresholds
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.01,           # Changed
    MAX_GENO=0.02,
    SNPS_ONLY=True,
    MAX_ALLELES=2,
)

# Change ML parameters
ML = MLConfig(
    TEST_SIZE=0.3,         # Changed from 0.2
    RANDOM_STATE=42,
    RF_N_ESTIMATORS=200,   # Changed from 100
    # ... other parameters
)
```

### Method 2: Override in Notebook

```python
# In notebook
from config import HARD_FILTERS, PathConfig
from dataclasses import replace

# Create modified config
custom_filters = replace(HARD_FILTERS, MIN_AF=0.01)

# Use in function calls
result = filter_maf(
    input_pfile="data",
    output_pfile="output",
    min_af=custom_filters.MIN_AF,  # Use custom value
)
```

### Method 3: Command Line Arguments

```bash
# Pass parameters as script arguments
python hard_filters.py \
  --input data.pgen \
  --output results \
  --maf 0.01 \
  --geno 0.01

python ml_training.py data.csv \
  --test-size 0.3 \
  --n-estimators 200
```

---

## Configuration Best Practices

1. **Version Control** - Keep config.py in git to track changes

   ```bash
   git commit scripts/config.py -m "Update filter thresholds"
   ```

2. **Document Changes** - Add comments explaining why

   ```python
   # Rare variant study: lower MAF threshold
   MIN_AF = 0.001
   ```

3. **Test Changes** - Validate on small dataset first

   ```bash
   # Test with 100 samples
   plink2 --pfile data --keep sample_100.txt ...
   ```

4. **Create Presets** - Define common configurations

   ```python
   # Strict QC preset
   STRICT_CONFIG = HardFilterThresholds(MIN_AF=0.05, MAX_GENO=0.01)
   
   # Lenient preset
   LENIENT_CONFIG = HardFilterThresholds(MIN_AF=0.001, MAX_GENO=0.05)
   ```

5. **Backup Original** - Keep copy of default config

   ```bash
   cp scripts/config.py scripts/config.py.bak
   ```

---

## Common Configuration Scenarios

### Scenario 1: Quick Prototyping

```python
# Fast, permissive settings
HARD_FILTERS = HardFilterThresholds(MIN_AF=0.05, MAX_GENO=0.05)
SITUATIONAL_FILTERS = SituationalFilterThresholds(HWE_P_THRESHOLD=0.05)
ML = MLConfig(CV_FOLDS=3)  # Faster CV
```

### Scenario 2: Publication-Ready Analysis

```python
# Stringent, well-documented settings
HARD_FILTERS = HardFilterThresholds(MIN_AF=0.0016, MAX_GENO=0.02)
SITUATIONAL_FILTERS = SituationalFilterThresholds(HWE_P_THRESHOLD=1e-6)
ML = MLConfig(CV_FOLDS=10, RANDOM_STATE=42)  # Robust validation
```

### Scenario 3: Rare Variant Study

```python
# Include rare variants
HARD_FILTERS = HardFilterThresholds(MIN_AF=0.001, MAX_GENO=0.05)
SITUATIONAL_FILTERS = SituationalFilterThresholds(HWE_P_THRESHOLD=0.05)
```

### Scenario 4: Large-Scale GWAS

```python
# Common variants only, aggressive LD pruning
HARD_FILTERS = HardFilterThresholds(MIN_AF=0.05, MAX_GENO=0.01)
SITUATIONAL_FILTERS = SituationalFilterThresholds(
    LD_R2_THRESHOLD=0.5,  # More aggressive
    HWE_P_THRESHOLD=1e-6,
)
```

---

## Validation

Verify configuration is correctly applied:

```python
# In notebook
from config import PATHS, HARD_FILTERS, ML

print("Paths:")
print(f"  VCF: {PATHS.VCF_FILE}")
print(f"  Output: {PATHS.PLINK_DIR}")

print("\nHard Filters:")
print(f"  MIN_AF: {HARD_FILTERS.MIN_AF}")
print(f"  MAX_GENO: {HARD_FILTERS.MAX_GENO}")

print("\nML Config:")
print(f"  TEST_SIZE: {ML.TEST_SIZE}")
print(f"  CV_FOLDS: {ML.CV_FOLDS}")
```

---

## Reference: All Config Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `VCF_FILE` | `1000genomes/snps_only.vcf` | str | Input VCF file |
| `MIN_AF` | 0.0016 | float | Minimum allele frequency |
| `MAX_GENO` | 0.02 | float | Maximum genotype missingness |
| `HWE_P_THRESHOLD` | 1e-6 | float | HWE p-value threshold |
| `LD_R2_THRESHOLD` | 0.1 | float | LD pruning r² threshold |
| `TEST_SIZE` | 0.2 | float | Train/test split ratio |
| `CV_FOLDS` | 5 | int | Cross-validation folds |
| `RF_N_ESTIMATORS` | 100 | int | Random Forest trees |
| `XGBOOST_N_ESTIMATORS` | 300 | int | XGBoost boosting rounds |
