# Quick Start Guide

Get up and running with the bioinformatics pipeline in 10 minutes.

---

## 1. Setup (2 minutes)

### Verify Prerequisites

```bash
# Activate conda environment
conda activate GeneAnalysisET4596E

# Verify PLINK2 and bcftools installed
which plink2
which bcftools

# Check Python packages
python -c "import pandas, numpy, sklearn, xgboost; print('✓ All packages installed')"
```

### Navigate to Project

```bash
cd /path/to/BioinfoMidterm
ls -la  # Verify: should see notebooks/, scripts/, 1000genomes/
```

---

## 2. Run Hard Filters (3 minutes)

**Option A: Jupyter Notebook** (Interactive, Recommended)

```bash
jupyter notebook scripts/notebooks/01_hard_filtering.ipynb
```

1. Click first cell (setup)
2. Press Shift+Enter to run
3. Run remaining cells sequentially

**Option B: Command Line** (Automated)

```bash
cd scripts
./run_hard_filters.sh
```

**What it does:**

- ✓ Keeps SNPs only
- ✓ Keeps biallelic variants
- ✓ Removes rare variants (MAF < 0.0016)
- ✓ Filters by call rate

**Outputs:** `1000genomes/output/EAS_AND_SNP_filtered_data*`

---

## 3. Apply Situational Filters (3 minutes)

**Jupyter Notebook:**

```bash
jupyter notebook scripts/notebooks/02_situational_filtering.ipynb
```

**Command Line:**

```bash
cd scripts
./run_situational_filters.sh
```

**What it does:**

- ✓ HWE filter
- ✓ LD pruning
- ✓ Standardize variant IDs

---

## 4. Quick Analysis (2 minutes)

Run analysis notebooks for visualizations:

```bash
# PCA and FST
jupyter notebook scripts/notebooks/03_fst_and_pca.ipynb

# ML Training
jupyter notebook scripts/notebooks/04_ml_training.ipynb

# Model Evaluation
jupyter notebook scripts/notebooks/05_model_evaluation.ipynb
```

---

## Configuration Quick Reference

### Adjust Filter Thresholds

Edit `scripts/config.py`:

```python
# For lenient filtering:
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.01,       # Higher = keep only common variants
    MAX_GENO=0.05,     # Higher = allow more missing data
)

# For strict filtering:
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.001,      # Lower = keep rare variants too
    MAX_GENO=0.01,     # Lower = require high quality
)
```

Then restart notebook kernel and re-run.

### Change ML Parameters

```python
# In scripts/config.py
ML = MLConfig(
    TEST_SIZE=0.3,              # Change train/test split
    RF_N_ESTIMATORS=200,        # More Random Forest trees
    XGBOOST_N_ESTIMATORS=500,   # More XGBoost rounds
)
```

### Change Populations

```python
# In scripts/config.py
POPULATIONS = PopulationConfig(
    EAS_SUBPOPS=["CHB", "JPT"],  # Study only CHB and JPT
    NUM_SAMPLES=207,
    # ... update counts
)
```

---

## Common Tasks

### Just run everything

```bash
cd scripts
./run_hard_filters.sh && ./run_situational_filters.sh
python ml_training.py ../vcf_numeric_transposed.csv -o ../models
```

### Train only Random Forest

Edit `scripts/ml_training.py`:

```python
# Line with models to train
MODELS_TO_TRAIN = ["rf"]  # Only Random Forest
```

### Use different data

Edit `scripts/config.py`:

```python
PATHS = PathConfig(
    VCF_FILE="my_data.vcf.gz",
    ML_DATA="my_genotypes.csv",
    # ... other paths
)
```

### Save detailed logs

Add to any script:

```bash
./run_hard_filters.sh 2>&1 | tee run_hard_filters.log
```

---

## Output Locations

| Step | Output | Location |
|------|--------|----------|
| Hard filters | PLINK files | `1000genomes/output/EAS_AND_SNP_filtered_data*` |
| Situational filters | HWE/LD filtered | `1000genomes/output/EAS_AND_SNP_filtered_data_*_filtered*` |
| FST | FST values | `1000genomes/output/*.fst` |
| PCA | Eigenvectors | `1000genomes/output/*.eigenvec` |
| ML | Models | `models/` |
| Reports | Text summaries | `reports/` |

---

## Troubleshooting Quick Fixes

| Problem | Solution |
|---------|----------|
| `ModuleNotFoundError` | Run cell 1 in notebook first |
| `FileNotFoundError` | Check paths in config.py |
| `plink2: command not found` | `conda install -c bioconda plink2` |
| Notebook slow | Reduce `CV_FOLDS` or `RF_N_ESTIMATORS` in config |
| Results differ each run | Set `RANDOM_STATE=42` consistently |

---

## Full Documentation

For detailed information, see:

- **[RUNNING_SCRIPTS.md](RUNNING_SCRIPTS.md)** - Complete script reference
- **[CONFIGURATION.md](CONFIGURATION.md)** - All configurable parameters
- **[NOTEBOOK_GUIDE.md](NOTEBOOK_GUIDE.md)** - Notebook execution guide

---

## Architecture Overview

```
Input VCF
    ↓
Hard Filters (SNP, biallelic, MAF, call rate)
    ↓
Situational Filters (HWE, LD pruning)
    ↓
Analysis (PCA, FST)
    ↓
ML Training (RF, XGBoost, LR)
    ↓
Evaluation & Visualization
```

---

## Next Steps

1. ✅ Read [CONFIGURATION.md](CONFIGURATION.md) to understand filter thresholds
2. ✅ Run notebooks 01-02 for filtering
3. ✅ Run notebooks 03-05 for analysis and ML
4. ✅ Check results in output directories
5. ✅ Adjust config and iterate if needed

Happy analyzing! 🧬
