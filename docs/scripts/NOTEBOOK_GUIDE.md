# Notebook Configuration and Execution Guide

Step-by-step guide to configuring and running Jupyter notebooks for the bioinformatics pipeline.

---

## Before You Start

### Prerequisites

1. **Python environment** with required packages:

   ```bash
   conda activate GeneAnalysisET4596E
   # or your environment name
   ```

2. **Verify dependencies**:

   ```python
   import pandas, numpy, sklearn, xgboost, matplotlib
   print("All dependencies installed ✓")
   ```

3. **Check file structure**:

   ```bash
   cd /path/to/BioinfoMidterm
   ls -la 1000genomes/
   ls -la docs/
   ls -la scripts/
   ```

---

## Working Directory Setup

**Critical:** All notebooks must run from the project root directory.

### Method 1: Automatic (Recommended)

The first cell in each notebook automatically sets the working directory:

```python
import sys, os
from pathlib import Path

# Auto-detect and change to project root
os.chdir(Path.cwd().parent.parent)
sys.path.insert(0, str(Path.cwd() / "scripts"))
print(f"Working directory: {os.getcwd()}")
```

Just execute cell 1 before running other cells.

### Method 2: Manual Setup

If automatic setup doesn't work:

```bash
# In terminal
cd /path/to/BioinfoMidterm

# Then open notebook from here
jupyter notebook scripts/notebooks/01_hard_filtering.ipynb
```

### Verify Setup

```python
# In notebook cell
import os
print(f"Current directory: {os.getcwd()}")
print(f"Has 1000genomes: {os.path.exists('1000genomes')}")
print(f"Has scripts: {os.path.exists('scripts')}")
```

All three should show `True` or match the project root.

---

## Notebook Execution Workflow

Execute notebooks in order. Each notebook depends on outputs from previous ones.

### 01_hard_filtering.ipynb

**Purpose:** Apply SNP-only, biallelic, MAF, and call rate filters

**Key Configuration:**

```python
# Cell 1: These come from config.py
HARD_FILTERS.MIN_AF          # Minimum allele frequency
HARD_FILTERS.MAX_GENO        # Maximum genotype missingness
POPULATIONS.EAS_SUBPOPS      # ["CHB", "JPT", "KHV"]
POPULATIONS.NUM_SAMPLES      # 306
```

**To Customize:**

```python
# In any cell (after imports)
from config import HARD_FILTERS
from dataclasses import replace

# Override default MIN_AF
custom_filters = replace(HARD_FILTERS, MIN_AF=0.01)

# Use in filter call
result = filter_maf(
    input_pfile=...,
    output_pfile=...,
    min_af=custom_filters.MIN_AF,
)
```

**Outputs:**

- `1000genomes/output/EAS_AND_SNP_filtered_data.*` (PLINK format)
- Frequency statistics and visualizations

**Expected Runtime:** 5-15 minutes

---

### 02_situational_filtering.ipynb

**Purpose:** Apply HWE, LD pruning, and variant ID standardization

**Key Configuration:**

```python
SITUATIONAL_FILTERS.HWE_P_THRESHOLD     # 1e-6 by default
SITUATIONAL_FILTERS.LD_WINDOW_KB        # 1000 kb
SITUATIONAL_FILTERS.LD_R2_THRESHOLD     # 0.1
```

**To Customize HWE:**

```python
from situational_filters import filter_hwe

# Lenient HWE filter
result = filter_hwe(
    input_pfile="data",
    output_pfile="output",
    hwe_p_threshold=0.05,  # Less stringent than default 1e-6
    verbose=True,
)
```

**To Customize LD Pruning:**

```python
from situational_filters import calculate_ld_prune_list

# Aggressive LD pruning
prune_list = calculate_ld_prune_list(
    input_pfile="data",
    window_kb=500,         # Smaller window = more aggressive
    step=50,
    r2_threshold=0.5,      # Higher threshold = less aggressive
)
```

**Outputs:**

- Filtered PLINK files with HWE and LD applied
- LD pruning lists (*.prune.in,*.prune.out)

**Expected Runtime:** 10-30 minutes

---

### 03_fst_and_pca.ipynb

**Purpose:** Calculate FST values and perform PCA

**Key Configuration:**

```python
# FST calculation compares populations
POPULATIONS.EAS_SUBPOPS  # Populations to compare

# PCA parameters are automatic
# (uses PLINK2 defaults)
```

**To Change PCA Settings:**

```python
from situational_filters import calculate_pca

# Generate more principal components
result = calculate_pca(
    input_pfile="data",
    output_prefix="results/pca",
    n_pcs=20,  # Instead of default 10
    ld_prune=True,
)
```

**To Customize FST Selection:**

```python
from fst_selection import select_top_fst_variants

# Select top 50 variants instead of 100
top_variants = select_top_fst_variants(
    fst_file="results/fst_values.fst.var",
    n_variants=50,  # Changed from default
    output_file="results/top_50_variants.txt",
)
```

**Outputs:**

- FST values between populations
- PCA eigenvectors and eigenvalues
- Visualization plots

**Expected Runtime:** 5-10 minutes

---

### 04_ml_training.ipynb

**Purpose:** Train Random Forest, XGBoost, and Logistic Regression models

**Key Configuration:**

```python
ML.TEST_SIZE              # 0.2 (20% test)
ML.RANDOM_STATE           # 42 (reproducibility)
ML.CV_FOLDS               # 5 (cross-validation)
ML.RF_N_ESTIMATORS        # 100 (Random Forest trees)
ML.XGBOOST_N_ESTIMATORS   # 300 (XGBoost rounds)
```

**To Change Train/Test Split:**

```python
from sklearn.model_selection import train_test_split

# Use 30% test instead of 20%
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.3,      # Changed from default 0.2
    stratify=y,
    random_state=42,
)
```

**To Tune Random Forest:**

```python
from sklearn.ensemble import RandomForestClassifier

# More trees, deeper trees = higher accuracy (risk of overfitting)
rf = RandomForestClassifier(
    n_estimators=500,      # Changed from 100
    max_depth=20,          # Allow deeper trees
    min_samples_split=2,
    n_jobs=-1,
)
```

**To Tune XGBoost:**

```python
from xgboost import XGBClassifier

# Conservative settings
xgb = XGBClassifier(
    n_estimators=100,           # Fewer rounds
    learning_rate=0.1,          # Higher learning rate = faster convergence
    max_depth=4,                # Shallower trees
    subsample=0.8,
    colsample_bytree=0.8,
)
```

**Outputs:**

- Trained models (*.pkl files)
- Feature importances (CSV)
- Classification reports (TXT)

**Expected Runtime:** 10-30 minutes

---

### 05_model_evaluation.ipynb

**Purpose:** Evaluate, interpret, and visualize trained models

**Configuration:** Uses pre-trained models from notebook 04

**To Load Specific Model:**

```python
from ml_training import load_model

# Load only Random Forest
rf, metadata = load_model("models/random_forest_full.pkl")
print(metadata)
```

**To Compare Models Selectively:**

```python
# Compare only specific models
models_to_compare = [
    "random_forest_full.pkl",
    "xgboost.pkl",
]

results = {}
for model_name in models_to_compare:
    model, metadata = load_model(f"models/{model_name}")
    results[model_name] = model.score(X_test, y_test)
```

**To Generate Custom Plots:**

```python
import matplotlib.pyplot as plt

# ROC curve for specific class
from sklearn.metrics import roc_curve, auc

fpr, tpr, _ = roc_curve(y_test_binary[:, 0], y_proba[:, 0])
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.3f})')
plt.show()
```

**Outputs:**

- Confusion matrices
- ROC/PR curves
- Feature importance plots
- Model comparison visualizations

**Expected Runtime:** 5-10 minutes

---

## Common Configuration Changes

### Change All Population Thresholds

Edit `scripts/config.py`:

```python
# Make filters more stringent
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.01,       # Instead of 0.0016
    MAX_GENO=0.01,     # Instead of 0.02
    SNPS_ONLY=True,
    MAX_ALLELES=2,
)
```

Then restart notebook kernel and re-run all cells.

### Change Train/Test Split Globally

```python
# In scripts/config.py
ML = MLConfig(
    TEST_SIZE=0.3,     # Changed from 0.2
    RANDOM_STATE=42,
    # ... other params
)
```

### Use Top Features Only

In notebook 04, modify the feature selection:

```python
# Use only top 50 features instead of all
if feature_importance is not None:
    top_features = feature_importance.head(50)['feature'].tolist()
    X_train_top = X_train[top_features]
    X_test_top = X_test[top_features]
    
    # Train on reduced feature set
    rf = RandomForestClassifier(...)
    rf.fit(X_train_top, y_train)
```

---

## Troubleshooting Notebooks

### Problem: "ModuleNotFoundError: No module named 'config'"

**Solution:**

1. Execute cell 1 first
2. Verify working directory:

   ```python
   import os
   print(os.getcwd())
   ```

3. Should show: `/path/to/BioinfoMidterm`

### Problem: "FileNotFoundError: No such file or directory"

**Solution:**

1. Check file paths in config.py
2. Verify file exists:

   ```python
   from pathlib import Path
   from config import PATHS
   print(Path(PATHS.VCF_FILE).exists())
   ```

3. Update path if needed

### Problem: Notebook runs slowly / crashes with memory error

**Solution:**

1. Use subset of data:

   ```python
   # Take first 100k variants
   X_subset = X.iloc[:, :100000]
   ```

2. Reduce cross-validation folds:

   ```python
   ML = replace(ML, CV_FOLDS=3)
   ```

3. Reduce model complexity:

   ```python
   RF_N_ESTIMATORS = 50  # Instead of 100
   ```

### Problem: Plots not showing

**Solution:**

```python
import matplotlib.pyplot as plt

# Add after imports
%matplotlib inline

# Or explicitly show
plt.show()
```

### Problem: Results differ between runs

**Solution:** Ensure reproducibility:

```python
# Set random seeds
import numpy as np
import random

RANDOM_STATE = 42
random.seed(RANDOM_STATE)
np.random.seed(RANDOM_STATE)
```

---

## Performance Tips

1. **Use ML test datasets first** - Test on 1000 samples before full dataset
2. **Skip cross-validation initially** - Use simple train/test first
3. **Reduce PCs for visualization** - Plot PC1 vs PC2 only initially
4. **Use sparse data formats** - For very large datasets, use scipy sparse matrices
5. **Cache intermediate results** - Save outputs between cells to avoid recomputation

---

## Output Organization

Notebooks save outputs to:

```
results/
├── hard_filtered/               # 01_hard_filtering.ipynb
│   ├── *.pgen, *.pvar, *.psam
│   ├── *_info.afreq
│   └── visualizations/
├── situational_filtered/        # 02_situational_filtering.ipynb
│   ├── *.pgen, *.pvar, *.psam
│   ├── *.prune.in, *.prune.out
│   └── hwe_results/
├── fst_and_pca/                 # 03_fst_and_pca.ipynb
│   ├── *.fst files
│   ├── eigenvec, eigenval
│   └── pca_plots/
├── ml_models/                   # 04_ml_training.ipynb
│   ├── random_forest.pkl
│   ├── xgboost.pkl
│   ├── logistic_regression.pkl
│   ├── feature_importances.csv
│   └── classification_reports/
└── evaluation/                  # 05_model_evaluation.ipynb
    ├── confusion_matrices/
    ├── roc_curves/
    ├── feature_importance_plots/
    └── summary_report.txt
```

---

## Quick Reference: Key Variables

### Notebook 1-2: Filtering

```python
# From config
PATHS.VCF_FILE                   # Input VCF
HARD_FILTERS.MIN_AF              # MAF threshold
SITUATIONAL_FILTERS.HWE_P_THRESHOLD   # HWE p-value

# Generated in notebook
step1_output                      # After SNP/biallelic filter
step2_output                      # After MAF filter
```

### Notebook 3: Analysis

```python
# Generated
fst_results                       # FST between populations
pca_eigenvec                      # PCA coordinates
pca_eigenval                      # PCA explained variance
```

### Notebook 4-5: ML

```python
# Generated
X, y                              # Features and labels
models                            # Trained model dict
feature_importance                # Feature rankings
results                           # Performance metrics
```

---

## Next Steps

1. ✅ Configure paths in `docs/scripts/CONFIGURATION.md`
2. ✅ Review filter thresholds for your use case
3. ✅ Run notebook 01_hard_filtering.ipynb
4. ✅ Check outputs in results/hard_filtered/
5. ✅ Proceed to notebook 02 with same approach
6. ✅ Iterate: adjust config → run notebook → evaluate outputs
