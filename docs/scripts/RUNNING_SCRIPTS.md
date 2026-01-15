# Running Scripts and Notebooks

A comprehensive guide to executing the bioinformatics pipeline for 1000 Genomes EAS population analysis.

## Quick Start

### Option 1: Jupyter Notebooks (Recommended for Exploration)

Run notebooks sequentially in VS Code or Jupyter:

```bash
cd /path/to/BioinfoMidterm
jupyter notebook scripts/notebooks/
```

**Execution order:**

1. `01_hard_filtering.ipynb` - Apply hard filters
2. `02_situational_filtering.ipynb` - Apply situational filters
3. `03_fst_and_pca.ipynb` - FST and PCA analysis
4. `04_ml_training.ipynb` - Train ML models
5. `05_model_evaluation.ipynb` - Evaluate and interpret models

### Option 2: Command Line Scripts (For Automation)

Run filtering workflows via bash scripts:

```bash
cd /path/to/BioinfoMidterm/scripts

# Hard filters
./run_hard_filters.sh

# Situational filters
./run_situational_filters.sh
```

### Option 3: Python CLI (For Specific Tasks)

Run individual Python modules:

```bash
cd /path/to/BioinfoMidterm/scripts

# ML training
python ml_training.py /path/to/data.csv -o /path/to/output

# FST analysis
python fst_selection.py --input data.pgen --output fst_results

# Panel EDA
python panel_eda.py 1000genomes/integrated_call_samples_v3.20130502.ALL.panel --super-pop EAS
```

---

## Detailed Script Usage

### Filtering Workflows

#### Hard Filters (SNP-only, Biallelic, MAF, Call Rate)

**Bash CLI:**

```bash
./run_hard_filters.sh
```

**Python CLI:**

```bash
python hard_filters.py \
  --input 1000genomes/output/EAS_AND_SNP_filtered_data.pgen \
  --output results/hard_filtered \
  --maf 0.0016 \
  --geno 0.02
```

**In Notebook:**

```python
from hard_filters import apply_all_hard_filters

result = apply_all_hard_filters(
    input_pfile="1000genomes/output/input",
    output_pfile="results/hard_filtered",
    verbose=True
)
```

#### Situational Filters (HWE, LD Pruning, FST, PCA)

**Bash CLI:**

```bash
./run_situational_filters.sh
```

**Python CLI:**

```bash
python situational_filters.py \
  --input results/hard_filtered.pgen \
  --output results/situational_filtered \
  --hwe \
  --ld-prune \
  --pca
```

**In Notebook:**

```python
from situational_filters import apply_all_situational_filters

result = apply_all_situational_filters(
    input_pfile="results/hard_filtered",
    output_pfile="results/situational_filtered",
    apply_hwe=True,
    apply_ld_prune=True,
    apply_pca=True,
    verbose=True
)
```

### Data Preparation

#### Add Population Labels to Genotypes

```bash
python add_panels.py genotypes/verogen_genotyped_samples.csv \
  --panel 1000genomes/integrated_call_samples_v3.20130502.ALL.panel \
  -o genotypes/verogen_with_panels.csv
```

#### Extract Genotypes from VCF

```bash
python genome_data.py 1000genomes/loci.bed \
  --vcf-template "1000genomes/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
  -o genotypes/extracted_genotypes.csv
```

#### Convert SNP Coordinates to BED Format

```bash
# Amplicon format (with start/end positions)
python csv_to_bed.py snps/ancestry_informative_snp.csv \
  -o snps/ancestry_snps.bed \
  --format amplicon

# rsID format (with single position)
python csv_to_bed.py snps/rsid_list.csv \
  -o snps/rsid_coords.bed \
  --format rsid
```

### Analysis

#### FST-Based Variant Selection

```bash
python fst_selection.py \
  --input 1000genomes/output/EAS_FINAL_DATA_FOR_FST.pgen \
  --populations CHB JPT KHV \
  --output results/fst_top_variants.txt \
  --n-variants 100
```

#### Panel EDA

```bash
python panel_eda.py 1000genomes/integrated_call_samples_v3.20130502.ALL.panel \
  --super-pop EAS \
  --populations CHB JPT KHV \
  -o results/eas_panel_summary.csv
```

#### PCA Visualization

```bash
python pca_plot.py 1000genomes/output/eigenvec \
  --panel 1000genomes/integrated_call_samples_v3.20130502.ALL.panel \
  --pop-col pop \
  -o results/pca_plot.png \
  --title "PCA: EAS Populations"

# Grid of PC pairs
python pca_plot.py 1000genomes/output/eigenvec \
  --grid --title "PCA Grid" \
  -o results/pca_grid.png
```

### Machine Learning

#### Full ML Pipeline

```bash
python ml_training.py genotypes/verogen_with_panels.csv \
  --output-dir results/ml_models \
  --test-size 0.2 \
  --random-state 42
```

**Options:**

- `--models`: Comma-separated models to train (rf, xgboost, lr)
- `--feature-selection`: Method (none, top_n, fst)
- `--n-features`: Number of top features to select
- `--cross-val-folds`: Cross-validation folds (default: 5)

#### In Notebook (Detailed)

```python
from ml_training import load_ml_data, train_random_forest, evaluate_model

# Load data
X, y, features = load_ml_data("genotypes/data.csv", target_column="pop")

# Train model
model = train_random_forest(X, y, n_estimators=100)

# Evaluate
metrics = evaluate_model(model, X_test, y_test)
print(metrics)
```

---

## VCF Preprocessing

For raw VCF files, normalize and convert to PLINK format:

```bash
./vcf_preprocessing.sh 1000genomes allchr
```

**What it does:**

1. Normalizes VCF files (splits multiallelics)
2. Converts to PLINK2 format
3. Merges all chromosomes
4. Sets standard variant IDs

---

## Full Pipeline Workflow

```
Raw VCF Files
    │
    ├─→ vcf_preprocessing.sh (normalize & convert)
    │
    ▼
1. 01_hard_filtering.ipynb
   └─→ SNP-only, Biallelic, MAF, Call Rate filters
    │
    ▼
2. 02_situational_filtering.ipynb
   └─→ HWE, LD Pruning, Variant ID updates
    │
    ├─────────────────────────────┬─────────────────────────┐
    │                             │                         │
    ▼                             ▼                         ▼
3. 03_fst_and_pca.ipynb    Population Structure      Ancestry Selection
   └─→ FST Calculation         Analysis                (FST variants)
    │
    ▼
4. 04_ml_training.ipynb
   └─→ Train RF, XGBoost, LR models
    │
    ▼
5. 05_model_evaluation.ipynb
   └─→ Evaluate, interpret, visualize
```

---

## Troubleshooting

### Import Errors in Notebooks

**Problem:** `ModuleNotFoundError: No module named 'config'`

**Solution:** Make sure cell 1 is executed first (sets working directory and path):

```python
import sys, os
from pathlib import Path
os.chdir(Path.cwd().parent.parent)
sys.path.insert(0, str(Path.cwd() / "scripts"))
```

### Missing Input Files

**Problem:** `FileNotFoundError: No such file or directory`

**Solution:** Check paths in `config.py`:

```python
# Verify file exists
from pathlib import Path
from config import PATHS

print(PATHS.VCF_FILE)
print(Path(PATHS.VCF_FILE).exists())
```

### PLINK2 Not Found

**Problem:** `plink2: command not found`

**Solution:** Install PLINK2:

```bash
conda install -c bioconda plink2
```

### Memory Issues with Large Datasets

**Solution:** Reduce dataset size or increase memory:

```bash
# Use subset of samples
plink2 --pfile input \
       --keep sample_subset.txt \
       --make-pgen --out output
```

---

## Advanced Usage

### Custom Configuration

Edit `config.py` to change defaults for all scripts/notebooks:

```python
# config.py
HARD_FILTERS = HardFilterThresholds(
    MIN_AF=0.0016,           # Minimum allele frequency
    MAX_GENO=0.02,           # Maximum missing genotype rate
    # ...
)
```

### Parallel Execution

Run hard and situational filters in parallel:

```bash
./run_hard_filters.sh &
PID1=$!

./run_situational_filters.sh &
PID2=$!

wait $PID1 $PID2
echo "Both filters complete"
```

### Output Organization

```
results/
├── hard_filtered/          # After hard filters
│   ├── *.pgen, *.pvar, *.psam
│   └── *_info.afreq
├── situational_filtered/   # After situational filters
├── fst_results/            # FST analysis
├── pca/                    # PCA plots & data
├── ml_models/              # Trained models
│   ├── random_forest.pkl
│   ├── xgboost.pkl
│   └── feature_importances.csv
└── reports/                # Classification reports
```

---

## Performance Tips

1. **Use PLINK2** for large datasets (faster than Python-based filtering)
2. **Parallel chromosome processing** in VCF preprocessing
3. **Feature selection** reduces ML training time
4. **Cross-validation** validates model generalization
5. **Sample stratification** prevents data leakage in train/test splits
