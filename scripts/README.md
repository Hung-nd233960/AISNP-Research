# 1000 Genomes EAS Population Classification Pipeline

A modular bioinformatics pipeline for population classification using SNP data from the 1000 Genomes Project.

## Overview

This pipeline analyzes East Asian (EAS) subpopulations (CHB, JPT, KHV) from the 1000 Genomes Project and builds machine learning models for ancestry inference.

## Project Structure

```
BioinfoMidterm/
├── scripts/                    # Python modules and bash scripts
│   ├── config.py              # Centralized configuration
│   ├── utils.py               # Shared utilities
│   ├── hard_filters.py        # Hard filter functions
│   ├── situational_filters.py # Situational filter functions
│   ├── ml_training.py         # ML training pipeline
│   ├── fst_selection.py       # FST-based variant selection
│   ├── run_hard_filters.sh    # Bash script for hard filtering
│   └── run_situational_filters.sh  # Bash script for situational filtering
├── notebooks/                  # Jupyter notebooks for analysis
│   ├── 01_hard_filtering.ipynb
│   ├── 02_situational_filtering.ipynb
│   ├── 03_fst_and_pca.ipynb
│   ├── 04_ml_training.ipynb
│   └── 05_model_evaluation.ipynb
├── 1000genomes/               # Input data directory
│   ├── main_vcf/              # Source VCF files
│   ├── output/                # Filtered PLINK files
│   └── ...
├── output/                    # Analysis outputs
└── reports/                   # Classification reports
```

## Filter Hierarchy

### Hard Filters (Always Applied)

These filters are **strict quality control** measures applied to all datasets. They remove low-quality or uninformative variants.

| Filter | Parameter | Rationale |
|--------|-----------|-----------|
| **SNP-only** | `--snps-only` | Remove indels; SNPs are more reliable for population genetics |
| **Biallelic** | `--max-alleles 2` | Remove multiallelic sites; simplifies analysis |
| **MAF** | `--maf 0.0016` | Remove rare variants (1/612 alleles ≈ singletons) |
| **Call Rate** | `--geno 0.02` | Remove variants missing in >2% of samples |

### Situational Filters (Context-Dependent)

These filters depend on the **study design** and analysis goals. Apply selectively.

| Filter | Parameter | When to Use |
|--------|-----------|-------------|
| **HWE** | `--hwe 1e-6` | Remove variants deviating from Hardy-Weinberg equilibrium (may indicate genotyping errors) |
| **LD Pruning** | `--indep-pairwise 1000 100 0.1` | For PCA, GWAS, population structure analysis (removes correlated variants) |
| **FST Selection** | Custom | For ancestry inference (select high FST variants between populations) |
| **Variant ID Update** | `--set-all-var-ids` | Standardize variant IDs for merging datasets |

## When to Use Each Filter Type

### Hard Filters - Always Apply

- ✅ Any population genetics analysis
- ✅ Machine learning on genetic data
- ✅ GWAS and association studies
- ✅ Before any downstream analysis

### Situational Filters - Apply Selectively

| Analysis Type | HWE | LD Pruning | FST Selection |
|---------------|-----|------------|---------------|
| **PCA** | ✅ | ✅ | ❌ |
| **Population Structure** | ✅ | ✅ | ❌ |
| **Ancestry Inference** | ✅ | ❌ | ✅ |
| **GWAS** | ✅ | ❌ | ❌ |
| **ML Classification** | ✅ | Optional | ✅ |

## Installation

### Prerequisites

```bash
# PLINK2
conda install -c bioconda plink2

# bcftools
conda install -c bioconda bcftools

# Python packages
pip install pandas numpy scikit-learn xgboost matplotlib seaborn
```

### Configuration

Edit `scripts/config.py` to adjust:

- Input/output paths
- Filter thresholds
- ML parameters

## Usage

### Option 1: Command Line (Bash Scripts)

```bash
# Run hard filters
cd scripts
./run_hard_filters.sh

# Run situational filters
./run_situational_filters.sh
```

### Option 2: Python CLI

```bash
# Hard filters
python scripts/hard_filters.py /path/to/input.pgen /path/to/output_prefix

# Situational filters
python scripts/situational_filters.py /path/to/input.pgen /path/to/output_prefix

# ML training
python scripts/ml_training.py /path/to/data.csv --output-dir /path/to/output
```

### Option 3: Jupyter Notebooks

Run notebooks in order:

1. `01_hard_filtering.ipynb` - Apply hard filters
2. `02_situational_filtering.ipynb` - Apply situational filters
3. `03_fst_and_pca.ipynb` - FST and PCA analysis
4. `04_ml_training.ipynb` - Train ML models
5. `05_model_evaluation.ipynb` - Evaluate models

## Pipeline Workflow

```
VCF Files
    │
    ▼
┌─────────────────────┐
│   HARD FILTERS      │ ◄── Always applied
├─────────────────────┤
│ • SNP-only          │
│ • Biallelic         │
│ • MAF ≥ 0.0016      │
│ • Call rate ≥ 98%   │
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│ SITUATIONAL FILTERS │ ◄── Applied based on analysis
├─────────────────────┤
│ • HWE (p > 1e-6)    │
│ • LD pruning        │
│ • Variant ID update │
└─────────────────────┘
    │
    ├────────────────────┐
    │                    │
    ▼                    ▼
┌──────────────┐   ┌──────────────┐
│     PCA      │   │     FST      │
│  Analysis    │   │  Selection   │
└──────────────┘   └──────────────┘
                         │
                         ▼
                  ┌──────────────┐
                  │      ML      │
                  │   Training   │
                  └──────────────┘
```

## Key Parameters

### MAF Threshold Calculation

For a sample of N individuals (2N alleles):

- Singleton: 1/(2N)
- For 306 samples: 1/612 ≈ 0.0016

### LD Pruning Parameters

```
--indep-pairwise 1000 100 0.1
                  │    │   │
                  │    │   └── r² threshold
                  │    └────── step size (variants)
                  └─────────── window size (kb)
```

### HWE Threshold

- `1e-6`: Stringent, removes severe deviations
- `0.05`: Liberal, standard for association studies

## ML Models

| Model | Description | Best For |
|-------|-------------|----------|
| Random Forest | Ensemble of decision trees | General classification |
| XGBoost | Gradient boosted trees | High accuracy |
| Logistic Regression | Linear classifier | Interpretability |

## Output Files

### Filtering Outputs

- `*_MAF_filtered.{pgen,pvar,psam}` - Hard-filtered data
- `*_HWE_filtered.{pgen,pvar,psam}` - HWE-filtered data
- `*.prune.in`, `*.prune.out` - LD pruning lists

### Analysis Outputs

- `*.fst` - FST values between populations
- `*.eigenvec`, `*.eigenval` - PCA results
- `*.csv` - Genotype matrices for ML

### Model Outputs

- `*.pkl` - Trained model files
- `*_report.txt` - Classification reports
- `feature_importances.csv` - Feature rankings

## References

- 1000 Genomes Project: <https://www.internationalgenome.org/>
- PLINK2: <https://www.cog-genomics.org/plink/2.0/>
- bcftools: <http://samtools.github.io/bcftools/>

## License

MIT License
