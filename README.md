# Ancestry-Informative SNP Selection Pipeline

A comprehensive bioinformatics pipeline for identifying and validating Ancestry-Informative Single Nucleotide Polymorphisms (AISNPs) from 1000 Genomes Project data.

## 🎯 Project Overview

This project implements a two-part analysis pipeline:

- **Part 1**: Statistical SNP selection from 1000 Genomes data using quality control, statistical tests, and machine learning
- **Part 2**: Comparison with known AISNP panels from research papers and commercial products

### Target Populations

- **CHB**: Han Chinese in Beijing
- **JPT**: Japanese in Tokyo  
- **KHV**: Kinh in Ho Chi Minh City, Vietnam

## 📁 Project Structure

```
BioinfoMidterm/
├── 1000genomes/              # Input data and intermediate outputs
│   ├── main_vcf/             # Source VCF files
│   ├── output/               # PLINK output files
│   └── *.csv                 # Sample lists
│
├── data/
│   └── known_aisnps/         # rsID lists from papers/products (Part 2)
│
├── scripts/
│   ├── config.py             # Centralized configuration
│   ├── utils.py              # Utility functions
│   ├── hard_filters.py       # Hard filtering functions
│   ├── situational_filters.py # Situational filtering
│   ├── fst_selection.py      # FST-based SNP selection
│   ├── ml_training.py        # ML training utilities
│   │
│   ├── notebooks/            # Jupyter notebooks
│   │   ├── 01_hard_filtering.ipynb
│   │   ├── 02_situational_filtering.ipynb
│   │   ├── 02b_statistical_snp_selection.ipynb
│   │   ├── 03_fst_and_pca.ipynb
│   │   ├── 03b_statistical_snp_analysis.ipynb
│   │   ├── 04_ml_training.ipynb
│   │   ├── 04b_ml_training_statistical.ipynb
│   │   ├── 04c_ml_consensus_snps.ipynb
│   │   ├── 05_model_evaluation.ipynb
│   │   └── part2/            # Part 2 notebooks
│   │       ├── 06_rsid_to_bed.ipynb
│   │       ├── 07_bed_to_ml_matrix.ipynb
│   │       └── 08_known_aisnps_ml.ipynb
│   │
│   └── part2/                # Part 2 helper modules
│       ├── rsid_utils.py
│       ├── bed_to_matrix.py
│       └── ml_comparison.py
│
├── output/
│   ├── ml_models/            # Trained models
│   └── part2/                # Part 2 outputs
│
├── graphs/                   # Visualization outputs
├── reports/                  # Text reports
└── docs/                     # Documentation
```

## 🚀 Quick Start

### Prerequisites

```bash
# System requirements
- Python 3.10+
- PLINK2 (https://www.cog-genomics.org/plink/2.0/)
- bcftools (optional, for VCF processing)

# Python packages
pip install pandas numpy scipy scikit-learn xgboost matplotlib seaborn tqdm statsmodels requests
```

### Conda Environment (spec-file)

- To recreate the exact environment used in this project, use the provided [spec-file.txt](spec-file.txt):

```bash
conda create -n aisnp --file spec-file.txt
conda activate aisnp
```

### Data Download

- Recommended: Use Globus to download 1000 Genomes Project data for fast, reliable transfers.
- Merged all vcf with bcftools

### Running the Pipeline

#### Part 1: Statistical SNP Selection

```bash
# Navigate to notebooks
cd scripts/notebooks

# Run in order:
# 1. Hard filtering (quality control)
jupyter notebook 01_hard_filtering.ipynb

# 2. Situational filtering (HWE, LD pruning)
jupyter notebook 02_situational_filtering.ipynb

# 3. Statistical SNP selection (χ², MI, IG, KL divergence)
jupyter notebook 02b_statistical_snp_selection.ipynb

# 4. PCA visualization
jupyter notebook 03_fst_and_pca.ipynb

# 5. ML training with consensus SNPs
jupyter notebook 04c_ml_consensus_snps.ipynb
```

#### Part 2: Known AISNP Comparison

```bash
# Add rsID files to data/known_aisnps/
# Then run:
jupyter notebook part2/06_rsid_to_bed.ipynb
jupyter notebook part2/07_bed_to_ml_matrix.ipynb
jupyter notebook part2/08_known_aisnps_ml.ipynb
```

## 📊 Pipeline Workflows

### Part 1: Statistical Selection Workflow

```
1000 Genomes VCF
       ↓
[01] Hard Filters (SNP-only, biallelic, MAF, call rate)
       ↓
[02] Situational Filters (HWE, LD pruning)
       ↓
[02b] Statistical Tests (χ², MI, IG, KL) → All-4-tests consensus SNPs
       ↓
[03] PCA Visualization
       ↓
[04c] ML Training (RF, XGBoost, SVM, LR, etc.)
       ↓
Performance metrics, confusion matrices, feature importance
```

### Part 2: Known AISNP Comparison Workflow

```
rsID lists (CSV/TXT)
       ↓
[06] rsID → BED conversion (Ensembl API lookup)
       ↓
[07] Extract genotypes from pfile → ML matrices
       ↓
[08] Compare performance across all sources
       ↓
Heatmaps, confusion matrices, comparison plots
```

## ⚙️ Configuration

All parameters are centralized in `scripts/config.py`:

### Hard Filters (Always Applied)

| Parameter | Default | Description |
|-----------|---------|-------------|
| SNP_ONLY | True | Remove indels, CNVs |
| MAX_ALLELES | 2 | Biallelic variants only |
| MIN_AF | 0.0016 | Minimum allele frequency |
| MIN_CALL_RATE | 0.95 | Minimum genotyping rate |

### Situational Filters (Context-Dependent)

| Parameter | Default | Description |
|-----------|---------|-------------|
| HWE_P_THRESHOLD | 1e-6 | HWE p-value threshold |
| LD_WINDOW_KB | 1000 | LD pruning window |
| LD_R2_THRESHOLD | 0.1 | LD R² cutoff |
| FST_TOP_N | 1000 | Top FST variants |

### Statistical Tests (02b)

| Test | Purpose |
|------|---------|
| Pearson χ² | Test genotype-population independence |
| Mutual Information | Measure genotype-population association |
| Information Gain | Population entropy reduction |
| KL Divergence | Genotype distribution divergence |

**Consensus SNPs**: Only SNPs passing ALL 4 tests are selected.

## 📈 Key Outputs

### Part 1 Outputs

| File | Description |
|------|-------------|
| `statistical_all4_snps_02b.csv` | SNPs passing all 4 statistical tests |
| `statistical_ml_data_02b.csv` | ML-ready matrix (consensus SNPs only) |
| `consensus_snps_cv_results.csv` | Cross-validation results |
| `consensus_snps_importance.csv` | Feature importance rankings |

### Part 2 Outputs

| File | Description |
|------|-------------|
| `{source}.bed` | BED file with genomic coordinates |
| `{source}_ml_matrix.csv` | ML-ready genotype matrix |
| `ml_comparison_results.csv` | K-fold CV results across sources |
| `performance_heatmaps.png` | Visual comparison |

## 🔬 Methods

### Statistical SNP Selection

1. **χ² Test**: Tests if genotype frequencies differ across populations
2. **Mutual Information**: Measures information shared between genotype and population
3. **Information Gain**: Quantifies population entropy reduction from knowing genotype
4. **KL Divergence**: Measures divergence between population-specific genotype distributions

### Machine Learning Models

| Model | Use Case |
|-------|----------|
| Random Forest | Robust, handles non-linear relationships |
| XGBoost | High accuracy, handles imbalanced data |
| Logistic Regression | Interpretable, good baseline |
| SVM (RBF/Linear) | Effective for high-dimensional data |
| K-Nearest Neighbors | Simple, distance-based |
| Gradient Boosting | Ensemble method, high accuracy |
| MLP Neural Network | Captures complex patterns |

## 📚 Documentation

### Core Documentation

- [docs/PIPELINE.md](docs/PIPELINE.md) - Detailed pipeline description
- [docs/CONFIGURATION.md](docs/CONFIGURATION.md) - Configuration reference
- [docs/STATISTICAL_TESTS.md](docs/STATISTICAL_TESTS.md) - Statistical methods
- [docs/ML_MODELS.md](docs/ML_MODELS.md) - ML model details
- [scripts/part2/README.md](scripts/part2/README.md) - Part 2 documentation

### Results

- [docs/RESULTS.md](docs/RESULTS.md) - Part 1 Results: Statistical AISNP Selection
- [docs/RESULTS_PART2.md](docs/RESULTS_PART2.md) - Part 2 Results: Known AISNP Panel Comparison

### Presentation & Diagrams

- [docs/slides/PRESENTATION_OUTLINE.md](docs/slides/PRESENTATION_OUTLINE.md) - Slide-by-slide presentation guide (40 slides)
- [docs/diagrams/SYSTEM_DIAGRAMS.md](docs/diagrams/SYSTEM_DIAGRAMS.md) - System block diagrams & data flow

## 🔗 References

- 1000 Genomes Project: <https://www.internationalgenome.org/>
- PLINK2: <https://www.cog-genomics.org/plink/2.0/>
- Ensembl REST API: <https://rest.ensembl.org/>

## 📝 License

This project is for educational/research purposes.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

*Last updated: January 2026*
