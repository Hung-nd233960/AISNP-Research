# Pipeline Documentation

Detailed description of the AISNP selection pipeline stages.

## Table of Contents

1. [Overview](#overview)
2. [Part 1: Statistical SNP Selection](#part-1-statistical-snp-selection)
3. [Part 2: Known AISNP Comparison](#part-2-known-aisnp-comparison)
4. [Data Flow](#data-flow)

---

## Overview

This pipeline identifies Ancestry-Informative SNPs (AISNPs) that can distinguish between East Asian subpopulations (CHB, JPT, KHV) using:

1. **Quality control filtering** to ensure data reliability
2. **Statistical tests** to identify population-discriminating SNPs
3. **Machine learning** to validate classification performance
4. **Comparison** with established AISNP panels

---

## Part 1: Statistical SNP Selection

### Stage 1: Hard Filtering (Notebook 01)

**Purpose**: Apply strict quality control to remove unreliable variants.

**Filters Applied**:

| Filter | Command | Rationale |
|--------|---------|-----------|
| SNP-only | `--snps-only` | Focus on single nucleotide variants |
| Biallelic | `--max-alleles 2` | Simplify analysis, improve reliability |
| MAF | `--maf 0.0016` | Remove singletons/rare variants |
| Call rate | `--geno 0.05` | Ensure ≥95% genotyping completeness |

**Input**: Raw VCF from 1000 Genomes  
**Output**: `EAS_AND_SNP_filtered_data` (PLINK2 pfile)

### Stage 2: Situational Filtering (Notebook 02)

**Purpose**: Apply population-genetics-specific filters.

**Filters Applied**:

| Filter | Parameters | Rationale |
|--------|------------|-----------|
| HWE | p < 1e-6, keep-fewhet | Remove genotyping errors while preserving selection signals |
| Unique IDs | Deduplicate | Ensure variant uniqueness |
| LD pruning | window=1000kb, r²=0.1 | Remove redundant variants for ancestry analysis |

**Input**: Hard-filtered pfile  
**Output**: `EAS_FINAL_DATA_FOR_FST` (LD-pruned pfile)

### Stage 3: Statistical SNP Selection (Notebook 02b)

**Purpose**: Identify SNPs that statistically differentiate populations.

**Statistical Tests**:

#### 1. Pearson Chi-Squared Test (χ²)

```
H₀: Genotype frequencies are equal across populations
H₁: Genotype frequencies differ between populations

- Creates 3×3 contingency table (3 genotypes × 3 populations)
- Higher χ² indicates stronger population-genotype association
- FDR correction applied for multiple testing
```

#### 2. Mutual Information (MI)

```
MI(Genotype; Population) = H(Population) - H(Population | Genotype)

- Measures information shared between genotype and population
- Higher MI = more informative SNP
- Top 500 SNPs selected
```

#### 3. Information Gain (IG)

```
IG = H(Population) - Σ P(g) × H(Population | Genotype=g)

- Quantifies entropy reduction in population labels given genotype
- Higher IG = better population discrimination
- Top 500 SNPs selected
```

#### 4. Kullback-Leibler Divergence (KL)

```
KL(P||Q) = Σ P(x) log(P(x)/Q(x))

- Measures divergence between population-specific genotype distributions
- Average pairwise KL across all population pairs
- Top 500 SNPs selected
```

**Consensus Selection**:

```
Final SNPs = χ²_significant ∩ MI_top500 ∩ IG_top500 ∩ KL_top500
```

Only SNPs passing **ALL 4 tests** are selected.

**Caching**: Results are cached to `statistical_snp_scores_02b.csv` to skip re-analysis on subsequent runs.

**Output Files**:

- `statistical_all4_snps_02b.csv` - Full results for consensus SNPs
- `statistical_ml_data_02b.csv` - ML-ready matrix with only consensus SNPs

### Stage 4: PCA Visualization (Notebook 03)

**Purpose**: Visualize population structure using selected SNPs.

**Analysis**:

- Principal Component Analysis on consensus SNPs
- 2D and 3D scatter plots colored by population
- Scree plot showing variance explained
- Evaluate population separability

### Stage 5: ML Training (Notebook 04c)

**Purpose**: Train and evaluate classifiers on consensus SNPs.

**Models Tested**:

1. Random Forest
2. XGBoost
3. Logistic Regression
4. SVM (RBF and Linear kernels)
5. K-Nearest Neighbors
6. Naive Bayes
7. Gradient Boosting
8. MLP Neural Network

**Evaluation**:

- 5-Fold Stratified Cross-Validation
- Metrics: Accuracy, F1, Precision, Recall
- Train vs Test accuracy (overfit detection)
- Feature importance analysis (RF + LR)
- Reduced feature set testing (5, 10, 15, 20, 25 SNPs)

---

## Part 2: Known AISNP Comparison

### Stage 6: rsID to BED Conversion (Notebook 06)

**Purpose**: Convert rsID lists from papers/products to genomic coordinates.

**Process**:

1. Load rsID list from CSV/TXT file
2. Query Ensembl REST API for coordinates (GRCh37)
3. Generate BED file (chr, start, end, rsID)
4. Cache results to avoid repeated API calls

**Supported Formats**:

- CSV with rsID column
- TXT with one rsID per line

### Stage 7: BED to ML Matrix (Notebook 07)

**Purpose**: Extract genotypes for known AISNP panels.

**Process**:

1. Load SNP IDs from BED file
2. Check overlap with our dataset
3. Extract matching SNPs using PLINK2
4. Convert to ML-ready matrix (sample, pop, genotypes)

**Output**: `{source}_ml_matrix.csv` for each AISNP source

### Stage 8: ML Comparison (Notebook 08)

**Purpose**: Compare classification performance across SNP sources.

**Analysis**:

- K-fold CV on all sources with all classifiers
- Confusion matrices for best models
- Performance heatmaps (Accuracy, F1 by source × model)
- SNP count vs accuracy scatter plot
- Summary report

---

## Data Flow

```
                    1000 Genomes VCF
                          │
                          ▼
              ┌───────────────────────┐
              │   01: Hard Filters    │
              │  (SNP, biallelic,     │
              │   MAF, call rate)     │
              └───────────────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │ 02: Situational       │
              │    Filters (HWE,      │
              │    LD pruning)        │
              └───────────────────────┘
                          │
          ┌───────────────┴───────────────┐
          │                               │
          ▼                               ▼
┌─────────────────────┐       ┌─────────────────────┐
│ 02b: Statistical    │       │ 03: FST Selection   │
│ Tests (χ², MI,      │       │ (pairwise FST,      │
│ IG, KL)             │       │ top variants)       │
└─────────────────────┘       └─────────────────────┘
          │                               │
          ▼                               ▼
┌─────────────────────┐       ┌─────────────────────┐
│ All-4-Tests         │       │ FST-based           │
│ Consensus SNPs      │       │ Top SNPs            │
└─────────────────────┘       └─────────────────────┘
          │                               │
          └───────────────┬───────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │  04c: ML Training     │
              │  (Compare methods)    │
              └───────────────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │  05: Evaluation       │
              │  (Final assessment)   │
              └───────────────────────┘


          PART 2: Known AISNP Comparison
          ==============================

    rsID Lists                Statistical SNPs
    (papers,                  (from Part 1)
    products)                       │
        │                           │
        ▼                           │
┌─────────────────┐                 │
│ 06: rsID → BED  │                 │
│ (Ensembl API)   │                 │
└─────────────────┘                 │
        │                           │
        ▼                           │
┌─────────────────┐                 │
│ 07: BED → ML    │                 │
│ Matrix          │                 │
└─────────────────┘                 │
        │                           │
        └─────────────┬─────────────┘
                      │
                      ▼
          ┌───────────────────────┐
          │ 08: ML Comparison     │
          │ (All sources)         │
          └───────────────────────┘
                      │
                      ▼
          ┌───────────────────────┐
          │ Performance Report    │
          │ Confusion Matrices    │
          │ Comparison Plots      │
          └───────────────────────┘
```

---

## File Dependencies

```
Input Files:
├── 1000genomes/main_vcf/ALL_merged.vcf.gz
├── 1000genomes/main_vcf/integrated_call_samples_v3.20130502.ALL.panel
└── 1000genomes/EAS_subpopulation_samples.csv

Intermediate Files (Part 1):
├── EAS_AND_SNP_filtered_data.{pgen,pvar,psam}
├── EAS_AND_SNP_filtered_data_MAF_filtered.{pgen,pvar,psam}
├── EAS_FINAL_DATA_FOR_FST.{pgen,pvar,psam}
└── statistical_snp_scores_02b.csv (cached results)

Output Files (Part 1):
├── statistical_all4_snps_02b.csv
├── statistical_ml_data_02b.csv
├── consensus_snps_cv_results.csv
└── consensus_snps_importance.csv

Output Files (Part 2):
├── {source}.bed
├── {source}_ml_matrix.csv
├── ml_comparison_results.csv
└── graphs/part2/*.png
```

---

*See also: [CONFIGURATION.md](CONFIGURATION.md), [STATISTICAL_TESTS.md](STATISTICAL_TESTS.md), [ML_MODELS.md](ML_MODELS.md)*
