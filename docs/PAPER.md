# Identification and Validation of Ancestry-Informative SNPs for East Asian Subpopulation Discrimination: A Statistical and Machine Learning Approach

## Paper Outline (IEEE Biomedical Statistics Format)

---

## Document Information

- **Type**: Research Article
- **Format**: IEEE Transactions on Biomedical Statistics
- **Keywords**: Ancestry-Informative SNPs (AISNPs), Population Genetics, Machine Learning, East Asian Subpopulations, Feature Selection, 1000 Genomes Project

---

## I. ABSTRACT

> **Background**: Ancestry-informative single nucleotide polymorphisms (AISNPs) are critical biomarkers for population stratification in biomedical research, forensic genetics, and personalized medicine. While panels exist for continental-level ancestry inference, distinguishing closely related subpopulations remains challenging due to minimal genetic differentiation.
>
> **Objective**: To develop and validate a minimal AISNP panel capable of distinguishing three East Asian subpopulations (Han Chinese, Japanese, and Vietnamese) using a rigorous statistical selection framework combined with machine learning validation.
>
> **Methods**: We analyzed 306 samples from the 1000 Genomes Project representing CHB (n=103), JPT (n=104), and KHV (n=99) populations. A multi-stage pipeline applied: (1) hard quality control filtering, (2) Hardy-Weinberg equilibrium and linkage disequilibrium pruning, (3) four-test statistical consensus selection (χ², mutual information, information gain, Kullback-Leibler divergence), and (4) machine learning cross-validation with 10 classifiers.
>
> **Results**: A consensus panel of 37 AISNPs achieved 89.55% classification accuracy (5-fold CV) using Support Vector Machine with RBF kernel. This panel demonstrated superior efficiency (2.42% accuracy per SNP) compared to published panels containing 52-125 markers. Pairwise Hudson F_ST values (0.006-0.013) confirmed minimal genetic differentiation, highlighting the challenge addressed.
>
> **Conclusions**: The four-test statistical consensus approach effectively identifies population-discriminating variants despite low F_ST values. The resulting 37-SNP panel offers a cost-effective solution for East Asian ancestry inference applicable to forensic genetics and pharmacogenomics studies.

---

## II. INTRODUCTION

### A. Background and Motivation

*Suggested content (2-3 paragraphs):*

1. **Population stratification in biomedical research**
   - Confounding effects in genome-wide association studies (GWAS)
   - Importance of ancestry adjustment in clinical trials
   - Relevance to pharmacogenomics and adverse drug reactions

2. **AISNPs in forensic and clinical applications**
   - Definition and characteristics of ancestry-informative markers
   - Current panels (Kidd lab, ForenSeq, Seldin et al.)
   - Limitation: Most panels target continental-level discrimination

3. **Challenge of subpopulation differentiation**
   - East Asian populations share recent common ancestry
   - Low F_ST values indicate minimal genetic differentiation
   - Need for specialized panels for within-region classification

### B. Research Objectives

> Primary Objective: Identify a minimal set of AISNPs that reliably distinguish CHB, JPT, and KHV populations from the 1000 Genomes Project.

> Secondary Objectives:
>
> 1. Compare statistical selection methods for AISNP identification
> 2. Validate selected markers using multiple machine learning classifiers
> 3. Benchmark against published AISNP panels

### C. Key Contributions

1. Novel four-test consensus statistical framework for AISNP selection
2. Demonstration that 37 SNPs achieve ~90% accuracy for East Asian subpopulation classification
3. Comparative analysis showing population-specific panels outperform generic continental panels
4. Open-source pipeline for reproducible ancestry marker discovery

---

## III. RELATED WORK

### A. AISNP Panel Development

| Reference | Year | SNPs | Target Populations | Accuracy |
|-----------|------|------|-------------------|----------|
| Kidd et al. [1] | 2014 | 55 | Global (7 regions) | ~95% continental |
| Seldin et al. [2] | 2006 | 128 | European substructure | Variable |
| Phillips et al. [3] | 2007 | 34 | Global ancestry | ~90% continental |
| Pakstis et al. [4] | 2012 | 128 | Worldwide | High continental |

### B. East Asian Population Genetics

*Cite and discuss:*

- 1000 Genomes Project Consortium publications [5]
- East Asian population structure studies [6-8]
- F_ST and genetic distance metrics [9]

### C. Machine Learning in Population Genetics

*Cite and discuss:*

- Random Forest for SNP selection [10]
- SVM applications in ancestry inference [11]
- Deep learning approaches [12]

---

## IV. MATERIALS AND METHODS

### A. Dataset Description

#### 1. 1000 Genomes Project Data

| Characteristic | Value |
|----------------|-------|
| **Data Source** | 1000 Genomes Project Phase 3 |
| **Genome Build** | GRCh37/hg19 |
| **Total Samples** | 2,504 |
| **Super Populations** | 5 (AFR, AMR, EAS, EUR, SAS) |
| **East Asian Samples** | 504 |

#### 2. Target Populations

| Population Code | Description | Sample Size (n) | Geographic Origin |
|-----------------|-------------|-----------------|-------------------|
| **CHB** | Han Chinese in Beijing | 103 | Northern China |
| **JPT** | Japanese in Tokyo | 104 | Japan |
| **KHV** | Kinh in Ho Chi Minh City | 99 | Southern Vietnam |
| **Total** | — | **306** | East Asia |

> **Figure 1** (Recommended): Map of East Asia showing sample collection locations with population sizes.

### B. Quality Control Pipeline

#### Stage 1: Hard Filtering

```
Input: Raw VCF (1000 Genomes)
       ↓
Filter: SNPs only (--snps-only)
       ↓
Filter: Biallelic variants (--max-alleles 2)
       ↓
Filter: MAF > 0.0016 (singleton removal)
       ↓
Filter: Call rate ≥ 95% (--geno 0.05)
       ↓
Output: EAS_AND_SNP_filtered_data
```

**Table 1**: Hard Filtering Parameters

| Filter | PLINK2 Parameter | Rationale | Variants Removed |
|--------|------------------|-----------|------------------|
| SNP-only | `--snps-only` | Focus on substitutions | Indels, CNVs |
| Biallelic | `--max-alleles 2` | Computational simplicity | Multi-allelic |
| MAF | `--maf 0.0016` | Remove rare variants | Singletons |
| Call rate | `--geno 0.05` | Data completeness | Poor quality |

#### Stage 2: Situational Filtering

| Filter | Parameters | Biological Rationale |
|--------|------------|---------------------|
| Hardy-Weinberg Equilibrium | p < 10⁻⁶, keep-fewhet | Remove genotyping errors; preserve selection signals |
| Linkage Disequilibrium | window=1000kb, r²<0.1 | Remove redundant markers; ensure independence |
| Unique IDs | Deduplicate | Data integrity |

> **Figure 2** (Recommended): Flowchart of the complete quality control pipeline with variant counts at each stage.

### C. Statistical SNP Selection Methods

We employed four complementary statistical tests to identify population-discriminating SNPs:

#### 1. Pearson Chi-Squared Test (χ²)

**Hypothesis Testing Framework:**

$$H_0: \text{Genotype frequencies are independent of population}$$
$$H_1: \text{Genotype frequencies differ between populations}$$

**Test Statistic:**

$$\chi^2 = \sum_{i=1}^{r} \sum_{j=1}^{c} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}$$

where $O_{ij}$ is the observed count in cell $(i,j)$ and $E_{ij} = \frac{n_{i\cdot} \times n_{\cdot j}}{n}$ is the expected count under independence.

**Selection Criterion:** FDR-corrected q-value < 0.05 (Benjamini-Hochberg procedure)

#### 2. Mutual Information (MI)

**Definition:**

$$MI(G; P) = \sum_{g \in \mathcal{G}} \sum_{p \in \mathcal{P}} P(g, p) \log_2 \frac{P(g, p)}{P(g)P(p)}$$

where $G$ represents genotype (0, 1, 2) and $P$ represents population (CHB, JPT, KHV).

**Properties:**

- Symmetric: $MI(G; P) = MI(P; G)$
- Non-negative: $MI \geq 0$
- Bounded: $MI \leq \min(H(G), H(P))$

**Selection Criterion:** Top 500 SNPs by MI score

#### 3. Information Gain (IG)

**Definition:**

$$IG(P, G) = H(P) - H(P|G) = H(P) - \sum_{g \in \mathcal{G}} P(G=g) \cdot H(P|G=g)$$

where entropy:

$$H(P) = -\sum_{p \in \mathcal{P}} P(p) \log_2 P(p)$$

**Selection Criterion:** Top 500 SNPs by IG score

#### 4. Kullback-Leibler Divergence (KL)

**Pairwise Definition:**

$$D_{KL}(P_i || P_j) = \sum_{g \in \mathcal{G}} P_i(g) \log_2 \frac{P_i(g)}{P_j(g)}$$

**Mean Pairwise KL Divergence:**

$$\bar{D}_{KL} = \frac{1}{|\mathcal{C}|} \sum_{(i,j) \in \mathcal{C}} D_{KL}(P_i || P_j)$$

where $\mathcal{C}$ = {(CHB,JPT), (CHB,KHV), (JPT,KHV)}

**Selection Criterion:** Top 500 SNPs by mean pairwise KL divergence

#### 5. Consensus Selection Strategy

Final AISNP set defined as intersection:

$$\text{Consensus SNPs} = S_{\chi^2} \cap S_{MI} \cap S_{IG} \cap S_{KL}$$

where $S_{\chi^2}$ = FDR-significant SNPs, and $S_{MI}, S_{IG}, S_{KL}$ = Top 500 SNPs by respective metrics.

> **Figure 3** (Recommended): Venn diagram showing overlap of SNPs selected by each statistical test.

### D. Population Genetic Analysis

#### Hudson F_ST Estimation

Population differentiation quantified using Hudson's estimator:

$$F_{ST} = \frac{H_T - H_S}{H_T}$$

Computed via PLINK2: `--fst hudson method=wc`

> **Table 2**: Pairwise F_ST Matrix

| | CHB | JPT | KHV |
|-----|------|------|------|
| CHB | — | 0.0062 | 0.0056 |
| JPT | 0.0062 | — | 0.0125 |
| KHV | 0.0056 | 0.0125 | — |

### E. Machine Learning Validation

#### 1. Classifiers Evaluated

| Model | Type | Key Hyperparameters |
|-------|------|---------------------|
| Random Forest | Ensemble | n_estimators=100, max_depth=10 |
| XGBoost | Gradient Boosting | n_estimators=100, max_depth=5, lr=0.1 |
| Logistic Regression | Linear | max_iter=1000, multi_class='multinomial' |
| SVM (RBF) | Kernel | kernel='rbf', probability=True |
| SVM (Linear) | Linear | kernel='linear' |
| K-Nearest Neighbors | Instance-based | n_neighbors=5 |
| Naive Bayes | Probabilistic | Gaussian assumption |
| Gradient Boosting | Ensemble | n_estimators=100, max_depth=5 |
| MLP | Neural Network | hidden_layers=(100, 50), early_stopping |
| AdaBoost | Boosting | Default parameters |

#### 2. Evaluation Protocol

- **Cross-Validation**: 5-fold stratified
- **Metrics**: Accuracy, F1-score (macro), Precision, Recall
- **Overfitting Assessment**: Train-Test accuracy gap

#### 3. Feature Importance Analysis

Combined ranking from:

1. **Random Forest**: Mean decrease in Gini impurity
2. **Logistic Regression**: Absolute coefficient magnitudes

$$\text{Composite Rank}_i = \frac{\text{Rank}_{RF}(i) + \text{Rank}_{LR}(i)}{2}$$

### F. Comparative Panel Analysis (Part 2)

Known AISNP panels converted to GRCh37 coordinates via Ensembl REST API and evaluated using identical ML protocol.

| Panel Source | Reference | SNPs | Original Purpose |
|--------------|-----------|------|------------------|
| cal_et_al | [13] | 52 | East Asian ancestry |
| Kidd_55 | [1] | 55 | Global biogeographic |
| Seldin_128 | [2] | 128 | European substructure |
| ForenSeq | Verogen | 55 | Forensic ancestry |
| Hsiao_Lin_Hwa | [14] | 125 | General ancestry |

### G. Software and Reproducibility

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK2 | 2.0 | Genetic data processing |
| Python | 3.10+ | Statistical analysis |
| scikit-learn | 1.x | Machine learning |
| XGBoost | 1.x | Gradient boosting |
| pandas | 2.x | Data manipulation |
| scipy | 1.x | Statistical tests |
| statsmodels | 0.14+ | FDR correction |

---

## V. RESULTS

### A. Quality Control Summary

> **Table 3**: Variant Counts Through QC Pipeline

| Stage | Variants Remaining | Percentage Retained |
|-------|-------------------|---------------------|
| Raw input | X,XXX,XXX | 100% |
| Hard filtering | XXX,XXX | XX% |
| HWE + LD pruning | XX,XXX | XX% |
| Final (consensus) | **37** | — |

### B. Statistical SNP Selection Results

> **Table 4**: Consensus SNP Statistics

| Metric | Value | Range |
|--------|-------|-------|
| Total consensus SNPs | 37 | — |
| χ² p-values | All < 0.05 (FDR) | 2.10×10⁻⁶ – 3.29×10⁻¹⁶ |
| Mutual Information | 0.063 – 0.200 | bits |
| Information Gain | 0.091 – 0.289 | — |
| KL Divergence | 1.33 – 4.39 | — |

> **Figure 4** (Recommended): Histogram distributions of statistical test scores for all SNPs with consensus SNP thresholds marked.

> **Table 5**: Top 10 Consensus SNPs by Composite Rank

| Rank | Variant ID | Chr | χ² p-value | MI | IG | KL |
|------|------------|-----|------------|-----|-----|-----|
| 1 | 3:130239945[b37]G,A | 3 | 2.55×10⁻⁹ | 0.103 | — | 4.22 |
| 2 | 4:101002868[b37]A,T | 4 | 3.29×10⁻¹⁶ | 0.140 | — | 3.99 |
| 3 | 12:128054516[b37]C,T | 12 | 2.32×10⁻¹³ | 0.113 | — | 3.04 |
| 4 | 20:9928437[b37]G,T | 20 | 3.23×10⁻⁹ | 0.067 | — | 1.80 |
| 5 | 11:130759223[b37]G,T | 11 | 2.10×10⁻⁶ | 0.063 | — | 1.81 |
| ... | ... | ... | ... | ... | ... | ... |

> **Figure 5** (Recommended): Manhattan-style plot showing χ² significance across chromosomes with consensus SNPs highlighted.

> **Figure 6** (Recommended): Chromosomal distribution bar chart of 37 consensus SNPs.

### C. Machine Learning Performance

> **Table 6**: Cross-Validation Results (37 Consensus SNPs)

| Rank | Model | Accuracy (%) | Std (%) | F1-Score | Overfit Gap (%) |
|------|-------|--------------|---------|----------|-----------------|
| 1 | **SVM (RBF)** | **89.55** | ±3.80 | 0.895 | 3.3 |
| 2 | Logistic Regression | 86.27 | ±2.47 | 0.861 | 7.4 |
| 3 | SVM (Linear) | 84.96 | ±4.34 | 0.846 | 8.3 |
| 4 | Naive Bayes | 82.68 | ±1.65 | 0.823 | 2.0 |
| 5 | Random Forest | 81.69 | ±4.10 | 0.814 | 10.5 |
| 6 | XGBoost | 80.72 | ±1.91 | 0.802 | 12.5 |
| 7 | MLP | 79.44 | ±6.23 | 0.783 | 7.5 |
| 8 | Gradient Boosting | 78.10 | ±3.98 | 0.779 | 16.6 |
| 9 | AdaBoost | 63.70 | ±5.85 | 0.592 | 5.3 |
| 10 | K-NN | 57.49 | ±8.01 | 0.580 | 9.8 |

> **Figure 7** (Recommended): Bar chart comparing accuracy ± std across all classifiers.

> **Figure 8** (Recommended): Confusion matrix heatmap for best model (SVM-RBF) showing per-class performance.

### D. Feature Reduction Analysis

> **Table 7**: Accuracy vs. Number of SNPs

| SNPs | Best Model | Accuracy (%) | Notes |
|------|------------|--------------|-------|
| 5 | AdaBoost | 67.64 | Insufficient |
| 10 | Logistic Regression | 74.85 | Moderate |
| 15 | Random Forest | 84.31 | Strong improvement |
| 20 | SVM (Linear) | 84.32 | Plateau begins |
| 25 | SVM (RBF) | 86.93 | Near-optimal |
| 37 | SVM (RBF) | 89.55 | Best overall |

> **Figure 9** (Recommended): Line plot showing accuracy curves for each model as function of SNP count (5, 10, 15, 20, 25, 37).

### E. Known AISNP Panel Comparison

> **Table 8**: Cross-Source Performance Comparison

| Source | SNPs | Best Accuracy (%) | Best Model | Efficiency (%/SNP) |
|--------|------|-------------------|------------|-------------------|
| **cal_et_al** | 52 | **91.50** | XGBoost | 1.76 |
| **statistical_all4** | 37 | **89.55** | SVM (RBF) | **2.42** |
| hsiao_lin_hwa | 125 | 64.38 | SVM (RBF) | 0.51 |
| seldin_128 | 124 | 64.05 | SVM (RBF) | 0.52 |
| forenseq | 55 | 60.80 | SVM (RBF) | 1.11 |
| kidd_55 | 53 | 60.78 | Logistic Reg. | 1.15 |

> **Figure 10** (Recommended): Performance heatmap (sources × models) showing accuracy and F1 scores.

> **Figure 11** (Recommended): Scatter plot of SNP count vs. accuracy with panel labels.

### F. PCA Visualization

> **Figure 12** (Recommended): 2D PCA scatter plot (PC1 vs PC2) colored by population using consensus SNPs.

> **Figure 13** (Recommended): 3D PCA scatter plot showing population clustering.

> **Figure 14** (Recommended): Scree plot showing variance explained by principal components.

---

## VI. DISCUSSION

### A. Principal Findings

1. **Effective subpopulation discrimination despite low F_ST**
   - F_ST values (0.006-0.013) indicate minimal differentiation
   - Statistical consensus approach captures subtle frequency differences
   - 89.55% accuracy demonstrates feasibility of within-region classification

2. **Statistical method comparison**
   - Chi-squared provides hypothesis testing framework
   - Information-theoretic measures (MI, IG, KL) capture non-linear associations
   - Consensus intersection ensures robust marker selection

3. **Panel design superiority over panel size**
   - 37-SNP statistical panel comparable to 52-SNP cal_et_al (89.5% vs 91.5%)
   - Generic 124-125 SNP panels achieve only 60-64%
   - Population-specific design critical for subpopulation inference

### B. Comparison with Literature

*Discuss in relation to:*

- Previous East Asian population studies
- F_ST thresholds for population differentiation
- Existing AISNP panel performances

### C. Model Selection Considerations

| Application | Recommended Model | Rationale |
|-------------|-------------------|-----------|
| Production deployment | SVM (RBF) | Best accuracy, moderate complexity |
| Interpretability required | Logistic Regression | Transparent coefficients, 86.27% |
| Quick screening | Naive Bayes | Fast, low overfit (2%), 82.68% |
| Feature importance | Random Forest | Gini importance ranking |

### D. Biological Implications

*Discuss potential biological significance of top-ranked SNPs:*

- Chromosomal distribution patterns
- Proximity to known population-stratified genes
- Potential functional relevance

### E. Limitations

1. **Sample size**: n=306 may limit generalizability
2. **Population scope**: Only three East Asian subpopulations included
3. **Genome build**: GRCh37 coordinates may require updating
4. **External validation**: Performance on independent cohorts unknown
5. **Admixed individuals**: Panel not validated for mixed ancestry

### F. Clinical and Forensic Implications

*Discuss applications in:*

- Forensic DNA phenotyping
- Pharmacogenomics stratification
- Clinical trial design
- Transplant matching

---

## VII. CONCLUSIONS

### A. Summary

This study demonstrates that a minimal panel of 37 statistically-selected AISNPs achieves 89.55% accuracy in distinguishing Han Chinese, Japanese, and Vietnamese populations. The four-test consensus framework (χ², MI, IG, KL divergence) effectively identifies population-discriminating markers despite low genetic differentiation (F_ST < 0.015).

### B. Key Contributions

1. Novel multi-test statistical framework for AISNP selection
2. Validated 37-SNP panel for East Asian subpopulation inference
3. Evidence that targeted panel design outperforms larger generic panels
4. Open-source, reproducible bioinformatics pipeline

### C. Future Directions

1. **Expand populations**: Include CDX (Dai Chinese), CHS (Southern Han Chinese)
2. **External validation**: Test on independent Asian cohorts
3. **Functional annotation**: Investigate biological relevance of top SNPs
4. **Panel optimization**: Reduce to optimal 15-25 SNP subset
5. **Clinical validation**: Prospective testing in forensic/clinical settings

---

## VIII. DATA AVAILABILITY

- **Source Data**: 1000 Genomes Project Phase 3 (publicly available)
- **Code Repository**: [GitHub Link]
- **Processed Data**: Available upon request

---

## REFERENCES

[1] K. K. Kidd et al., "Progress toward an efficient panel of SNPs for ancestry inference," *Forensic Sci. Int. Genet.*, vol. 10, pp. 23–32, 2014.

[2] M. F. Seldin et al., "Application of ancestry informative markers to association studies in European Americans," *PLoS Genet.*, vol. 4, no. 1, e1000004, 2008.

[3] C. Phillips et al., "Inferring ancestral origin using a single multiplex assay of ancestry-informative marker SNPs," *Forensic Sci. Int. Genet.*, vol. 1, no. 3–4, pp. 273–280, 2007.

[4] A. J. Pakstis et al., "SNPs for a universal individual identification panel," *Hum. Genet.*, vol. 127, no. 3, pp. 315–324, 2010.

[5] The 1000 Genomes Project Consortium, "A global reference for human genetic variation," *Nature*, vol. 526, pp. 68–74, 2015.

[6] HUGO Pan-Asian SNP Consortium, "Mapping human genetic diversity in Asia," *Science*, vol. 326, no. 5959, pp. 1541–1545, 2009.

[7] J. Z. Li et al., "Worldwide human relationships inferred from genome-wide patterns of variation," *Science*, vol. 319, no. 5866, pp. 1100–1104, 2008.

[8] S. Xu et al., "Genomic dissection of population substructure of Han Chinese and its implication in association studies," *Am. J. Hum. Genet.*, vol. 85, no. 6, pp. 762–774, 2009.

[9] B. S. Weir and C. C. Cockerham, "Estimating F-statistics for the analysis of population structure," *Evolution*, vol. 38, no. 6, pp. 1358–1370, 1984.

[10] Y. S. Chen et al., "Random forests for genomic data analysis," *Genomics*, vol. 99, no. 6, pp. 323–329, 2012.

[11] A. L. Price et al., "Principal components analysis corrects for stratification in genome-wide association studies," *Nat. Genet.*, vol. 38, no. 8, pp. 904–909, 2006.

[12] A. Auton et al., "Genetic recombination is targeted towards gene promoter regions in dogs," *PLoS Genet.*, vol. 9, no. 12, e1003984, 2013.

[13] Cal et al., "East Asian ancestry informative SNP panel," *[Journal]*, 2025. [To be updated with actual citation]

[14] Hsiao, Lin, and Hwa, "Ancestry informative markers for East Asian populations," *[Journal]*, 2017. [To be updated with actual citation]

---

## SUPPLEMENTARY MATERIALS

### Supplementary Table S1

Complete list of 37 consensus SNPs with full statistical metrics

### Supplementary Table S2

Cross-validation results for all models across all SNP counts

### Supplementary Table S3

Known AISNP panel details and source information

### Supplementary Figure S1

Complete confusion matrices for all classifiers

### Supplementary Figure S2

ROC curves for multi-class classification

### Supplementary Code

Pipeline scripts available at: [Repository URL]

---

## RECOMMENDED FIGURES CHECKLIST

| Figure | Type | Description | Recommended Tool | Priority |
|--------|------|-------------|------------------|----------|
| Fig 1 | Map | Sample collection locations | matplotlib/folium | Medium |
| Fig 2 | Flowchart | QC pipeline with variant counts | draw.io/mermaid | High |
| Fig 3 | Venn | Statistical test overlap | matplotlib-venn | High |
| Fig 4 | Histogram | Statistical score distributions | seaborn | Medium |
| Fig 5 | Manhattan | χ² significance by chromosome | matplotlib | Medium |
| Fig 6 | Bar | Chromosomal distribution of SNPs | matplotlib | Low |
| Fig 7 | Bar | Classifier accuracy comparison | matplotlib | **High** |
| Fig 8 | Heatmap | Best model confusion matrix | seaborn | **High** |
| Fig 9 | Line | Accuracy vs SNP count curves | matplotlib | **High** |
| Fig 10 | Heatmap | Source × Model performance | seaborn | **High** |
| Fig 11 | Scatter | SNP count vs accuracy | matplotlib | Medium |
| Fig 12 | Scatter | PCA 2D by population | matplotlib | **High** |
| Fig 13 | 3D Scatter | PCA 3D by population | plotly | Medium |
| Fig 14 | Line | PCA scree plot | matplotlib | Low |

---

## AUTHOR CONTRIBUTIONS

- **Conceptualization**: [Authors]
- **Methodology**: [Authors]
- **Software**: [Authors]
- **Validation**: [Authors]
- **Formal Analysis**: [Authors]
- **Data Curation**: [Authors]
- **Writing – Original Draft**: [Authors]
- **Writing – Review & Editing**: [Authors]
- **Visualization**: [Authors]

---

## ACKNOWLEDGMENTS

The authors thank the 1000 Genomes Project Consortium for making data publicly available.

---

## CONFLICT OF INTEREST

The authors declare no conflict of interest.

---

*Last Updated: January 2026*
