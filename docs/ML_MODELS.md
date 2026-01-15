# Machine Learning Models Reference

Detailed documentation of ML classifiers used for ancestry inference.

## Table of Contents

1. [Overview](#overview)
2. [Model Descriptions](#model-descriptions)
3. [Hyperparameters](#hyperparameters)
4. [Model Selection Guidelines](#model-selection-guidelines)
5. [Evaluation Metrics](#evaluation-metrics)
6. [Feature Importance](#feature-importance)

---

## Overview

The pipeline uses multiple classifiers to:

1. Validate that selected SNPs can distinguish populations
2. Compare performance across different algorithms
3. Identify the most robust model for deployment

### Models Used

| Model | Type | Strengths | Weaknesses |
|-------|------|-----------|------------|
| Random Forest | Ensemble | Robust, handles non-linearity | Memory-intensive |
| XGBoost | Gradient Boosting | High accuracy, handles imbalance | Prone to overfitting |
| Logistic Regression | Linear | Interpretable, fast | Assumes linearity |
| SVM (RBF) | Kernel | Effective in high dimensions | Slow with large datasets |
| SVM (Linear) | Linear | Fast, interpretable | Assumes linear separability |
| K-Nearest Neighbors | Instance-based | Simple, no training | Slow prediction, curse of dimensionality |
| Naive Bayes | Probabilistic | Fast, works with small data | Assumes feature independence |
| Gradient Boosting | Ensemble | High accuracy | Slow training |
| MLP Neural Network | Deep Learning | Captures complex patterns | Requires more data, prone to overfit |

---

## Model Descriptions

### 1. Random Forest

**Algorithm**: Ensemble of decision trees with bootstrap aggregating (bagging).

```python
RandomForestClassifier(
    n_estimators=100,      # Number of trees
    max_depth=10,          # Maximum tree depth
    random_state=42,
    n_jobs=-1              # Use all CPUs
)
```

**How it works**:

1. Create n_estimators bootstrap samples
2. Train a decision tree on each sample
3. Each tree uses random subset of features at each split
4. Final prediction = majority vote

**Best for**: General-purpose classification, feature importance analysis.

### 2. XGBoost

**Algorithm**: Gradient boosting with regularization.

```python
XGBClassifier(
    n_estimators=100,
    max_depth=5,
    learning_rate=0.1,
    random_state=42,
    n_jobs=-1,
    verbosity=0
)
```

**How it works**:

1. Start with base prediction
2. Iteratively add trees that correct residual errors
3. Regularization prevents overfitting
4. Uses gradient descent for optimization

**Best for**: High-accuracy classification, handling imbalanced data.

### 3. Logistic Regression

**Algorithm**: Linear model with softmax for multi-class.

```python
LogisticRegression(
    max_iter=1000,
    random_state=42,
    multi_class='multinomial',  # Multi-class softmax
    n_jobs=-1
)
```

**How it works**:

1. Compute linear combination: z = Σ wᵢxᵢ + b
2. Apply softmax: P(class k) = exp(zₖ) / Σⱼ exp(zⱼ)
3. Optimize with maximum likelihood

**Best for**: Interpretable models, baseline comparison.

### 4. Support Vector Machine (SVM)

#### RBF Kernel

```python
SVC(kernel='rbf', random_state=42, probability=True)
```

**How it works**:

1. Map data to higher dimension using RBF kernel: K(x, y) = exp(-γ||x-y||²)
2. Find hyperplane that maximizes margin
3. Classify based on which side of hyperplane

**Best for**: Non-linear boundaries, high-dimensional data.

#### Linear Kernel

```python
SVC(kernel='linear', random_state=42, probability=True)
```

**Best for**: Linearly separable data, faster than RBF.

### 5. K-Nearest Neighbors (KNN)

```python
KNeighborsClassifier(n_neighbors=5, n_jobs=-1)
```

**How it works**:

1. For new sample, find k nearest training samples
2. Majority vote determines class
3. Distance metric: Euclidean (default)

**Best for**: Simple problems, local patterns.

### 6. Naive Bayes (Gaussian)

```python
GaussianNB()
```

**How it works**:

1. Assumes features are conditionally independent given class
2. Models each feature as Gaussian distribution
3. Uses Bayes' theorem: P(class|features) ∝ P(features|class) × P(class)

**Best for**: Small datasets, fast training/prediction.

### 7. Gradient Boosting

```python
GradientBoostingClassifier(
    n_estimators=100,
    max_depth=5,
    random_state=42
)
```

**How it works**:

1. Similar to XGBoost but without regularization tricks
2. Sequentially adds trees to minimize loss function
3. Each tree fits the negative gradient of loss

**Best for**: High accuracy, interpretable feature importance.

### 8. MLP Neural Network

```python
MLPClassifier(
    hidden_layer_sizes=(100, 50),  # Two hidden layers
    max_iter=500,
    random_state=42,
    early_stopping=True
)
```

**How it works**:

1. Input layer → hidden layers → output layer
2. Each neuron: activation(Σ wᵢxᵢ + b)
3. Backpropagation for learning
4. Early stopping prevents overfitting

**Best for**: Complex non-linear patterns, large datasets.

---

## Hyperparameters

### Default Settings (in config.py)

```python
@dataclass(frozen=True)
class MLConfig:
    # Data splits
    TEST_SIZE: float = 0.2
    RANDOM_STATE: int = 42
    
    # Cross-validation
    CV_FOLDS: int = 5
    
    # Random Forest
    RF_N_ESTIMATORS: int = 200
    RF_MAX_DEPTH: Optional[int] = None
    RF_MIN_SAMPLES_SPLIT: int = 2
    RF_MIN_SAMPLES_LEAF: int = 1
    
    # XGBoost
    XGB_N_ESTIMATORS: int = 200
    XGB_MAX_DEPTH: int = 6
    XGB_LEARNING_RATE: float = 0.1
    XGB_SUBSAMPLE: float = 0.8
    
    # Logistic Regression
    LR_MAX_ITER: int = 1000
    LR_SOLVER: str = "lbfgs"
    
    # Feature selection
    TOP_N_FEATURES: int = 25
```

### Tuning Recommendations

| Model | Parameter | Low Value Effect | High Value Effect |
|-------|-----------|------------------|-------------------|
| RF | n_estimators | Underfit, high variance | More stable, slower |
| RF | max_depth | Underfit | Overfit |
| XGB | learning_rate | Need more trees | Faster convergence, may overfit |
| XGB | max_depth | Underfit | Overfit |
| SVM | C (regularization) | Underfit | Overfit |
| SVM | gamma (RBF) | Smooth boundary | Complex boundary |
| KNN | n_neighbors | Sensitive to noise | Over-smooth |
| MLP | hidden_layers | Underfit | Overfit, slow |

---

## Model Selection Guidelines

### For This Pipeline (Ancestry Classification)

**Recommended**: Random Forest or XGBoost

Reasons:

1. Handle non-linear relationships in genotype data
2. Robust to outliers and missing values
3. Provide feature importance for SNP ranking
4. Good performance with moderate sample sizes

### When to Use Each Model

| Scenario | Recommended Model |
|----------|-------------------|
| Quick baseline | Logistic Regression |
| Interpretability needed | Logistic Regression |
| Maximum accuracy | XGBoost |
| Feature importance | Random Forest |
| Small dataset | Naive Bayes |
| Non-linear boundaries | SVM (RBF) |
| Very large dataset | Logistic Regression, Linear SVM |

### Data Preprocessing

| Model | Scaling Required | Encoding |
|-------|------------------|----------|
| Random Forest | No | Numeric (0, 1, 2) |
| XGBoost | No | Numeric |
| Logistic Regression | **Yes** | Numeric |
| SVM | **Yes** | Numeric |
| KNN | **Yes** | Numeric |
| Naive Bayes | Optional | Numeric |
| MLP | **Yes** | Numeric |

**Scaling**: Use `StandardScaler` for SVM, LR, MLP, KNN.

---

## Evaluation Metrics

### Cross-Validation

```python
from sklearn.model_selection import StratifiedKFold, cross_validate

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scoring = {
    'accuracy': 'accuracy',
    'f1_weighted': 'f1_weighted',
    'precision_weighted': 'precision_weighted',
    'recall_weighted': 'recall_weighted'
}
scores = cross_validate(clf, X, y, cv=cv, scoring=scoring)
```

### Metrics Explained

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **Accuracy** | (TP + TN) / Total | Overall correctness |
| **Precision** | TP / (TP + FP) | Of predicted positives, how many correct |
| **Recall** | TP / (TP + FN) | Of actual positives, how many found |
| **F1 Score** | 2 × (P × R) / (P + R) | Harmonic mean of precision and recall |

### Multi-Class Averaging

| Type | Description |
|------|-------------|
| **macro** | Unweighted mean across classes |
| **weighted** | Weighted by class support (used in pipeline) |
| **micro** | Global TP, FP, FN counts |

### Overfit Detection

```python
overfit_gap = train_accuracy - test_accuracy
```

| Gap | Interpretation |
|-----|----------------|
| < 0.05 | Good generalization |
| 0.05 - 0.10 | Mild overfitting |
| > 0.10 | Significant overfitting |

---

## Feature Importance

### Random Forest Importance

Based on mean decrease in impurity (Gini importance):

```python
rf = RandomForestClassifier(n_estimators=100)
rf.fit(X, y)
importance = rf.feature_importances_
```

**Interpretation**: How much each feature reduces impurity across all trees.

### Logistic Regression Coefficients

```python
lr = LogisticRegression(multi_class='multinomial')
lr.fit(X_scaled, y)

# Mean absolute coefficient across classes
coef_importance = np.abs(lr.coef_).mean(axis=0)
```

**Interpretation**: Magnitude of effect on log-odds per unit change in feature.

### Combined Ranking

```python
# Rank by both RF and LR
importance_df['rank_rf'] = importance_df['rf_importance'].rank(ascending=False)
importance_df['rank_lr'] = importance_df['lr_importance'].rank(ascending=False)
importance_df['avg_rank'] = (importance_df['rank_rf'] + importance_df['rank_lr']) / 2
```

### Reduced Feature Set Testing

Test performance with fewer SNPs:

```python
feature_sizes = [5, 10, 15, 20, 25]

for n in feature_sizes:
    top_snps = importance_df.nsmallest(n, 'avg_rank')['snp'].tolist()
    X_reduced = X[:, snp_indices]
    scores = cross_val_score(clf, X_reduced, y, cv=cv)
```

This helps identify the minimum number of SNPs needed for good classification.

---

## Code Examples

### Full Training Pipeline

```python
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

# Load data
df = pd.read_csv('statistical_ml_data_02b.csv')
snp_cols = [c for c in df.columns if c not in ['sample', 'pop']]
X = df[snp_cols].values
y = df['pop'].values

# Encode labels
le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Scale for some models
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Cross-validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Train Random Forest
rf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
rf_scores = cross_validate(rf, X, y_encoded, cv=cv, 
                           scoring=['accuracy', 'f1_weighted'])

print(f"RF Accuracy: {rf_scores['test_accuracy'].mean():.4f}")
print(f"RF F1: {rf_scores['test_f1_weighted'].mean():.4f}")
```

---

*See also: [PIPELINE.md](PIPELINE.md), [STATISTICAL_TESTS.md](STATISTICAL_TESTS.md)*
