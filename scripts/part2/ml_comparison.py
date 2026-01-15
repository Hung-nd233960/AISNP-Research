"""
ml_comparison.py - ML model comparison utilities for Part 2

This module provides functions to:
1. Run multiple classifiers on different SNP sets
2. Generate confusion matrices
3. Create performance comparison plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Optional, Tuple
from pathlib import Path

# ML imports
from sklearn.model_selection import cross_val_score, StratifiedKFold, cross_validate
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
)

# Classifiers
from sklearn.ensemble import (
    RandomForestClassifier,
    GradientBoostingClassifier,
    AdaBoostClassifier,
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier


def get_default_classifiers() -> Dict:
    """Get dictionary of default classifiers for comparison."""
    return {
        "Random Forest": RandomForestClassifier(
            n_estimators=100, max_depth=10, random_state=42, n_jobs=-1
        ),
        "Logistic Regression": LogisticRegression(
            max_iter=1000, random_state=42, multi_class="multinomial", n_jobs=-1
        ),
        "SVM (RBF)": SVC(kernel="rbf", random_state=42, probability=True),
        "SVM (Linear)": SVC(kernel="linear", random_state=42, probability=True),
        "K-Nearest Neighbors": KNeighborsClassifier(n_neighbors=5, n_jobs=-1),
        "Naive Bayes": GaussianNB(),
        "Gradient Boosting": GradientBoostingClassifier(
            n_estimators=100, max_depth=5, random_state=42
        ),
        "MLP Neural Network": MLPClassifier(
            hidden_layer_sizes=(100, 50),
            max_iter=500,
            random_state=42,
            early_stopping=True,
        ),
    }


def prepare_data(
    df: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, LabelEncoder]:
    """
    Prepare data for ML training.

    Args:
        df: DataFrame with sample, pop, and SNP columns

    Returns:
        X, X_scaled, y_encoded, label_encoder
    """
    snp_cols = [c for c in df.columns if c not in ["sample", "pop", "source"]]

    X = df[snp_cols].values
    y = df["pop"].values

    # Handle missing values
    X = np.nan_to_num(X, nan=0)

    # Encode labels
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    return X, X_scaled, y_encoded, le


def run_model_comparison(
    df: pd.DataFrame,
    classifiers: Optional[Dict] = None,
    n_folds: int = 5,
    source_name: str = "unknown",
) -> pd.DataFrame:
    """
    Run K-fold cross-validation on multiple classifiers.

    Args:
        df: ML-ready DataFrame
        classifiers: Dict of classifier name -> classifier object
        n_folds: Number of CV folds
        source_name: Name of SNP source for labeling

    Returns:
        DataFrame with CV results
    """
    if classifiers is None:
        classifiers = get_default_classifiers()

    X, X_scaled, y_encoded, le = prepare_data(df)
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)

    n_snps = X.shape[1]
    print(f"\n{'='*70}")
    print(f"Running {n_folds}-Fold CV on {source_name}")
    print(f"Samples: {len(df)}, SNPs: {n_snps}, Classes: {list(le.classes_)}")
    print(f"{'='*70}")

    results = []

    for name, clf in classifiers.items():
        # Select appropriate data
        if any(x in name for x in ["SVM", "Logistic", "MLP"]):
            X_use = X_scaled
        else:
            X_use = X

        try:
            scoring = {
                "accuracy": "accuracy",
                "f1_weighted": "f1_weighted",
                "precision_weighted": "precision_weighted",
                "recall_weighted": "recall_weighted",
            }

            scores = cross_validate(
                clf, X_use, y_encoded, cv=cv, scoring=scoring, return_train_score=True
            )

            result = {
                "Source": source_name,
                "N_SNPs": n_snps,
                "Model": name,
                "Accuracy_Mean": scores["test_accuracy"].mean(),
                "Accuracy_Std": scores["test_accuracy"].std(),
                "F1_Mean": scores["test_f1_weighted"].mean(),
                "F1_Std": scores["test_f1_weighted"].std(),
                "Precision_Mean": scores["test_precision_weighted"].mean(),
                "Recall_Mean": scores["test_recall_weighted"].mean(),
                "Train_Accuracy": scores["train_accuracy"].mean(),
                "Overfit_Gap": scores["train_accuracy"].mean()
                - scores["test_accuracy"].mean(),
            }
            results.append(result)

            print(
                f"  {name}: Acc={result['Accuracy_Mean']:.4f}±{result['Accuracy_Std']:.4f}, "
                f"F1={result['F1_Mean']:.4f}"
            )

        except Exception as e:
            print(f"  {name}: Error - {e}")
            results.append(
                {
                    "Source": source_name,
                    "N_SNPs": n_snps,
                    "Model": name,
                    "Accuracy_Mean": np.nan,
                    "Accuracy_Std": np.nan,
                    "F1_Mean": np.nan,
                    "F1_Std": np.nan,
                    "Precision_Mean": np.nan,
                    "Recall_Mean": np.nan,
                    "Train_Accuracy": np.nan,
                    "Overfit_Gap": np.nan,
                }
            )

    return pd.DataFrame(results)


def generate_confusion_matrices(
    df: pd.DataFrame, classifiers: Optional[Dict] = None, test_size: float = 0.2
) -> Dict[str, np.ndarray]:
    """
    Generate confusion matrices for each classifier.

    Args:
        df: ML-ready DataFrame
        classifiers: Dict of classifiers
        test_size: Fraction for test split

    Returns:
        Dict of classifier name -> confusion matrix
    """
    from sklearn.model_selection import train_test_split

    if classifiers is None:
        classifiers = get_default_classifiers()

    X, X_scaled, y_encoded, le = prepare_data(df)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, test_size=test_size, random_state=42, stratify=y_encoded
    )
    X_train_scaled, X_test_scaled, _, _ = train_test_split(
        X_scaled, y_encoded, test_size=test_size, random_state=42, stratify=y_encoded
    )

    conf_matrices = {}

    for name, clf in classifiers.items():
        if any(x in name for x in ["SVM", "Logistic", "MLP"]):
            clf.fit(X_train_scaled, y_train)
            y_pred = clf.predict(X_test_scaled)
        else:
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)

        cm = confusion_matrix(y_test, y_pred)
        conf_matrices[name] = cm

    return conf_matrices, le.classes_


def plot_confusion_matrices(
    conf_matrices: Dict[str, np.ndarray],
    class_names: List[str],
    figsize: Tuple[int, int] = (16, 12),
    save_path: Optional[str] = None,
):
    """
    Plot confusion matrices for all classifiers.

    Args:
        conf_matrices: Dict of classifier name -> confusion matrix
        class_names: List of class labels
        figsize: Figure size
        save_path: Optional path to save figure
    """
    n_classifiers = len(conf_matrices)
    n_cols = 3
    n_rows = (n_classifiers + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()

    for idx, (name, cm) in enumerate(conf_matrices.items()):
        ax = axes[idx]
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap="Blues",
            ax=ax,
            xticklabels=class_names,
            yticklabels=class_names,
        )
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("Actual")

    # Hide empty subplots
    for idx in range(n_classifiers, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Confusion matrices saved: {save_path}")

    plt.show()


def plot_performance_comparison(
    results_df: pd.DataFrame,
    metric: str = "Accuracy_Mean",
    figsize: Tuple[int, int] = (14, 8),
    save_path: Optional[str] = None,
):
    """
    Plot performance comparison across sources and models.

    Args:
        results_df: DataFrame from run_model_comparison
        metric: Metric to plot
        figsize: Figure size
        save_path: Optional path to save figure
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Bar chart by source
    ax = axes[0]
    sources = results_df["Source"].unique()
    x = np.arange(len(results_df["Model"].unique()))
    width = 0.8 / len(sources)

    for i, source in enumerate(sources):
        source_data = results_df[results_df["Source"] == source]
        source_data = source_data.sort_values("Model")
        ax.bar(x + i * width, source_data[metric], width, label=source, alpha=0.8)

    ax.set_xlabel("Model")
    ax.set_ylabel(metric)
    ax.set_title(f"{metric} by Model and Source")
    ax.set_xticks(x + width * (len(sources) - 1) / 2)
    ax.set_xticklabels(source_data["Model"].values, rotation=45, ha="right", fontsize=8)
    ax.legend(title="Source")
    ax.grid(True, alpha=0.3, axis="y")

    # Plot 2: Best model per source
    ax = axes[1]
    best_per_source = results_df.loc[results_df.groupby("Source")[metric].idxmax()]
    colors = plt.cm.Set2(np.linspace(0, 1, len(best_per_source)))

    bars = ax.bar(best_per_source["Source"], best_per_source[metric], color=colors)
    ax.set_xlabel("Source")
    ax.set_ylabel(metric)
    ax.set_title(f"Best {metric} per Source")
    ax.set_xticklabels(best_per_source["Source"], rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    # Add best model name on bars
    for bar, model in zip(bars, best_per_source["Model"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            model,
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=45,
        )

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Performance comparison saved: {save_path}")

    plt.show()


def create_summary_report(results_df: pd.DataFrame, output_path: str) -> str:
    """
    Create text summary report of ML comparison.

    Args:
        results_df: DataFrame from run_model_comparison
        output_path: Path to save report

    Returns:
        Report text
    """
    lines = []
    lines.append("=" * 70)
    lines.append("ML MODEL COMPARISON REPORT")
    lines.append("=" * 70)
    lines.append("")

    for source in results_df["Source"].unique():
        source_data = results_df[results_df["Source"] == source]
        n_snps = source_data["N_SNPs"].iloc[0]

        lines.append(f"\n{'='*50}")
        lines.append(f"Source: {source}")
        lines.append(f"Number of SNPs: {n_snps}")
        lines.append(f"{'='*50}")

        best = source_data.loc[source_data["Accuracy_Mean"].idxmax()]
        lines.append(f"\nBest Model: {best['Model']}")
        lines.append(
            f"  Accuracy: {best['Accuracy_Mean']:.4f} ± {best['Accuracy_Std']:.4f}"
        )
        lines.append(f"  F1 Score: {best['F1_Mean']:.4f}")
        lines.append(f"  Overfit Gap: {best['Overfit_Gap']:.4f}")

        lines.append(f"\nAll Models:")
        for _, row in source_data.sort_values(
            "Accuracy_Mean", ascending=False
        ).iterrows():
            lines.append(
                f"  {row['Model']}: {row['Accuracy_Mean']:.4f} ± {row['Accuracy_Std']:.4f}"
            )

    lines.append("\n" + "=" * 70)
    lines.append("END OF REPORT")
    lines.append("=" * 70)

    report = "\n".join(lines)

    with open(output_path, "w") as f:
        f.write(report)

    print(f"Report saved: {output_path}")
    return report


if __name__ == "__main__":
    print("ml_comparison.py - ML comparison utilities for Part 2")
    print("Import this module to use the functions.")
