"""
Machine Learning Training Pipeline for Population Classification.

Supports multiple classifiers:
- Random Forest
- XGBoost
- Logistic Regression

Features:
- Train/test splitting with stratification
- Cross-validation
- Feature importance analysis
- Model evaluation and reporting
"""

import os
import pickle
from pathlib import Path
from typing import Optional, List, Tuple, Dict, Any, Union
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    f1_score,
)

# Import XGBoost if available
try:
    import xgboost as xgb

    HAS_XGBOOST = True
except ImportError:
    HAS_XGBOOST = False

# Import configuration
try:
    from config import PATHS, ML
    from utils import ensure_dir, save_report
except ImportError:
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from config import PATHS, ML
    from utils import ensure_dir, save_report


# =============================================================================
# Data Loading
# =============================================================================


def load_ml_data(
    data_path: Union[str, Path] = None,
    target_column: str = "pop",
    sample_column: str = "sample",
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.Series, List[str]]:
    """
    Load ML-ready data from CSV.

    Args:
        data_path: Path to CSV with genotype matrix and population labels
        target_column: Column name for target labels
        sample_column: Column name for sample IDs
        verbose: Print progress

    Returns:
        X: Feature matrix (genotypes)
        y: Target labels (populations)
        feature_names: List of feature (variant) names
    """
    if data_path is None:
        data_path = PATHS.ML_DATA

    if verbose:
        print(f"Loading data from: {data_path}")

    df = pd.read_csv(data_path)

    # Separate features and target
    drop_cols = [sample_column, target_column]
    if "super_pop" in df.columns:
        drop_cols.append("super_pop")

    X = df.drop(columns=drop_cols, errors="ignore")
    y = df[target_column]
    feature_names = X.columns.tolist()

    if verbose:
        print(f"  Samples: {len(df)}")
        print(f"  Features: {len(feature_names)}")
        print(f"  Classes: {y.nunique()} ({y.value_counts().to_dict()})")

    return X, y, feature_names


# =============================================================================
# Model Training
# =============================================================================


def train_random_forest(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_estimators: int = ML.RF_N_ESTIMATORS,
    max_depth: Optional[int] = ML.RF_MAX_DEPTH,
    random_state: int = ML.RANDOM_STATE,
    verbose: bool = True,
) -> RandomForestClassifier:
    """
    Train a Random Forest classifier.

    Args:
        X_train: Training features
        y_train: Training labels
        n_estimators: Number of trees
        max_depth: Maximum tree depth (None for unlimited)
        random_state: Random seed
        verbose: Print progress

    Returns:
        Trained RandomForestClassifier
    """
    if verbose:
        print(f"\nTraining Random Forest...")
        print(f"  n_estimators: {n_estimators}")
        print(f"  max_depth: {max_depth}")

    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=ML.RF_MIN_SAMPLES_SPLIT,
        min_samples_leaf=ML.RF_MIN_SAMPLES_LEAF,
        random_state=random_state,
        n_jobs=-1,
    )

    clf.fit(X_train, y_train)

    if verbose:
        print("  Training complete!")

    return clf


def train_xgboost(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    n_estimators: int = ML.XGB_N_ESTIMATORS,
    max_depth: int = ML.XGB_MAX_DEPTH,
    learning_rate: float = ML.XGB_LEARNING_RATE,
    random_state: int = ML.RANDOM_STATE,
    verbose: bool = True,
) -> Any:
    """
    Train an XGBoost classifier.

    Args:
        X_train: Training features
        y_train: Training labels (must be encoded as integers)
        n_estimators: Number of boosting rounds
        max_depth: Maximum tree depth
        learning_rate: Learning rate
        random_state: Random seed
        verbose: Print progress

    Returns:
        Trained XGBClassifier
    """
    if not HAS_XGBOOST:
        raise ImportError("XGBoost not installed. Install with: pip install xgboost")

    if verbose:
        print(f"\nTraining XGBoost...")
        print(f"  n_estimators: {n_estimators}")
        print(f"  max_depth: {max_depth}")
        print(f"  learning_rate: {learning_rate}")

    # Encode labels if needed
    if y_train.dtype == object:
        le = LabelEncoder()
        y_encoded = le.fit_transform(y_train)
    else:
        y_encoded = y_train
        le = None

    clf = xgb.XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        subsample=ML.XGB_SUBSAMPLE,
        use_label_encoder=False,
        eval_metric="mlogloss",
        random_state=random_state,
        n_jobs=-1,
    )

    clf.fit(X_train, y_encoded)

    # Store label encoder for later use
    clf.label_encoder_ = le

    if verbose:
        print("  Training complete!")

    return clf


def train_logistic_regression(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    max_iter: int = ML.LR_MAX_ITER,
    solver: str = ML.LR_SOLVER,
    random_state: int = ML.RANDOM_STATE,
    verbose: bool = True,
) -> LogisticRegression:
    """
    Train a Logistic Regression classifier.

    Args:
        X_train: Training features
        y_train: Training labels
        max_iter: Maximum iterations
        solver: Optimization algorithm
        random_state: Random seed
        verbose: Print progress

    Returns:
        Trained LogisticRegression
    """
    if verbose:
        print(f"\nTraining Logistic Regression...")
        print(f"  max_iter: {max_iter}")
        print(f"  solver: {solver}")

    clf = LogisticRegression(
        multi_class="multinomial",
        solver=solver,
        max_iter=max_iter,
        random_state=random_state,
        n_jobs=-1,
    )

    clf.fit(X_train, y_train)

    if verbose:
        print("  Training complete!")

    return clf


# =============================================================================
# Model Evaluation
# =============================================================================


def evaluate_model(
    clf: Any,
    X_test: pd.DataFrame,
    y_test: pd.Series,
    model_name: str = "Model",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Evaluate a trained classifier.

    Args:
        clf: Trained classifier
        X_test: Test features
        y_test: Test labels
        model_name: Name for reporting
        verbose: Print results

    Returns:
        Dictionary with evaluation metrics
    """
    # Handle XGBoost label encoding
    if hasattr(clf, "label_encoder_") and clf.label_encoder_ is not None:
        y_test_encoded = clf.label_encoder_.transform(y_test)
        y_pred = clf.predict(X_test)
        y_test_for_report = y_test_encoded
        target_names = clf.label_encoder_.classes_
    else:
        y_pred = clf.predict(X_test)
        y_test_for_report = y_test
        target_names = None

    # Calculate metrics
    accuracy = accuracy_score(y_test_for_report, y_pred)
    f1_macro = f1_score(y_test_for_report, y_pred, average="macro")
    f1_weighted = f1_score(y_test_for_report, y_pred, average="weighted")

    report = classification_report(
        y_test_for_report,
        y_pred,
        target_names=target_names,
    )

    conf_matrix = confusion_matrix(y_test_for_report, y_pred)

    results = {
        "model_name": model_name,
        "accuracy": accuracy,
        "f1_macro": f1_macro,
        "f1_weighted": f1_weighted,
        "classification_report": report,
        "confusion_matrix": conf_matrix,
    }

    if verbose:
        print(f"\n{'='*60}")
        print(f"{model_name} EVALUATION")
        print(f"{'='*60}")
        print(f"Accuracy: {accuracy:.4f}")
        print(f"F1 (macro): {f1_macro:.4f}")
        print(f"F1 (weighted): {f1_weighted:.4f}")
        print(f"\nClassification Report:")
        print(report)

    return results


def cross_validate_model(
    clf: Any,
    X: pd.DataFrame,
    y: pd.Series,
    cv: int = ML.CV_FOLDS,
    model_name: str = "Model",
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Perform cross-validation on a classifier.

    Args:
        clf: Classifier (will be cloned for each fold)
        X: Full feature matrix
        y: Full labels
        cv: Number of folds
        model_name: Name for reporting
        verbose: Print results

    Returns:
        Dictionary with CV scores
    """
    if verbose:
        print(f"\nCross-validating {model_name} ({cv}-fold)...")

    # Handle label encoding for XGBoost
    y_cv = y
    if hasattr(clf, "label_encoder_") and clf.label_encoder_ is not None:
        le = LabelEncoder()
        y_cv = le.fit_transform(y)

    scores = cross_val_score(clf, X, y_cv, cv=cv, scoring="accuracy")

    results = {
        "model_name": model_name,
        "cv_folds": cv,
        "cv_scores": scores,
        "cv_mean": scores.mean(),
        "cv_std": scores.std(),
    }

    if verbose:
        print(f"  CV Accuracy: {scores.mean():.4f} (+/- {scores.std()*2:.4f})")

    return results


# =============================================================================
# Feature Importance
# =============================================================================


def get_feature_importances(
    clf: Any,
    feature_names: List[str],
    top_n: int = ML.TOP_N_FEATURES,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Extract feature importances from a trained classifier.

    Args:
        clf: Trained classifier (must have feature_importances_ attribute)
        feature_names: List of feature names
        top_n: Number of top features to return
        verbose: Print results

    Returns:
        DataFrame with feature importances sorted by importance
    """
    if not hasattr(clf, "feature_importances_"):
        raise ValueError("Classifier does not have feature_importances_ attribute")

    importances = pd.DataFrame(
        {
            "feature": feature_names,
            "importance": clf.feature_importances_,
        }
    )

    importances = importances.sort_values("importance", ascending=False)
    top_features = importances.head(top_n)

    if verbose:
        print(f"\nTop {top_n} Features by Importance:")
        for i, row in top_features.iterrows():
            print(f"  {row['feature']}: {row['importance']:.6f}")

    return importances


def train_on_top_features(
    clf_class: type,
    X: pd.DataFrame,
    y: pd.Series,
    feature_importances: pd.DataFrame,
    top_n: int = ML.TOP_N_FEATURES,
    test_size: float = ML.TEST_SIZE,
    random_state: int = ML.RANDOM_STATE,
    verbose: bool = True,
    **clf_kwargs,
) -> Tuple[Any, Dict[str, Any]]:
    """
    Train a new classifier using only top N features.

    Args:
        clf_class: Classifier class to instantiate
        X: Full feature matrix
        y: Labels
        feature_importances: DataFrame from get_feature_importances
        top_n: Number of top features to use
        test_size: Test split ratio
        random_state: Random seed
        verbose: Print results
        **clf_kwargs: Additional arguments for classifier

    Returns:
        Trained classifier and evaluation results
    """
    # Get top features
    top_feature_names = feature_importances.head(top_n)["feature"].tolist()
    X_top = X[top_feature_names]

    if verbose:
        print(f"\nTraining on top {top_n} features...")

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X_top,
        y,
        test_size=test_size,
        stratify=y,
        random_state=random_state,
    )

    # Train
    clf = clf_class(**clf_kwargs)
    clf.fit(X_train, y_train)

    # Evaluate
    results = evaluate_model(
        clf,
        X_test,
        y_test,
        model_name=f"{clf_class.__name__} (top {top_n})",
        verbose=verbose,
    )

    return clf, results


# =============================================================================
# Model Persistence
# =============================================================================


def save_model(
    clf: Any,
    output_path: Union[str, Path],
    metadata: Optional[Dict] = None,
    verbose: bool = True,
) -> None:
    """
    Save a trained model to disk.

    Args:
        clf: Trained classifier
        output_path: Path for saving
        metadata: Optional metadata dict to save with model
        verbose: Print progress
    """
    output_path = Path(output_path)
    ensure_dir(output_path.parent)

    save_obj = {
        "model": clf,
        "metadata": metadata or {},
        "saved_at": datetime.now().isoformat(),
    }

    with open(output_path, "wb") as f:
        pickle.dump(save_obj, f)

    if verbose:
        print(f"Model saved to: {output_path}")


def load_model(
    model_path: Union[str, Path],
    verbose: bool = True,
) -> Tuple[Any, Dict]:
    """
    Load a trained model from disk.

    Args:
        model_path: Path to saved model
        verbose: Print progress

    Returns:
        Trained classifier and metadata dict
    """
    with open(model_path, "rb") as f:
        save_obj = pickle.load(f)

    if verbose:
        print(f"Model loaded from: {model_path}")
        print(f"  Saved at: {save_obj.get('saved_at', 'unknown')}")

    return save_obj["model"], save_obj.get("metadata", {})


# =============================================================================
# Full Training Pipeline
# =============================================================================


def run_full_pipeline(
    data_path: Union[str, Path] = None,
    target_column: str = "pop",
    test_size: float = ML.TEST_SIZE,
    random_state: int = ML.RANDOM_STATE,
    train_rf: bool = True,
    train_xgb: bool = True,
    train_lr: bool = True,
    save_models: bool = True,
    output_dir: Union[str, Path] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run complete ML training pipeline.

    Args:
        data_path: Path to ML data CSV
        target_column: Target column name
        test_size: Test split ratio
        random_state: Random seed
        train_rf: Train Random Forest
        train_xgb: Train XGBoost
        train_lr: Train Logistic Regression
        save_models: Save trained models
        output_dir: Directory for outputs
        verbose: Print progress

    Returns:
        Dictionary with all results
    """
    if output_dir is None:
        output_dir = PATHS.ML_MODELS_DIR
    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    print("=" * 60)
    print("ML TRAINING PIPELINE")
    print("=" * 60)

    # Load data
    X, y, feature_names = load_ml_data(data_path, target_column, verbose=verbose)

    # Encode labels
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=test_size,
        stratify=y,
        random_state=random_state,
    )

    y_train_encoded = le.transform(y_train)
    y_test_encoded = le.transform(y_test)

    if verbose:
        print(f"\nData split:")
        print(f"  Train: {len(X_train)}")
        print(f"  Test: {len(X_test)}")

    results = {
        "feature_names": feature_names,
        "label_encoder": le,
        "models": {},
        "evaluations": {},
    }

    # Train Random Forest
    if train_rf:
        rf = train_random_forest(X_train, y_train, verbose=verbose)
        rf_eval = evaluate_model(rf, X_test, y_test, "Random Forest", verbose=verbose)
        results["models"]["random_forest"] = rf
        results["evaluations"]["random_forest"] = rf_eval

        # Feature importance
        rf_importance = get_feature_importances(rf, feature_names, verbose=verbose)
        results["feature_importances"] = rf_importance

        if save_models:
            save_model(rf, output_dir / "random_forest.pkl", verbose=verbose)

    # Train XGBoost
    if train_xgb and HAS_XGBOOST:
        xgb_clf = train_xgboost(X_train, y_train_encoded, verbose=verbose)
        xgb_eval = evaluate_model(
            xgb_clf, X_test, y_test_encoded, "XGBoost", verbose=verbose
        )
        results["models"]["xgboost"] = xgb_clf
        results["evaluations"]["xgboost"] = xgb_eval

        if save_models:
            save_model(xgb_clf, output_dir / "xgboost.pkl", verbose=verbose)

    # Train Logistic Regression
    if train_lr:
        lr = train_logistic_regression(X_train, y_train, verbose=verbose)
        lr_eval = evaluate_model(
            lr, X_test, y_test, "Logistic Regression", verbose=verbose
        )
        results["models"]["logistic_regression"] = lr
        results["evaluations"]["logistic_regression"] = lr_eval

        if save_models:
            save_model(lr, output_dir / "logistic_regression.pkl", verbose=verbose)

    # Generate summary report
    report_lines = [
        "=" * 60,
        "ML TRAINING SUMMARY",
        "=" * 60,
        f"Date: {datetime.now().isoformat()}",
        f"Samples: {len(X)}",
        f"Features: {len(feature_names)}",
        f"Classes: {len(le.classes_)} ({list(le.classes_)})",
        "",
        "Model Performance:",
    ]

    for model_name, eval_result in results["evaluations"].items():
        report_lines.append(f"  {model_name}:")
        report_lines.append(f"    Accuracy: {eval_result['accuracy']:.4f}")
        report_lines.append(f"    F1 (macro): {eval_result['f1_macro']:.4f}")

    report_text = "\n".join(report_lines)
    save_report(report_text, output_dir / "training_summary.txt", verbose=verbose)

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)

    return results


# =============================================================================
# CLI Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Train ML models for population classification"
    )

    parser.add_argument("--data", help="Path to ML data CSV")
    parser.add_argument("--target", default="pop", help="Target column name")
    parser.add_argument(
        "--test-size", type=float, default=ML.TEST_SIZE, help="Test split ratio"
    )
    parser.add_argument("--output", "-o", help="Output directory")
    parser.add_argument("--no-rf", action="store_true", help="Skip Random Forest")
    parser.add_argument("--no-xgb", action="store_true", help="Skip XGBoost")
    parser.add_argument("--no-lr", action="store_true", help="Skip Logistic Regression")
    parser.add_argument("--no-save", action="store_true", help="Don't save models")

    args = parser.parse_args()

    run_full_pipeline(
        data_path=args.data,
        target_column=args.target,
        test_size=args.test_size,
        output_dir=args.output,
        train_rf=not args.no_rf,
        train_xgb=not args.no_xgb,
        train_lr=not args.no_lr,
        save_models=not args.no_save,
    )
