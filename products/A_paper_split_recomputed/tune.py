#!/usr/bin/env python
"""Product A — Optuna hyperparameter tuning.

Tunes XGBoost for each of the four tasks (SaPS, PdPS, hSaPS, hPdPS) on the
paper S2 training partition only, using **recomputed** feature values
(`data/interim/recomputed_features_full.csv`), with a stratified 5-fold
ROC-AUC objective. The held-out S3 partition is sealed until train.py runs.

Outputs:
  products/A_paper_split_recomputed/tuned_params.json
  products/A_paper_split_recomputed/train_accessions.tsv
  products/A_paper_split_recomputed/test_accessions.tsv

Anti-leakage:
  - Tuning never sees S3.
  - Inner CV folds split S2 by the task's stratification.
  - Train/test accession files are written and the train script refuses to
    proceed if their intersection is non-empty.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS  # noqa: E402
from phasepred.paper import load_paper_split  # noqa: E402

PRODUCT_DIR = REPO_ROOT / "products" / "A_paper_split_recomputed"
RECOMPUTED_PATH = REPO_ROOT / "data" / "interim" / "recomputed_features_full.csv"
S2_PATH = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd02.xlsx"
S3_PATH = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd03.xlsx"

SEED = 42
N_TRIALS = 100
CV_FOLDS = 5

TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")

optuna.logging.set_verbosity(optuna.logging.WARNING)


def _feature_columns(task: str) -> list[str]:
    return list(HUMAN_FEATURE_COLUMNS) if task.startswith("h") else list(BASE_FEATURE_COLUMNS)


def _join_features(split_df: pd.DataFrame, recomputed: pd.DataFrame, task: str) -> pd.DataFrame:
    """Take the paper split rows for one task and attach recomputed features.

    Rows where the recomputed feature table is missing the protein entirely
    are dropped — they would contribute NaN to every column. Rows where only
    some features are NaN are kept; SimpleImputer handles them in the pipeline.
    """
    task_rows = split_df.query("task == @task").copy()
    merged = task_rows.merge(recomputed, on="UniprotEntry", how="left")
    cols = _feature_columns(task)
    missing_all = merged[cols].isna().all(axis=1)
    dropped = merged[missing_all]
    kept = merged[~missing_all]
    if len(dropped):
        print(f"  [{task}] dropped {len(dropped)} rows missing all features (no recomputed match)")
    return kept


def main() -> None:
    print(f"Loading paper split from {S2_PATH.name} / {S3_PATH.name}...")
    split = load_paper_split(S2_PATH, S3_PATH)
    print(f"  S2 (train): {len(split.train)} rows; S3 (test): {len(split.test)} rows")

    print(f"Loading recomputed features from {RECOMPUTED_PATH.name}...")
    recomputed = pd.read_csv(RECOMPUTED_PATH)
    print(f"  {len(recomputed)} proteins, columns: {list(recomputed.columns)}")

    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    tuned: dict[str, dict] = {}
    train_accessions: dict[str, list[str]] = {}
    test_accessions: dict[str, list[str]] = {}

    for task in TASKS:
        cols = _feature_columns(task)
        train_df = _join_features(split.train, recomputed, task)
        test_df = _join_features(split.test, recomputed, task)
        # Drop within-task duplicate accessions, keep first occurrence
        train_df = train_df.drop_duplicates("UniprotEntry").reset_index(drop=True)
        test_df = test_df.drop_duplicates("UniprotEntry").reset_index(drop=True)

        X_train = train_df[cols].astype(float)
        y_train = train_df["label"].astype(int).to_numpy()

        # Sealing: record accessions before any model touches the data.
        train_accessions[task] = train_df["UniprotEntry"].tolist()
        test_accessions[task] = test_df["UniprotEntry"].tolist()
        overlap = set(train_accessions[task]) & set(test_accessions[task])
        if overlap:
            raise RuntimeError(
                f"Leakage detected in paper split for {task}: "
                f"{len(overlap)} accessions are in BOTH train and test"
            )

        pos = int((y_train == 1).sum())
        neg = int((y_train == 0).sum())
        pos_weight = float(neg) / float(pos) if pos else 1.0

        print()
        print(f"=== Tuning {task} ===")
        print(f"  Train rows: {len(X_train)}  (pos={pos}, neg={neg}, pos_weight={pos_weight:.2f})")
        print(f"  Features ({len(cols)}): {cols}")

        cv = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=SEED)

        def objective(trial: optuna.Trial, X=X_train, y=y_train, w=pos_weight, cv=cv) -> float:
            pipe = Pipeline([
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("clf", XGBClassifier(
                    random_state=SEED,
                    eval_metric="logloss",
                    scale_pos_weight=w,
                    n_estimators=trial.suggest_int("n_estimators", 100, 600, step=50),
                    max_depth=trial.suggest_int("max_depth", 3, 8),
                    learning_rate=trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
                    subsample=trial.suggest_float("subsample", 0.5, 1.0),
                    colsample_bytree=trial.suggest_float("colsample_bytree", 0.5, 1.0),
                    min_child_weight=trial.suggest_int("min_child_weight", 1, 20),
                    gamma=trial.suggest_float("gamma", 0.0, 2.0),
                    reg_alpha=trial.suggest_float("reg_alpha", 0.0, 2.0),
                    reg_lambda=trial.suggest_float("reg_lambda", 0.1, 10.0, log=True),
                )),
            ])
            scores = cross_val_score(pipe, X, y, cv=cv, scoring="roc_auc", n_jobs=1)
            return float(scores.mean())

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=SEED))
        study.optimize(objective, n_trials=N_TRIALS, n_jobs=8, show_progress_bar=False)

        best = study.best_params
        print(f"  Best 5-fold CV AUC: {study.best_value:.4f}")
        tuned[task] = {
            "best_cv_auc": float(study.best_value),
            "best_params": best,
            "feature_columns": cols,
            "pos_weight": pos_weight,
            "n_trials": N_TRIALS,
            "cv_folds": CV_FOLDS,
            "seed": SEED,
        }

    out = PRODUCT_DIR / "tuned_params.json"
    out.write_text(json.dumps(tuned, indent=2))
    print(f"\nWrote {out}")

    # Persist sealed accession lists for the train.py leakage check.
    for task in TASKS:
        (PRODUCT_DIR / f"train_accessions_{task}.tsv").write_text(
            "\n".join(train_accessions[task]) + "\n"
        )
        (PRODUCT_DIR / f"test_accessions_{task}.tsv").write_text(
            "\n".join(test_accessions[task]) + "\n"
        )
    print(f"Wrote sealed train/test accession files for {len(TASKS)} tasks.")


if __name__ == "__main__":
    main()
