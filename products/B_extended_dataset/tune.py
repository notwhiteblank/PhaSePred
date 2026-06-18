#!/usr/bin/env python
"""Product B — Optuna hyperparameter tuning on the updated dataset.

Builds the updated label table (PhaSepDB3 + paper's 10-organism negatives
+ LLPSDB2 confirmed negatives) and a stratified 80/20 train/test split
per task, then tunes XGBoost on the training partition only with 5-fold
stratified ROC-AUC.

Outputs:
  products/B_extended_dataset/tuned_params.json
  products/B_extended_dataset/train_accessions_<task>.tsv
  products/B_extended_dataset/test_accessions_<task>.tsv

Anti-leakage:
  - The test partition is sealed before any tuning starts.
  - Optuna inner CV runs on the training partition only.
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
from phasepred.updated_data import build_label_table, build_task_split  # noqa: E402

PRODUCT_DIR = REPO_ROOT / "products" / "B_extended_dataset"
RECOMPUTED_PATH = REPO_ROOT / "data" / "interim" / "recomputed_features_full.csv"

SEED = 42
N_TRIALS = 100
CV_FOLDS = 5
TEST_SIZE = 0.20

TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")

optuna.logging.set_verbosity(optuna.logging.WARNING)


def _feature_columns(task: str) -> list[str]:
    return list(HUMAN_FEATURE_COLUMNS) if task.startswith("h") else list(BASE_FEATURE_COLUMNS)


def _join(split_df: pd.DataFrame, recomputed: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    merged = split_df.merge(recomputed, on="UniprotEntry", how="left")
    missing_all = merged[cols].isna().all(axis=1)
    return merged[~missing_all].drop_duplicates("UniprotEntry").reset_index(drop=True)


def main() -> None:
    print(f"Building updated label table...")
    label_table = build_label_table(REPO_ROOT)
    print(f"  {len(label_table)} label rows")

    print(f"Loading recomputed features from {RECOMPUTED_PATH.name}...")
    recomputed = pd.read_csv(RECOMPUTED_PATH)
    print(f"  {len(recomputed)} proteins, columns: {list(recomputed.columns)}")

    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    tuned: dict[str, dict] = {}

    for task in TASKS:
        cols = _feature_columns(task)
        train_split, test_split = build_task_split(
            label_table, RECOMPUTED_PATH,
            task=task, test_size=TEST_SIZE, random_state=SEED,
        )
        train_df = _join(train_split, recomputed, cols)
        test_df = _join(test_split, recomputed, cols)

        # Sealing first.
        (PRODUCT_DIR / f"train_accessions_{task}.tsv").write_text(
            "\n".join(train_df["UniprotEntry"].tolist()) + "\n"
        )
        (PRODUCT_DIR / f"test_accessions_{task}.tsv").write_text(
            "\n".join(test_df["UniprotEntry"].tolist()) + "\n"
        )
        overlap = set(train_df["UniprotEntry"]) & set(test_df["UniprotEntry"])
        if overlap:
            raise RuntimeError(
                f"Leakage detected for {task}: {len(overlap)} accessions in "
                f"both train and test"
            )

        X_train = train_df[cols].astype(float)
        y_train = train_df["label"].astype(int).to_numpy()
        pos = int((y_train == 1).sum())
        neg = int((y_train == 0).sum())
        if pos == 0 or neg == 0:
            raise RuntimeError(f"Degenerate training set for {task}: pos={pos}, neg={neg}")
        pos_weight = float(neg) / float(pos)

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

        print(f"  Best 5-fold CV AUC: {study.best_value:.4f}")
        tuned[task] = {
            "best_cv_auc": float(study.best_value),
            "best_params": study.best_params,
            "feature_columns": cols,
            "pos_weight": pos_weight,
            "n_trials": N_TRIALS,
            "cv_folds": CV_FOLDS,
            "seed": SEED,
            "test_size": TEST_SIZE,
        }

    out = PRODUCT_DIR / "tuned_params.json"
    out.write_text(json.dumps(tuned, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
