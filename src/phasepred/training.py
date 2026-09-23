"""Paper-protocol training for the four PhaSePred models (v2022).

These functions implement the published models' protocol (Chen et al. 2022
PNAS; ``docs/RETRAIN_PROTOCOL.md``) over the committed training tables written
by :mod:`phasepred.training_data`. Byte-for-byte reproduction of the shipped
``8f_model_*.joblib`` artifacts depends on :func:`make_xgb`, :func:`fit_model`,
:func:`ensemble_score`, :func:`_fold_auc`, :func:`neg_sets`,
:func:`run_cross_validation`, :func:`train_final_ensemble` and
:func:`feature_columns` staying character-for-character identical to the
revision that produced them. ``feature_columns`` fixes the column order the
matrices are built in and is part of the frozen set. Any change to those eight
functions (or to :data:`SEED`, :data:`N_MODELS`, :data:`NEG_RATIO`,
:data:`N_FOLDS`) invalidates the reproduction and must be accompanied by a fresh
run of the byte-identical gate plus a ``docs/RETRAIN_PROTOCOL.md`` errata entry.

The XGBoost classifier keeps every library default except ``random_state``,
``n_jobs`` and ``eval_metric``; there is no imputation, no scaling and no
``scale_pos_weight``.

The public names drop the leading underscores of the originating script. The
private aliases ``_make_xgb``, ``_fit_model`` and ``_neg_sets`` below are the
original call tokens: the protocol bodies use them so that they remain
verbatim, which matters because ``run_cross_validation`` and
``train_final_ensemble`` already bind a local variable named ``neg_sets``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier

from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS
from phasepred.training_data import load_task_frame

SEED = 42
N_MODELS = 10
NEG_RATIO = 2
N_FOLDS = 5
TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")


def feature_columns(task: str) -> list[str]:
    return list(HUMAN_FEATURE_COLUMNS) if task.startswith("h") else list(BASE_FEATURE_COLUMNS)


def make_xgb(seed: int) -> XGBClassifier:
    """Paper-protocol XGBoost: library defaults, seed fixed."""
    return XGBClassifier(random_state=seed, n_jobs=4, eval_metric="logloss")


def fit_model(X: pd.DataFrame, y: np.ndarray, seed: int) -> XGBClassifier:
    model = _make_xgb(seed)
    model.fit(X, y)
    return model


def ensemble_score(models: list[XGBClassifier], X: pd.DataFrame) -> np.ndarray:
    probs = np.column_stack([m.predict_proba(X)[:, 1] for m in models])
    return probs.mean(axis=1)


def _fold_auc(
    models: list[XGBClassifier], X_val: pd.DataFrame, y_val: np.ndarray
) -> tuple[float, float]:
    """Return (ensemble_auc, per_model_mean_auc) for one fold.

    The paper's Table 1 "Average AUC" matches the *per-model* mean AUC
    (10 models scored individually, then averaged), not the ensemble AUC
    (predictions averaged before one AUC). Both are computed here; the
    per-model mean is the paper-layer comparator.
    """
    if len(set(y_val.tolist())) < 2:
        return float("nan"), float("nan")
    probs = np.column_stack([m.predict_proba(X_val)[:, 1] for m in models])
    ensemble_auc = float(roc_auc_score(y_val, probs.mean(axis=1)))
    per_model_auc = float(
        np.mean([roc_auc_score(y_val, probs[:, i]) for i in range(probs.shape[1])])
    )
    return ensemble_auc, per_model_auc


def neg_sets(
    neg_pool: pd.DataFrame,
    pos_train_size: int,
    seed: int,
    *,
    exclude: pd.Index | None = None,
) -> list[pd.DataFrame]:
    """Draw the 10 paper-protocol negative training sets (2:1, seed-separated)."""
    base = neg_pool
    if exclude is not None:
        base = base[~base.index.isin(exclude)]
    pool = base.reset_index(drop=True)
    n_neg = min(len(pool), pos_train_size * NEG_RATIO)
    sets: list[pd.DataFrame] = []
    for i in range(N_MODELS):
        rng = np.random.default_rng(seed + i)
        idx = rng.choice(len(pool), size=n_neg, replace=False)
        sets.append(pool.iloc[idx])
    return sets


def run_cross_validation(
    X_pos: pd.DataFrame,
    X_neg_pool: pd.DataFrame,
    cols: list[str],
    seed: int,
) -> tuple[float, list[float], float, list[float]]:
    """RETRAIN_PROTOCOL 3.1: 5-fold x 5 rounds x 10 neg sets.

    Returns ``(mean_per_model_auc, per_model_fold_aucs, mean_ensemble_auc,
    ensemble_fold_aucs)``. ``mean_per_model_auc`` is the paper-comparable
    Table 1 metric; ``mean_ensemble_auc`` is what the shipped ensemble
    actually achieves at inference.
    """
    rng = np.random.default_rng(seed)
    indices = rng.permutation(len(X_pos))
    fold_splits = np.array_split(indices, N_FOLDS)

    fold_aucs: list[float] = []  # per-model mean (paper comparator)
    fold_ensembles: list[float] = []  # ensemble AUC
    for fold, val_idx in enumerate(fold_splits):
        train_idx = np.concatenate([s for f, s in enumerate(fold_splits) if f != fold])
        pos_train = X_pos.iloc[train_idx].reset_index(drop=True)
        pos_val = X_pos.iloc[val_idx].reset_index(drop=True)

        rng_val = np.random.default_rng(seed + 1000 + fold)
        n_val = min(len(X_neg_pool), len(pos_val) * NEG_RATIO)
        val_neg_idx = rng_val.choice(len(X_neg_pool), size=n_val, replace=False)
        neg_val = X_neg_pool.iloc[val_neg_idx].reset_index(drop=True)

        neg_sets = _neg_sets(
            X_neg_pool,
            len(pos_train),
            seed,
            exclude=X_neg_pool.index[val_neg_idx],
        )

        models: list[XGBClassifier] = []
        for m_idx, neg_set in enumerate(neg_sets):
            X = pd.concat([pos_train, neg_set], ignore_index=True)[cols]
            y = np.concatenate(
                [np.ones(len(pos_train), dtype=int), np.zeros(len(neg_set), dtype=int)]
            )
            models.append(_fit_model(X, y, seed + m_idx))

        X_val = pd.concat([pos_val, neg_val], ignore_index=True)[cols]
        y_val = np.concatenate(
            [np.ones(len(pos_val), dtype=int), np.zeros(len(neg_val), dtype=int)]
        )
        ens_auc, per_auc = _fold_auc(models, X_val, y_val)
        fold_ensembles.append(ens_auc)
        fold_aucs.append(per_auc)
        print(
            f"    fold {fold + 1}/{N_FOLDS}: pos_val={len(pos_val)} neg_val={len(neg_val)} "
            f"per_model_auc={per_auc:.4f} ensemble_auc={ens_auc:.4f}"
        )

    finite = [a for a in fold_aucs if a == a]
    mean_auc = float(np.mean(finite)) if finite else float("nan")
    finite_ens = [a for a in fold_ensembles if a == a]
    mean_ensemble = float(np.mean(finite_ens)) if finite_ens else float("nan")
    return mean_auc, fold_aucs, mean_ensemble, fold_ensembles


def train_final_ensemble(
    X_pos: pd.DataFrame,
    X_neg_pool: pd.DataFrame,
    cols: list[str],
    seed: int,
) -> list[XGBClassifier]:
    """RETRAIN_PROTOCOL 3.2: all positives + 10 neg sets -> 10 models."""
    neg_sets = _neg_sets(X_neg_pool, len(X_pos), seed)
    models: list[XGBClassifier] = []
    for m_idx, neg_set in enumerate(neg_sets):
        X = pd.concat([X_pos, neg_set], ignore_index=True)[cols]
        y = np.concatenate(
            [np.ones(len(X_pos), dtype=int), np.zeros(len(neg_set), dtype=int)]
        )
        models.append(_fit_model(X, y, seed + m_idx))
        print(f"    final model {m_idx + 1}/10: n_neg={len(neg_set)}")
    return models


def task_frame(task: str, split: str, *, root: Path | None = None) -> pd.DataFrame:
    """One task's rows from the committed training tables (E1)."""
    scope = "human" if task.startswith("h") else "base"
    return load_task_frame(scope, split, task, root=root)


def split_pos_neg(frame: pd.DataFrame, cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    pos = frame[frame["label"] == 1].reset_index(drop=True)
    neg_pool = frame[frame["label"] == 0].reset_index(drop=True)
    return pos[cols].astype(float), neg_pool[cols].astype(float)


def check_leakage(train_frame: pd.DataFrame, test_frame: pd.DataFrame) -> int:
    train_accs = set(train_frame["UniprotEntry"].astype(str))
    test_accs = set(test_frame["UniprotEntry"].astype(str))
    return len(train_accs & test_accs)


# Original call tokens, retained so the protocol bodies above stay verbatim.
_make_xgb = make_xgb
_fit_model = fit_model
_neg_sets = neg_sets
