#!/usr/bin/env python
"""Product A — train the 10-model XGBoost ensemble using tuned params.

Loads tuned_params.json (produced by tune.py), trains 10 models per task
with 2:1 negative subsampling (the paper protocol), evaluates the
ensemble on the sealed S3 test partition, and also reports leakage-safe
external AUCs on PhaSePro and LLPSDB2.

Anti-leakage:
  - Refuses to start if train_accessions_*.tsv and test_accessions_*.tsv
    intersect.
  - External validation sets (PhaSePro, LLPSDB2) are reduced to the
    accessions that do NOT appear in the training set for the task.
    Removed counts are recorded in leakage_report.json.

Outputs:
  products/A_paper_split_recomputed/models/<task>/8f_model_{0..9}.joblib
  products/A_paper_split_recomputed/metrics.json
  products/A_paper_split_recomputed/leakage_report.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from phasepred.features import BASE_FEATURE_COLUMNS, HUMAN_FEATURE_COLUMNS  # noqa: E402
from phasepred.paper import load_paper_split  # noqa: E402
from phasepred.updated_data import (  # noqa: E402
    load_llpsdb2_labels,
    load_phasepro_labels,
)

PRODUCT_DIR = REPO_ROOT / "products" / "A_paper_split_recomputed"
RECOMPUTED_PATH = REPO_ROOT / "data" / "interim" / "recomputed_features_full.csv"
S2_PATH = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd02.xlsx"
S3_PATH = REPO_ROOT / "PhaSePred_article&data" / "pnas.2115369119.sd03.xlsx"

SEED = 42
N_MODELS = 10
NEG_RATIO = 2
TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")


def _feature_columns(task: str) -> list[str]:
    return list(HUMAN_FEATURE_COLUMNS) if task.startswith("h") else list(BASE_FEATURE_COLUMNS)


def _join_features(split_df: pd.DataFrame, recomputed: pd.DataFrame, task: str) -> pd.DataFrame:
    task_rows = split_df.query("task == @task").copy()
    merged = task_rows.merge(recomputed, on="UniprotEntry", how="left")
    cols = _feature_columns(task)
    missing_all = merged[cols].isna().all(axis=1)
    kept = merged[~missing_all].drop_duplicates("UniprotEntry").reset_index(drop=True)
    return kept


def _train_one_subsample(
    X_pos: pd.DataFrame,
    X_neg: pd.DataFrame,
    params: dict,
    seed: int,
) -> Pipeline:
    """Train one XGB pipeline on positives + a 2:1 sample of negatives."""
    rng = np.random.default_rng(seed)
    n_neg = min(len(X_pos) * NEG_RATIO, len(X_neg))
    idx = rng.choice(len(X_neg), size=n_neg, replace=False)
    X = pd.concat([X_pos, X_neg.iloc[idx]], ignore_index=True)
    y = np.concatenate([np.ones(len(X_pos), dtype=int), np.zeros(n_neg, dtype=int)])
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", XGBClassifier(
            random_state=seed,
            eval_metric="logloss",
            n_jobs=4,  # cap so we share cores with concurrent B tune
            **params,
        )),
    ])
    pipe.fit(X, y)
    return pipe


def _ensemble_score(models: list[Pipeline], X: pd.DataFrame) -> np.ndarray:
    """Mean probability across the ensemble."""
    probs = np.column_stack([m.predict_proba(X)[:, 1] for m in models])
    return probs.mean(axis=1)


def _external_labels() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return positive PhaSePro / LLPSDB2 positives / LLPSDB2 negatives."""
    pp = load_phasepro_labels(REPO_ROOT)
    llps_pos, llps_neg = load_llpsdb2_labels(REPO_ROOT)
    return pp, llps_pos, llps_neg


def _filter_no_leak(ext_df: pd.DataFrame, train_accs: set[str]) -> tuple[pd.DataFrame, int]:
    """Drop accessions present in train; return filtered + removed count."""
    mask = ~ext_df["UniprotEntry"].isin(train_accs)
    return ext_df[mask].reset_index(drop=True), int((~mask).sum())


def main() -> None:
    tuned_path = PRODUCT_DIR / "tuned_params.json"
    if not tuned_path.exists():
        raise FileNotFoundError(
            f"tuned_params.json not found. Run tune.py first:\n"
            f"  uv run python {PRODUCT_DIR.relative_to(REPO_ROOT)}/tune.py"
        )
    tuned = json.loads(tuned_path.read_text())

    print(f"Loading paper split + recomputed features...")
    split = load_paper_split(S2_PATH, S3_PATH)
    recomputed = pd.read_csv(RECOMPUTED_PATH)

    pp, llps_pos, llps_neg = _external_labels()
    pp_acc = set(pp["UniprotEntry"]) if "UniprotEntry" in pp.columns else set()
    llps_pos_acc = set(llps_pos["UniprotEntry"]) if "UniprotEntry" in llps_pos.columns else set()
    llps_neg_acc = set(llps_neg["UniprotEntry"]) if "UniprotEntry" in llps_neg.columns else set()

    metrics: dict[str, dict] = {}
    leakage: dict[str, dict] = {}

    for task in TASKS:
        if task not in tuned:
            raise KeyError(f"Tuned params missing for {task}")
        cols = _feature_columns(task)
        params = dict(tuned[task]["best_params"])
        pos_weight = tuned[task]["pos_weight"]
        params.setdefault("scale_pos_weight", pos_weight)

        train_df = _join_features(split.train, recomputed, task)
        test_df = _join_features(split.test, recomputed, task)

        # Sealed-set leakage check.
        train_accs = set(train_df["UniprotEntry"])
        test_accs = set(test_df["UniprotEntry"])
        if train_accs & test_accs:
            raise RuntimeError(
                f"Leakage: {len(train_accs & test_accs)} accessions in both "
                f"train and test for {task}. Refusing to train."
            )

        sealed_train = (PRODUCT_DIR / f"train_accessions_{task}.tsv").read_text().split()
        if set(sealed_train) != train_accs:
            print(
                f"  WARNING [{task}]: sealed train accessions disagree with "
                f"current join — sealed={len(sealed_train)} now={len(train_accs)}. "
                f"Continuing with current join, but re-running tune.py is "
                f"recommended."
            )

        X_train = train_df[cols].astype(float)
        y_train = train_df["label"].astype(int).to_numpy()
        X_test = test_df[cols].astype(float)
        y_test = test_df["label"].astype(int).to_numpy()

        pos_mask = y_train == 1
        X_pos = X_train[pos_mask].reset_index(drop=True)
        X_neg = X_train[~pos_mask].reset_index(drop=True)
        n_pos = len(X_pos)
        n_neg = len(X_neg)
        print()
        print(f"=== Training {task} ===")
        print(f"  train pos/neg: {n_pos}/{n_neg}  test rows: {len(X_test)}")

        models_out = PRODUCT_DIR / "models" / task
        models_out.mkdir(parents=True, exist_ok=True)
        ensemble: list[Pipeline] = []
        for i in range(N_MODELS):
            pipe = _train_one_subsample(X_pos, X_neg, params, seed=SEED + i)
            joblib.dump(pipe, models_out / f"8f_model_{i}.joblib")
            ensemble.append(pipe)

        # CV AUC was already computed by tune.py; report it here too for clarity.
        cv_auc = float(tuned[task]["best_cv_auc"])
        # Held-out test AUC (S3).
        test_scores = _ensemble_score(ensemble, X_test)
        test_auc = float(roc_auc_score(y_test, test_scores)) if len(set(y_test)) > 1 else float("nan")
        print(f"  CV AUC (S2, 5-fold): {cv_auc:.4f}")
        print(f"  TEST AUC (S3):       {test_auc:.4f}")

        # External validation with leakage filter.
        is_human = task.startswith("h")
        ext_auc: dict[str, float] = {}
        ext_leakage: dict[str, dict] = {}

        # PhaSePro positives (vs S3 negatives as ad-hoc background).
        pp_match = pp[pp["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        pp_filtered, pp_removed = _filter_no_leak(pp_match, train_accs)
        ext_leakage["phasepro_positives_removed"] = {
            "total_phasepro_with_recomputed": int(len(pp_match)),
            "removed_due_to_train_overlap": int(pp_removed),
            "kept": int(len(pp_filtered)),
        }
        if len(pp_filtered) > 0:
            neg_pool = test_df[test_df["label"] == 0]
            joined_pos = pp_filtered.merge(recomputed, on="UniprotEntry", how="left").dropna(subset=cols, how="all")
            X_ext = pd.concat([joined_pos[cols], neg_pool[cols]], ignore_index=True).astype(float)
            y_ext = np.concatenate([np.ones(len(joined_pos), dtype=int), np.zeros(len(neg_pool), dtype=int)])
            if len(set(y_ext)) > 1:
                ext_scores = _ensemble_score(ensemble, X_ext)
                ext_auc["phasepro_vs_S3_negatives"] = float(roc_auc_score(y_ext, ext_scores))

        # LLPSDB2: positives + confirmed negatives, both leakage-filtered.
        llps_pos_match = llps_pos[llps_pos["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        llps_neg_match = llps_neg[llps_neg["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        llps_pos_filtered, lp_removed = _filter_no_leak(llps_pos_match, train_accs)
        llps_neg_filtered, ln_removed = _filter_no_leak(llps_neg_match, train_accs)
        ext_leakage["llpsdb2_positives_removed"] = {
            "total": int(len(llps_pos_match)),
            "removed_due_to_train_overlap": int(lp_removed),
            "kept": int(len(llps_pos_filtered)),
        }
        ext_leakage["llpsdb2_negatives_removed"] = {
            "total": int(len(llps_neg_match)),
            "removed_due_to_train_overlap": int(ln_removed),
            "kept": int(len(llps_neg_filtered)),
        }
        if len(llps_pos_filtered) > 0 and len(llps_neg_filtered) > 0:
            joined_pos = llps_pos_filtered.merge(recomputed, on="UniprotEntry", how="left")
            joined_neg = llps_neg_filtered.merge(recomputed, on="UniprotEntry", how="left")
            X_ext = pd.concat([joined_pos[cols], joined_neg[cols]], ignore_index=True).astype(float)
            y_ext = np.concatenate([
                np.ones(len(joined_pos), dtype=int),
                np.zeros(len(joined_neg), dtype=int),
            ])
            if len(set(y_ext)) > 1:
                ext_scores = _ensemble_score(ensemble, X_ext)
                ext_auc["llpsdb2_pos_vs_neg"] = float(roc_auc_score(y_ext, ext_scores))

        metrics[task] = {
            "train_rows": int(len(X_train)),
            "test_rows": int(len(X_test)),
            "n_pos_train": int(n_pos),
            "n_neg_train": int(n_neg),
            "cv_auc": cv_auc,
            "test_auc_S3": test_auc,
            "external_auc": ext_auc,
        }
        leakage[task] = ext_leakage

    (PRODUCT_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2))
    (PRODUCT_DIR / "leakage_report.json").write_text(json.dumps(leakage, indent=2))

    print()
    print("=== Product A summary ===")
    print(f"{'task':<7} {'cv_auc':<8} {'test_auc_S3':<14} {'external':<60}")
    for task in TASKS:
        m = metrics[task]
        ext = ", ".join(f"{k}={v:.3f}" for k, v in m["external_auc"].items())
        print(f"{task:<7} {m['cv_auc']:<8.4f} {m['test_auc_S3']:<14.4f} {ext}")
    print()
    print(f"Wrote {PRODUCT_DIR / 'metrics.json'}")
    print(f"Wrote {PRODUCT_DIR / 'leakage_report.json'}")
    print(f"Models in {PRODUCT_DIR / 'models'}/<task>/8f_model_{{0..9}}.joblib")


if __name__ == "__main__":
    main()
