#!/usr/bin/env python
"""Product B — train the 10-model XGBoost ensemble on the updated dataset.

Uses tuned hyperparameters from tune.py and the same train/test split that
tune.py sealed (recomputed deterministically from the random_state). Class
imbalance is handled with `scale_pos_weight` rather than 2:1 negative
subsampling, because the negative pool here is ~60K and subsampling
discards too much signal.

Anti-leakage:
  - Recompute the deterministic split with the same seed; verify it matches
    the sealed train_accessions_*.tsv (otherwise: abort).
  - External validation (PhaSePro / LLPSDB2) is filtered to drop any
    accession present in the task's training set. Removed counts are in
    leakage_report.json.

Outputs:
  products/B_extended_dataset/models/<task>/8f_model_{0..9}.joblib
  products/B_extended_dataset/metrics.json
  products/B_extended_dataset/leakage_report.json
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
from phasepred.updated_data import (  # noqa: E402
    build_label_table,
    build_task_split,
    load_llpsdb2_labels,
    load_phasepro_labels,
)

PRODUCT_DIR = REPO_ROOT / "products" / "B_extended_dataset"
RECOMPUTED_PATH = REPO_ROOT / "data" / "interim" / "recomputed_features_full.csv"

SEED = 42
N_MODELS = 10
TEST_SIZE = 0.20
TASKS = ("SaPS", "PdPS", "hSaPS", "hPdPS")


def _feature_columns(task: str) -> list[str]:
    return list(HUMAN_FEATURE_COLUMNS) if task.startswith("h") else list(BASE_FEATURE_COLUMNS)


def _join(split_df: pd.DataFrame, recomputed: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    merged = split_df.merge(recomputed, on="UniprotEntry", how="left")
    missing_all = merged[cols].isna().all(axis=1)
    return merged[~missing_all].drop_duplicates("UniprotEntry").reset_index(drop=True)


def _train_one(X: pd.DataFrame, y: np.ndarray, params: dict, seed: int) -> Pipeline:
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", XGBClassifier(
            eval_metric="logloss",
            n_jobs=8,
            **params,
        )),
    ])
    pipe.fit(X, y)
    return pipe


def _ensemble_score(models: list[Pipeline], X: pd.DataFrame) -> np.ndarray:
    probs = np.column_stack([m.predict_proba(X)[:, 1] for m in models])
    return probs.mean(axis=1)


def _filter_no_leak(ext_df: pd.DataFrame, train_accs: set[str]) -> tuple[pd.DataFrame, int]:
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

    print(f"Building updated label table...")
    label_table = build_label_table(REPO_ROOT)
    print(f"Loading recomputed features...")
    recomputed = pd.read_csv(RECOMPUTED_PATH)

    pp = load_phasepro_labels(REPO_ROOT)
    llps_pos, llps_neg = load_llpsdb2_labels(REPO_ROOT)

    metrics: dict[str, dict] = {}
    leakage: dict[str, dict] = {}

    for task in TASKS:
        if task not in tuned:
            raise KeyError(f"Tuned params missing for {task}")
        cols = _feature_columns(task)
        params = dict(tuned[task]["best_params"])
        params.setdefault("scale_pos_weight", tuned[task]["pos_weight"])

        train_split, test_split = build_task_split(
            label_table, RECOMPUTED_PATH,
            task=task, test_size=TEST_SIZE, random_state=SEED,
        )
        train_df = _join(train_split, recomputed, cols)
        test_df = _join(test_split, recomputed, cols)

        # Confirm the sealed accession files agree.
        sealed_train = set((PRODUCT_DIR / f"train_accessions_{task}.tsv").read_text().split())
        sealed_test = set((PRODUCT_DIR / f"test_accessions_{task}.tsv").read_text().split())
        now_train = set(train_df["UniprotEntry"])
        now_test = set(test_df["UniprotEntry"])
        if sealed_train != now_train or sealed_test != now_test:
            raise RuntimeError(
                f"[{task}] sealed accessions disagree with current split. "
                f"Re-run tune.py before train.py."
            )
        if now_train & now_test:
            raise RuntimeError(f"[{task}] leakage: {len(now_train & now_test)} accessions in both partitions")

        X_train = train_df[cols].astype(float)
        y_train = train_df["label"].astype(int).to_numpy()
        X_test = test_df[cols].astype(float)
        y_test = test_df["label"].astype(int).to_numpy()
        n_pos = int((y_train == 1).sum())
        n_neg = int((y_train == 0).sum())

        print()
        print(f"=== Training {task} ===")
        print(f"  Train: {len(X_train)} ({n_pos} pos / {n_neg} neg);  Test: {len(X_test)}")

        models_out = PRODUCT_DIR / "models" / task
        models_out.mkdir(parents=True, exist_ok=True)
        ensemble: list[Pipeline] = []
        for i in range(N_MODELS):
            # vary the random seed across models so the ensemble has diversity
            iter_params = dict(params)
            iter_params["random_state"] = SEED + i
            pipe = _train_one(X_train, y_train, iter_params, seed=SEED + i)
            joblib.dump(pipe, models_out / f"8f_model_{i}.joblib")
            ensemble.append(pipe)

        cv_auc = float(tuned[task]["best_cv_auc"])
        test_scores = _ensemble_score(ensemble, X_test)
        test_auc = float(roc_auc_score(y_test, test_scores)) if len(set(y_test)) > 1 else float("nan")
        print(f"  CV AUC (train, 5-fold): {cv_auc:.4f}")
        print(f"  TEST AUC (sealed 20%):  {test_auc:.4f}")

        # External validation with leakage filter.
        ext_auc: dict[str, float] = {}
        ext_leakage: dict[str, dict] = {}

        pp_match = pp[pp["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        pp_filtered, pp_removed = _filter_no_leak(pp_match, now_train)
        ext_leakage["phasepro_positives"] = {
            "total_with_recomputed": int(len(pp_match)),
            "removed_due_to_train_overlap": int(pp_removed),
            "kept": int(len(pp_filtered)),
        }
        if len(pp_filtered) > 0:
            neg_pool = test_df[test_df["label"] == 0]
            jp = pp_filtered.merge(recomputed, on="UniprotEntry", how="left")
            X_ext = pd.concat([jp[cols], neg_pool[cols]], ignore_index=True).astype(float)
            y_ext = np.concatenate([np.ones(len(jp), dtype=int), np.zeros(len(neg_pool), dtype=int)])
            if len(set(y_ext)) > 1:
                ext_scores = _ensemble_score(ensemble, X_ext)
                ext_auc["phasepro_vs_test_negatives"] = float(roc_auc_score(y_ext, ext_scores))

        llps_pos_match = llps_pos[llps_pos["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        llps_neg_match = llps_neg[llps_neg["UniprotEntry"].isin(recomputed["UniprotEntry"])].drop_duplicates("UniprotEntry")
        llps_pos_filtered, lp_removed = _filter_no_leak(llps_pos_match, now_train)
        llps_neg_filtered, ln_removed = _filter_no_leak(llps_neg_match, now_train)
        ext_leakage["llpsdb2_positives"] = {
            "total": int(len(llps_pos_match)),
            "removed_due_to_train_overlap": int(lp_removed),
            "kept": int(len(llps_pos_filtered)),
        }
        ext_leakage["llpsdb2_negatives"] = {
            "total": int(len(llps_neg_match)),
            "removed_due_to_train_overlap": int(ln_removed),
            "kept": int(len(llps_neg_filtered)),
        }
        if len(llps_pos_filtered) > 0 and len(llps_neg_filtered) > 0:
            jp = llps_pos_filtered.merge(recomputed, on="UniprotEntry", how="left")
            jn = llps_neg_filtered.merge(recomputed, on="UniprotEntry", how="left")
            X_ext = pd.concat([jp[cols], jn[cols]], ignore_index=True).astype(float)
            y_ext = np.concatenate([
                np.ones(len(jp), dtype=int),
                np.zeros(len(jn), dtype=int),
            ])
            if len(set(y_ext)) > 1:
                ext_scores = _ensemble_score(ensemble, X_ext)
                ext_auc["llpsdb2_pos_vs_neg"] = float(roc_auc_score(y_ext, ext_scores))

        metrics[task] = {
            "train_rows": int(len(X_train)),
            "test_rows": int(len(X_test)),
            "n_pos_train": n_pos,
            "n_neg_train": n_neg,
            "cv_auc": cv_auc,
            "test_auc_holdout20": test_auc,
            "external_auc": ext_auc,
        }
        leakage[task] = ext_leakage

    (PRODUCT_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2))
    (PRODUCT_DIR / "leakage_report.json").write_text(json.dumps(leakage, indent=2))

    print()
    print("=== Product B summary ===")
    print(f"{'task':<7} {'cv_auc':<8} {'test_auc':<10} {'external':<60}")
    for task in TASKS:
        m = metrics[task]
        ext = ", ".join(f"{k}={v:.3f}" for k, v in m["external_auc"].items())
        print(f"{task:<7} {m['cv_auc']:<8.4f} {m['test_auc_holdout20']:<10.4f} {ext}")
    print()
    print(f"Wrote {PRODUCT_DIR / 'metrics.json'}")
    print(f"Wrote {PRODUCT_DIR / 'leakage_report.json'}")
    print(f"Models in {PRODUCT_DIR / 'models'}/<task>/")


if __name__ == "__main__":
    main()
